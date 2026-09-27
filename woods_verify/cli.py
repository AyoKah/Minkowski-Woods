"""Command-line interface.  Run `python -m woods_verify --help`.

This lives in its own module (not in __main__.py) because worker processes
started with the "spawn" method -- the default on macOS and Windows -- cannot
import functions defined in a package's __main__.py.
"""
import argparse
import json
import multiprocessing
import os
import sys
import time

from . import core

CERT_KEYS = ("nodes_explored", "max_depth", "pruned_infeasible", "pruned_bound",
             "pruned_propagation", "surviving_count")


# ----------------------------------------------------------------- run
def _parse_roots(spec, total):
    if spec is None:
        return list(range(total))
    out = []
    for part in spec.split(","):
        if ":" in part:
            a, b = part.split(":")
            out.extend(range(int(a) if a else 0, int(b) if b else total))
        else:
            out.append(int(part))
    bad = [k for k in out if not 0 <= k < total]
    if bad:
        sys.exit(f"root indices out of range 0..{total - 1}: {bad[:5]}")
    return sorted(set(out))


def _worker(args):
    n, idx, L, U, cert_dir, surv_dir = args
    cert, surv = core.run_subtree(n, L, U)
    cert = {"box_idx": idx, **cert}
    # Write survivors first, then the certificate: a certificate on disk
    # therefore always means the subtree finished (used by --resume).
    if surv:
        with open(os.path.join(surv_dir, f"subtree_{idx}.json"), "w") as f:
            json.dump(surv, f)
    tmp = os.path.join(cert_dir, f".subtree_{idx}.json.tmp")
    with open(tmp, "w") as f:
        json.dump(cert, f)
    os.replace(tmp, os.path.join(cert_dir, f"subtree_{idx}.json"))
    return cert


def cmd_run(a):
    n = a.n
    if not 2 <= n <= 11:
        print(f"warning: the analytic part of the proof and Lemma 3.1's initial box are only "
              f"established for n <= 11 (n = {n} requested)", file=sys.stderr)
    base = os.path.join(a.out, f"n{n}")
    cert_dir = os.path.join(base, "certificates")
    surv_dir = os.path.join(base, "survivors")
    os.makedirs(cert_dir, exist_ok=True)
    os.makedirs(surv_dir, exist_ok=True)

    roots = core.root_boxes(n)
    todo = _parse_roots(a.roots, len(roots))
    done = {int(f[8:-5]) for f in os.listdir(cert_dir) if f.startswith("subtree_")}
    todo = [k for k in todo if k not in done]
    print(f"n={n}: {len(roots)} root boxes; {len(todo)} to run "
          f"({len(done)} already have certificates); workers={a.workers}", flush=True)

    tasks = [(n, k, roots[k][1], roots[k][2], cert_dir, surv_dir) for k in todo]
    t0 = time.time()
    finished = 0
    # Use "spawn" on every platform so that Linux, macOS and Windows behave identically.
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(a.workers) as pool:
        for cert in pool.imap_unordered(_worker, tasks):
            finished += 1
            if a.verbose or finished == len(tasks) or finished % max(1, len(tasks) // 20) == 0:
                print(f"  [{finished}/{len(tasks)}] subtree {cert['box_idx']}: "
                      f"{cert['nodes_explored']} nodes, {cert['surviving_count']} survivors "
                      f"({time.time() - t0:.0f}s elapsed)", flush=True)
    with open(os.path.join(base, "run_log.jsonl"), "a") as f:
        f.write(json.dumps({"n": n, "roots": a.roots, "subtrees_run": len(tasks),
                            "workers": a.workers, "wall_clock_seconds": time.time() - t0,
                            "mpmath": core.mpmath.__version__,
                            "python": sys.version.split()[0]}) + "\n")
    print(f"done in {time.time() - t0:.1f}s. Next: python -m woods_verify summarize --n {n} --out {a.out}")


# ----------------------------------------------------------- summarize
def _load_certs(cert_dir):
    certs = {}
    for f in os.listdir(cert_dir):
        if f.startswith("subtree_") and f.endswith(".json"):
            with open(os.path.join(cert_dir, f)) as fh:
                c = json.load(fh)
            certs[int(c["box_idx"])] = c
    return certs


def _load_survivors(base):
    surv = []
    d = os.path.join(base, "survivors")
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.endswith(".json"):
                with open(os.path.join(d, f)) as fh:
                    surv.extend(json.load(fh))
    legacy = os.path.join(base, "datasets", "surviving_boxes.json")   # version-1 layout
    if not surv and os.path.exists(legacy):
        with open(legacy) as fh:
            surv = json.load(fh)
    return surv


def cmd_summarize(a):
    n = a.n
    base = os.path.join(a.out, f"n{n}")
    certs = _load_certs(os.path.join(base, "certificates"))
    n_roots = len(core.root_boxes(n))
    tot = {k: sum(c[k] for c in certs.values()) for k in CERT_KEYS if k != "max_depth"}
    depth = max((c["max_depth"] for c in certs.values()), default=0)
    bisections = tot["nodes_explored"] - tot["pruned_bound"] - tot["pruned_infeasible"] - tot["surviving_count"]
    audit = tot["nodes_explored"] == len(certs) + 2 * bisections - tot["pruned_propagation"]

    surv = _load_survivors(base)
    lo = min((min(L) for L, U in surv), default=None)
    hi = max((max(U) for L, U in surv), default=None)
    inside = all(min(L) >= float(core.C_LO) and max(U) <= float(core.C_HI) and
                 all(l <= u for l, u in zip(L, U)) for L, U in surv)

    complete = sorted(certs) == list(range(n_roots))
    ok = complete and audit and inside and len(surv) == tot["surviving_count"]
    summary = {
        "n": n, "root_boxes": n_roots, "subtrees_with_certificate": len(certs),
        "complete": complete, "nodes": tot["nodes_explored"], "max_depth": depth,
        "pruned_by_bound": tot["pruned_bound"], "empty_after_reduction": tot["pruned_propagation"],
        "infeasible": tot["pruned_infeasible"], "survivors": tot["surviving_count"],
        "survivor_files_count": len(surv), "survivor_hull": [lo, hi],
        "audit_identity_holds": audit, "all_survivors_inside_C": inside,
        "worker_cpu_seconds": sum(c.get("worker_cpu_seconds", 0) for c in certs.values()),
        "verified": ok,
    }
    with open(os.path.join(base, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    for k, v in summary.items():
        print(f"  {k:28s} {v}")
    if ok:
        print(f"\nVERIFIED: sup B_DP <= {n} on every box of Omega_{n} outside C = [0.95,1.05]^{n}.\n"
              f"Together with Proposition 5.1 (inside C) this gives B_DP <= {n} on Omega_{n}.")
    elif not complete:
        print(f"\nPARTIAL: {len(certs)}/{n_roots} subtrees present; missing e.g. "
              f"{sorted(set(range(n_roots)) - set(certs))[:10]}")
    else:
        print("\nFAILED: see fields above.")
        sys.exit(1)


# ------------------------------------------------------------- compare
def cmd_compare(a):
    mine = _load_certs(a.certs)
    ref = _load_certs(a.ref)
    common = sorted(set(mine) & set(ref))
    if not common:
        sys.exit("no subtrees in common")
    bad = [k for k in common if any(mine[k][f] != ref[k][f] for f in CERT_KEYS)]
    for k in bad[:20]:
        print(f"  MISMATCH subtree {k}: mine={[mine[k][f] for f in CERT_KEYS]} "
              f"ref={[ref[k][f] for f in CERT_KEYS]}")
    print(f"{len(common) - len(bad)}/{len(common)} subtrees identical in "
          f"{', '.join(CERT_KEYS)}")
    sys.exit(1 if bad else 0)


# -------------------------------------------------------- other checks
def cmd_selftest(a):
    from . import checks
    sys.exit(0 if checks.selftest() else 1)


def cmd_witnesses(a):
    from . import checks
    sys.exit(0 if checks.witnesses() else 1)


def cmd_local(a):
    from . import checks
    sys.exit(0 if checks.local_lemma() else 1)


def main():
    p = argparse.ArgumentParser(prog="python -m woods_verify",
                                description="Certified verification of B_DP(X) <= n on Omega_n "
                                            "(Woods' conjecture, blocks of size 1 and 2).")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run Algorithm 1 on some or all root boxes")
    r.add_argument("--n", type=int, required=True)
    r.add_argument("--roots", default=None,
                   help="root indices, e.g. '0:100' or '3,17,900:910' (default: all). "
                        "Subtrees that already have a certificate are skipped, so runs can be resumed.")
    r.add_argument("--workers", type=int, default=os.cpu_count())
    r.add_argument("--out", default="results")
    r.add_argument("--verbose", action="store_true")
    r.set_defaults(f=cmd_run)

    s = sub.add_parser("summarize", help="aggregate certificates and check the run")
    s.add_argument("--n", type=int, required=True)
    s.add_argument("--out", default="results")
    s.set_defaults(f=cmd_summarize)

    c = sub.add_parser("compare", help="compare two certificate directories subtree by subtree")
    c.add_argument("certs")
    c.add_argument("ref")
    c.set_defaults(f=cmd_compare)

    sub.add_parser("selftest", help="unit checks + full runs for n=4,5 (about a minute)").set_defaults(f=cmd_selftest)
    sub.add_parser("witnesses", help="check the n=9 and n=12 witnesses of Remark 6.2 and Section 7").set_defaults(f=cmd_witnesses)
    sub.add_parser("local-lemma", help="check the numerical condition in Proposition 5.1").set_defaults(f=cmd_local)

    a = p.parse_args()
    a.f(a)

