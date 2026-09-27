# woods-verify — certified verification of Woods' conjecture for n ≤ 11

Companion code for

> A. Kahlon, *A certified computational proof of the conjectures of Woods and Minkowski in dimension 11* (submitted).

The program proves, for a given dimension `n`, the inequality

  **B_DP(X) ≤ n for every X ∈ Ω_n outside the cube C = [0.95, 1.05]^n**   (paper, Proposition 4.2)

using interval branch and bound with outward rounding. Together with the analytic Proposition 5.1 for the cube `C`, this gives `B_DP ≤ n` on all of `Ω_n`. By Theorem 3.4 of the paper, that implies Woods' conjecture in dimension `n`.

## Requirements

* Python ≥ 3.8
* `mpmath` 1.3.0 (`pip install -r requirements.txt`)

Nothing else is needed. Run all commands from this directory. Worker processes are always started with Python's `spawn` method, so the program behaves the same on Linux, macOS and Windows.

## Five-minute check

```bash
pip install -r requirements.txt
python -m woods_verify selftest      # ~10 s: unit checks of Prop. 3.2 and Lemma 4.1; full runs n=4,5
python -m woods_verify witnesses     # ~1 s: the n=12 obstruction (Sec. 7) and Remark 6.2 (n=9)
python -m woods_verify local-lemma   # ~5 s: the numerical condition in Prop. 5.1 for 2 ≤ n ≤ 15
python -m woods_verify summarize --n 11   # audits the archived n = 11 run in results/n11
```

`summarize` checks the following:

* every root box has a certificate;
* the tree-count identity `nodes = roots + 2·bisections − empty` holds;
* no box was discarded as infeasible without being counted;
* every surviving box lies inside `C`.

It prints `VERIFIED` only if all of these hold.

## Re-proving a dimension from scratch

```bash
python -m woods_verify run --n 8 --out myrun        # all 2^(n-1) root boxes, all cores
python -m woods_verify summarize --n 8 --out myrun
python -m woods_verify compare myrun/n8/certificates results/n8/certificates
```

Measured on an otherwise idle 2-core cloud VM (Python 3.11), all cores:

| n | root boxes | nodes | time |
|---|---:|---:|---:|
| 4 | 8 | 131 | 0.2 s |
| 5 | 16 | 573 | 0.6 s |
| 6 | 32 | 2,559 | 3 s |
| 7 | 64 | 13,644 | 19 s |
| 8 | 128 | 84,768 | 2.5 min |
| 9 | 256 | 595,547 | 19 min |
| 10 | 512 | 4,536,259 | ≈ 2.2 h |
| 11 | 1,024 | 39,991,484 | 268,715 s wall-clock on a 1.6 GHz dual-core laptop (961,275 s total worker CPU) |

## Replicating n = 11

The full `n = 11` run is long, but each of the 1,024 root boxes is an independent subtree with its own certificate in `results/n11/certificates/`. Any subset can therefore be re-run and compared exactly against the archive:

```bash
python -m woods_verify run --n 11 --roots $(paste -sd, data/n11_sample_roots.txt) --out myrun
python -m woods_verify compare myrun/n11/certificates results/n11/certificates
# -> 8/8 subtrees identical in nodes_explored, max_depth, pruned_infeasible, ...
```

The sample in `data/n11_sample_roots.txt` takes under a minute on a current machine. You can choose any other indices with `--roots`, e.g. `--roots 0:100` or `--roots 961` (the largest subtree, 3.3 million nodes).

A full run can be split across machines with disjoint `--roots` ranges. It can also be resumed: subtrees that already have a certificate are skipped. After copying all certificates and survivors into one directory, run `summarize`.

## What is trusted

* **Interval arithmetic.** `mpmath`'s interval type (`mpmath.iv`) is used at 53-bit precision. Every `+ − × ÷ √` returns an interval with outward-rounded endpoints.
* **The Python code in `woods_verify/core.py`** (≈ 280 lines). It maps to the paper as follows:
  * `B2_interval` — Lemma 4.1
  * `ub_dp` — the recursion for the upper bound of `B_DP` in §4.1
  * `propagate_bounds` and `is_feasible` — the domain reduction of §4.2
  * `run_subtree` — Algorithm 1
* **The analytic results of the paper.** These are Theorem 3.4 (4ρ² ≤ B_DP), Proposition 4.2 (soundness) and Proposition 5.1 (the cube `C`).

Floating-point values are used in only two places:

* **Bisection midpoints.** Both children share the same midpoint, so the two halves always cover the parent box.
* **Stored survivor coordinates.** These are 53-bit values that convert to `float` exactly.

## Layout

```
woods_verify/core.py        numerical core (frozen; identical to v1 except for the dimension parameter)
woods_verify/cli.py         command line: run / summarize / compare / selftest / witnesses / local-lemma
woods_verify/__main__.py    entry point for `python -m woods_verify`
woods_verify/checks.py      independent checks of auxiliary claims in the paper
data/n11_sample_roots.txt   a sample of n = 11 root boxes for quick replication
results/n{4..10}/           complete runs with this version
results/n11/                the archived version-1 run (certificates, surviving boxes, telemetry, log)
```

Output of `run --n N --out DIR`:

* `DIR/nN/certificates/subtree_K.json` — the certificate for root box `K`: nodes, depth, pruning counts, CPU time.
* `DIR/nN/survivors/subtree_K.json` — the boxes of subtree `K` that reached `C`.
* `DIR/nN/run_log.jsonl` — one line per invocation.

`summarize` writes `DIR/nN/summary.json`.


## License

MIT (see `LICENSE`).
