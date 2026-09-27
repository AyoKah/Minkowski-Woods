"""Independent checks of the auxiliary claims in the paper."""
import math
import random
from fractions import Fraction as F

import mpmath
from mpmath import iv

from . import core


def _ok(flag, msg):
    print(("  PASS  " if flag else "  FAIL  ") + msg)
    return flag


# ------------------------------------------------------------------ helpers
def b2_float(u, v):
    return u + v if u <= v else 2 * u / (1 + math.sqrt(1 - v / u))


def covering_radius_sq_2d(u, v, mu, grid=160):
    """4*rho^2 of the lattice with basis (sqrt u, 0), (mu sqrt u, sqrt v), by brute force."""
    b1 = (math.sqrt(u), 0.0)
    b2 = (mu * math.sqrt(u), math.sqrt(v))
    pts = [(i * b1[0] + j * b2[0], i * b1[1] + j * b2[1]) for i in range(-2, 4) for j in range(-2, 4)]
    best = 0.0
    for a in range(grid + 1):
        for b in range(grid + 1):
            s, t = a / grid, b / grid
            x, y = s * b1[0] + t * b2[0], s * b1[1] + t * b2[1]
            best = max(best, min((x - p) ** 2 + (y - q) ** 2 for p, q in pts))
    return 4 * best


def bdp_float(X):
    D = [0.0, X[0]]
    for i in range(2, len(X) + 1):
        D.append(min(D[i - 1] + X[i - 1], D[i - 2] + b2_float(X[i - 2], X[i - 1])))
    return D[-1]


# ------------------------------------------------------------------ selftest
def selftest():
    rng = random.Random(1)
    ok = True
    print("Proposition 3.2: B_2 bounds the covering radius of reduced 2D lattices (brute force)")
    worst = -1.0
    for _ in range(40):
        v = rng.uniform(0.75, 1.5)
        mu = rng.uniform(0, 0.5)
        if mu * mu + v < 1:
            mu = math.sqrt(1 - v)          # KZ reduced: |b2| >= |b1|
        worst = max(worst, covering_radius_sq_2d(1.0, v, mu) - b2_float(1.0, v))
    ok &= _ok(worst <= 1e-9, f"max(4rho^2 - B_2) over 40 random lattices = {worst:.2e} (<= 0)")
    gap = b2_float(1.0, 0.91) - covering_radius_sq_2d(1.0, 0.91, 0.3, grid=400)
    ok &= _ok(0 <= gap < 5e-3, f"bound attained at mu = sqrt(1-v/u): gap {gap:.1e} (grid resolution)")

    print("Lemma 4.1: interval upper bound of B_2 dominates dense sampling on random boxes")
    viol = 0
    for _ in range(300):
        ua = rng.uniform(0.3, 3.0); ub = ua + rng.uniform(0, 1.0)
        va = rng.uniform(0.3, 3.0); vb = va + rng.uniform(0, 1.0)
        sup = core.B2_interval(mpmath.mpi(ua, ub), mpmath.mpi(va, vb)).b
        m = max(b2_float(ua + (ub - ua) * i / 40, va + (vb - va) * j / 40)
                for i in range(41) for j in range(41))
        viol += m > float(sup) * (1 + 1e-12)
    ok &= _ok(viol == 0, f"{viol} violations in 300 boxes")

    print("Algorithm 1: full runs for n = 4 and 5 reproduce the published counts")
    expected = {4: (8, 131, 13, 63, 1, 6), 5: (16, 573, 17, 275, 5, 17)}
    for n, exp in expected.items():
        roots = core.root_boxes(n)
        tot = [0, 0, 0, 0, 0]
        for k, L, U in roots:
            c, _ = core.run_subtree(n, L, U)
            tot[0] += c["nodes_explored"]; tot[1] = max(tot[1], c["max_depth"])
            tot[2] += c["pruned_bound"]; tot[3] += c["pruned_propagation"]; tot[4] += c["surviving_count"]
        got = (len(roots), *tot)
        ok &= _ok(got == exp, f"n={n}: (roots, nodes, depth, pruned, empty, survivors) = {got}")
    print("ALL PASSED" if ok else "SOME CHECKS FAILED")
    return ok


# ------------------------------------------------------------------ witnesses
def _in_omega(X, pvz=False):
    n = len(X)
    ok = all(x <= X[0] for x in X)
    ok &= all(X[i + 1] >= F(3, 4) * X[i] for i in range(n - 1))
    ok &= all(X[i + 2] >= F(2, 3) * X[i] for i in range(n - 2))
    ok &= all(X[i + 3] >= F(1, 2) * X[i] for i in range(n - 3))
    if pvz:
        ok &= all(X[i + 4] >= F(46873, 100000) * X[i] for i in range(n - 4))
    return ok     # prod X_i = 1 is imposed by scaling, which preserves all of the above


def _bdp_interval_normalized(Xs):
    """Enclosure of B_DP(s X) with s = (prod X)^(-1/n), in 200-bit interval arithmetic."""
    X = [F(x) for x in Xs]
    n = len(X)
    iv.prec = 200
    P = F(1)
    for x in X:
        P *= x
    s = (iv.mpf(P.numerator) / iv.mpf(P.denominator)) ** (iv.mpf(-1) / n)
    Y = [iv.mpf(x.numerator) / iv.mpf(x.denominator) * s for x in X]
    D = [iv.mpf(0), Y[0]]
    for i in range(2, n + 1):
        u, v = Y[i - 2], Y[i - 1]
        b2 = u + v if X[i - 2] <= X[i - 1] else 2 * u / (1 + iv.sqrt(1 - v / u))
        a, b = D[i - 1] + Y[i - 1], D[i - 2] + b2
        D.append(iv.mpf([min(a.a, b.a), min(a.b, b.b)]))
    iv.prec = 53
    return D[n]


def witnesses():
    ok = True
    print("Section 7: points of Omega_12 with B_DP > 12")
    w1 = ["2.3634", "2.3634", "2.1532", "1.6149", "1.444", "1.0836", "0.9627", "0.7247",
          "0.6442", "0.5007", "0.4357", "0.3356"]
    w2 = ["2.3143", "2.3143", "2.0561", "1.5438", "1.4271", "1.0848", "0.9643", "0.7295",
          "0.6702", "0.5085", "0.452", "0.3506"]
    for w, pvz in ((w1, False), (w2, True)):
        inside = _in_omega([F(x) for x in w], pvz)
        enc = _bdp_interval_normalized(w)
        ok &= _ok(inside and enc.a > 12,
                  f"{'with' if pvz else 'without'} X_(i+4) >= 0.46873 X_i: in domain = {inside}, "
                  f"B_DP in {iv.nstr(enc, 10)}")

    print("Remark 6.2: with Woods' 2D bound 4u - 2u^2/v, the size-1/size-2 bound exceeds 9 on Omega_9")
    w = ["922523/500000", "918483/500000", "77567/50000", "31043/25000", "517541/500000",
         "165697/200000", "690679/1000000", "138193/250000", "467821/1000000"]
    X = [F(x) for x in w]
    n = len(X)

    def W(u, v):
        return u + v if u <= v else min(u + v, 4 * u - 2 * u * u / v)
    D = [F(0), X[0]]
    for i in range(2, n + 1):
        D.append(min(D[i - 1] + X[i - 1], D[i - 2] + W(X[i - 2], X[i - 1])))
    P = F(1)
    for x in X:
        P *= x
    # the bound is homogeneous of degree 1, so B(sX) > n  <=>  B(X)^n > n^n * prod(X)
    exceeds = D[n] ** n > F(n) ** n * P
    ok &= _ok(_in_omega(X) and exceeds,
              f"in domain = {_in_omega(X)}, exact test B > 9: {exceeds} "
              f"(value {float(D[n]) * float(P) ** (-1 / n):.5f})")
    print("ALL PASSED" if ok else "SOME CHECKS FAILED")
    return ok


# ------------------------------------------------------------------ Proposition 5.1
def local_lemma():
    """The final condition of Proposition 5.1: (2 eta)^(3/2) < 2 (1-eta)^(3/2) / (n sqrt(n-1))."""
    ok = True
    eta = iv.mpf(1) / 20
    lhs = (2 * eta) ** (iv.mpf(3) / 2)
    print("Proposition 5.1, eta = 0.05")
    for n in range(2, 16):
        rhs = 2 * (1 - eta) ** (iv.mpf(3) / 2) / (n * iv.sqrt(iv.mpf(n - 1)))
        ok &= _ok(lhs.b < rhs.a, f"n={n:2d}: 0.1^(3/2) = {float(lhs.b):.6f} < {float(rhs.a):.6f}")
    print("Sanity check by sampling Omega_11 near X = 1 (not part of the proof)")
    rng = random.Random(7)
    worst, count = -1e9, 0
    while count < 20000:
        z = [rng.uniform(-0.05, 0.05) for _ in range(11)]
        m = sum(math.log1p(t) for t in z) / 11
        X = [(1 + t) / math.exp(m) for t in z]
        if not all(0.95 <= x <= 1.05 for x in X) or max(X) != X[0]:
            continue
        if any(X[i + 1] < 0.75 * X[i] for i in range(10)):
            continue
        count += 1
        worst = max(worst, bdp_float(X) - 11)
    ok &= _ok(worst < 0, f"max(B_DP - 11) over {count} samples = {worst:.4f}")
    print("ALL PASSED" if ok else "SOME CHECKS FAILED")
    return ok
