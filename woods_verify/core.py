"""Certified interval branch-and-bound for the inequality B_DP(X) <= n on Omega_n.

The numerical routines in this module (_eval_b2_point, B2_interval, ub_dp,
propagate_bounds, is_feasible, and the loop in run_subtree) are the ones used
for the n = 11 computation reported in the paper (archived as verify_m11.py in
version 1 of this repository). Only the dimension has been made a parameter and
input/output has been moved out of the numerical code. Any change here
invalidates the comparison against the archived per-subtree certificates, so
treat this file as frozen.

Paper references: B2_interval implements Lemma 4.1; ub_dp the enclosure of
B_DP (Section 4.1); propagate_bounds the domain reduction (Section 4.2);
run_subtree Algorithm 1 on one root box.
"""
import itertools
import time

import mpmath
from mpmath import iv

# 53-bit (double precision) mantissas; mpmath intervals round outward.
mpmath.mp.dps = 15
iv.dps = 15

# Exclusion cube C = [0.95, 1.05]^n (Proposition 5.1 covers it analytically).
C_LO = mpmath.mpf('0.95')
C_HI = mpmath.mpf('1.05')


def _eval_b2_point(u_val, v_val):
    """Interval enclosure of 2u/(1+sqrt(1-v/u)) at a point (u > v)."""
    u_iv = mpmath.mpi(u_val, u_val)
    v_iv = mpmath.mpi(v_val, v_val)
    one = mpmath.mpi('1.0', '1.0')
    two = mpmath.mpi('2.0', '2.0')
    inner = one - (v_iv / u_iv)
    if inner.a < 0:
        inner = mpmath.mpi('0.0', inner.b)
    return (two * u_iv) / (one + iv.sqrt(inner))


def B2_interval(X1, X2):
    """Enclosure of B_2 on the box X1 x X2 (Lemma 4.1).

    Only the upper endpoint is used by the algorithm; it equals the maximum of
    B_2(., v_b) over {u_a, u_b} and u = v_b (if v_b lies in X1), rounded up.
    """
    u_a, u_b = X1.a, X1.b
    v_a, v_b = X2.a, X2.b
    res_intervals = []

    if u_a <= v_b:
        min_u, max_u = u_a, min(u_b, v_b)
        min_v, max_v = max(u_a, v_a), v_b
        min_val = (mpmath.mpi(min_u, min_u) + mpmath.mpi(min_v, min_v)).a
        max_val = (mpmath.mpi(max_u, max_u) + mpmath.mpi(max_v, max_v)).b
        res_intervals.append(mpmath.mpi(min_val, max_val))

    if u_b > v_a:
        min_cands = []
        u_start_min = max(u_a, v_a)
        if u_start_min > v_a:
            min_cands.append(_eval_b2_point(u_start_min, v_a).a)
        elif u_start_min == v_a:
            min_cands.append((mpmath.mpi(u_start_min, u_start_min) + mpmath.mpi(v_a, v_a)).a)

        if u_b > v_a:
            min_cands.append(_eval_b2_point(u_b, v_a).a)

        nine_eighths = mpmath.mpi('9.0', '9.0') / mpmath.mpi('8.0', '8.0')
        crit_u_iv = nine_eighths * mpmath.mpi(v_a, v_a)

        if u_a <= crit_u_iv.b and crit_u_iv.a < u_b:
            twenty_seven_sixteenths = mpmath.mpi('27.0', '27.0') / mpmath.mpi('16.0', '16.0')
            slice_min = (twenty_seven_sixteenths * mpmath.mpi(v_a, v_a)).a
            min_cands.append(slice_min)

        max_cands = []
        u_start_max = max(u_a, v_b)
        if u_start_max > v_b:
            max_cands.append(_eval_b2_point(u_start_max, v_b).b)
        elif u_start_max == v_b:
            max_cands.append((mpmath.mpi(u_start_max, u_start_max) + mpmath.mpi(v_b, v_b)).b)

        if u_b > v_b:
            max_cands.append(_eval_b2_point(u_b, v_b).b)

        if min_cands and max_cands:
            res_intervals.append(mpmath.mpi(min(min_cands), max(max_cands)))

    if not res_intervals:
        return X1 + X2

    global_min = min(i.a for i in res_intervals)
    global_max = max(i.b for i in res_intervals)
    return mpmath.mpi(global_min, global_max)


def ub_dp(L, U, n):
    """Rigorous upper bound for sup of B_DP over the box prod [L_i, U_i]."""
    dp = [mpmath.mpi(0, 0)] * (n + 1)
    decisions = [0] * (n + 1)
    X = [mpmath.mpi(L[i], U[i]) for i in range(n)]
    dp[1] = dp[0] + X[0]
    decisions[1] = 1
    for i in range(2, n + 1):
        c1 = dp[i - 1] + X[i - 1]
        c2 = dp[i - 2] + B2_interval(X[i - 2], X[i - 1])
        if c2.b <= c1.b:
            dp[i] = c2
            decisions[i] = 2
        else:
            dp[i] = c1
            decisions[i] = 1

    path = []
    curr = n
    while curr > 0:
        d = decisions[curr]
        path.append(d)
        curr -= d
    path.reverse()
    return dp[n].b, path


def propagate_bounds(L, U, n):
    """Domain reduction (Section 4.2). Returns (None, None) if the box is empty."""
    L = list(L)
    U = list(U)
    changed = True
    iters = 0
    frac_3_4 = mpmath.mpi('3.0', '3.0') / mpmath.mpi('4.0', '4.0')
    frac_2_3 = mpmath.mpi('2.0', '2.0') / mpmath.mpi('3.0', '3.0')
    frac_1_2 = mpmath.mpi('1.0', '1.0') / mpmath.mpi('2.0', '2.0')
    one = mpmath.mpi('1.0', '1.0')
    tol = mpmath.mpf('1e-9')

    while changed and iters < 10:
        changed = False
        iters += 1
        for i in range(1, n):
            if U[i] > U[0]: U[i] = U[0]; changed = True
            if L[0] < L[i]: L[0] = L[i]; changed = True

        for i in range(n - 1):
            L_b = mpmath.mpf((mpmath.mpi(L[i], L[i]) * frac_3_4).a)
            if L[i + 1] < L_b: L[i + 1] = L_b; changed = True
            U_b = mpmath.mpf((mpmath.mpi(U[i + 1], U[i + 1]) / frac_3_4).b)
            if U[i] > U_b: U[i] = U_b; changed = True

        for i in range(n - 2):
            L_b = mpmath.mpf((mpmath.mpi(L[i], L[i]) * frac_2_3).a)
            if L[i + 2] < L_b: L[i + 2] = L_b; changed = True
            U_b = mpmath.mpf((mpmath.mpi(U[i + 2], U[i + 2]) / frac_2_3).b)
            if U[i] > U_b: U[i] = U_b; changed = True

        for i in range(n - 3):
            L_b = mpmath.mpf((mpmath.mpi(L[i], L[i]) * frac_1_2).a)
            if L[i + 3] < L_b: L[i + 3] = L_b; changed = True
            U_b = mpmath.mpf((mpmath.mpi(U[i + 3], U[i + 3]) / frac_1_2).b)
            if U[i] > U_b: U[i] = U_b; changed = True

        for i in range(n):
            prod_L_exc = one
            for j in range(n):
                if j != i: prod_L_exc *= mpmath.mpi(L[j], L[j])
            if prod_L_exc.a > 0:
                max_xi = mpmath.mpf((one / prod_L_exc).b)
                if max_xi < U[i]: U[i] = max_xi; changed = True

            prod_U_exc = one
            for j in range(n):
                if j != i: prod_U_exc *= mpmath.mpi(U[j], U[j])
            if prod_U_exc.a > 0:
                min_xi = mpmath.mpf((one / prod_U_exc).a)
                if min_xi > L[i]: L[i] = min_xi; changed = True

        for i in range(n):
            if L[i] > U[i] + tol: return None, None

    return L, U


def is_feasible(L, U):
    """Cheap infeasibility test: X_1 >= 1 and prod X_i = 1 must be possible."""
    if U[0] < mpmath.mpf('1.0'): return False
    prod_L = mpmath.mpi('1.0', '1.0')
    for x in L: prod_L *= mpmath.mpi(x, x)
    if prod_L.a > mpmath.mpf('1.0'): return False
    prod_U = mpmath.mpi('1.0', '1.0')
    for x in U: prod_U *= mpmath.mpi(x, x)
    if prod_U.b < mpmath.mpf('1.0'): return False
    return True


def root_boxes(n):
    """The 2^(n-1) root boxes of Algorithm 1, reduced, with their indices.

    Index k corresponds to the k-th element of itertools.product over
    X_2..X_n in {[0.05,1], [1,4.25]} (all roots are nonempty for 4 <= n <= 11).
    """
    L0, U0 = mpmath.mpf('1.0'), mpmath.mpf('4.25')
    intervals = [('0.05', '1.0'), ('1.0', '4.25')]
    boxes = []
    for p in itertools.product(intervals, repeat=n - 1):
        L_init = [L0] + [mpmath.mpf(i[0]) for i in p]
        U_init = [U0] + [mpmath.mpf(i[1]) for i in p]
        L_prop, U_prop = propagate_bounds(L_init, U_init, n)
        if L_prop is not None:
            boxes.append((len(boxes), L_prop, U_prop))
    return boxes


def run_subtree(n, L_init, U_init):
    """Algorithm 1 on one root box. Returns (certificate dict, survivors)."""
    queue = [(L_init, U_init, 0)]
    nodes_explored = 0
    max_depth = 0
    pruned_infeasible = 0
    pruned_bound = 0
    pruned_propagation = 0
    surviving = []

    threshold = mpmath.mpf(f'{n}.0')
    split_divisor = mpmath.mpf('2.0')
    cpu0 = time.process_time()

    while queue:
        L, U, depth = queue.pop()
        nodes_explored += 1
        if depth > max_depth:
            max_depth = depth

        if not is_feasible(L, U):
            pruned_infeasible += 1
            continue

        ub_mpf, path = ub_dp(L, U, n)
        if ub_mpf <= threshold:
            pruned_bound += 1
            continue

        if max(U) <= C_HI and min(L) >= C_LO:
            # mpf values at 53-bit precision convert to float exactly.
            surviving.append(([float(x) for x in L], [float(x) for x in U]))
            continue

        widths = [U[i] - L[i] for i in range(n)]
        max_w = max(widths)
        split_dim = widths.index(max_w)
        mid = L[split_dim] + max_w / split_divisor

        L1, U1 = list(L), list(U)
        U1[split_dim] = mid
        L1_prop, U1_prop = propagate_bounds(L1, U1, n)
        if L1_prop is not None:
            queue.append((L1_prop, U1_prop, depth + 1))
        else:
            pruned_propagation += 1

        L2, U2 = list(L), list(U)
        L2[split_dim] = mid
        L2_prop, U2_prop = propagate_bounds(L2, U2, n)
        if L2_prop is not None:
            queue.append((L2_prop, U2_prop, depth + 1))
        else:
            pruned_propagation += 1

    cert = {
        "nodes_explored": nodes_explored,
        "max_depth": max_depth,
        "pruned_infeasible": pruned_infeasible,
        "pruned_bound": pruned_bound,
        "pruned_propagation": pruned_propagation,
        "surviving_count": len(surviving),
        "worker_cpu_seconds": time.process_time() - cpu0,
    }
    return cert, surviving
