#!/usr/bin/env python3
"""
Experiment 13: Differential Column Attack on A
==============================================

Goal:
  Test whether chosen-plaintext differential queries along basis directions
  recover columns of A much more cheaply than exhaustive search over GL(d, F_p).

Idea:
  For plaintexts on a single axis s = t * e_j, we have

      c(t e_j) = t a_j + B sigma(t a_j)

  where a_j is column j of A.

  Taking consecutive differences removes the need to reason about the whole A:

      Δc_t = c((t+1)e_j) - c(t e_j)
           = a_j + B ( sigma((t+1)a_j) - sigma(t a_j) )

  For a guessed column a_j, B appears linearly. This means:
    1) enumerate candidate columns a_j != 0
    2) solve B linearly from the first fit_diffs (= d) differentials
    3) validate against longer differential prefixes

Why this matters:
  If columns can be recovered independently, the attack surface on A shrinks from
      search over GL(d, F_p)
  to roughly
      d * (p^d - 1)
  candidate columns, plus a consistency check on B.

Usage:
    python3 experiments/13_differential_column_attack.py
"""

import sys
import os
import random
from itertools import product
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


def solve_linear_Fp(equations, n_unknowns, p):
    n_eq = len(equations)
    aug = np.zeros((n_eq, n_unknowns + 1), dtype=object)
    for i, (coeffs, rhs) in enumerate(equations):
        aug[i, :n_unknowns] = np.array(coeffs, dtype=object) % p
        aug[i, n_unknowns] = int(rhs) % p

    pivots, row = [], 0
    for col in range(n_unknowns):
        found = False
        for r in range(row, n_eq):
            if aug[r][col] % p != 0:
                aug[[row, r]] = aug[[r, row]]
                found = True
                break
        if not found:
            continue
        pivots.append(col)
        iv = pow(int(aug[row][col]), p - 2, p)
        aug[row] = (aug[row] * iv) % p
        for r in range(n_eq):
            if r != row and aug[r][col] % p != 0:
                aug[r] = (aug[r] - aug[r][col] * aug[row]) % p
        row += 1

    for r in range(len(pivots), n_eq):
        if aug[r, n_unknowns] % p != 0:
            return None

    sol = np.zeros(n_unknowns, dtype=int)
    for i, pc in enumerate(pivots):
        sol[pc] = int(aug[i, n_unknowns]) % p
    return sol


def iter_nonzero_columns(Fp, d):
    for vals in product(range(Fp.p), repeat=d):
        v = np.array(vals, dtype=int)
        if np.any(v % Fp.p != 0):
            yield v % Fp.p


def make_axis_plaintext(d, p, axis, scalar):
    s = np.zeros(d, dtype=int)
    s[axis] = scalar % p
    return s


def collect_axis_differentials(proto, key, axis, n_scalars):
    """
    Collect consecutive line differentials on axis j:
      Δc_t = c((t+1)e_j) - c(t e_j),  t = 0..n_scalars-2
    Returns list of tuples (t, delta_c).
    """
    p = proto.Fp.p
    d = proto.dv
    node = 0
    line_ct = []
    for t in range(n_scalars):
        sv = make_axis_plaintext(d, p, axis, t)
        cv = proto.encrypt_node(sv, key.A[node], key.B[node]) % p
        line_ct.append(cv)

    diffs = []
    for t in range(n_scalars - 1):
        dc = (line_ct[t + 1] - line_ct[t]) % p
        diffs.append((t, dc))
    return diffs


def sigma_diff(proto, col, t):
    p = proto.Fp.p
    y0 = (t * col) % p
    y1 = ((t + 1) * col) % p
    return (proto.sigma(y1) - proto.sigma(y0)) % p


def solve_B_for_column_from_diffs(proto, col, diffs_fit):
    p = proto.Fp.p
    d = col.shape[0]
    equations = []
    for t, delta_c in diffs_fit:
        dz = sigma_diff(proto, col, t)
        rhs = (delta_c - col) % p
        for i in range(d):
            coeffs = [0] * (d * d)
            for k in range(d):
                coeffs[i * d + k] = int(dz[k]) % p
            equations.append((coeffs, int(rhs[i]) % p))
    sol = solve_linear_Fp(equations, d * d, p)
    if sol is None:
        return None
    return sol.reshape((d, d)) % p


def first_failure_index_diff(proto, col, B, diffs):
    p = proto.Fp.p
    for idx, (t, delta_c) in enumerate(diffs, start=1):
        dz = sigma_diff(proto, col, t)
        pred = (col + proto.Fp.mat_vec(B, dz)) % p
        if not np.array_equal(pred, delta_c % p):
            return idx
    return len(diffs) + 1


def analyze_axis(proto, key, axis, n_scalars):
    p = proto.Fp.p
    d = proto.dv
    diffs = collect_axis_differentials(proto, key, axis, n_scalars)
    fit_diffs = d
    fit_set = diffs[:fit_diffs]

    true_col = key.A[0][:, axis] % p
    true_B = key.B[0] % p

    tested_cols = 0
    survivors_by_cap = {k: 0 for k in range(fit_diffs, len(diffs) + 1)}
    true_survives_by_cap = {k: False for k in range(fit_diffs, len(diffs) + 1)}
    false_pos_by_cap = {k: 0 for k in range(fit_diffs, len(diffs) + 1)}
    true_B_recovered_by_cap = {k: False for k in range(fit_diffs, len(diffs) + 1)}

    for col in iter_nonzero_columns(proto.Fp, d):
        tested_cols += 1
        B = solve_B_for_column_from_diffs(proto, col, fit_set)
        if B is None:
            continue

        fail_idx = first_failure_index_diff(proto, col, B, diffs)
        is_true_col = np.array_equal(col % p, true_col)
        is_true_B = np.array_equal(B % p, true_B)

        for k in range(fit_diffs, len(diffs) + 1):
            if fail_idx > k:
                survivors_by_cap[k] += 1
                if is_true_col:
                    true_survives_by_cap[k] = True
                else:
                    false_pos_by_cap[k] += 1
                if is_true_col and is_true_B:
                    true_B_recovered_by_cap[k] = True

    min_cap_unique = None
    for k in range(fit_diffs, len(diffs) + 1):
        if survivors_by_cap[k] == 1 and true_survives_by_cap[k]:
            min_cap_unique = k
            break

    return {
        "axis": axis,
        "d": d,
        "p": p,
        "sigma": proto.params.sigma.name,
        "tested_cols": tested_cols,
        "fit_diffs": fit_diffs,
        "n_diffs": len(diffs),
        "min_cap_unique": min_cap_unique,
        "survivors_by_cap": survivors_by_cap,
        "true_survives_by_cap": true_survives_by_cap,
        "false_pos_by_cap": false_pos_by_cap,
        "true_B_recovered_by_cap": true_B_recovered_by_cap,
    }


def print_axis_table(res):
    print("\n" + "=" * 108)
    print(f"CONFIG: d={res['d']} p={res['p']} sigma={res['sigma']} axis={res['axis']} tested_cols={res['tested_cols']}")
    print(f"fit_diffs={res['fit_diffs']} min_cap_unique={res['min_cap_unique']}")
    print("=" * 108)
    print(f"{'diffs':>5} {'survivors':>12} {'survival_rate':>14} {'true_col':>10} {'false_pos':>10} {'true_B':>10}")
    print("-" * 108)
    tested = max(res['tested_cols'], 1)
    for k in range(res['fit_diffs'], res['n_diffs'] + 1):
        surv = res['survivors_by_cap'][k]
        rate = surv / tested
        tc = res['true_survives_by_cap'][k]
        fp = res['false_pos_by_cap'][k]
        tb = res['true_B_recovered_by_cap'][k]
        print(f"{k:>5} {surv:>12} {100*rate:>13.4f}% {str(tc):>10} {fp:>10} {str(tb):>10}")


def run_config(d, p, sigma_name, n_scalars, seed):
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, sigma_name))
    key = proto.keygen(rng=rng)

    results = []
    for axis in range(d):
        results.append(analyze_axis(proto, key, axis, n_scalars))
    return results


def print_summary(rows):
    print("\n" + "=" * 120)
    print("SUMMARY")
    print("=" * 120)
    print(f"{'d':>3} {'p':>4} {'sigma':>8} {'axis':>4} {'cols_tested':>12} {'min_unique_diff':>16} {'surv@final':>12} {'fp@final':>10}")
    print("-" * 120)
    for r in rows:
        final_surv = r['survivors_by_cap'][r['n_diffs']]
        final_fp = r['false_pos_by_cap'][r['n_diffs']]
        print(f"{r['d']:>3} {r['p']:>4} {r['sigma']:>8} {r['axis']:>4} {r['tested_cols']:>12} {str(r['min_cap_unique']):>16} {final_surv:>12} {final_fp:>10}")


def main():
    print("=" * 108)
    print("EXPERIMENT 13: DIFFERENTIAL COLUMN ATTACK ON A")
    print("=" * 108)

    all_rows = []

    campaigns = [
        (2, 17, 'inverse', 8, 5100),
        (2, 17, 'cube',    8, 5200),
        (3, 17, 'inverse', 8, 5300),
        (3, 17, 'cube',    8, 5400),
    ]

    for d, p, sigma_name, n_scalars, seed in campaigns:
        rows = run_config(d=d, p=p, sigma_name=sigma_name, n_scalars=n_scalars, seed=seed)
        for r in rows:
            all_rows.append(r)
            print_axis_table(r)

    print_summary(all_rows)


if __name__ == '__main__':
    main()
