#!/usr/bin/env python3
"""
Experiment 12: Key-Recovery Diagnostics
=======================================
Implements four audit experiments:

1) Minimal captures needed for unique recovery (cap = 2..20).
2) Extension to d=3 (sampled A-search, since exhaustive GL(3,p) is huge).
3) Survival rate of A candidates as captures increase.
4) False positives: distinct (A,B) pairs matching all captures.

Method:
  - Node-wise known-plaintext setting.
  - For each A candidate, solve B linearly from the first fit_captures
    (fit_captures = d), then validate against capture prefixes.

Usage:
    python3 experiments/12_key_recovery_diagnostics.py
"""

import sys
import os
import random
import numpy as np
from itertools import product

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


def all_gl_matrices_d2(Fp):
    p = Fp.p
    for a, b, c, d in product(range(p), repeat=4):
        if (a * d - b * c) % p == 0:
            continue
        yield np.array([[a, b], [c, d]], dtype=int)


def sampled_gl_matrices(Fp, d, budget, rng):
    for _ in range(budget):
        yield Fp.random_gl(d, rng=rng)


def collect_captures(proto, key, n_caps, seed):
    rng = random.Random(seed)
    sheaf = proto.sheaf
    node = 0
    caps = []
    for _ in range(n_caps):
        s = sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)
        s_v = sheaf.get_node_value(s, node)
        c_v = sheaf.get_node_value(c, node)
        caps.append((s_v, c_v))
    return caps


def solve_B_for_A(proto, A, captures_fit):
    Fp = proto.Fp
    p = Fp.p
    d = A.shape[0]
    equations = []
    for s_v, c_v in captures_fit:
        y = Fp.mat_vec(A, s_v)
        z = proto.sigma(y)
        rhs = (c_v - y) % p
        for i in range(d):
            coeffs = [0] * (d * d)
            for k in range(d):
                coeffs[i * d + k] = int(z[k]) % p
            equations.append((coeffs, int(rhs[i])))
    sol = solve_linear_Fp(equations, d * d, p)
    if sol is None:
        return None
    return sol.reshape((d, d)) % p


def first_failure_index(proto, A, B, captures):
    """Return 1-based index of first failing capture, or len(captures)+1 if all pass."""
    Fp = proto.Fp
    p = Fp.p
    for idx, (s_v, c_v) in enumerate(captures, start=1):
        y = Fp.mat_vec(A, s_v)
        z = proto.sigma(y)
        c_hat = (y + Fp.mat_vec(B, z)) % p
        if not np.array_equal(c_hat, c_v % p):
            return idx
    return len(captures) + 1


def analyze_config(d, p, sigma_name, max_caps, mode, sample_budget, seed):
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, sigma_name))
    key = proto.keygen(rng=rng)

    caps = collect_captures(proto, key, max_caps, seed + 999)
    fit_captures = d
    fit_set = caps[:fit_captures]

    if mode == "exhaustive" and d == 2:
        A_stream = all_gl_matrices_d2(Fp)
    elif mode == "sampled":
        A_stream = sampled_gl_matrices(Fp, d, sample_budget, rng)
    else:
        raise ValueError("Unsupported mode/d combination")

    A_true = key.A[0] % p
    B_true = key.B[0] % p

    tested_A = 0
    survivors_by_cap = {k: 0 for k in range(fit_captures, max_caps + 1)}
    true_survives_by_cap = {k: False for k in range(fit_captures, max_caps + 1)}
    false_pos_by_cap = {k: 0 for k in range(fit_captures, max_caps + 1)}

    for A in A_stream:
        tested_A += 1
        B = solve_B_for_A(proto, A, fit_set)
        if B is None:
            continue

        fail_idx = first_failure_index(proto, A, B, caps)
        is_true = np.array_equal(A % p, A_true) and np.array_equal(B % p, B_true)

        for k in range(fit_captures, max_caps + 1):
            if fail_idx > k:
                survivors_by_cap[k] += 1
                if is_true:
                    true_survives_by_cap[k] = True
                else:
                    false_pos_by_cap[k] += 1

    min_cap_unique = None
    for k in range(fit_captures, max_caps + 1):
        if survivors_by_cap[k] == 1 and true_survives_by_cap[k]:
            min_cap_unique = k
            break

    return {
        "d": d,
        "p": p,
        "sigma": sigma_name,
        "mode": mode,
        "tested_A": tested_A,
        "fit_captures": fit_captures,
        "max_caps": max_caps,
        "min_cap_unique": min_cap_unique,
        "survivors_by_cap": survivors_by_cap,
        "true_survives_by_cap": true_survives_by_cap,
        "false_pos_by_cap": false_pos_by_cap,
    }


def print_config_tables(res):
    print("\n" + "=" * 108)
    print(f"CONFIG: d={res['d']} p={res['p']} sigma={res['sigma']} mode={res['mode']} tested_A={res['tested_A']}")
    print(f"fit_captures={res['fit_captures']} min_cap_unique={res['min_cap_unique']}")
    print("=" * 108)
    print(f"{'cap':>4} {'survivors':>12} {'survival_rate':>14} {'true_in':>8} {'false_pos':>10}")
    print("-" * 108)
    tested = max(res["tested_A"], 1)
    for k in range(res["fit_captures"], res["max_caps"] + 1):
        surv = res["survivors_by_cap"][k]
        rate = surv / tested
        true_in = res["true_survives_by_cap"][k]
        fp = res["false_pos_by_cap"][k]
        print(f"{k:>4} {surv:>12} {100*rate:>13.4f}% {str(true_in):>8} {fp:>10}")


def print_summary(rows):
    print("\n" + "=" * 108)
    print("SUMMARY")
    print("=" * 108)
    print(f"{'d':>3} {'p':>4} {'sigma':>8} {'mode':>10} {'A_tested':>10} {'min_unique_cap':>14} "
          f"{'surv@20':>10} {'fp@20':>10}")
    print("-" * 108)
    for r in rows:
        surv20 = r["survivors_by_cap"][r["max_caps"]]
        fp20 = r["false_pos_by_cap"][r["max_caps"]]
        print(f"{r['d']:>3} {r['p']:>4} {r['sigma']:>8} {r['mode']:>10} {r['tested_A']:>10} "
              f"{str(r['min_cap_unique']):>14} {surv20:>10} {fp20:>10}")


def main():
    print("=" * 108)
    print("EXPERIMENT 12: KEY-RECOVERY DIAGNOSTICS")
    print("=" * 108)

    all_results = []

    # A/B/C/D metrics for d=2 with exhaustive scan (p=17).
    for sigma in ["inverse", "cube"]:
        r = analyze_config(
            d=2, p=17, sigma_name=sigma, max_caps=20,
            mode="exhaustive", sample_budget=0, seed=3300 + (0 if sigma == "inverse" else 100)
        )
        all_results.append(r)
        print_config_tables(r)

    # d=3 sampled diagnostics for two primes.
    for sigma in ["inverse", "cube"]:
        for p in [17, 23]:
            r = analyze_config(
                d=3, p=p, sigma_name=sigma, max_caps=20,
                mode="sampled", sample_budget=60000, seed=4400 + p + (0 if sigma == "inverse" else 200)
            )
            all_results.append(r)
            print_config_tables(r)

    print_summary(all_results)


if __name__ == "__main__":
    main()
