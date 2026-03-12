#!/usr/bin/env python3
"""
Experiment 9: Known-Plaintext Key Recovery via A Search
=======================================================
Attack vector:
  Given captured pairs (s, c) at one node with fixed key, recover (A, B)
  by:
    1) Search over A in GL(d, F_p) (exhaustive or sampled)
    2) For each A, solve B linearly from:
         c - A*s = B*sigma(A*s)
    3) Validate candidate (A, B) against all captures

This is practical for small toy parameters (e.g., d=2, p up to ~29).

Usage:
    python3 experiments/09_exhaustive_A_known_plaintext_attack.py
"""

import sys
import os
import random
import time
from itertools import product
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


def solve_linear_Fp(equations, n_unknowns, p):
    """Solve/rank linear system over F_p."""
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

    inconsistent = False
    for r in range(len(pivots), n_eq):
        if aug[r, n_unknowns] % p != 0:
            inconsistent = True
            break

    solution = None
    if not inconsistent:
        # One solution with free vars = 0 if underdetermined.
        solution = np.zeros(n_unknowns, dtype=int)
        for i, pc in enumerate(pivots):
            solution[pc] = int(aug[i, n_unknowns]) % p

    return len(pivots), inconsistent, solution


def gl_size(d: int, p: int) -> int:
    """|GL(d, p)| = prod_{i=0}^{d-1} (p^d - p^i)."""
    total = 1
    pd = p ** d
    for i in range(d):
        total *= (pd - p ** i)
    return total


def all_gl_matrices(Fp: FiniteField, d: int):
    """Enumerate GL(d, F_p) by brute force (feasible for small d)."""
    p = Fp.p
    if d == 2:
        for a, b, c, d_ in product(range(p), repeat=4):
            det = (a * d_ - b * c) % p
            if det == 0:
                continue
            yield np.array([[a, b], [c, d_]], dtype=int)
        return

    for vals in product(range(p), repeat=d * d):
        M = np.array(vals, dtype=int).reshape((d, d))
        try:
            Fp.mat_inv(M)
            yield M
        except ValueError:
            continue


def sampled_gl_matrices(Fp: FiniteField, d: int, budget: int, rng):
    """Sample random invertible matrices from GL(d, F_p)."""
    for _ in range(budget):
        yield Fp.random_gl(d, rng=rng)


def recover_key_from_captures(proto: Protocol, captures, d: int, A_stream,
                              stop_when_found=None):
    """
    Recover (A, B) from node captures:
      captures = [(s_v, c_v), ...]
    """
    Fp = proto.Fp
    p = Fp.p
    tested_A = 0
    candidates = []
    fit_captures = captures[: max(2, d)]  # enough to solve B quickly in many cases

    for A in A_stream:
        tested_A += 1

        equations = []
        for s_v, c_v in fit_captures:
            y = Fp.mat_vec(A, s_v)
            z = proto.sigma(y)
            rhs = (c_v - y) % p
            # rhs_i = sum_k b_{i,k} z_k
            for i in range(d):
                coeffs = [0] * (d * d)
                for k in range(d):
                    coeffs[i * d + k] = int(z[k]) % p
                equations.append((coeffs, int(rhs[i]) % p))

        rank, inconsistent, sol = solve_linear_Fp(equations, d * d, p)
        if inconsistent or sol is None:
            continue

        B = sol.reshape((d, d)) % p

        # Validate exact fit across all captures.
        ok = True
        for s_v, c_v in captures:
            y = Fp.mat_vec(A, s_v)
            z = proto.sigma(y)
            c_hat = (y + Fp.mat_vec(B, z)) % p
            if not np.array_equal(c_hat, c_v):
                ok = False
                break
        if ok:
            candidates.append((A.copy(), B.copy(), rank))
            if stop_when_found is not None:
                A_target, B_target = stop_when_found
                if np.array_equal(A % p, A_target % p) and np.array_equal(B % p, B_target % p):
                    break

    return tested_A, candidates


def run_attack(d: int, p: int, sigma_name: str, n_captures: int, seed: int,
               mode: str = "exhaustive", sampled_budget: int = 0):
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, sigma_name))
    key = proto.keygen(rng=rng)

    node = 0
    captures = []
    for _ in range(n_captures):
        s = sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)
        captures.append(
            (sheaf.get_node_value(s, node), sheaf.get_node_value(c, node))
        )

    A_true = key.A[node] % p
    B_true = key.B[node] % p

    if mode == "exhaustive":
        A_stream = all_gl_matrices(Fp, d)
    elif mode == "sampled":
        A_stream = sampled_gl_matrices(Fp, d, sampled_budget, rng)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    t0 = time.time()
    tested_A, candidates = recover_key_from_captures(
        proto, captures, d, A_stream,
        stop_when_found=(A_true, B_true) if mode == "exhaustive" else None
    )
    elapsed = time.time() - t0

    true_recovered = False
    for A, B, _ in candidates:
        if np.array_equal(A % p, A_true) and np.array_equal(B % p, B_true):
            true_recovered = True
            break

    return {
        "d": d,
        "p": p,
        "sigma": sigma_name,
        "captures": n_captures,
        "mode": mode,
        "tested_A": tested_A,
        "n_candidates": len(candidates),
        "true_recovered": true_recovered,
        "elapsed_s": elapsed,
    }


def main():
    print("=" * 92)
    print("EXPERIMENT 9: KNOWN-PLAINTEXT KEY RECOVERY (EXHAUSTIVE + SAMPLED)")
    print("=" * 92)

    # Multi-prime campaign. Exhaustive mode is bounded by GL size.
    primes = [17, 23, 29, 47, 59, 89, 131, 251]
    sigmas = ["inverse", "cube"]
    d = 2
    capture_grid = [4, 8, 12]
    max_gl_to_run = 100_000  # exhaustive threshold (keeps runtime practical)
    sampled_budget = 5000    # random A trials for large-prime rows

    for sigma_name in sigmas:
        print("\n" + "=" * 92)
        print(f"SIGMA={sigma_name}")
        print("=" * 92)
        print(f"{'p':>5} {'cap':>5} {'|GL|':>12} {'mode':>12} {'A_tested':>12} {'cands':>8} "
              f"{'true_key':>10} {'time(s)':>10}")
        print("-" * 92)

        for p in primes:
            Fp = FiniteField(p)
            if sigma_name == "cube" and not Fp.cube_invertible():
                for cap in capture_grid:
                    print(f"{p:>5} {cap:>5} {'-':>12} {'skip-cube':>12} {'-':>12} {'-':>8} {'-':>10} {'-':>10}")
                continue

            gl = gl_size(d, p)
            mode = "exhaustive" if gl <= max_gl_to_run else "sampled"
            for cap in capture_grid:
                r = run_attack(
                    d=d,
                    p=p,
                    sigma_name=sigma_name,
                    n_captures=cap,
                    seed=7000 + p + 100 * cap + (0 if sigma_name == "inverse" else 10000),
                    mode=mode,
                    sampled_budget=sampled_budget,
                )
                print(f"{p:>5} {cap:>5} {gl:>12} {mode:>12} {r['tested_A']:>12} "
                      f"{r['n_candidates']:>8} {str(r['true_recovered']):>10} {r['elapsed_s']:>10.2f}")

    return


if __name__ == "__main__":
    main()
