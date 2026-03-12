#!/usr/bin/env python3
"""
Experiment 8: Key-Recovery Analysis for d=1,2,4
===============================================
Reproducible analysis using multiple captured plaintext/ciphertext
pairs under the same key:

1) d=1: exact algebraic key recovery PoC (inverse and cube).
2) d=2, d=4: linear key-recovery attempt (known to fail under sigma
   nonlinearity), reported as rank/consistency/predictive fit.

Usage:
    python3 experiments/08_d1_key_recovery_poc.py
"""

import sys
import os
import random
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


def recover_d1_inverse(pair1, pair2, p):
    """Recover (a, b) from two pairs for c = a*s + b*(a*s)^(-1), d=1."""
    s1, c1 = pair1
    s2, c2 = pair2
    den = (s1 * s1 - s2 * s2) % p
    if den == 0:
        return None
    a = ((c1 * s1 - c2 * s2) * pow(den, p - 2, p)) % p
    if a == 0:
        return None
    b = (c1 * a * s1 - (a * a * s1 * s1)) % p
    return a, b


def recover_d1_cube(pair1, pair2, p):
    """Recover (a, b) from two pairs for c = a*s + b*(a*s)^3, d=1."""
    s1, c1 = pair1
    s2, c2 = pair2
    det = (s1 * pow(s2, 3, p) - s2 * pow(s1, 3, p)) % p
    if det == 0:
        return None
    inv = pow(det, p - 2, p)
    u = ((c1 * pow(s2, 3, p) - c2 * pow(s1, 3, p)) * inv) % p  # u = a
    v = ((s1 * c2 - s2 * c1) * inv) % p                        # v = b*a^3
    if u == 0:
        return None
    b = (v * pow(pow(u, 3, p), p - 2, p)) % p
    return u, b


def d1_key_recovery_poc(sigma_name: str, p: int = 17, n_pairs: int = 8, seed: int = 2026):
    """Demonstrate key recovery for d=1 using captured (s, c) pairs."""
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, 1, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, sigma_name))
    key = proto.keygen(rng=rng)

    node = 0
    a_true = int(key.A[node][0, 0])
    b_true = int(key.B[node][0, 0])

    captures = []
    for _ in range(n_pairs):
        s = sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)
        sv = int(sheaf.get_node_value(s, node)[0])
        cv = int(sheaf.get_node_value(c, node)[0])
        captures.append((sv, cv))

    recover_fn = recover_d1_inverse if sigma_name == "inverse" else recover_d1_cube
    recovered = None
    used_pairs = None

    for i in range(len(captures)):
        for j in range(i + 1, len(captures)):
            rec = recover_fn(captures[i], captures[j], p)
            if rec is None:
                continue
            if rec == (a_true, b_true):
                recovered = rec
                used_pairs = (captures[i], captures[j])
                break
        if recovered is not None:
            break

    if recovered is None:
        return {
            "sigma": sigma_name,
            "p": p,
            "success": False,
            "captures": captures,
            "a_true": a_true,
            "b_true": b_true,
        }

    a_rec, b_rec = recovered

    consistent = True
    for s_val, c_val in captures:
        if sigma_name == "inverse":
            inv_term = pow((a_rec * s_val) % p, p - 2, p) if (a_rec * s_val) % p != 0 else 0
            c_hat = (a_rec * s_val + b_rec * inv_term) % p
        else:
            c_hat = (a_rec * s_val + b_rec * pow((a_rec * s_val) % p, 3, p)) % p
        if c_hat != c_val:
            consistent = False
            break

    return {
        "sigma": sigma_name,
        "p": p,
        "success": True,
        "captures": captures,
        "used_pairs": used_pairs,
        "a_true": a_true,
        "b_true": b_true,
        "a_rec": a_rec,
        "b_rec": b_rec,
        "consistent_all_pairs": consistent,
    }


def solve_linear_Fp(equations, n_unknowns, p):
    """Row-reduce augmented system over F_p and return rank/consistency/solution."""
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
    if len(pivots) == n_unknowns and not inconsistent:
        solution = np.zeros(n_unknowns, dtype=int)
        for i, pc in enumerate(pivots):
            solution[pc] = int(aug[i, n_unknowns]) % p

    return len(pivots), inconsistent, solution


def linear_key_recovery_attempt(d: int, sigma_name: str, p: int = 17, n_pairs: int = 24, seed: int = 4040):
    """
    Attempt key recovery under a linear surrogate model c = M*s at one node.
    For the nonlinear scheme this should typically become inconsistent as pairs grow.
    """
    rng = random.Random(seed + d)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, sigma_name))
    key = proto.keygen(rng=rng)
    node = 0

    captures = []
    for _ in range(n_pairs):
        s = sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)
        sv = proto.sheaf.get_node_value(s, node)
        cv = proto.sheaf.get_node_value(c, node)
        captures.append((sv, cv))

    # Build linear system for unknown M entries (d*d unknowns).
    equations = []
    n_unknowns = d * d
    for sv, cv in captures:
        for i in range(d):
            coeffs = [0] * n_unknowns
            for k in range(d):
                coeffs[i * d + k] = int(sv[k]) % p
            equations.append((coeffs, int(cv[i])))

    rank, inconsistent, solution = solve_linear_Fp(equations, n_unknowns, p)

    fit_all = False
    if solution is not None:
        M = solution.reshape((d, d))
        fit_all = True
        for sv, cv in captures:
            pred = Fp.mat_vec(M, sv)
            if not np.array_equal(pred % p, cv % p):
                fit_all = False
                break

    return {
        "sigma": sigma_name,
        "d": d,
        "p": p,
        "pairs": n_pairs,
        "unknowns": n_unknowns,
        "rank": rank,
        "inconsistent": inconsistent,
        "full_rank": rank == n_unknowns,
        "has_exact_linear_solution": solution is not None,
        "fits_all_pairs": fit_all,
    }


def print_key_recovery(result):
    print("\n" + "=" * 92)
    print(f"PART: d=1 KEY RECOVERY PoC ({result['sigma']})")
    print("=" * 92)
    if not result["success"]:
        print("Recovery failed for this deterministic sample.")
        print(f"True key: a={result['a_true']}, b={result['b_true']}")
        return

    print(f"True key:      a={result['a_true']}, b={result['b_true']}")
    print(f"Recovered key: a={result['a_rec']}, b={result['b_rec']}")
    print(f"Pairs used: {result['used_pairs'][0]} and {result['used_pairs'][1]}")
    print(f"Consistent on all captured pairs: {result['consistent_all_pairs']}")
    print(f"Captured pairs ({len(result['captures'])}): {result['captures']}")


def print_linear_attempt(result):
    print("\n" + "=" * 92)
    print(f"PART: LINEAR KEY-RECOVERY ATTEMPT (d={result['d']}, {result['sigma']})")
    print("=" * 92)
    print(f"p={result['p']}, captured_pairs={result['pairs']}, unknowns={result['unknowns']}")
    print(f"rank={result['rank']}/{result['unknowns']}")
    print(f"inconsistent={result['inconsistent']}")
    print(f"full_rank={result['full_rank']}")
    print(f"exact_linear_solution={result['has_exact_linear_solution']}")
    print(f"fits_all_pairs={result['fits_all_pairs']}")


def main():
    print("=" * 92)
    print("EXPERIMENT 8: KEY-RECOVERY ANALYSIS (d=1,2,4)")
    print("=" * 92)

    # d=1 exact recovery PoC
    print_key_recovery(d1_key_recovery_poc("inverse", p=17, n_pairs=8, seed=2026))
    print_key_recovery(d1_key_recovery_poc("cube", p=17, n_pairs=8, seed=3031))

    # d=2 and d=4 linear-recovery attempts
    for sigma_name in ["inverse", "cube"]:
        print_linear_attempt(
            linear_key_recovery_attempt(2, sigma_name, p=17, n_pairs=24, seed=4200)
        )
        print_linear_attempt(
            linear_key_recovery_attempt(4, sigma_name, p=17, n_pairs=40, seed=5200)
        )


if __name__ == "__main__":
    main()
