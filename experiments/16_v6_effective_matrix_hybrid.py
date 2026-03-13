#!/usr/bin/env python3
"""
Experiment 16: Hybrid attack on the effective matrix M_eff=(I+B)A
=================================================================

Goal:
  Test whether learning M_eff alone drastically reduces the key space.

Method:
  - Recover M_eff by interpolation (Experiment 14 primitive).
  - Enumerate all A in GL(2,p), derive B = M_eff A^{-1} - I.
  - Keep only pairs with B in GL(2,p).
  - Measure how many random local captures are needed to isolate the true pair.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.v6_audit_utils import (
    all_gl_matrices_d2,
    build_v6_instance,
    collect_local_captures,
    coeff_to_node_map,
    gl_size_d2,
    recover_linear_coefficient_matrix,
)


def candidate_pairs_from_Meff(proto, M_eff):
    Fp = proto.Fp
    p = Fp.p
    I = np.eye(2, dtype=int)
    pairs = []
    for A in all_gl_matrices_d2(Fp):
        A_inv = Fp.mat_inv(A)
        B = (M_eff @ A_inv - I) % p
        det_B = int((B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]) % p)
        if det_B == 0:
            continue
        pairs.append((A.copy(), B.copy()))
    return pairs


def first_failure_index(proto, A, B, captures):
    for idx, (s_v, c_v) in enumerate(captures, start=1):
        y = proto.Fp.mat_vec(A, s_v)
        z = proto.sigma(y)
        c_hat = (y + proto.Fp.mat_vec(B, z)) % proto.Fp.p
        if not np.array_equal(c_hat, c_v):
            return idx
    return len(captures) + 1


def run_config(p, seed):
    proto, key, _ = build_v6_instance(d=2, p=p, seed=seed)
    M_coeff, _ = recover_linear_coefficient_matrix(proto, key, node=0)
    P = coeff_to_node_map(proto.sheaf, 0)
    P_inv = proto.Fp.mat_inv(P)
    M_eff = (M_coeff @ P_inv) % p
    pairs = candidate_pairs_from_Meff(proto, M_eff)
    captures = collect_local_captures(proto, key, n_caps=12, node=0, seed=seed + 99)

    true_A = key.A[0] % p
    true_B = key.B[0] % p
    fail_rows = []
    for A, B in pairs:
        fail_rows.append((
            first_failure_index(proto, A, B, captures),
            np.array_equal(A, true_A) and np.array_equal(B, true_B),
        ))
    survivors_by_cap = {}
    min_cap_unique = None

    for k in range(1, len(captures) + 1):
        surv = 0
        true_in = False
        for fail_idx, is_true in fail_rows:
            if fail_idx > k:
                surv += 1
                if is_true:
                    true_in = True
        survivors_by_cap[k] = surv
        if min_cap_unique is None and surv == 1 and true_in:
            min_cap_unique = k

    total_key_pairs = gl_size_d2(p) ** 2
    return {
        "p": p,
        "P": P,
        "candidate_pairs": len(pairs),
        "total_key_pairs": total_key_pairs,
        "reduction_factor": total_key_pairs / max(len(pairs), 1),
        "min_cap_unique": min_cap_unique,
        "survivors_by_cap": survivors_by_cap,
    }


def print_result(r):
    print("\n" + "=" * 108)
    print(f"CONFIG: d=2 p={r['p']} sigma=id_spn")
    print("=" * 108)
    print(f"P_node =\n{r['P']}")
    print(f"Total key pairs |GL|^2: {r['total_key_pairs']}")
    print(f"Candidates after M_eff known: {r['candidate_pairs']}")
    print(f"Reduction factor: {r['reduction_factor']:.2f}x")
    print(f"min_cap_unique: {r['min_cap_unique']}")
    print("survivors_by_cap:")
    for k, v in r["survivors_by_cap"].items():
        print(f"  cap={k:>2} -> {v}")


def main():
    print("=" * 108)
    print("EXPERIMENT 16: V6 HYBRID ATTACK ON M_eff=(I+B)A")
    print("=" * 108)
    print("Measures whether interpolation-derived M_eff collapses the key space.")

    rows = [
        run_config(p=17, seed=16017),
        run_config(p=23, seed=16023),
    ]
    for r in rows:
        print_result(r)

    print("\nFindings:")
    print("- If M_eff leaves only a small number of candidate pairs, V6 may still leak")
    print("  enough structure for a practical hybrid CPA attack.")


if __name__ == "__main__":
    main()
