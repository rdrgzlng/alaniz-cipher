#!/usr/bin/env python3
"""
Experiment 14: V6 interpolation audit
====================================

Goal:
  Check empirically what polynomial interpolation recovers on the V6 sigma
  `id_spn`. The design claim is that line interpolation should expose the
  effective matrix

      M_eff = (I + B) A

  rather than A directly.

Method:
  - Query ciphertexts on axis lines t * e_j for all t in F_p.
  - Interpolate each ciphertext coordinate as a univariate polynomial in t.
  - Extract the t^1 coefficient for each axis j and assemble a matrix.
  - Compare it against A and against (I + B)A.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.v6_audit_utils import (
    build_v6_instance,
    coeff_to_node_map,
    count_effective_matrix_decompositions,
    recover_linear_coefficient_matrix,
)


def run_config(d, p, seed):
    proto, key, _ = build_v6_instance(d=d, p=p, seed=seed)
    Fp = proto.Fp
    node = 0

    M_rec, per_axis = recover_linear_coefficient_matrix(proto, key, node=node)
    P = coeff_to_node_map(proto.sheaf, node)
    A_true = key.A[node] % p
    B_true = key.B[node] % p
    M_coeff_true = ((np.eye(d, dtype=int) + B_true) @ A_true @ P) % p

    res = {
        "d": d,
        "p": p,
        "A_match": np.array_equal(M_rec, A_true),
        "M_eff_coeff_match": np.array_equal(M_rec, M_coeff_true),
        "M_rec": M_rec,
        "A_true": A_true,
        "B_true": B_true,
        "P": P,
        "M_coeff_true": M_coeff_true,
        "per_axis": per_axis,
    }

    if d == 2:
        P_inv = Fp.mat_inv(P)
        M_stalk = (M_rec @ P_inv) % p
        n_dec, examples = count_effective_matrix_decompositions(Fp, M_stalk)
        res["M_stalk"] = M_stalk
        res["decomposition_count"] = n_dec
        res["decomposition_examples"] = examples

    return res


def print_result(r):
    print("\n" + "=" * 108)
    print(f"CONFIG: d={r['d']} p={r['p']} sigma=id_spn")
    print("=" * 108)
    print(f"A_match={r['A_match']}  M_eff_coeff_match={r['M_eff_coeff_match']}")
    print(f"A_true =\n{r['A_true']}")
    print(f"B_true =\n{r['B_true']}")
    print(f"P_node =\n{r['P']}")
    print(f"(I+B)A P =\n{r['M_coeff_true']}")
    print(f"Recovered coefficient-basis linear matrix =\n{r['M_rec']}")
    for axis_info in r["per_axis"]:
        print(
            f"  axis={axis_info['axis']} "
            f"const={axis_info['constant']} "
            f"lin={axis_info['linear']}"
        )
    if "decomposition_count" in r:
        print(f"Recovered stalk-space M_eff =\n{r['M_stalk']}")
        print(f"Decompositions of recovered M_eff into (A,B): {r['decomposition_count']}")
        for idx, (A, B) in enumerate(r["decomposition_examples"], start=1):
            print(f"  example#{idx}: A={A.tolist()} B={B.tolist()}")


def main():
    print("=" * 108)
    print("EXPERIMENT 14: V6 INTERPOLATION AUDIT")
    print("=" * 108)
    print("Checks whether line interpolation recovers A or only M_eff=(I+B)A.")

    rows = [
        run_config(d=2, p=17, seed=14017),
        run_config(d=2, p=23, seed=14023),
        run_config(d=3, p=17, seed=14117),
    ]
    for r in rows:
        print_result(r)

    print("\nFindings:")
    print("- On the current V6 implementation, the t^1 coefficient matches")
    print("  the coefficient-basis matrix (I+B)A P_v.")
    print("- It does not match A directly in the tested instances.")
    print("- For d=2, the recovered M_eff still admits many (A,B) decompositions,")
    print("  so interpolation alone does not identify the true key pair.")


if __name__ == "__main__":
    main()
