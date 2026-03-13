#!/usr/bin/env python3
"""
Experiment 18: V6 node decoupling audit
======================================

Goal:
  Test whether CPA access still lets the attacker recover node-local effective
  matrices independently, despite the sheaf coupling.

Method:
  - Use the same global coefficient queries t * e_j for all nodes.
  - For each node v, project the ciphertext to c_v and interpolate the t^1 term.
  - Recover M_eff(v) = (I + B_v) A_v independently for every node.
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


def run_config(d, p, n_nodes, seed):
    proto, key, _ = build_v6_instance(d=d, p=p, n_nodes=n_nodes, seed=seed)
    rows = []
    for node in proto.sheaf.graph.nodes:
        M_rec, _ = recover_linear_coefficient_matrix(proto, key, node=node)
        P = coeff_to_node_map(proto.sheaf, node)
        A = key.A[node] % p
        B = key.B[node] % p
        M_true = ((np.eye(d, dtype=int) + B) @ A @ P) % p
        row = {
            "node": node,
            "match": np.array_equal(M_rec, M_true),
            "M_rec": M_rec,
            "M_true": M_true,
        }
        if d == 2:
            P_inv = proto.Fp.mat_inv(P)
            count, _ = count_effective_matrix_decompositions(proto.Fp, (M_rec @ P_inv) % p)
            row["decompositions"] = count
        rows.append(row)
    return {
        "d": d,
        "p": p,
        "n_nodes": n_nodes,
        "global_queries": d * p,
        "rows": rows,
    }


def print_result(r):
    print("\n" + "=" * 108)
    print(f"CONFIG: d={r['d']} p={r['p']} nodes={r['n_nodes']} sigma=id_spn")
    print("=" * 108)
    print(f"Global CPA queries used for all nodes together: {r['global_queries']}")
    for row in r["rows"]:
        extra = ""
        if "decompositions" in row:
            extra = f"  decompositions={row['decompositions']}"
        print(f"node={row['node']:>2}  match={row['match']}{extra}")


def main():
    print("=" * 108)
    print("EXPERIMENT 18: V6 NODE DECOUPLING AUDIT")
    print("=" * 108)
    print("Checks whether node-local effective matrices are still recoverable independently.")

    rows = [
        run_config(d=2, p=17, n_nodes=6, seed=18017),
        run_config(d=2, p=23, n_nodes=8, seed=18023),
    ]
    for r in rows:
        print_result(r)

    print("\nFindings:")
    print("- If every node-local M_eff(v) is recoverable from the same global queries,")
    print("  the CPA surface remains naturally decomposable by node.")


if __name__ == "__main__":
    main()
