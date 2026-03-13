#!/usr/bin/env python3
"""
Experiment 17: V6 degeneracy audit
=================================

Goal:
  Look for structural degenerations introduced by the constant and mixing layer.

Metrics:
  - sigma(0)
  - number of fixed points of sigma
  - number of invariant 1-D linear subspaces
  - fraction of points with singular discrete Jacobian
  - fraction of points where some output component equals the input component
"""

from __future__ import annotations

import os
import sys
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.crypto.sigma import Sigma


def all_vectors(p, d):
    for vals in product(range(p), repeat=d):
        yield np.array(vals, dtype=int)


def normalized_lines(p):
    seen = set()
    lines = []
    for a, b in product(range(p), repeat=2):
        if a == 0 and b == 0:
            continue
        vec = np.array([a, b], dtype=int)
        nz = next(i for i, v in enumerate(vec) if v % p != 0)
        inv = pow(int(vec[nz]), p - 2, p)
        key = tuple(int((inv * v) % p) for v in vec)
        if key not in seen:
            seen.add(key)
            lines.append(np.array(key, dtype=int))
    return lines


def discrete_jacobian(sigma, x, p):
    d = len(x)
    cols = []
    for j in range(d):
        ej = np.zeros(d, dtype=int)
        ej[j] = 1
        cols.append((sigma((x + ej) % p) - sigma(x)) % p)
    return np.column_stack(cols) % p


def run_config(p):
    Fp = FiniteField(p)
    sigma = Sigma("id_spn", Fp, d=2)
    sigma_zero = sigma(np.zeros(2, dtype=int)) % p

    fixed_points = 0
    singular_j = 0
    component_fixed = 0
    total = p * p

    for x in all_vectors(p, 2):
        y = sigma(x) % p
        if np.array_equal(y, x % p):
            fixed_points += 1
        if np.any(y == (x % p)):
            component_fixed += 1
        J = discrete_jacobian(sigma, x, p)
        det = int((J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]) % p)
        if det == 0:
            singular_j += 1

    invariant_lines = 0
    for v in normalized_lines(p):
        closed = True
        for t in range(p):
            x = (t * v) % p
            y = sigma(x) % p
            if not any(np.array_equal(y, (u * v) % p) for u in range(p)):
                closed = False
                break
        if closed:
            invariant_lines += 1

    return {
        "p": p,
        "sigma_zero": sigma_zero,
        "fixed_points": fixed_points,
        "invariant_lines": invariant_lines,
        "singular_j_fraction": singular_j / total,
        "component_fixed_fraction": component_fixed / total,
    }


def print_result(r):
    print("\n" + "=" * 108)
    print(f"CONFIG: d=2 p={r['p']} sigma=id_spn")
    print("=" * 108)
    print(f"sigma(0) = {r['sigma_zero']}")
    print(f"fixed points = {r['fixed_points']}")
    print(f"invariant 1-D linear subspaces = {r['invariant_lines']}")
    print(f"singular discrete Jacobian fraction = {r['singular_j_fraction']:.6f}")
    print(f"component-fixed fraction = {r['component_fixed_fraction']:.6f}")


def main():
    print("=" * 108)
    print("EXPERIMENT 17: V6 DEGENERACY AUDIT")
    print("=" * 108)
    print("Searches for low-rank, invariant-subspace, and fixed-point pathologies.")

    rows = [run_config(p) for p in [5, 7, 11, 17, 23]]
    for r in rows:
        print_result(r)

    print("\nFindings:")
    print("- Any nontrivial invariant line, unusually large fixed-point set, or high")
    print("  singular-Jacobian rate is a concrete structural warning for V6.")


if __name__ == "__main__":
    main()
