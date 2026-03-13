#!/usr/bin/env python3
"""
Experiment 15: V6 SPN differential audit
=======================================

Goal:
  Measure exact differential biases of sigma=id_spn for small parameters.

Method:
  - For each nonzero input difference Delta, enumerate all x in F_p^d.
  - Count the distribution of sigma(x + Delta) - sigma(x).
  - Report differential uniformity, support size, and the strongest trails.
"""

from __future__ import annotations

import os
import sys
from collections import defaultdict
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.crypto.sigma import Sigma


def all_vectors(p, d):
    for vals in product(range(p), repeat=d):
        yield np.array(vals, dtype=int)


def differential_profile(p, d):
    Fp = FiniteField(p)
    sigma = Sigma("id_spn", Fp, d=d)
    best = None
    axis_rows = []
    total_nonzero = p ** d - 1

    for delta in all_vectors(p, d):
        if not np.any(delta):
            continue
        counts = defaultdict(int)
        for x in all_vectors(p, d):
            dy = tuple(int(v) for v in (sigma((x + delta) % p) - sigma(x)) % p)
            counts[dy] += 1

        max_count = max(counts.values())
        support = len(counts)
        row = {
            "delta": delta.copy(),
            "max_count": max_count,
            "max_prob": max_count / (p ** d),
            "support": support,
            "is_axis": int(np.count_nonzero(delta) == 1),
        }
        if best is None or row["max_count"] > best["max_count"]:
            most_likely = max(counts.items(), key=lambda kv: kv[1])
            best = {
                **row,
                "best_output_diff": np.array(most_likely[0], dtype=int),
            }
        if row["is_axis"]:
            axis_rows.append(row)

    axis_max = max(r["max_count"] for r in axis_rows)
    axis_support_min = min(r["support"] for r in axis_rows)
    return {
        "p": p,
        "d": d,
        "uniformity": best["max_count"],
        "uniformity_prob": best["max_prob"],
        "best_delta": best["delta"],
        "best_output_diff": best["best_output_diff"],
        "axis_max_count": axis_max,
        "axis_max_prob": axis_max / (p ** d),
        "axis_support_min": axis_support_min,
        "nonzero_deltas": total_nonzero,
    }


def print_result(r):
    print("\n" + "=" * 108)
    print(f"CONFIG: d={r['d']} p={r['p']} sigma=id_spn")
    print("=" * 108)
    print(f"Differential uniformity: {r['uniformity']} / {r['p'] ** r['d']} "
          f"= {r['uniformity_prob']:.6f}")
    print(f"Best input diff:  {r['best_delta']}")
    print(f"Best output diff: {r['best_output_diff']}")
    print(f"Axis-aligned max count: {r['axis_max_count']} "
          f"(prob={r['axis_max_prob']:.6f})")
    print(f"Minimum support size for axis-aligned differences: {r['axis_support_min']}")


def main():
    print("=" * 108)
    print("EXPERIMENT 15: V6 SPN DIFFERENTIAL AUDIT")
    print("=" * 108)
    print("Exact DDT-style audit for sigma=id_spn on small dimensions.")

    rows = [
        differential_profile(p=17, d=2),
        differential_profile(p=23, d=2),
    ]
    for r in rows:
        print_result(r)

    print("\nFindings:")
    print("- Any large-probability differential here is a concrete attack lead on V6.")
    print("- Axis-aligned input differences are reported separately because they are")
    print("  the natural continuation of the old column-wise CPA analysis.")


if __name__ == "__main__":
    main()
