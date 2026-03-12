#!/usr/bin/env python3
"""
Experiment 7: Collision Table for g(y)=y+B*sigma(y)
====================================================
Reproducible collision analysis for the local encryption map:

    g(y) = y + B * sigma(y),  y in F_p^d

Usage:
    python3 experiments/07_collision_table.py
"""

import sys
import os
import random
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.crypto.sigma import Sigma


def _int_to_vec(n: int, p: int, d: int) -> np.ndarray:
    vals = []
    x = n
    for _ in range(d):
        vals.append(x % p)
        x //= p
    return np.array(vals, dtype=int)


def local_collision_stats(p: int, d: int, sigma_name: str, seed: int):
    """Compute exact collision stats of g over the full input space F_p^d."""
    Fp = FiniteField(p)
    sigma = Sigma(sigma_name, Fp)
    rng = random.Random(seed)
    B = Fp.random_gl(d, rng=rng)

    preimages = {}
    total = p ** d
    for n in range(total):
        y = _int_to_vec(n, p, d)
        c = (y + Fp.mat_vec(B, sigma(y))) % p
        ct = tuple(int(v) for v in c)
        preimages.setdefault(ct, []).append(tuple(int(v) for v in y))

    unique_ciphertexts = len(preimages)
    collisions = total - unique_ciphertexts
    max_mult = max(len(v) for v in preimages.values())
    avg_mult = total / unique_ciphertexts

    example = None
    for ct, ys in preimages.items():
        if len(ys) > 1:
            example = (ys[0], ys[1], ct)
            break

    return {
        "p": p,
        "d": d,
        "sigma": sigma_name,
        "total_inputs": total,
        "unique_outputs": unique_ciphertexts,
        "collisions": collisions,
        "collision_rate": collisions / total,
        "avg_preimages": avg_mult,
        "max_preimages": max_mult,
        "example": example,
        "B00": int(B[0, 0]),
    }


def random_collision_probe(p: int, d: int, sigma_name: str, seed: int, n_samples: int):
    """
    Random collision search for larger dimensions.

    Samples random y values and reports the first collision found.
    """
    Fp = FiniteField(p)
    sigma = Sigma(sigma_name, Fp)
    rng = random.Random(seed)
    B = Fp.random_gl(d, rng=rng)

    seen = {}
    found = None

    for i in range(n_samples):
        y = np.array([rng.randint(0, p - 1) for _ in range(d)], dtype=int)
        c = tuple(int(v) for v in ((y + Fp.mat_vec(B, sigma(y))) % p))
        y_tuple = tuple(int(v) for v in y)

        prev = seen.get(c)
        if prev is not None and prev != y_tuple:
            found = (prev, y_tuple, c, i + 1)
            break
        seen[c] = y_tuple

    return {
        "p": p,
        "d": d,
        "sigma": sigma_name,
        "samples": n_samples,
        "unique_outputs_seen": len(seen),
        "collision_found": found is not None,
        "collision_example": found,
        "B00": int(B[0, 0]),
    }


def print_collision_table(rows):
    print("\n" + "=" * 92)
    print("PART 1: COLLISION TABLE FOR g(y)=y+B*sigma(y)")
    print("=" * 92)
    print(f"{'sigma':<8} {'p':>3} {'d':>3} {'inputs':>8} {'unique':>8} "
          f"{'coll':>8} {'coll%':>8} {'avg_pre':>8} {'max_pre':>8}")
    print("-" * 92)
    for r in rows:
        print(f"{r['sigma']:<8} {r['p']:>3} {r['d']:>3} {r['total_inputs']:>8} "
              f"{r['unique_outputs']:>8} {r['collisions']:>8} "
              f"{100 * r['collision_rate']:>7.2f}% {r['avg_preimages']:>8.3f} "
              f"{r['max_preimages']:>8}")

    print("\nExamples of collisions (first found per row):")
    for r in rows:
        ex = r["example"]
        if ex is None:
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: no collision found")
        else:
            y1, y2, c = ex
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: y1={y1}, y2={y2}, c={c}")


def print_random_probe(rows):
    print("\n" + "=" * 92)
    print("PART 2: RANDOM COLLISION PROBE (d=4)")
    print("=" * 92)
    print(f"{'sigma':<8} {'p':>3} {'d':>3} {'samples':>8} {'unique':>8} {'found':>8} {'iter':>8}")
    print("-" * 92)
    for r in rows:
        found = "yes" if r["collision_found"] else "no"
        at_iter = r["collision_example"][3] if r["collision_example"] else "-"
        print(f"{r['sigma']:<8} {r['p']:>3} {r['d']:>3} {r['samples']:>8} "
              f"{r['unique_outputs_seen']:>8} {found:>8} {str(at_iter):>8}")

    print("\nExamples (if collision found):")
    for r in rows:
        ex = r["collision_example"]
        if ex is None:
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: no collision in {r['samples']} samples")
        else:
            y1, y2, c, i = ex
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: y1={y1}, y2={y2}, c={c}, iter={i}")


def main():
    print("=" * 92)
    print("EXPERIMENT 7: COLLISION TABLE")
    print("=" * 92)

    primes = [17, 23, 29, 47, 59, 89, 131, 251]
    dimensions = [1, 2]
    sigmas = ["inverse", "cube"]

    rows = []
    for sigma_name in sigmas:
        for p in primes:
            Fp = FiniteField(p)
            if sigma_name == "cube" and not Fp.cube_invertible():
                continue
            for d in dimensions:
                seed = 1000 + 100 * p + 10 * d + (0 if sigma_name == "inverse" else 1)
                rows.append(local_collision_stats(p, d, sigma_name, seed))

    print_collision_table(rows)

    d4_rows = []
    d4_samples = 80000
    for sigma_name in sigmas:
        for p in primes:
            Fp = FiniteField(p)
            if sigma_name == "cube" and not Fp.cube_invertible():
                continue
            seed = 900000 + 100 * p + (0 if sigma_name == "inverse" else 1)
            d4_rows.append(
                random_collision_probe(p, 4, sigma_name, seed=seed, n_samples=d4_samples)
            )

    print_random_probe(d4_rows)


if __name__ == "__main__":
    main()
