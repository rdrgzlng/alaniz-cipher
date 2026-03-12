#!/usr/bin/env python3
"""
Experiment 7: Collision Table for g(y)=y+B*sigma(y)
====================================================
Reproducible collision analysis for the local encryption map:

    g(y) = y + B * sigma(y),  y in F_p^d

Metrics computed per (p, d, sigma):
  - image_fraction:        |Im(g_B)| / p^d  (surjectivity measure)
  - non_surjective_fraction: 1 - image_fraction  (complement)
  - collision_prob:        Pr[g(y1)=g(y2)] for uniform independent y1,y2
                           = sum_c (|g^{-1}(c)| / p^d)^2
  - avg_preimages:         p^d / |Im(g_B)|  (mean multiplicity)
  - max_preimages:         maximum preimage count over all outputs

Each configuration is averaged over b_trials random B matrices so that
results are not skewed by a single atypical B.

Part 1 (d=1, d=2): exact enumeration of F_p^d.
Part 2 (d=4):      random birthday-probe with sample count scaled to
                   3 * ceil(sqrt(p^d)) to guarantee high collision
                   probability regardless of prime size.

Seed convention: seed = 7000 + 100*p + 10*d + sigma_offset
                 (consistent across parts and dimensions).

Note: these collision properties affect decryption ambiguity but do not
determine CPA security. Experiment 10 shows a direct CPA break in O(d)
queries via the scaling homogeneity of sigma, independently of g_B's
collision structure.

Usage:
    python3 experiments/07_collision_table.py
"""

import math
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


def _collision_prob(preimages: dict, total: int) -> float:
    """
    True pairwise collision probability: Pr[g(y1) = g(y2)] for uniform y1, y2.

    = sum_c (|g^{-1}(c)| / total)^2

    Distinct from (total - |Im|) / total, which is the non-surjective fraction.
    """
    return sum((len(ys) / total) ** 2 for ys in preimages.values())


def local_collision_stats(p: int, d: int, sigma_name: str, seed: int,
                          b_trials: int = 8):
    """
    Exact collision stats of g_B over F_p^d, averaged over b_trials random B.

    Returns aggregate metrics and the collision example from the B with the
    lowest image fraction (worst-case surjectivity).
    """
    Fp = FiniteField(p)
    sigma = Sigma(sigma_name, Fp)
    rng = random.Random(seed)
    total = p ** d

    image_fractions = []
    collision_probs = []
    avg_mults = []
    max_mults = []

    worst_image_frac = float("inf")
    worst_example = None  # collision example from worst B

    for _ in range(b_trials):
        B = Fp.random_gl(d, rng=rng)

        preimages = {}
        for n in range(total):
            y = _int_to_vec(n, p, d)
            c = (y + Fp.mat_vec(B, sigma(y))) % p
            ct = tuple(int(v) for v in c)
            preimages.setdefault(ct, []).append(tuple(int(v) for v in y))

        unique = len(preimages)
        image_frac = unique / total

        image_fractions.append(image_frac)
        collision_probs.append(_collision_prob(preimages, total))
        avg_mults.append(total / unique)
        max_mults.append(max(len(v) for v in preimages.values()))

        # Track the B with the lowest image fraction for the collision example.
        if image_frac < worst_image_frac:
            worst_image_frac = image_frac
            example = None
            for ct, ys in preimages.items():
                if len(ys) > 1:
                    # Verify algebraic structure: Δy = B * (sigma(y2) - sigma(y1))
                    y1 = np.array(ys[0], dtype=int)
                    y2 = np.array(ys[1], dtype=int)
                    delta_y = (y1 - y2) % p
                    delta_sigma = (sigma(y2) - sigma(y1)) % p
                    predicted = Fp.mat_vec(B, delta_sigma) % p
                    structure_ok = np.array_equal(delta_y, predicted)
                    example = (ys[0], ys[1], ct, structure_ok)
                    break
            worst_example = example

    return {
        "p": p,
        "d": d,
        "sigma": sigma_name,
        "b_trials": b_trials,
        "total_inputs": total,
        "mean_image_fraction": float(np.mean(image_fractions)),
        "std_image_fraction": float(np.std(image_fractions)),
        "mean_non_surjective_fraction": float(np.mean([1 - f for f in image_fractions])),
        "mean_collision_prob": float(np.mean(collision_probs)),
        "mean_avg_preimages": float(np.mean(avg_mults)),
        "mean_max_preimages": float(np.mean(max_mults)),
        "worst_collision_example": worst_example,
    }


def random_collision_probe(p: int, d: int, sigma_name: str, seed: int):
    """
    Random birthday-probe collision search for d=4.

    Sample count is scaled to 3 * ceil(sqrt(p^d)) so that the expected
    number of collisions exceeds 1 regardless of p, even when |Im(g_B)| ~ p^d.

    Also verifies the algebraic condition Δy = B*(sigma(y2)-sigma(y1)) on
    any collision found.
    """
    Fp = FiniteField(p)
    sigma = Sigma(sigma_name, Fp)
    rng = random.Random(seed)
    B = Fp.random_gl(d, rng=rng)

    domain_size = p ** d
    # Birthday bound: collisions expected after ~sqrt(|Im|) samples.
    # Use 3x sqrt(domain_size) as a conservative upper bound (works even if
    # image fraction is close to 1).
    n_samples = 3 * math.isqrt(domain_size)

    seen = {}
    found = None

    for i in range(n_samples):
        y = np.array([rng.randint(0, p - 1) for _ in range(d)], dtype=int)
        c = tuple(int(v) for v in ((y + Fp.mat_vec(B, sigma(y))) % p))
        y_tuple = tuple(int(v) for v in y)

        prev = seen.get(c)
        if prev is not None and prev != y_tuple:
            y1 = np.array(prev, dtype=int)
            y2 = np.array(y_tuple, dtype=int)
            delta_y = (y1 - y2) % p
            delta_sigma = (sigma(y2) - sigma(y1)) % p
            predicted = Fp.mat_vec(B, delta_sigma) % p
            structure_ok = np.array_equal(delta_y, predicted)
            found = (prev, y_tuple, c, i + 1, structure_ok)
            break
        seen[c] = y_tuple

    return {
        "p": p,
        "d": d,
        "sigma": sigma_name,
        "birthday_bound": n_samples,
        "unique_outputs_seen": len(seen),
        "collision_found": found is not None,
        "collision_example": found,
    }


def print_collision_table(rows):
    print("\n" + "=" * 108)
    print("PART 1: COLLISION TABLE FOR g(y)=y+B*sigma(y)  [averaged over b_trials B matrices]")
    print("=" * 108)
    print(f"{'sigma':<8} {'p':>3} {'d':>3} {'B_trials':>8} {'inputs':>8} "
          f"{'img_frac':>10} {'img_std':>8} {'coll_prob':>10} {'avg_pre':>8} {'max_pre':>8}")
    print("-" * 108)
    for r in rows:
        print(f"{r['sigma']:<8} {r['p']:>3} {r['d']:>3} {r['b_trials']:>8} "
              f"{r['total_inputs']:>8} "
              f"{r['mean_image_fraction']:>10.4f} {r['std_image_fraction']:>8.4f} "
              f"{r['mean_collision_prob']:>10.6f} "
              f"{r['mean_avg_preimages']:>8.3f} {r['mean_max_preimages']:>8.1f}")

    print("\nWorst-B collision examples (lowest image fraction among b_trials):")
    print("  Algebraic condition: Δy = B*(sigma(y2)-sigma(y1)) must hold for every collision.")
    for r in rows:
        ex = r["worst_collision_example"]
        if ex is None:
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: no collision found in any trial")
        else:
            y1, y2, c, structure_ok = ex
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: "
                  f"y1={y1}, y2={y2}, c={c}  [Δy=B·Δσ: {structure_ok}]")


def print_random_probe(rows):
    print("\n" + "=" * 108)
    print("PART 2: RANDOM BIRTHDAY-PROBE (d=4, n_samples=3*ceil(sqrt(p^d)))")
    print("=" * 108)
    print(f"{'sigma':<8} {'p':>3} {'d':>3} {'bday_bound':>12} {'unique':>8} "
          f"{'found':>8} {'at_iter':>8} {'Δy=B·Δσ':>8}")
    print("-" * 108)
    for r in rows:
        found_str = "yes" if r["collision_found"] else "no"
        ex = r["collision_example"]
        at_iter = ex[3] if ex else "-"
        structure = str(ex[4]) if ex else "-"
        print(f"{r['sigma']:<8} {r['p']:>3} {r['d']:>3} {r['birthday_bound']:>12} "
              f"{r['unique_outputs_seen']:>8} {found_str:>8} {str(at_iter):>8} {structure:>8}")

    print("\nCollision examples (if found):")
    for r in rows:
        ex = r["collision_example"]
        if ex is None:
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: "
                  f"no collision in {r['birthday_bound']} samples")
        else:
            y1, y2, c, i, structure_ok = ex
            print(f"  {r['sigma']}, p={r['p']}, d={r['d']}: "
                  f"y1={y1}, y2={y2}, c={c}, iter={i}  [Δy=B·Δσ: {structure_ok}]")


def main():
    print("=" * 108)
    print("EXPERIMENT 7: COLLISION TABLE")
    print("=" * 108)

    primes = [17, 23, 29, 47, 59, 89, 131, 251]
    sigmas = ["inverse", "cube"]
    # Seed formula: 7000 + 100*p + 10*d + sigma_offset (consistent across parts).
    sigma_offset = {"inverse": 0, "cube": 1}

    # Part 1: exact stats for d=1, d=2, averaged over b_trials B matrices.
    rows = []
    for sigma_name in sigmas:
        for p in primes:
            if sigma_name == "cube" and not FiniteField(p).cube_invertible():
                continue
            for d in [1, 2]:
                seed = 7000 + 100 * p + 10 * d + sigma_offset[sigma_name]
                rows.append(local_collision_stats(p, d, sigma_name, seed, b_trials=8))

    print_collision_table(rows)

    # Part 2: birthday-scaled random probe for d=4.
    d4_rows = []
    for sigma_name in sigmas:
        for p in primes:
            if sigma_name == "cube" and not FiniteField(p).cube_invertible():
                continue
            seed = 7000 + 100 * p + 10 * 4 + sigma_offset[sigma_name]
            d4_rows.append(random_collision_probe(p, 4, sigma_name, seed))

    print_random_probe(d4_rows)

    print("\n" + "=" * 108)
    print("NOTE ON INTERPRETATION")
    print("=" * 108)
    print("- collision_prob is Pr[g(y1)=g(y2)] = sum_c (|g^{-1}(c)|/p^d)^2,")
    print("  the true pairwise collision probability. It differs from")
    print("  non_surjective_fraction = 1 - image_fraction = (p^d - |Im|)/p^d.")
    print("- Collisions in g_B cause decryption ambiguity (multiple valid plaintexts)")
    print("  but do not determine CPA security.")
    print("- Experiment 10 shows a direct CPA key recovery in O(d) queries via the")
    print("  scaling homogeneity sigma(lambda*y) = lambda^e * sigma(y), independently")
    print("  of g_B's collision structure.")


if __name__ == "__main__":
    main()
