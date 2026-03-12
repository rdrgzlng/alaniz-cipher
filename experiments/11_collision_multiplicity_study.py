#!/usr/bin/env python3
"""
Experiment 11: Deep Collision/Multi-Preimage Study of g_B(y)
=============================================================
For g_B(y) = y + B*sigma(y), this experiment measures:

1) Mean image fraction |Im(g_B)| / p^d over random B.
2) Mean preimage multiplicity p^d / |Im(g_B)|.
3) Distribution of output multiplicities (1,2,3,... preimages).
4) "Bad B" instances (lowest image fraction), per (p,d,sigma).

Usage:
    python3 experiments/11_collision_multiplicity_study.py
"""

import sys
import os
import random
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField


def sigma_apply(y: np.ndarray, p: int, sigma_name: str) -> np.ndarray:
    if sigma_name == "inverse":
        return np.array([pow(int(v), p - 2, p) if int(v) % p != 0 else 0 for v in y], dtype=int)
    if sigma_name == "cube":
        return np.array([pow(int(v), 3, p) for v in y], dtype=int)
    raise ValueError(f"Unknown sigma {sigma_name}")


def int_to_vec(n: int, p: int, d: int) -> np.ndarray:
    vals = [0] * d
    x = n
    for i in range(d):
        vals[i] = x % p
        x //= p
    return np.array(vals, dtype=int)


def map_stats_exact(B: np.ndarray, p: int, d: int, sigma_name: str):
    """Exact image/multiplicity stats by full enumeration of F_p^d."""
    counts = defaultdict(int)
    total = p ** d
    for n in range(total):
        y = int_to_vec(n, p, d)
        z = sigma_apply(y, p, sigma_name)
        c = (y + (B @ z) % p) % p
        counts[tuple(int(v) for v in c)] += 1

    unique_outputs = len(counts)
    multiplicity_hist = defaultdict(int)
    for m in counts.values():
        multiplicity_hist[m] += 1

    max_mult = max(counts.values())
    image_fraction = unique_outputs / total
    avg_preimages = total / unique_outputs
    return {
        "unique_outputs": unique_outputs,
        "image_fraction": image_fraction,
        "avg_preimages": avg_preimages,
        "max_preimages": max_mult,
        "multiplicity_hist": dict(multiplicity_hist),
    }


def map_stats_sampled(B: np.ndarray, p: int, d: int, sigma_name: str, y_samples: int, rng):
    """
    Approximate stats for larger spaces via random y sampling.
    Tracks repeated outputs among sampled points.
    """
    counts = defaultdict(int)
    for _ in range(y_samples):
        y = np.array([rng.randint(0, p - 1) for _ in range(d)], dtype=int)
        z = sigma_apply(y, p, sigma_name)
        c = (y + (B @ z) % p) % p
        counts[tuple(int(v) for v in c)] += 1

    unique_outputs = len(counts)
    multiplicity_hist = defaultdict(int)
    for m in counts.values():
        multiplicity_hist[m] += 1
    max_mult = max(counts.values()) if counts else 0

    # In sampled mode these are sample-space metrics, not full-space exact values.
    unique_fraction_sample = unique_outputs / y_samples if y_samples > 0 else 0.0
    avg_preimages_sample = y_samples / unique_outputs if unique_outputs > 0 else float("inf")
    return {
        "unique_outputs": unique_outputs,
        "image_fraction": unique_fraction_sample,
        "avg_preimages": avg_preimages_sample,
        "max_preimages": max_mult,
        "multiplicity_hist": dict(multiplicity_hist),
    }


def summarize_hist(hists: list[dict[int, int]]):
    agg = defaultdict(int)
    for h in hists:
        for m, cnt in h.items():
            agg[m] += cnt
    total_outputs = sum(agg.values())
    if total_outputs == 0:
        return {}
    return {m: agg[m] / total_outputs for m in sorted(agg)}


def run_config(p: int, d: int, sigma_name: str, b_samples: int, exact: bool, y_samples: int, seed: int):
    Fp = FiniteField(p)
    if sigma_name == "cube" and not Fp.cube_invertible():
        return None

    rng = random.Random(seed)
    image_fractions = []
    avg_preimages = []
    max_preimages = []
    hists = []
    per_B = []

    for i in range(b_samples):
        B = Fp.random_gl(d, rng=rng)
        if exact:
            s = map_stats_exact(B, p, d, sigma_name)
        else:
            s = map_stats_sampled(B, p, d, sigma_name, y_samples, rng)

        image_fractions.append(s["image_fraction"])
        avg_preimages.append(s["avg_preimages"])
        max_preimages.append(s["max_preimages"])
        hists.append(s["multiplicity_hist"])
        per_B.append((s["image_fraction"], B.copy(), s["max_preimages"]))

    fr = np.array(image_fractions, dtype=float)
    ap = np.array(avg_preimages, dtype=float)
    mp = np.array(max_preimages, dtype=float)

    q10 = float(np.quantile(fr, 0.10))
    bad_idx = [i for i, v in enumerate(fr) if v <= q10]
    bad_rate = len(bad_idx) / len(fr) if len(fr) else 0.0

    per_B_sorted = sorted(per_B, key=lambda t: t[0])
    worst = per_B_sorted[: min(3, len(per_B_sorted))]
    worst_serialized = [
        {
            "image_fraction": float(v),
            "max_preimages": int(m),
            "B": np.array(B, dtype=int).tolist(),
        }
        for v, B, m in worst
    ]

    return {
        "p": p,
        "d": d,
        "sigma": sigma_name,
        "exact": exact,
        "b_samples": b_samples,
        "y_samples": y_samples if not exact else p ** d,
        "mean_image_fraction": float(np.mean(fr)),
        "std_image_fraction": float(np.std(fr)),
        "mean_avg_preimages": float(np.mean(ap)),
        "std_avg_preimages": float(np.std(ap)),
        "mean_max_preimages": float(np.mean(mp)),
        "q10_image_fraction": q10,
        "bad_rate_q10": bad_rate,
        "agg_multiplicity_hist": summarize_hist(hists),
        "worst_B_examples": worst_serialized,
    }


def print_summary_table(rows):
    print("\n" + "=" * 120)
    print("PART 1: SUMMARY TABLE")
    print("=" * 120)
    print(f"{'sigma':<8} {'p':>4} {'d':>3} {'mode':>8} {'B_samp':>7} {'Y_samp':>8} "
          f"{'img_mean':>10} {'img_std':>9} {'pre_mean':>10} {'maxpre':>8} {'bad_q10':>9}")
    print("-" * 120)
    for r in rows:
        mode = "exact" if r["exact"] else "sample"
        print(f"{r['sigma']:<8} {r['p']:>4} {r['d']:>3} {mode:>8} {r['b_samples']:>7} {r['y_samples']:>8} "
              f"{r['mean_image_fraction']:>10.4f} {r['std_image_fraction']:>9.4f} "
              f"{r['mean_avg_preimages']:>10.4f} {r['mean_max_preimages']:>8.2f} {100*r['bad_rate_q10']:>8.2f}%")


def print_histograms(rows):
    print("\n" + "=" * 120)
    print("PART 2: AGGREGATED MULTIPLICITY HISTOGRAMS (TOP MULTIPLICITIES)")
    print("=" * 120)
    for r in rows:
        hist = r["agg_multiplicity_hist"]
        top = sorted(hist.items(), key=lambda kv: kv[1], reverse=True)[:6]
        mode = "exact" if r["exact"] else "sample"
        top_s = ", ".join([f"m={m}:{100*v:.2f}%" for m, v in top]) if top else "-"
        print(f"{r['sigma']}, p={r['p']}, d={r['d']}, mode={mode} -> {top_s}")


def print_worst_B(rows):
    print("\n" + "=" * 120)
    print("PART 3: WORST-B EXAMPLES (LOWEST IMAGE FRACTION)")
    print("=" * 120)
    for r in rows:
        mode = "exact" if r["exact"] else "sample"
        print(f"\n{r['sigma']} p={r['p']} d={r['d']} mode={mode}")
        for w in r["worst_B_examples"]:
            print(f"  img={w['image_fraction']:.4f}, max_pre={w['max_preimages']}, B={w['B']}")


def main():
    print("=" * 120)
    print("EXPERIMENT 11: DEEP COLLISION/MULTIPREIMAGE STUDY OF g_B")
    print("=" * 120)

    primes = [17, 23, 29, 47, 59, 89, 131, 251]
    sigmas = ["inverse", "cube"]

    rows = []

    # d=2 exact campaign across all primes.
    for sigma_name in sigmas:
        for p in primes:
            exact = True
            b_samples = 24 if p <= 89 else 12
            cfg = run_config(
                p=p,
                d=2,
                sigma_name=sigma_name,
                b_samples=b_samples,
                exact=exact,
                y_samples=0,
                seed=11000 + p + (0 if sigma_name == "inverse" else 7000),
            )
            if cfg is not None:
                rows.append(cfg)

    # d=4 sampled campaign (exact would be too large).
    for sigma_name in sigmas:
        for p in [17, 23, 29, 47, 89, 251]:
            cfg = run_config(
                p=p,
                d=4,
                sigma_name=sigma_name,
                b_samples=12,
                exact=False,
                y_samples=30000,
                seed=21000 + p + (0 if sigma_name == "inverse" else 9000),
            )
            if cfg is not None:
                rows.append(cfg)

    print_summary_table(rows)
    print_histograms(rows)
    print_worst_B(rows)


if __name__ == "__main__":
    main()
