#!/usr/bin/env python3
"""
Experiment 10: CPA Key Recovery via Scaling Closed Form
=======================================================
Attack vector:
  For the final local encryption law
      c = A*s + B*sigma(A*s)
  and any sigma satisfying a scaling law
      sigma(lambda * y) = lambda^e * sigma(y),
  a chosen-plaintext attacker can query the pair (s, lambda*s) and solve
  directly for y = A*s and w = B*sigma(A*s):

      c1 = y + w
      c2 = lambda*y + lambda^e*w

  Hence, if Delta = lambda^e - lambda != 0 mod p,

      y = (lambda^e * c1 - c2) / Delta
      w = (-lambda * c1 + c2) / Delta

  Repeating this on the basis vectors e_j recovers all columns A*e_j and thus A.
  Once A is known, B is recovered linearly from

      c - A*s = B*sigma(A*s)

  This is a direct CPA break for sigma=cube and sigma=inverse, avoiding any
  exhaustive search over GL(d, F_p).

Usage:
    python3 experiments/10_cpa_scaling_closed_form_attack.py
"""

import sys
import os
import random
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


# -----------------------------------------------------------------------------
# Small linear-algebra helpers over F_p
# -----------------------------------------------------------------------------

def solve_linear_Fp(equations, n_unknowns, p):
    """Solve/rank a linear system over F_p via Gauss-Jordan."""
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
        solution = np.zeros(n_unknowns, dtype=int)
        for i, pc in enumerate(pivots):
            solution[pc] = int(aug[i, n_unknowns]) % p

    return len(pivots), inconsistent, solution


def mat_from_columns(cols, p):
    """Build matrix from a list of column vectors."""
    return np.column_stack([np.array(c, dtype=int) % p for c in cols]) % p


def coeff_to_node_map(sheaf, node: int):
    """
    Matrix P such that s_node = P * coeffs, where coeffs are H^0 basis coords.
    """
    p = sheaf.Fp.p
    d = sheaf.dv
    cols = []
    for bv in sheaf.H0_basis:
        cols.append(sheaf.get_node_value(bv.astype(int), node) % p)
    return mat_from_columns(cols, p)


# -----------------------------------------------------------------------------
# Sigma scaling model
# -----------------------------------------------------------------------------

def sigma_scaling_exponent(p: int, sigma_name: str) -> int:
    """
    Return the exponent e in sigma(lambda*y) = lambda^e sigma(y), componentwise.

    For the two public sigma maps discussed in the paper and scripts:
      - cube:    e = 3
      - inverse: e = -1  (mod p-1 in multiplicative notation)
    """
    if sigma_name == "cube":
        return 3
    if sigma_name == "inverse":
        return -1
    raise ValueError(f"Unsupported sigma for scaling attack: {sigma_name}")


def scalar_power_mod(lam: int, e: int, p: int) -> int:
    """Compute lam^e mod p for integer e, allowing e < 0."""
    lam %= p
    if lam == 0:
        raise ValueError("lambda must be nonzero in F_p")
    if e >= 0:
        return pow(lam, e, p)
    inv = pow(lam, p - 2, p)
    return pow(inv, -e, p)


def choose_lambda(p: int, sigma_name: str, preferred=(2, 3, 5, 7, 11, 13)):
    """
    Pick a nonzero scalar lambda such that Delta = lambda^e - lambda != 0 mod p.
    This ensures the 2x2 system can be inverted.
    """
    e = sigma_scaling_exponent(p, sigma_name)
    for lam in preferred:
        if lam % p == 0:
            continue
        lam_e = scalar_power_mod(lam, e, p)
        delta = (lam_e - lam) % p
        if delta != 0:
            return lam % p, lam_e, delta
    for lam in range(1, p):
        lam_e = scalar_power_mod(lam, e, p)
        delta = (lam_e - lam) % p
        if delta != 0:
            return lam, lam_e, delta
    raise RuntimeError(f"No usable lambda found for p={p}, sigma={sigma_name}")


# -----------------------------------------------------------------------------
# Core CPA oracle helpers at one node
# -----------------------------------------------------------------------------

def encrypt_root_value(proto: Protocol, key, root_value, node: int = 0):
    """
    Build a global section from a d-vector of coefficients, encrypt it, and
    return the ciphertext value at the selected node.

    IMPORTANT:
      - In the current repo, Sheaf exposes ``section_from_coeffs`` (not
        ``section_from_root``).
      - On the path/tree constructions used in these experiments, H^0 has
        dimension d and a coefficient vector gives a valid CPA plaintext that
        we can scale as needed for the attack.

    This helper is written defensively so it works both on the current repo and
    on older variants that may still expose ``section_from_root``.
    """
    coeffs = np.array(root_value, dtype=int) % proto.Fp.p

    if hasattr(proto.sheaf, "section_from_coeffs"):
        s = proto.sheaf.section_from_coeffs(coeffs)
    elif hasattr(proto.sheaf, "section_from_root"):
        s = proto.sheaf.section_from_root(coeffs)
    else:
        raise AttributeError(
            "Sheaf has neither section_from_coeffs nor section_from_root; "
            "adapt encrypt_root_value() to the repo API"
        )

    c = proto.encrypt(s, key)
    return proto.sheaf.get_node_value(c, node) % proto.Fp.p


def recover_A_closed_form(proto: Protocol, key, d: int, sigma_name: str,
                          node: int = 0, lambda_hint=None):
    """
    Recover A column-by-column using CPA queries on e_j and lambda*e_j.

    Returns:
      A_rec, details_dict
    """
    p = proto.Fp.p
    lam, lam_e, delta = choose_lambda(p, sigma_name) if lambda_hint is None else lambda_hint
    delta_inv = pow(int(delta), p - 2, p)

    cols = []
    per_col = []
    for j in range(d):
        basis = np.zeros(d, dtype=int)
        basis[j] = 1
        s1 = basis.copy()
        s2 = (lam * basis) % p

        c1 = encrypt_root_value(proto, key, s1, node=node)
        c2 = encrypt_root_value(proto, key, s2, node=node)

        # y = (lam^e * c1 - c2) / (lam^e - lam)
        y = ((lam_e * c1 - c2) * delta_inv) % p
        # optional recovered nonlinear part for diagnostics
        w = ((-lam * c1 + c2) * delta_inv) % p

        cols.append(y)
        per_col.append({
            "j": j,
            "s": s1,
            "lambda_s": s2,
            "c1": c1,
            "c2": c2,
            "A_col": y,
            "w_col": w,
        })

    A_rec = mat_from_columns(cols, p)
    return A_rec, {
        "lambda": lam,
        "lambda_e": lam_e,
        "delta": delta,
        "queries_for_A": 2 * d,
        "per_column": per_col,
    }


def recover_B_from_A(proto: Protocol, key, A_rec, d: int, node: int = 0,
                     max_trials: int = 64, rng=None):
    """
    Recover B once A is known, by collecting d CPA queries whose z_j = sigma(A*s_j)
    form an invertible matrix, then solving W = B Z.
    """
    if rng is None:
        rng = random.Random(0)

    Fp = proto.Fp
    p = Fp.p
    trials = []

    # Try basis first (deterministic and cheap), then random vectors if needed.
    candidates = []
    for j in range(d):
        e = np.zeros(d, dtype=int)
        e[j] = 1
        candidates.append(e)
    for _ in range(max(0, max_trials - d)):
        candidates.append(np.array([rng.randrange(p) for _ in range(d)], dtype=int))

    for s in candidates:
        c = encrypt_root_value(proto, key, s, node=node)
        y = Fp.mat_vec(A_rec, s) % p
        z = proto.sigma(y) % p
        w = (c - y) % p
        trials.append((s, c, y, z, w))

        if len(trials) < d:
            continue

        selected = trials[-d:]
        Z = mat_from_columns([z for (_, _, _, z, _) in selected], p)
        try:
            Z_inv = Fp.mat_inv(Z)
        except ValueError:
            continue

        W = mat_from_columns([w for (_, _, _, _, w) in selected], p)
        B_rec = (W @ Z_inv) % p
        return B_rec, {
            "queries_for_B": len(trials),
            "selected": selected,
            "Z": Z,
        }

    raise RuntimeError("Failed to find d queries with invertible Z while recovering B")


def validate_key_on_random_queries(proto: Protocol, key, A_rec, B_rec,
                                   n_tests: int = 16, seed: int = 0, node: int = 0):
    """Validate recovered (A, B) on fresh random CPA-style local queries."""
    rng = random.Random(seed)
    Fp = proto.Fp
    p = Fp.p

    for _ in range(n_tests):
        s = np.array([rng.randrange(p) for _ in range(A_rec.shape[1])], dtype=int)
        c = encrypt_root_value(proto, key, s, node=node)
        y = Fp.mat_vec(A_rec, s) % p
        z = proto.sigma(y) % p
        c_hat = (y + Fp.mat_vec(B_rec, z)) % p
        if not np.array_equal(c_hat, c % p):
            return False
    return True


# -----------------------------------------------------------------------------
# One full experiment instance
# -----------------------------------------------------------------------------

def run_attack(d: int, p: int, sigma_name: str, seed: int):
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, sigma_name))
    key = proto.keygen(rng=rng)

    node = 0
    A_true = key.A[node] % p
    B_true = key.B[node] % p
    P = coeff_to_node_map(sheaf, node) % p

    t0 = time.time()
    A_rec, a_info = recover_A_closed_form(proto, key, d, sigma_name, node=node)
    tA = time.time() - t0

    t1 = time.time()
    B_rec, b_info = recover_B_from_A(proto, key, A_rec, d, node=node, rng=rng)
    tB = time.time() - t1

    total = time.time() - t0

    # A_rec is in coefficient coordinates: y = A_rec * coeffs.
    # True local key acts on stalk coords: y = A_true * s_node with s_node = P*coeffs.
    # Therefore expected relation is A_rec = A_true * P.
    A_eff_true = (A_true @ P) % p
    A_eff_match = np.array_equal(A_rec % p, A_eff_true % p)

    A_stalk_match = False
    A_rec_stalk = None
    try:
        P_inv = Fp.mat_inv(P)
        A_rec_stalk = (A_rec @ P_inv) % p
        A_stalk_match = np.array_equal(A_rec_stalk % p, A_true % p)
    except ValueError:
        # Non-invertible P should not happen on the tree configs used here.
        pass

    B_match = np.array_equal(B_rec % p, B_true % p)
    valid = validate_key_on_random_queries(proto, key, A_rec, B_rec,
                                           n_tests=20, seed=seed + 999, node=node)

    return {
        "d": d,
        "p": p,
        "sigma": sigma_name,
        "lambda": a_info["lambda"],
        "delta": a_info["delta"],
        "qA": a_info["queries_for_A"],
        "qB": b_info["queries_for_B"],
        "qTotal": a_info["queries_for_A"] + b_info["queries_for_B"],
        "A_eff_match": A_eff_match,
        "A_stalk_match": A_stalk_match,
        "B_match": B_match,
        "valid": valid,
        "time_A_s": tA,
        "time_B_s": tB,
        "time_total_s": total,
        "A_true": A_true,
        "A_eff_true": A_eff_true,
        "A_rec": A_rec,
        "A_rec_stalk": A_rec_stalk,
        "B_true": B_true,
        "B_rec": B_rec,
    }


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def main():
    print("=" * 108)
    print("EXPERIMENT 10: CPA KEY RECOVERY VIA SCALING CLOSED FORM")
    print("=" * 108)
    print("Breaks c = A*s + B*sigma(A*s) for homogeneous sigma using paired CPA queries (s, lambda*s)")

    primes = [17, 23, 29, 47, 59, 89, 131, 251]
    sigmas = ["inverse", "cube"]
    dimensions = [2, 4]

    for d in dimensions:
        for sigma_name in sigmas:
            print("\n" + "=" * 108)
            print(f"SIGMA={sigma_name} | d={d}")
            print("=" * 108)
            print(f"{'p':>5} {'lambda':>8} {'Delta':>8} {'qA':>6} {'qB':>6} {'qTot':>6} "
                  f"{'Aeff_ok':>8} {'Astalk':>8} {'B_ok':>8} {'valid':>8} "
                  f"{'time_A':>10} {'time_B':>10} {'time_tot':>10}")
            print("-" * 108)

            for p in primes:
                Fp = FiniteField(p)
                if sigma_name == "cube" and not Fp.cube_invertible():
                    print(f"{p:>5} {'-':>8} {'-':>8} {'-':>6} {'-':>6} {'-':>6} "
                          f"{'skip':>8} {'skip':>8} {'skip':>8} {'skip':>8} "
                          f"{'-':>10} {'-':>10} {'-':>10}")
                    continue

                r = run_attack(
                    d=d,
                    p=p,
                    sigma_name=sigma_name,
                    seed=12000 + 100 * d + p + (0 if sigma_name == "inverse" else 5000),
                )
                print(f"{p:>5} {r['lambda']:>8} {r['delta']:>8} {r['qA']:>6} {r['qB']:>6} {r['qTotal']:>6} "
                      f"{str(r['A_eff_match']):>8} {str(r['A_stalk_match']):>8} "
                      f"{str(r['B_match']):>8} {str(r['valid']):>8} "
                      f"{r['time_A_s']:>10.4f} {r['time_B_s']:>10.4f} {r['time_total_s']:>10.4f}")

    print("\n" + "=" * 108)
    print("INTERPRETATION")
    print("=" * 108)
    print("- qA = 2d queries recover A in the coefficient query basis (A_eff).")
    print("- If P (coeff->stalk map) is invertible, Astalk=True confirms exact stalk-space A.")
    print("- qB is the number of extra CPA queries needed to find an invertible Z and solve for B.")
    print("- No search over GL(d, F_p) is performed.")
    print("- For sigma=inverse and sigma=cube, successful recovery falsifies the claim that the final design")
    print("  resists CPA merely because sigma is applied after the secret map A.")


if __name__ == "__main__":
    main()
