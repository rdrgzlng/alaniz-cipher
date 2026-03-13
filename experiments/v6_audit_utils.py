#!/usr/bin/env python3
"""
Shared helpers for V6 audit experiments.
"""

from __future__ import annotations

import os
import random
import sys
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


def solve_linear_Fp(equations, n_unknowns, p):
    """Solve a linear system over F_p. Return None if inconsistent."""
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
        inv = pow(int(aug[row][col]), p - 2, p)
        aug[row] = (aug[row] * inv) % p
        for r in range(n_eq):
            if r != row and aug[r][col] % p != 0:
                aug[r] = (aug[r] - aug[r][col] * aug[row]) % p
        row += 1

    for r in range(len(pivots), n_eq):
        if aug[r, n_unknowns] % p != 0:
            return None

    sol = np.zeros(n_unknowns, dtype=int)
    for i, pc in enumerate(pivots):
        sol[pc] = int(aug[i, n_unknowns]) % p
    return sol


def interpolate_univariate(values, p):
    """
    Recover degree < p coefficients from f(t), t = 0..p-1.
    """
    equations = []
    for t, val in enumerate(values):
        coeffs = [pow(t, k, p) for k in range(p)]
        equations.append((coeffs, int(val) % p))
    sol = solve_linear_Fp(equations, p, p)
    if sol is None:
        raise RuntimeError("Interpolation system is inconsistent")
    return sol % p


def build_v6_instance(d, p, n_nodes=6, seed=0):
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(n_nodes)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    proto = Protocol(PublicParams.generate(sheaf, "id_spn"))
    key = proto.keygen(rng=rng)
    return proto, key, rng


def make_coeff_plaintext(proto: Protocol, coeffs) -> np.ndarray:
    return proto.sheaf.section_from_coeffs(np.array(coeffs, dtype=int) % proto.Fp.p)


def encrypt_coeffs_at_node(proto: Protocol, key, coeffs, node=0) -> np.ndarray:
    s = make_coeff_plaintext(proto, coeffs)
    c = proto.encrypt(s, key)
    return proto.sheaf.get_node_value(c, node) % proto.Fp.p


def coeff_to_node_map(sheaf, node: int) -> np.ndarray:
    """
    Matrix P_v such that s_v = P_v * coeffs for coeff coordinates in the H^0 basis.
    """
    p = sheaf.Fp.p
    cols = []
    for basis_vec in sheaf.H0_basis:
        cols.append(sheaf.get_node_value(basis_vec.astype(int), node) % p)
    return np.column_stack(cols) % p


def line_values_for_axis(proto: Protocol, key, axis: int, node=0):
    """
    Return c_v(t e_axis) for all t in F_p.
    """
    p = proto.Fp.p
    d = proto.sheaf.H0_dim
    vals = []
    for t in range(p):
        coeffs = np.zeros(d, dtype=int)
        coeffs[axis] = t
        vals.append(encrypt_coeffs_at_node(proto, key, coeffs, node=node))
    return vals


def recover_linear_coefficient_matrix(proto: Protocol, key, node=0):
    """
    Recover the t^1 coefficient column-by-column from axis line interpolation.
    """
    p = proto.Fp.p
    d = proto.sheaf.H0_dim
    cols = []
    per_axis = []
    for axis in range(d):
        vals = line_values_for_axis(proto, key, axis, node=node)
        coeff_by_coord = [interpolate_univariate([int(v[i]) for v in vals], p) for i in range(proto.dv)]
        col = np.array([int(coeffs[1]) % p for coeffs in coeff_by_coord], dtype=int)
        cols.append(col)
        per_axis.append({
            "axis": axis,
            "constant": np.array([int(coeffs[0]) % p for coeffs in coeff_by_coord], dtype=int),
            "linear": col.copy(),
        })
    return np.column_stack(cols) % p, per_axis


def gl_size_d2(p):
    return (p * p - 1) * (p * p - p)


def all_gl_matrices_d2(Fp):
    p = Fp.p
    for a, b, c, d in product(range(p), repeat=4):
        if (a * d - b * c) % p == 0:
            continue
        yield np.array([[a, b], [c, d]], dtype=int)


def count_effective_matrix_decompositions(Fp, M_eff):
    """
    Count (A,B) such that (I + B) A = M_eff with A,B in GL(2,p).
    Only for d=2, exhaustive.
    """
    p = Fp.p
    count = 0
    examples = []
    I = np.eye(2, dtype=int)
    for A in all_gl_matrices_d2(Fp):
        A_inv = Fp.mat_inv(A)
        B = (M_eff @ A_inv - I) % p
        det_B = int((B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]) % p)
        det_I_plus_B = int(((I + B)[0, 0] * (I + B)[1, 1] - (I + B)[0, 1] * (I + B)[1, 0]) % p)
        if det_B == 0 or det_I_plus_B == 0:
            continue
        count += 1
        if len(examples) < 3:
            examples.append((A.copy(), B.copy()))
    return count, examples


def collect_local_captures(proto: Protocol, key, n_caps, node=0, seed=0):
    rng = random.Random(seed)
    caps = []
    for _ in range(n_caps):
        s = proto.sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)
        s_v = proto.sheaf.get_node_value(s, node) % proto.Fp.p
        c_v = proto.sheaf.get_node_value(c, node) % proto.Fp.p
        caps.append((s_v, c_v))
    return caps
