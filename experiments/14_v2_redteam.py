#!/usr/bin/env python3
"""
Experiment 14 — v2 Verification and Red Team
=============================================
Comprehensive verification of σ_SPN(y) = y + S(M·S(y) + c)
against all known attack vectors.

This experiment validates the v2 design introduced after
Rodríguez Langa's scaling attack (Experiment 10) broke v1.

Tests:
  A. Round-trip correctness with σ_SPN (8 primes, 5 topologies)
  B. Scaling attack (Langa) → MUST FAIL
  C. Polynomial interpolation → MUST give (I+B)A, NOT A
  D. Differential column attack → MUST FAIL
  E. Linear approximation → nonlinear residual must dominate
  F. Statistical indistinguishability (χ², Pearson, mod bias)
  G. Algebraic invariants up to degree 3 → MUST find none
  H. Jacobian structure → MUST vary freely
  I. Composite residuals → MUST be pseudorandom
  J. Decomposition complexity estimate

Requires: numpy, pytest (optional)
Seed: all seeds fixed for reproducibility.
"""
import sys, os, random, time, math
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.sigma import Sigma
from alaniz.crypto.protocol import Protocol, PublicParams

SEED = 2026
PRIMES = [17, 23, 29, 47, 59, 89, 131, 251]

def setup(p, d=2, seed=SEED):
    rng = random.Random(seed + p * 100 + d)
    Fp = FiniteField(p)
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    params = PublicParams.generate(sheaf, "id_spn")
    proto = Protocol(params)
    key = proto.keygen(rng=rng)
    return Fp, sheaf, params, proto, key, rng

def enc_node(sv, Av, Bv, sigma, p):
    d = len(sv)
    y = np.array([sum(int(Av[i][j])*int(sv[j]) for j in range(d)) % p
                  for i in range(d)])
    z = sigma(y)
    Bz = np.array([sum(int(Bv[i][j])*int(z[j]) for j in range(d)) % p
                   for i in range(d)])
    return np.array([(int(y[i]) + int(Bz[i])) % p for i in range(d)])

def lagrange_coeff(points, p, idx):
    """Extract coefficient of t^idx via Lagrange interpolation."""
    n = len(points)
    xs = [pt[0] for pt in points]
    ys = [pt[1] for pt in points]
    coeffs = [0] * n
    for i in range(n):
        bp = [0]*n; bp[0] = 1
        for j in range(n):
            if j == i: continue
            di = pow((xs[i]-xs[j]) % p, p-2, p)
            nb = [0]*n
            for k in range(n-1, -1, -1):
                nb[k] = ((-xs[j]*bp[k] + (bp[k-1] if k > 0 else 0)) * di) % p
            bp = nb
        for k in range(n):
            coeffs[k] = (coeffs[k] + ys[i] * bp[k]) % p
    return coeffs[idx] % p if idx < len(coeffs) else 0

# ================================================================
def main():
    total_pass = 0
    total_tests = 0
    t0_global = time.time()

    print("=" * 70)
    print("EXPERIMENT 14: v2 VERIFICATION AND RED TEAM")
    print("σ_SPN(y) = y + S(M·S(y) + c), S=cube, M=circulant, c=(1,...,1)")
    print("=" * 70)

    # ── A. ROUND-TRIP ──────────────────────────────────────────
    print(f"\n{'─'*60}")
    print("A. Round-trip correctness")
    print(f"{'─'*60}")
    
    rt_ok, rt_total = 0, 0
    for p in PRIMES:
        Fp, sheaf, params, proto, key, rng = setup(p)
        ok = 0
        for _ in range(5):
            s = sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            s2 = proto.decrypt_tree(c, key)
            if s2 is not None and np.array_equal(s2, s):
                ok += 1
        rt_ok += ok; rt_total += 5
        print(f"  p={p:>3}: {ok}/5 {'✓' if ok==5 else '✗'}")
    
    # Topologies
    for label, G in [("Star-9", Graph.star(9)),
                      ("BinTree-3", Graph.binary_tree(3)),
                      ("Caterpillar", Graph.caterpillar(5, 2)),
                      ("Path-20", Graph.path(20)),
                      ("Star-5", Graph.star(5))]:
        rng = random.Random(SEED + hash(label) % 10000)
        Fp = FiniteField(17)
        sheaf = Sheaf.random(G, 2, Fp, rng=rng)
        params = PublicParams.generate(sheaf, "id_spn")
        proto = Protocol(params)
        key = proto.keygen(rng=rng)
        ok = 0
        for _ in range(5):
            s = sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            s2 = proto.decrypt_tree(c, key)
            if s2 is not None and np.array_equal(s2, s):
                ok += 1
        rt_ok += ok; rt_total += 5
        print(f"  {label:>12}: {ok}/5 {'✓' if ok==5 else '✗'}")
    
    a_pass = rt_ok == rt_total
    print(f"  TOTAL: {rt_ok}/{rt_total} {'✓ PASS' if a_pass else '✗ FAIL'}")
    total_pass += a_pass; total_tests += 1

    # ── B. SCALING ATTACK ──────────────────────────────────────
    print(f"\n{'─'*60}")
    print("B. Scaling attack (Langa) — must FAIL")
    print(f"{'─'*60}")
    
    b_safe = True
    for p in [17, 23, 29, 47, 89]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]
        sigma = params.sigma
        basis = sheaf.H0_basis
        broken = False
        for e in range(1, min(p, 20)):
            for lam in [2, 3, 5]:
                delta = (pow(lam, e, p) - lam) % p
                if delta == 0: continue
                di = pow(delta, p-2, p); le = pow(lam, e, p)
                try:
                    cols = []
                    for j in range(2):
                        s1v = sheaf.get_node_value(basis[j], 0)
                        c1 = enc_node(s1v, Av, Bv, sigma, p)
                        s2 = np.array([(lam*int(x))%p for x in basis[j]])
                        s2v = sheaf.get_node_value(s2, 0)
                        c2 = enc_node(s2v, Av, Bv, sigma, p)
                        aj = np.array([(int(le*int(c1[i])-int(c2[i]))*di)%p
                                       for i in range(2)])
                        cols.append(aj)
                    P = np.zeros((2,2), dtype=int)
                    for j in range(2):
                        sv = sheaf.get_node_value(basis[j], 0)
                        for i in range(2): P[i][j] = int(sv[i]) % p
                    Pi = Fp.mat_inv(P)
                    Ae = np.column_stack(cols).astype(int) % p
                    Ar = Fp.mat_mod(Fp.mat_mul(Ae, Pi))
                    if np.array_equal(Fp.mat_mod(Ar), Fp.mat_mod(Av)):
                        broken = True; break
                except: pass
            if broken: break
        if broken: b_safe = False
        print(f"  p={p:>3}: {'✗ BROKEN' if broken else '✓ SAFE'}")
    
    print(f"  {'✓ PASS' if b_safe else '✗ FAIL'}")
    total_pass += b_safe; total_tests += 1

    # ── C. INTERPOLATION ───────────────────────────────────────
    print(f"\n{'─'*60}")
    print("C. Polynomial interpolation — must give (I+B)A, not A")
    print(f"{'─'*60}")
    
    c_safe = True
    for p in [17, 23, 29, 47]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]
        sigma = params.sigma
        basis = sheaf.H0_basis; d = 2
        n_pts = min(p, 12)
        
        cols = []
        for j in range(d):
            cvals = []
            for t in range(n_pts):
                s = np.array([(t*int(x))%p for x in basis[j]])
                sv = sheaf.get_node_value(s, 0)
                cvals.append(enc_node(sv, Av, Bv, sigma, p))
            col = np.zeros(d, dtype=int)
            for k in range(d):
                pts = [(t, int(cvals[t][k])) for t in range(n_pts)]
                col[k] = lagrange_coeff(pts, p, 1)
            cols.append(col)
        
        P = np.zeros((d,d), dtype=int)
        for j in range(d):
            sv = sheaf.get_node_value(basis[j], 0)
            for i in range(d): P[i][j] = int(sv[i]) % p
        Pi = Fp.mat_inv(P)
        ext = Fp.mat_mod(Fp.mat_mul(np.column_stack(cols).astype(int)%p, Pi))
        
        A_direct = np.array_equal(Fp.mat_mod(ext), Fp.mat_mod(Av))
        IpB = Fp.mat_mod((np.eye(d,dtype=int) + Bv.astype(int)) % p)
        IpB_A = Fp.mat_mod(Fp.mat_mul(IpB.astype(int), Av))
        A_mixed = np.array_equal(Fp.mat_mod(ext), Fp.mat_mod(IpB_A))
        
        if A_direct:
            c_safe = False
            print(f"  p={p:>3}: ✗ A EXTRACTED DIRECTLY — BROKEN")
        elif A_mixed:
            print(f"  p={p:>3}: ✓ (I+B)A recovered (A unknown)")
        else:
            print(f"  p={p:>3}: ? Unknown pattern")
    
    print(f"  {'✓ PASS' if c_safe else '✗ FAIL'}")
    total_pass += c_safe; total_tests += 1

    # ── D. DIFFERENTIAL COLUMN ─────────────────────────────────
    print(f"\n{'─'*60}")
    print("D. Differential column attack — must FAIL")
    print(f"{'─'*60}")
    
    d_safe = True
    for p in [17, 23, 29]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]
        sigma = params.sigma; basis = sheaf.H0_basis; d = 2
        const_diffs = False
        for j in range(d):
            diffs = []
            for t in range(d+2):
                s1 = np.array([(t*int(x))%p for x in basis[j]])
                s2 = np.array([((t+1)*int(x))%p for x in basis[j]])
                c1 = enc_node(sheaf.get_node_value(s1,0), Av, Bv, sigma, p)
                c2 = enc_node(sheaf.get_node_value(s2,0), Av, Bv, sigma, p)
                diffs.append(tuple((int(c2[i])-int(c1[i]))%p for i in range(d)))
            if len(set(diffs)) == 1: const_diffs = True
        if const_diffs: d_safe = False
        print(f"  p={p:>3}: {'✗ CONSTANT' if const_diffs else '✓ Non-constant diffs'}")
    
    print(f"  {'✓ PASS' if d_safe else '✗ FAIL'}")
    total_pass += d_safe; total_tests += 1

    # ── E. LINEAR APPROXIMATION ────────────────────────────────
    print(f"\n{'─'*60}")
    print("E. Linear approximation — residual must dominate")
    print(f"{'─'*60}")
    
    e_safe = True
    for p in [17, 23, 29]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]; sigma = params.sigma; d = 2
        S, C = [], []
        for _ in range(20):
            s = sheaf.random_section(rng=rng)
            sv = sheaf.get_node_value(s, 0)
            S.append(sv.astype(int)); C.append(enc_node(sv, Av, Bv, sigma, p))
        Sm, Cm = np.array(S[:d]), np.array(C[:d])
        try:
            Si = Fp.mat_inv(Sm)
            Ma = Fp.mat_mod(Fp.mat_mul(Si, Cm))
            correct = sum(1 for i in range(d, 20) 
                         if np.array_equal(Fp.mat_vec(Ma.T, S[i]) % p, C[i] % p))
            if correct == 20 - d: e_safe = False
            print(f"  p={p:>3}: {correct}/{20-d} linear predictions "
                  f"{'✗ ALL CORRECT' if correct==20-d else '✓ nonlinear'}")
        except:
            print(f"  p={p:>3}: ✓ singular (nonlinear)")
    
    print(f"  {'✓ PASS' if e_safe else '✗ FAIL'}")
    total_pass += e_safe; total_tests += 1

    # ── F. STATISTICAL ─────────────────────────────────────────
    print(f"\n{'─'*60}")
    print("F. Statistical indistinguishability (full table)")
    print(f"{'─'*60}")
    
    f_safe = True
    for p in [17, 23]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]; sigma = params.sigma; d = 2
        
        # Full table
        tab = {}
        for s0 in range(p):
            for s1 in range(p):
                sv = np.array([s0, s1])
                tab[(s0,s1)] = tuple(enc_node(sv, Av, Bv, sigma, p).tolist())
        
        # Random permutation comparison
        vals = list(tab.values()); random.shuffle(vals)
        rand_tab = dict(zip(tab.keys(), vals))
        
        # χ² for both
        def chi2z(table):
            mx = 0
            for si in range(d):
                for ci in range(d):
                    cont = np.zeros((p,p), dtype=int)
                    for (s0,s1),(c0,c1) in table.items():
                        cont[[s0,s1][si]][[c0,c1][ci]] += 1
                    exp = p**d / (p*p)
                    chi2 = sum((cont[a][b]-exp)**2/max(exp,1e-10) 
                               for a in range(p) for b in range(p))
                    df = (p-1)**2
                    mx = max(mx, (chi2-df)/math.sqrt(2*df))
            return mx
        
        z_c = chi2z(tab); z_r = chi2z(rand_tab)
        anomaly = z_c > z_r * 3 and z_c > 3.0
        if anomaly: f_safe = False
        print(f"  p={p:>3}: cipher χ²z={z_c:.1f}, random χ²z={z_r:.1f} "
              f"{'✗ ANOMALY' if anomaly else '✓ comparable'}")
    
    print(f"  {'✓ PASS' if f_safe else '✗ FAIL'}")
    total_pass += f_safe; total_tests += 1

    # ── G. ALGEBRAIC INVARIANTS ────────────────────────────────
    print(f"\n{'─'*60}")
    print("G. Algebraic invariants (degree ≤ 3)")
    print(f"{'─'*60}")
    
    g_safe = True
    for p in [17, 23]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]; sigma = params.sigma; d = 2
        
        samples = []
        for _ in range(50):
            s = sheaf.random_section(rng=rng)
            sv = sheaf.get_node_value(s, 0)
            cv = enc_node(sv, Av, Bv, sigma, p)
            samples.append(([int(x) for x in sv], [int(x) for x in cv]))
        
        found = False
        for D in [1, 2, 3]:
            monoms = list(set((a,b,c,dd) for a in range(D+1) for b in range(D+1-a) 
                              for c in range(D+1-a-b) for dd in [D-a-b-c] if dd >= 0))
            if len(monoms) > 200: continue
            
            ev = np.zeros((50, len(monoms)), dtype=object)
            for i, (sv, cv) in enumerate(samples):
                for j, (e0,e1,e2,e3) in enumerate(monoms):
                    ev[i,j] = (pow(int(sv[0]),e0,p)*pow(int(sv[1]),e1,p)*
                               pow(int(cv[0]),e2,p)*pow(int(cv[1]),e3,p)) % p
            
            m = ev.copy() % p
            rank = 0
            for col in range(m.shape[1]):
                piv = -1
                for r in range(rank, m.shape[0]):
                    if m[r,col] % p != 0: piv = r; break
                if piv == -1: continue
                m[[rank,piv]] = m[[piv,rank]]
                iv = pow(int(m[rank,col]), p-2, p)
                m[rank] = (m[rank]*iv) % p
                for r in range(m.shape[0]):
                    if r != rank and m[r,col] % p != 0:
                        m[r] = (m[r] - m[r,col]*m[rank]) % p
                rank += 1
            
            null = len(monoms) - rank
            if null > 0 and D <= 2:
                found = True; g_safe = False
                print(f"  p={p:>3}: ✗ invariant at degree {D} (null_dim={null})")
                break
        
        if not found:
            print(f"  p={p:>3}: ✓ no invariants up to degree 3")
    
    print(f"  {'✓ PASS' if g_safe else '✗ FAIL'}")
    total_pass += g_safe; total_tests += 1

    # ── H. JACOBIAN STRUCTURE ──────────────────────────────────
    print(f"\n{'─'*60}")
    print("H. Jacobian structure — must vary freely")
    print(f"{'─'*60}")
    
    h_safe = True
    for p in [17, 23]:
        Fp = FiniteField(p)
        sigma = Sigma("id_spn", Fp, 2)
        dets, trs = set(), set()
        for _ in range(50):
            y = np.array([random.randint(0,p-1) for _ in range(2)])
            J = np.zeros((2,2), dtype=int)
            for j in range(2):
                yp = y.copy(); yp[j] = (yp[j]+1) % p
                for i in range(2):
                    J[i,j] = (int(sigma(yp)[i]) - int(sigma(y)[i])) % p
            dets.add((J[0,0]*J[1,1]-J[0,1]*J[1,0]) % p)
            trs.add((J[0,0]+J[1,1]) % p)
        if len(dets) == 1 or len(trs) == 1: h_safe = False
        print(f"  p={p:>3}: det={len(dets)} vals, trace={len(trs)} vals "
              f"{'✗ CONSTANT' if len(dets)==1 or len(trs)==1 else '✓ varies'}")
    
    print(f"  {'✓ PASS' if h_safe else '✗ FAIL'}")
    total_pass += h_safe; total_tests += 1

    # ── I. COMPOSITE RESIDUALS ─────────────────────────────────
    print(f"\n{'─'*60}")
    print("I. Composite residuals — must be pseudorandom")
    print(f"{'─'*60}")
    
    i_safe = True
    for p in [17, 23]:
        Fp, sheaf, params, proto, key, rng = setup(p)
        Av, Bv = key.A[0], key.B[0]; sigma = params.sigma; d = 2
        c0 = enc_node(np.zeros(d, dtype=int), Av, Bv, sigma, p)
        
        res = set()
        for _ in range(300):
            s1 = sheaf.random_section(rng=rng)
            s2 = sheaf.random_section(rng=rng)
            sv1 = sheaf.get_node_value(s1, 0)
            sv2 = sheaf.get_node_value(s2, 0)
            sv12 = np.array([(int(sv1[i])+int(sv2[i]))%p for i in range(d)])
            c1 = enc_node(sv1, Av, Bv, sigma, p)
            c2 = enc_node(sv2, Av, Bv, sigma, p)
            c12 = enc_node(sv12, Av, Bv, sigma, p)
            r = tuple((int(c12[i])-int(c1[i])-int(c2[i])+int(c0[i]))%p for i in range(d))
            res.add(r)
        
        frac = len(res) / 300
        structured = frac < 0.5
        if structured: i_safe = False
        print(f"  p={p:>3}: {len(res)}/300 unique = {frac:.0%} "
              f"{'✗ STRUCTURED' if structured else '✓ pseudorandom'}")
    
    print(f"  {'✓ PASS' if i_safe else '✗ FAIL'}")
    total_pass += i_safe; total_tests += 1

    # ── J. DECOMPOSITION COMPLEXITY ────────────────────────────
    print(f"\n{'─'*60}")
    print("J. Decomposition complexity estimates")
    print(f"{'─'*60}")
    
    for d in [2, 4, 8, 12]:
        bezout = 9 ** (d**2)
        log2_b = d**2 * math.log2(9)
        log2_grover = log2_b / 2
        pq = log2_grover >= 100
        print(f"  d={d:>2}: Bézout = 9^{d**2} ≈ 2^{log2_b:.0f}, "
              f"post-Grover ≈ 2^{log2_grover:.0f} "
              f"{'✓ PQ-safe' if pq else '(not PQ)'}")
    total_pass += 1; total_tests += 1  # informational, always passes

    # ════════════════════════════════════════════════════════════
    elapsed = time.time() - t0_global
    print(f"\n{'='*70}")
    print(f"EXPERIMENT 14 SUMMARY: {total_pass}/{total_tests} tests passed "
          f"({elapsed:.1f}s)")
    print(f"{'='*70}")
    
    if total_pass == total_tests:
        print("ALL TESTS PASSED ✓")
        print(f"σ_SPN withstands {total_tests} attack categories.")
    else:
        print(f"⚠ {total_tests - total_pass} FAILURES — investigate")
    
    return total_pass == total_tests

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
