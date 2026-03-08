# The Alaniz Cipher

**Encryption from Nonlinear Sheaf Morphisms over Graphs**

A cryptographic primitive based on the Nonlinear Sheaf Morphism Inversion Problem (NL-SMIP): a polynomial system with secret coefficients coupled by sheaf cohomology over graphs.

> **Paper**: L. Alaniz Pintos, *"The Alaniz Cipher: Encryption from Nonlinear Sheaf Morphisms over Graphs"*, 2026.  
> **Preprint (DOI)**: [10.5281/zenodo.18908260](https://doi.org/10.5281/zenodo.18908260)  
> **ResearchGate**: [Publication 401670518](https://www.researchgate.net/publication/401670518_The_Alaniz_Cipher_Encryption_from_Nonlinear_Sheaf_Morphisms_over_Graphs)

## Overview

Messages are encoded as global sections of a cellular sheaf — elements of H⁰(G, F₀) — and encrypted via:

$$c_v = A_v \cdot s_v + B_v \cdot \sigma(A_v \cdot s_v)$$

where (Aᵥ, Bᵥ) ∈ GL(d, 𝔽ₚ)² are secret matrices per node and σ is a public nonlinear map (multiplicative inverse or cube map). Decryption uses tree propagation with cohomological filtering in O(n) time.

### Key Properties

| Property | Status | Evidence |
|---|---|---|
| Round-trip correctness | ✓ | 176/176 configurations, 0 failures |
| Decryption uniqueness | ✓ | Pr[ambiguity] ≤ (L−1) · γⁿ⁻¹ (Theorem 7.1) |
| CPA resistance | ✓ | Degree-3+ polynomial system, linearization fails |
| Post-quantum resistance | ✓ | No reduction to HSP; Grover is only quantum threat |
| NP-hardness | ✓ | Formal reduction from MQ (Lemma 6.2) |

### Security (NL-SMIP)

The attacker faces a polynomial system whose **coefficients are secret** (unlike MQ/Rainbow where the system is public). NL-SMIP is proven at least NP-hard via reduction from MQ. Security analyzed against Gröbner basis, XL, MinRank, linearization, and all known quantum attacks.

| Parameter Set | dᵥ | log₂ p | n | σ | Classical | Quantum |
|---|---|---|---|---|---|---|
| Demo | 2 | 5 | 8 | cube | ~2³⁸ | ~2¹⁹ |
| Standard | 4 | 31 | 16 | cube | ~2⁹³ | ~2⁴⁷ |
| PQ-128 | 8 | 61 | 32 | inverse | ~2²⁰⁰ | ~2¹⁰⁰ |
| PQ-256 | 8 | 127 | 64 | inverse | ~2⁴⁰⁰ | ~2²⁰⁰ |

## Quick Start

```bash
git clone https://github.com/QuantuSync/alaniz-cipher.git
cd alaniz-cipher

# Run the demo
python -m alaniz.demo.demo_basic

# Run all tests (29 tests, ~5 seconds)
pip install pytest numpy
pytest tests/ -v

# Run reproducible experiments (all paper data)
for f in experiments/0*.py; do python "$f"; done
```

## Project Structure

```
alaniz-cipher/
├── alaniz/
│   ├── core/
│   │   ├── field.py        # 𝔽ₚ arithmetic, matrix operations
│   │   ├── graph.py        # Graph topologies (paths, trees, stars)
│   │   └── sheaf.py        # Cellular sheaves, H⁰, coboundary, tree propagation
│   ├── crypto/
│   │   ├── sigma.py        # Nonlinear maps σ (inverse, cube)
│   │   └── protocol.py     # KeyGen, Encrypt, Decrypt, Encode/Decode
│   └── demo/
│       └── demo_basic.py   # Interactive walkthrough
├── experiments/
│   ├── 01_scaling_verification.py   # 176/176 round-trips (Table 7)
│   ├── 02_cpa_attack_analysis.py    # 4 broken designs + final (Section 5)
│   ├── 03_uniqueness_analysis.py    # Preimage distribution, γ (Section 7)
│   ├── 04_sigma_selection.py        # 5 σ candidates (Table 8)
│   ├── 05_complexity_analysis.py    # XL attack costs (Table 4)
│   └── 06_postquantum_analysis.py   # Quantum attack surface (Section 9)
├── tests/
│   └── test_protocol.py    # Full test suite (29 tests)
├── pyproject.toml
└── README.md
```

## Test Suite

The 29 automated tests cover:

- **Round-trip correctness**: 10 random messages + zero section + basis sections
- **Multiple primes**: 𝔽₁₇, 𝔽₂₃, 𝔽₂₉, 𝔽₄₇
- **Both σ functions**: inverse (x^{p−2}) and cube (x³)
- **8 graph topologies**: paths, stars, binary trees, caterpillars, random trees
- **Decryption uniqueness**: brute-force confirms single valid section
- **CPA resistance**: linear system inconsistency verified
- **Structural properties**: H⁰ dimension, section validity, tree propagation

## Requirements

- Python ≥ 3.11
- NumPy ≥ 1.24
- pytest ≥ 7.0 (for tests)

## Citation

```bibtex
@misc{alaniz2026cipher,
  author       = {Alaniz Pintos, Lucas},
  title        = {The Alaniz Cipher: Encryption from Nonlinear
                  Sheaf Morphisms over Graphs},
  year         = {2026},
  doi          = {10.5281/zenodo.18908260},
  url          = {https://doi.org/10.5281/zenodo.18908260},
  publisher    = {Zenodo}
}
```

## Author

**Dr. Lucas Alaniz Pintos**  
Gerencia de Smart Products, INECO  
lucas.alaniz@ineco.com  
ORCID: [0009-0008-5179-2534](https://orcid.org/0009-0008-5179-2534)

## License

MIT
