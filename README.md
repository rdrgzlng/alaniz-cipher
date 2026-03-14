# The Alaniz Cipher

**Encryption from Nonlinear Sheaf Morphisms over Graphs**

A cryptographic primitive based on the Nonlinear Sheaf Morphism Inversion Problem (NL-SMIP): a polynomial system with secret coefficients coupled by sheaf cohomology over graphs.

> **Paper**: L. Alaniz Pintos, *"The Alaniz Cipher: Encryption from Nonlinear Sheaf Morphisms over Graphs"*, 2026.  
> **Preprint (DOI)**: [10.5281/zenodo.18998497](https://doi.org/10.5281/zenodo.18998497)  
> **ResearchGate**: [Publication 401670518](https://www.researchgate.net/publication/401670518_The_Alaniz_Cipher_Encryption_from_Nonlinear_Sheaf_Morphisms_over_Graphs)

## Overview

Messages are encoded as global sections of a cellular sheaf — elements of H⁰(G, F₀) — and encrypted via:

$$c_v = A_v \cdot s_v + B_v \cdot \sigma(A_v \cdot s_v)$$

where (Aᵥ, Bᵥ) ∈ GL(d, 𝔽ₚ)² are secret matrices per node and σ is a public nonlinear map. The current implementation defaults to `id_spn`, a V6 sigma with identity linear part, cross-component mixing, and an internal constant term. Legacy componentwise maps `inverse` and `cube` are retained only for comparison and historical experiments.

### Key Properties

| Property | Status | Evidence |
|---|---|---|
| Round-trip correctness | ✓ | Code path validated for the V6 default `id_spn` |
| Decryption uniqueness | ✓ | Root solve + tree propagation + global verification |
| CPA audit status | In progress | Legacy CPA breaks target `inverse` and `cube`, not `id_spn` |
| Legacy sigma support | Historical | `inverse` and `cube` kept for comparison and archived analysis |

### Security (NL-SMIP)

The attacker faces a polynomial system whose **coefficients are secret** (unlike MQ/Rainbow where the system is public). The repository currently contains both the legacy componentwise sigma family and the V6 `id_spn` repair candidate. Security claims and experiments must therefore be read in that context: several archived analyses target the legacy family and should not be interpreted as a break of V6.

| Parameter Set | dᵥ | log₂ p | n | σ | Classical | Quantum |
|---|---|---|---|---|---|---|
| Demo | 2 | 5 | 8 | id_spn | pending | pending |
| Standard | 4 | 31 | 16 | id_spn | pending | pending |
| Legacy comparison | 2 | 5 | 8 | cube | archived | archived |

## Quick Start

```bash
git clone https://github.com/QuantuSync/alaniz-cipher.git
cd alaniz-cipher

# Run the demo
python -m alaniz.demo.demo_basic

# Run the test suite
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
│   │   ├── sigma.py        # Nonlinear maps σ (default: id_spn; legacy: inverse, cube)
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
│   └── test_protocol.py    # Full test suite
├── pyproject.toml
└── README.md
```

## Test Suite

The automated tests cover:

- **Round-trip correctness**: 10 random messages + zero section + basis sections
- **Multiple primes**: 𝔽₁₇, 𝔽₂₃, 𝔽₂₉, 𝔽₄₇
- **V6 default sigma**: `id_spn`
- **Legacy sigma compatibility**: inverse (x^{p−2}) and cube (x³)
- **8 graph topologies**: paths, stars, binary trees, caterpillars, random trees
- **Decryption uniqueness**: brute-force confirms single valid section
- **Basic CPA sanity check**: linear model inconsistency on the V6 default path
- **Structural properties**: H⁰ dimension, section validity, tree propagation

## Versioning Notes

- The implementation in `alaniz/crypto` is V6-oriented and defaults to `id_spn`.
- Several experiments and analysis files still study the legacy `inverse` and `cube` designs because those attacks motivated the V6 change.
- When an experiment or document discusses scaling homogeneity or componentwise sigma, it refers to the legacy family unless it explicitly names `id_spn`.

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
