# Experiments

Reproducible scripts for every claim in the paper.

| Script | Paper Section | What it reproduces |
|---|---|---|
| `01_scaling_verification.py` | Table 5 | 106/106 round-trip correctness across 22 configurations |
| `02_cpa_attack_analysis.py` | Section 5 | Four broken designs + final design CPA resistance |
| `03_uniqueness_analysis.py` | Section 7 | Preimage distribution, collision probability γ, tree propagation |
| `04_sigma_selection.py` | Table 6 | Bijectivity, differential uniformity, round-trip for 5 σ candidates |
| `05_complexity_analysis.py` | Table 4 | XL attack degree, matrix size, cost per parameter set |
| `06_postquantum_analysis.py` | Section 9 | Shor inapplicability, Grover quantification, scheme comparison |

## Running

```bash
# Run all experiments
python experiments/01_scaling_verification.py
python experiments/02_cpa_attack_analysis.py
python experiments/03_uniqueness_analysis.py
python experiments/04_sigma_selection.py
python experiments/05_complexity_analysis.py
python experiments/06_postquantum_analysis.py

# Or run them all at once
for f in experiments/0*.py; do echo "=== $f ===" && python "$f" && echo; done
```

Experiments 1, 3, and 4 involve brute-force search over F_p^d and may take 1-5 minutes each. Experiments 2, 5, and 6 are fast (<5 seconds).
