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
| `07_collision_table.py` | Security stress test | Collision table for g(y)=y+Bσ(y) across (p, d, σ) |
| `08_d1_key_recovery_poc.py` | Security stress test | Key-recovery analysis: exact PoC for d=1 + linear recovery attempts for d=2 and d=4 |
| `09_exhaustive_A_known_plaintext_attack.py` | Security stress test | Node-wise key recovery campaign for d=2 across primes/captures via A-search (exhaustive+sampled) + linear solve for B |
| `11_collision_multiplicity_study.py` | Security stress test | Deep study of g_B(y): image fraction, multiplicity histograms, and worst-B families across (p,d,σ) |
| `12_key_recovery_diagnostics.py` | Security stress test | cap=2..20 diagnostics: minimum captures, A-survival rates, false positives, and d=3 sampled extension |

## Running

```bash
# Run all experiments
python experiments/01_scaling_verification.py
python experiments/02_cpa_attack_analysis.py
python experiments/03_uniqueness_analysis.py
python experiments/04_sigma_selection.py
python experiments/05_complexity_analysis.py
python experiments/06_postquantum_analysis.py
python experiments/07_collision_table.py
python experiments/08_d1_key_recovery_poc.py
python experiments/09_exhaustive_A_known_plaintext_attack.py
python experiments/11_collision_multiplicity_study.py
python experiments/12_key_recovery_diagnostics.py

# Or run them all at once
for f in experiments/0*.py; do echo "=== $f ===" && python "$f" && echo; done
```

Experiments 1, 3, and 4 involve brute-force search over F_p^d and may take 1-5 minutes each. Experiments 2, 5, and 6 are fast (<5 seconds).
