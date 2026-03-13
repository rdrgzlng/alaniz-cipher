# Experiments

Reproducible scripts for the paper and the subsequent security analysis.

## Status

The repository now defaults to the V6 sigma `id_spn` in `alaniz/crypto`.
Many experiments below were written for the legacy componentwise sigma family
(`inverse`, `cube`) and are retained as archival attack evidence that motivated
the V6 redesign. Unless a script explicitly names `id_spn`, do not read its
results as claims about the current default implementation.

| Script | Paper Section | What it reproduces |
|---|---|---|
| `01_scaling_verification.py` | Table 5 | 106/106 round-trip correctness across 22 configurations |
| `02_cpa_attack_analysis.py` | Section 5 | Legacy design study: four broken designs + pre-V6 final design |
| `03_uniqueness_analysis.py` | Section 7 | Preimage distribution, collision probability γ, tree propagation |
| `04_sigma_selection.py` | Table 6 | Bijectivity, differential uniformity, round-trip for 5 σ candidates |
| `05_complexity_analysis.py` | Table 4 | XL attack degree, matrix size, cost per parameter set |
| `06_postquantum_analysis.py` | Section 9 | Shor inapplicability, Grover quantification, scheme comparison |
| `07_collision_table.py` | Security stress test | Legacy sigma collision table for g(y)=y+Bσ(y) across (p, d, σ) |
| `08_d1_key_recovery_poc.py` | Security stress test | Legacy key-recovery analysis: exact PoC for d=1 + linear recovery attempts for d=2 and d=4 |
| `09_exhaustive_A_known_plaintext_attack.py` | Security stress test | Legacy node-wise key recovery campaign for d=2 across primes/captures via A-search |
| `10_cpa_scaling_closed_form_attack.py` | Security stress test | Legacy direct CPA key recovery via componentwise σ homogeneity |
| `11_collision_multiplicity_study.py` | Security stress test | Legacy deep study of g_B(y): image fraction, multiplicity histograms, and worst-B families |
| `12_key_recovery_diagnostics.py` | Security stress test | Legacy cap=2..20 diagnostics for exhaustive/sampled A-search |
| `13_differential_column_attack.py` | Security stress test | Legacy chosen-plaintext differential attack on column recovery |
| `14_v6_interpolation_matrix_audit.py` | V6 audit | Interpolation audit: verifies recovery of `(I+B)A` rather than `A` and measures decomposition ambiguity |
| `15_v6_spn_differential_audit.py` | V6 audit | Exact differential profile of `id_spn` on small parameters |
| `16_v6_effective_matrix_hybrid.py` | V6 audit | Hybrid attack audit: candidate reduction and extra captures once `(I+B)A` is known |
| `17_v6_degeneracy_audit.py` | V6 audit | Degeneracy audit: fixed points, invariant lines, singular discrete Jacobians, component fixation |
| `18_v6_node_decoupling_audit.py` | V6 audit | Node-decoupling audit: recovers node-local `(I+B_v)A_v` from shared global CPA queries |

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
python experiments/10_cpa_scaling_closed_form_attack.py
python experiments/11_collision_multiplicity_study.py
python experiments/12_key_recovery_diagnostics.py
python experiments/13_differential_column_attack.py
python experiments/14_v6_interpolation_matrix_audit.py
python experiments/15_v6_spn_differential_audit.py
python experiments/16_v6_effective_matrix_hybrid.py
python experiments/17_v6_degeneracy_audit.py
python experiments/18_v6_node_decoupling_audit.py

# Or run them all at once
for f in experiments/0*.py; do echo "=== $f ===" && python "$f" && echo; done
```

Experiments 1, 3, and 4 involve brute-force search over F_p^d and may take 1-5 minutes each. Experiments 2, 5, and 6 are fast (<5 seconds). Experiment 10 is fast (<5 seconds, O(d) queries). Experiments 11 and 13 may take several minutes for larger parameters.
