# Extended Security Analysis of the Alaniz Cipher
## Chosen-Plaintext Attacks, Key Recovery, and Structural Properties of the Local Encryption Map

---

**Companion to**: L. Alaniz Pintos, *"The Alaniz Cipher: Encryption from Nonlinear Sheaf Morphisms over Graphs"*, 2026.
DOI: [10.5281/zenodo.18908260](https://doi.org/10.5281/zenodo.18908260)

---

## Abstract

This document reports seven reproducible security experiments (Experiments 07–13) that complement and stress-test the security claims of the Alaniz Cipher (SHEC). The experiments probe the local encryption map $g_B(y) = y + B\,\sigma(y)$ for collision properties, attempt key recovery under known- and chosen-plaintext models, and analyze the algebraic structure exploitable by an adversary.

The central finding is contained in **Experiment 10**: both supported nonlinear maps ($\sigma = x^3$ and $\sigma = x^{p-2}$) satisfy a componentwise scaling law $\sigma(\lambda y) = \lambda^e\,\sigma(y)$, which allows a chosen-plaintext attacker to recover the per-node key matrices $(A_v, B_v)$ using only $O(d)$ encryption queries — independently of the field size $p$ and the size of $\mathrm{GL}(d, \mathbb{F}_p)$. This result directly contradicts the CPA-resistance claim in the original paper and is not addressed by the NL-SMIP hardness reduction, which concerns a distinct adversarial model. The remaining experiments characterize the collision structure of $g_B$ (Experiments 07 and 11), measure the concrete complexity of exhaustive and differential key recovery under the known-plaintext model (Experiments 08, 09, 12, and 13), and provide the algebraic ground-truth necessary for evaluating any proposed fix.

---

## 1. Introduction

The Alaniz Cipher encrypts a message encoded as a global section $s \in H^0(G, \mathcal{F})$ of a cellular sheaf over a graph $G$. At each node $v$ the local ciphertext is

$$c_v = A_v \cdot s_v + B_v \cdot \sigma(A_v \cdot s_v), \qquad (A_v, B_v) \in \mathrm{GL}(d, \mathbb{F}_p)^2,$$

where $\sigma : \mathbb{F}_p^d \to \mathbb{F}_p^d$ is a public componentwise nonlinear map. The original paper argues that this construction resists chosen-plaintext attacks (CPA) because $\sigma$ defeats linearization: any attempt to model $c_v$ as a linear function of $s_v$ produces an inconsistent system as soon as enough ciphertext pairs are collected.

The experiments in this document examine whether non-linear attacks succeed. They are organized in order of increasing algebraic sophistication:

- **Experiments 07 and 11** characterize the map $g_B(y) = y + B\,\sigma(y)$, the core function the attacker must invert.
- **Experiment 08** establishes the exact threshold at which the key becomes trivially recoverable ($d = 1$) and confirms algebraically that linear attacks fail for $d \geq 2$.
- **Experiments 09 and 12** measure the empirical cost of key recovery under the known-plaintext model via exhaustive A-search.
- **Experiment 10** presents a direct CPA break requiring only $O(d)$ queries.
- **Experiment 13** proposes a differential column attack that recovers the columns of $A$ independently, reducing the search space from $|\mathrm{GL}(d, p)|$ to $d\,(p^d - 1)$ candidates.

All experiments use the production implementation in `alaniz/` without modification. Seeds are fixed for full reproducibility.

---

## 2. Protocol Notation

For reference, the relevant quantities are:

| Symbol | Meaning |
|--------|---------|
| $p$ | Field characteristic (prime) |
| $d$ | Fiber dimension |
| $G$ | Graph with $n$ nodes |
| $H^0(G, \mathcal{F})$ | Space of global sections (plaintext space) |
| $s_v \in \mathbb{F}_p^d$ | Local plaintext at node $v$ |
| $(A_v, B_v)$ | Secret key at node $v$; both in $\mathrm{GL}(d, \mathbb{F}_p)$ |
| $\sigma$ | Public componentwise map: `inverse` ($x_i \mapsto x_i^{p-2}$) or `cube` ($x_i \mapsto x_i^3$) |
| $y_v = A_v s_v$ | Secret intermediate vector at node $v$ |
| $g_{B}(y) = y + B\,\sigma(y)$ | Local encryption function (with $A = I$) |

Decryption solves $y + B_r\,\sigma(y) = c_r$ at the root, propagates via restriction maps, and uses cohomological filtering to select the valid global section.

---

## 3. Experiment 07 — Collision Table for the Local Encryption Map

### 3.1 Objective

Quantify the surjectivity and preimage multiplicity of $g_B$ over the full domain $\mathbb{F}_p^d$ for $d \in \{1, 2, 4\}$ and $p \in \{17, 23, 29, 47, 59, 89, 131, 251\}$. This provides concrete evidence of whether $g_B$ is injective (one-to-one), which determines decryption ambiguity.

### 3.2 Methodology

For $d \in \{1, 2\}$, all $p^d$ inputs are enumerated exactly. The preimage map $\{c \mapsto g_B^{-1}(c)\}$ is built and the following metrics computed per random $B \in \mathrm{GL}(d, \mathbb{F}_p)$ (averaged over 8 independent draws):

$$\text{image fraction} = \frac{|\mathrm{Im}(g_B)|}{p^d}, \qquad
\text{collision probability} = \sum_c \!\left(\frac{|g_B^{-1}(c)|}{p^d}\right)^{\!2}.$$

Note that the *image fraction* (the proportion of outputs that are reachable) and the *pairwise collision probability* $\Pr[g_B(y_1) = g_B(y_2)]$ are distinct metrics; the experiment reports both. For $d = 4$, where exhaustive enumeration is infeasible, a birthday probe with $3\lceil\sqrt{p^d}\rceil$ random samples guarantees a collision is found with high probability by the birthday bound. Every collision example is verified against the algebraic condition $\Delta y = B(\sigma(y_2) - \sigma(y_1))$, which holds with probability 1 by construction.

### 3.3 Findings

Across all tested primes and both $\sigma$ functions, the image fraction stabilizes near **0.53 for $\sigma = \text{inverse}$** and **0.66 for $\sigma = \text{cube}$** in dimension $d = 1$, and near **0.62–0.64** for $d = 2$. The map $g_B$ is therefore consistently non-surjective: roughly one third of ciphertexts at each node are unreachable under a fixed key. Collisions are confirmed at all parameter sizes, including $d = 4$, $p = 251$.

### 3.4 Security implications

Non-surjectivity implies two consequences. First, decryption at the root node admits multiple preimages (at most the maximum multiplicity, observed at 2–5 for $d = 2$), introducing ambiguity that the cohomological filter resolves. Second, an eavesdropper who observes enough ciphertexts can statistically distinguish reachable from unreachable values, potentially narrowing the effective plaintext space. The collision probability also provides a direct estimate of the parameter $\gamma$ from Theorem 7.1 of the original paper.

---

## 4. Experiment 08 — Key Recovery at $d = 1$ and Linear Attack Diagnostics

### 4.1 Objective

Establish the exact threshold below which the key is algebraically recoverable, and confirm that the linearization defense holds for $d \geq 2$.

### 4.2 Methodology — Exact Recovery for $d = 1$

For $d = 1$ the encryption reduces to a scalar equation. Given two plaintext–ciphertext pairs $(s_1, c_1)$ and $(s_2, c_2)$:

**Case $\sigma = \text{inverse}$:** Multiplying $c_i = a s_i + b(as_i)^{-1}$ through by $as_i$ and subtracting yields
$$a = \frac{c_1 s_1 - c_2 s_2}{s_1^2 - s_2^2}, \qquad b = c_1(as_1) - (as_1)^2.$$

**Case $\sigma = \text{cube}$:** The two equations form a $2\times 2$ Vandermonde-type system in $(a,\, ba^3)$, solved directly by Cramer's rule with determinant $s_1 s_2^3 - s_2 s_1^3$.

Both formulas are exact over $\mathbb{F}_p$ and require a single pair of equations. The experiment demonstrates exact recovery of $(a, b) = (A_v[0,0],\, B_v[0,0])$ from 2 captures using deterministic arithmetic.

### 4.3 Methodology — Linear Surrogate Failure for $d \geq 2$

For $d \in \{2, 4\}$, the experiment builds the linear system $c_v = M s_v$ (treating the key as a single $d\times d$ matrix $M$) and checks whether it becomes inconsistent as captures accumulate. As expected, the nonlinearity of $\sigma$ makes the system inconsistent after enough pairs, confirming that pure linearization fails.

### 4.4 Findings and Security Implications

Dimension $d = 1$ provides **no computational security**: the key is recoverable with 2 plaintext-ciphertext pairs using only modular arithmetic. The original paper does not claim security at $d = 1$, but this result establishes the baseline: the minimum secure dimension is $d \geq 2$. The linear attack failure at $d \geq 2$ is consistent with the paper's CPA resistance argument, but Experiment 10 shows that a non-linear CPA attack succeeds at all dimensions.

---

## 5. Experiment 09 — Exhaustive Known-Plaintext Key Recovery

### 5.1 Objective

Measure the concrete cost of recovering the per-node key $(A_v, B_v)$ under the known-plaintext model via exhaustive search over $A \in \mathrm{GL}(d, \mathbb{F}_p)$ with linear solve for $B$.

### 5.2 Methodology

For each candidate $A$, compute $y_i = A s_i$ and $z_i = \sigma(y_i)$ for $d$ fit captures, then solve the linear system

$$c_i - y_i = B\,z_i, \quad i = 1, \ldots, d$$

for the $d^2$ entries of $B$. Validate the recovered $(A, B)$ against all $k$ captures. The campaign covers $p \in \{17, 23, 29, 47, 59, 89, 131, 251\}$, $d = 2$, with $k \in \{4, 8, 12\}$ captures. For $|\mathrm{GL}(2, p)| \leq 10^5$ exhaustive enumeration is used; for larger primes 5,000 random samples are drawn.

### 5.3 Findings

For $p \leq 29$, $|\mathrm{GL}(2, p)| \leq 7\times 10^5$, and the true key is recovered exactly in all trials. For $p \geq 47$ with exhaustive enumeration disabled, the true key is not found in 5,000 random samples (probability $\approx 5000/|\mathrm{GL}(2,p)|$). The number of false positives drops sharply with $k$: at $k = 12$ captures, no false positive survives in any exhaustive run.

### 5.4 Security Implications

The exhaustive A-search is intractable for the proposed standard parameters ($d = 4$, $p \approx 2^{31}$), confirming that known-plaintext brute force is not a practical threat at those sizes. However, Experiment 10 renders this bound irrelevant by showing that the key can be recovered with $O(d)$ chosen-plaintext queries regardless of $p$.

---

## 6. Experiment 10 — Direct CPA Key Recovery via Scaling Homogeneity

*This is the main security result of this analysis.*

### 6.1 Objective

Demonstrate a chosen-plaintext attack that recovers $(A_v, B_v)$ in $O(d)$ encryption queries, without any exhaustive search over $\mathrm{GL}(d, \mathbb{F}_p)$.

### 6.2 Mathematical Basis

Both supported nonlinear maps satisfy a **componentwise scaling law**:

$$\sigma(\lambda\, y)_i = \lambda^e\, \sigma(y)_i \quad \forall\, y \in \mathbb{F}_p^d,\; \lambda \in \mathbb{F}_p^*,$$

with $e = 3$ for $\sigma = \text{cube}$ and $e = -1$ (i.e., $e \equiv p-2 \pmod{p-1}$) for $\sigma = \text{inverse}$. Because $B$ acts linearly on $\sigma(y)$, this extends to the full ciphertext:

$$c(t\,e_j) = t\,a_j + B\,\sigma(t\,a_j) = t\,a_j + t^e\,B\,\sigma(a_j),$$

where $a_j = A\,e_j$ is the $j$-th column of $A$ and $e_j$ is the $j$-th canonical basis vector of the H⁰ coefficient space.

Querying $s_1 = e_j$ and $s_2 = \lambda\,e_j$ and writing $y_j = a_j$, $w_j = B\,\sigma(a_j)$:

$$\begin{pmatrix} 1 & 1 \\ \lambda & \lambda^e \end{pmatrix} \begin{pmatrix} y_j \\ w_j \end{pmatrix} = \begin{pmatrix} c_1 \\ c_2 \end{pmatrix}.$$

The system has determinant $\Delta = \lambda^e - \lambda$. Choosing any $\lambda$ with $\Delta \neq 0$ (which exists for all $p > 3$) gives a closed-form solution:

$$y_j = \frac{\lambda^e c_1 - c_2}{\Delta}, \qquad w_j = \frac{-\lambda c_1 + c_2}{\Delta}.$$

Since $y_j$ is the $j$-th column of $A\,P$ (where $P$ is the H⁰ coefficient-to-stalk map at the root node), repeating for $j = 0, \ldots, d-1$ recovers $A_{\mathrm{eff}} := A\,P$ in $2d$ queries. Since $A_{\mathrm{eff}}$ is sufficient for decryption in the coefficient basis, and since $P$ is publicly computable from the sheaf, $A$ itself is recovered as $A_{\mathrm{eff}}\,P^{-1}$.

With $A_{\mathrm{eff}}$ known, $B$ is recovered by collecting $d$ additional queries $s^{(i)}$ with linearly independent $z^{(i)} = \sigma(A_{\mathrm{eff}}\,s^{(i)})$ and solving the linear system $W = B\,Z$, where $W_i = c^{(i)} - A_{\mathrm{eff}}\,s^{(i)}$ and $Z = [z^{(1)} \cdots z^{(d)}]$:

$$B = W\,Z^{-1}.$$

The total query count is $2d + d = 3d$ in the worst case, reduced to $2d + 1$ when the basis queries already provide an invertible $Z$.

### 6.3 Experimental Validation

The attack is executed against all primes $p \in \{17, 23, 29, 47, 59, 89, 131, 251\}$ and both $\sigma$ functions at dimensions $d \in \{2, 4\}$. In every case:

- $A_{\mathrm{eff}}$ is recovered exactly (`A_eff_match = True`).
- $A$ in stalk coordinates is recovered exactly (`A_stalk_match = True`).
- $B$ is recovered exactly (`B_match = True`).
- The recovered key passes 20 fresh validation queries (`valid = True`).

For $d = 4$ the total number of encryption queries used is $2d + d = 12$, independently of $p$.

### 6.4 Security Implications

This result has the following consequences for the security claims of the original paper:

1. **The CPA-resistance claim is falsified for all parameter sets.** The attack does not require linearization and is not blocked by the nonlinearity of $\sigma$. On the contrary, the scaling homogeneity of $\sigma$ is the algebraic property that enables the attack. The statement "CPA resistant: linearization fails because $B_v\,\sigma(A_v s_v)$ mixes key and message nonlinearly" remains true as stated but does not capture this attack vector.

2. **The NL-SMIP hardness reduction is not affected.** The reduction from MQ establishes that recovering the key from a *polynomial system with secret coefficients* is NP-hard. The scaling attack does not solve NL-SMIP; instead it avoids it entirely by exploiting a structural property of $\sigma$ that is not captured by the polynomial formulation. These are complementary results addressing different adversarial models.

3. **Post-quantum claims are unaffected in substance.** The scaling attack is purely algebraic and applies with the same $O(d)$ complexity classically and quantum-mechanically. It does not exploit a quantum speedup; it is simply faster than any exhaustive search.

4. **The attack extends to all nodes simultaneously.** Because each node's key $(A_v, B_v)$ is independent and the attack is node-local, querying all $n$ nodes requires $3dn$ encryption queries. For a 6-node path with $d = 4$ this is 72 queries, still negligible.

5. **The attack is repairable if $\sigma$ is replaced by a map that lacks componentwise scaling homogeneity.** Any $\sigma$ satisfying $\sigma(\lambda y) \neq \lambda^e \sigma(y)$ for all $e$ would defeat this specific attack. The analysis of Experiment 04 of the original paper evaluates five candidate $\sigma$ functions on bijectivity and differential uniformity; extending that evaluation to scaling homogeneity is the natural direction for a fix.

---

## 7. Experiment 11 — Statistical Collision Study of $g_B$

### 7.1 Objective

Obtain statistically robust estimates of the collision properties of $g_B$ across the full parameter space, extending the spot-check of Experiment 07 to large samples of $B$ matrices.

### 7.2 Methodology

For $d = 2$ and all primes, exact enumeration is performed over $p^d$ inputs and repeated for 12–24 independent draws of $B \in \mathrm{GL}(d, \mathbb{F}_p)$. For $d = 4$, 30,000 random inputs are sampled per $B$. Per configuration the following are reported:

- **Mean image fraction** $\pm$ standard deviation over the $B$ sample.
- **Mean preimage multiplicity** $p^d / |\mathrm{Im}(g_B)|$.
- **Worst-$B$ examples**: the three $B$ matrices with the lowest image fraction.
- **Aggregated multiplicity histogram**: fraction of outputs with $1, 2, 3, \ldots$ preimages.

### 7.3 Findings

The image fraction for $\sigma = \text{inverse}$ (approximately 0.53 at $d = 1$, rising to 0.63 at $d = 2$) shows very low variance ($\sigma_{\text{std}} < 0.01$ for $p \geq 59$), indicating that the non-surjectivity of $g_B$ is a stable property of $\sigma$, not an artifact of unlucky $B$ choices. The maximum preimage multiplicity grows slowly with $p$ but remains bounded (below 7 for $d = 2$, $p \leq 251$), consistent with Theorem 7.1. The worst-$B$ families show image fractions up to 15% lower than the mean, suggesting that certain $B$ matrices could be avoided by design to reduce decryption ambiguity.

### 7.4 Relation to Experiment 07

Experiment 11 is the rigorous statistical extension of the spot-check in Experiment 07. Together they establish: (i) collisions are universal for all tested parameters; (ii) their frequency is stable and characterizable; (iii) the $\gamma$ parameter in Theorem 7.1 can be estimated empirically from the mean image fraction.

---

## 8. Experiment 12 — Key Recovery Diagnostics and Minimum Capture Analysis

### 8.1 Objective

Measure the minimum number of known-plaintext pairs required for unique key recovery at $d = 2$ (exhaustive) and $d = 3$ (sampled), and track how quickly false positives are eliminated as captures accumulate.

### 8.2 Methodology

For each $A$ candidate (exhaustive for $d = 2$, $p = 17$; sampled with 60,000 draws for $d = 3$), $B$ is solved from the first $d$ captures and the resulting candidate is tested against an increasing prefix of up to 20 captures. The *first failure index* — the smallest $k$ at which the candidate fails — determines how many captures it survives. The minimum $k$ at which exactly one candidate (the true key) remains is recorded as `min_cap_unique`.

### 8.3 Findings

For $d = 2$, $p = 17$: the number of surviving candidates drops from thousands at $k = d = 2$ to single digits by $k = 6$, and to exactly 1 (the true key, no false positives) by $k \approx 10$–$15$, depending on $\sigma$. For $d = 3$ in sampled mode, when the true key falls within the 60,000-sample budget, the survival curve is qualitatively similar, with unique recovery achieved by $k \approx 12$–$18$ captures. These results quantify the empirical "information content" of each capture: roughly $d^2 \log_2 p$ bits of key entropy are resolved per capture, consistent with the polynomial system structure.

### 8.4 Security Implications

Under the known-plaintext model, the key is recoverable in practice for small parameters ($d = 2$, $p = 17$) with 10–15 captures. For the proposed standard parameters ($d = 4$, $p \approx 2^{31}$), the exhaustive A-search is infeasible (see Experiment 09), so this analysis is primarily relevant for toy configurations. The data provide concrete lower bounds on the number of plaintext-ciphertext pairs an attacker would need to observe.

---

## 9. Experiment 13 — Differential Column Attack

### 9.1 Objective

Test whether the columns of $A$ can be recovered *independently* via chosen-plaintext differential queries along coordinate axes, and measure how much this reduces the attack complexity compared to a full search over $\mathrm{GL}(d, \mathbb{F}_p)$.

### 9.2 Mathematical Basis

For plaintexts $s = t\,e_j$ (scalar multiples of the $j$-th basis vector), the local ciphertext is $c(t\,e_j) = t\,a_j + B\,\sigma(t\,a_j)$. Taking consecutive differences:

$$\Delta c_t = c\bigl((t+1)\,e_j\bigr) - c\bigl(t\,e_j\bigr) = a_j + B\bigl(\sigma\bigl((t+1)a_j\bigr) - \sigma(t\,a_j)\bigr).$$

For a guessed column $\hat{a}_j$, this is a linear equation in $B$. Collecting $d$ such differentials and denoting $\Delta\sigma_t = \sigma((t+1)\hat{a}_j) - \sigma(t\hat{a}_j)$:

$$\Delta c_t - \hat{a}_j = B\,\Delta\sigma_t.$$

If the $d \times d$ matrix $[\Delta\sigma_0 \cdots \Delta\sigma_{d-1}]$ is invertible, $B$ is determined uniquely. Each column of $A$ is then recoverable by searching the $p^d - 1$ nonzero candidates for $\hat{a}_j$, reducing the total search space from $|\mathrm{GL}(d, \mathbb{F}_p)| \approx p^{d^2}$ to $d \cdot (p^d - 1)$.

### 9.3 Findings

For $d = 2$, $p = 17$: $p^2 - 1 = 288$ column candidates per axis, versus $|\mathrm{GL}(2, 17)| = 69{,}360$ matrices. The true column is uniquely identified with 4–6 differentials in all tested configurations, and the recovered $B$ matches the true key in every case where the true column has no zero component (see limitation below).

For $d = 3$, $p = 17$: $p^3 - 1 = 4{,}912$ candidates versus $|\mathrm{GL}(3, 17)| \approx 1.4 \times 10^7$. Unique recovery is achieved with 4–7 differentials.

**Limitation.** If a column $a_j$ has a zero in component $i$, then $\sigma(t\,a_j)_i = 0$ for all $t$ (since $\sigma(0) = 0$), making the $i$-th row of $\Delta\sigma$ identically zero for all differentials. The corresponding entries of $B$ remain undetermined and the attack produces false negatives for those configurations.

### 9.4 Security Implications

The differential column attack reduces the per-column search space from $p^{d^2}$ to $p^d - 1$ candidates — an exponential improvement for $d \geq 2$. However, it still requires $O(d\,(p^d - 1))$ linear solves, which is substantially more expensive than the $O(d)$ query attack of Experiment 10. As a standalone attack it is more relevant under scenarios where only *some* columns of $A$ need to be recovered, or as a diagnostic tool for understanding the algebraic structure of the key space.

---

## 10. Discussion

### 10.1 The Scaling Homogeneity as a Structural Vulnerability

The root cause of the CPA break in Experiment 10 is the interplay between three design choices:

1. **$\sigma$ is applied componentwise** — enabling the scaling law $\sigma(\lambda y)_i = \lambda^e \sigma(y)_i$.
2. **$B$ acts linearly on $\sigma(y)$** — so that $B\,\sigma(\lambda y) = \lambda^e\,B\,\sigma(y)$.
3. **The plaintext space includes scalar multiples of any vector** — allowing the attacker to query $(s, \lambda s)$ and decouple $y = As$ from $w = B\sigma(As)$ via a $2\times 2$ linear system.

Any two of these three properties in isolation would not be sufficient for the attack. A fix must therefore break at least one of them. The most natural option — replacing $\sigma$ with a map that lacks componentwise scaling — would also need to preserve bijectivity, CPA-resistance against linearization (as verified by Experiment 02), and compatibility with the decryption algorithm.

### 10.2 Relationship to the Original Security Proof

The original paper's CPA resistance argument (Section 5) considers designs where $\sigma$ is applied *before* or *independently of* the secret map, showing that in those cases the key can be recovered linearly. The final design places $\sigma$ *after* the secret linear map $A_v$, which does defeat linearization. Experiment 10 shows that defeating linearization is necessary but not sufficient for CPA security: the scaling homogeneity provides an alternative non-linear attack that does not require the attacker to linearize the system.

The NL-SMIP hardness proof (Lemma 6.2) establishes that recovering the key from a polynomial system with secret coefficients is NP-hard. This reduction models a different scenario — one in which the attacker has access to a fixed polynomial system and tries to find its solution. In the CPA model, the attacker controls the plaintexts and can structure their queries to exploit algebraic structure in $\sigma$. These two models are not equivalent, and the NP-hardness of NL-SMIP does not imply CPA security.

### 10.3 Collision Properties and Decryption Correctness

Experiments 07 and 11 together establish that $g_B$ is non-surjective with image fraction $\approx 0.63$ for $d = 2$ across all tested primes. The preimage multiplicity distribution is stable and bounded (maximum $\approx 5$–$7$ for $d = 2$). These facts are consistent with the correctness bound in Theorem 7.1 ($\Pr[\text{ambiguity}] \leq (L-1)\,\gamma^{n-1}$) and provide empirical estimates of $\gamma \approx 1 - 0.63 = 0.37$ and $L \approx 1/0.63 \approx 1.59$ for this setting. The bound decays exponentially in $n$, so for graphs with $n \geq 6$ nodes, decryption ambiguity is negligible in practice.

---

## 11. Conclusion

This analysis presents seven reproducible security experiments that probe the Alaniz Cipher from multiple adversarial perspectives. The key conclusions are:

1. **Dimension $d = 1$ is trivially insecure**: the key is algebraically recoverable from 2 plaintext-ciphertext pairs (Experiment 08).

2. **The local encryption map $g_B$ is consistently non-surjective** (image fraction $\approx 0.63$) but bounded in preimage multiplicity, validating the decryption uniqueness analysis (Experiments 07 and 11).

3. **Known-plaintext exhaustive key recovery** is feasible for toy parameters ($d = 2$, $p \leq 29$) and requires approximately 10–15 captures for unique recovery; it is infeasible at proposed standard parameters (Experiments 09 and 12).

4. **A direct CPA break exists for all parameters**: the scaling homogeneity $\sigma(\lambda y) = \lambda^e\,\sigma(y)$ allows recovery of the per-node key $(A_v, B_v)$ in $O(d)$ chosen-plaintext queries, independently of $p$ (Experiment 10). This is the principal security finding of this analysis.

5. **A differential column attack** reduces the per-column search space from $|\mathrm{GL}(d, p)|$ to $d\,(p^d - 1)$ candidates, though it is dominated by the scaling attack for practical purposes (Experiment 13).

The scaling attack of Experiment 10 is a significant vulnerability that requires attention before the scheme can be considered secure under standard CPA definitions. The natural direction for a repair is to replace the componentwise scaling-homogeneous $\sigma$ with a map that mixes components across the $d$-dimensional fiber while preserving the algebraic properties required for efficient decryption.

---

## Appendix: Experiment Summary Table

| Exp | Model | Attack / Analysis | Complexity | Main Result |
|-----|-------|------------------|------------|-------------|
| 07 | — | $g_B$ collision characterization | $O(p^d \cdot b)$ exact | Image fraction $\approx 0.63$, collisions universal |
| 08 | KPA | Algebraic key recovery | $O(1)$ for $d=1$ | $d=1$ trivially broken; linearization fails at $d\geq 2$ |
| 09 | KPA | Exhaustive A-search + linear B | $O(|\mathrm{GL}(d,p)|)$ | Feasible for $p\leq 29$; intractable at standard params |
| 10 | **CPA** | **Scaling homogeneity attack** | **$O(d)$ queries** | **Full key recovery at all params; CPA claim falsified** |
| 11 | — | $g_B$ statistical multiplicity | $O(p^d \cdot b)$ exact | Stable image fraction; worst-$B$ families identified |
| 12 | KPA | Capture count diagnostics | $O(|\mathrm{GL}(d,p)| \cdot k)$ | Unique recovery in $\approx 10$–$15$ captures for $d=2$ |
| 13 | CPA | Differential column attack | $O(d \cdot (p^d-1))$ | Search space reduced from $p^{d^2}$ to $p^d - 1$ per column |

KPA = Known-Plaintext Attack, CPA = Chosen-Plaintext Attack.

---

*Experiments 07–13 are implemented in `experiments/07_collision_table.py` through `experiments/13_differential_column_attack.py`. All seeds are fixed; results are fully reproducible.*
