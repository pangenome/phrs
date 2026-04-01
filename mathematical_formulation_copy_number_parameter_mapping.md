# Mathematical Formulation of Copy-Number Parameter Mapping for Hypergeometric Testing

## Abstract

This document provides a rigorous mathematical foundation for transforming hypergeometric test parameters from gene counts to copy-number weighted gene instance counts. We establish formal notation, prove mathematical equivalence with instance expansion methods, and analyze parameter validity constraints.

## 1. Mathematical Framework and Notation

### 1.1 Basic Definitions

Let **G** be the set of all genes in the genome, and let **C**: G → ℕ be a copy-number function that maps each gene to its copy number.

**Definition 1.1 (Gene Instance Space):** For a gene set S ⊆ G, the instance space I(S) is defined as:
```
I(S) = {(g, i) : g ∈ S, i ∈ {1, 2, ..., C(g)}}
```
where each element (g, i) represents the i-th instance of gene g.

**Definition 1.2 (Instance Count):** The instance count of a gene set S is:
```
|I(S)| = Σ_{g∈S} C(g)
```

### 1.2 Hypergeometric Test Framework

**Definition 1.3 (Standard Hypergeometric Model):** 
Let:
- B ⊆ G be the background gene set (universe)
- P ⊆ B be the pathway gene set (items of interest)
- Q ⊆ B be the query gene set (sample)
- O = Q ∩ P be the observed overlap

The standard hypergeometric test uses parameters:
```
k₀ = |Q|           (sample size)
m₀ = |P|           (items of interest in population)
n₀ = |B| - |P|     (items not of interest in population)
q₀ = |O|           (observed overlap)
```

**Definition 1.4 (Copy-Number Weighted Model):**
The weighted hypergeometric test transforms to instance space with parameters:
```
k_w = |I(Q)| = Σ_{g∈Q} C(g)                    (weighted sample size)
m_w = |I(P∩B)| = Σ_{g∈P∩B} C(g)               (weighted items of interest)
n_w = |I(B)| - |I(P∩B)|                        (weighted items not of interest)
q_w = |I(Q∩P)| = Σ_{g∈Q∩P} C(g)               (weighted observed overlap)
```

## 2. Parameter Transformation Theorems

### 2.1 Fundamental Transformation Relations

**Theorem 2.1 (Parameter Transformation):** 
For any gene sets Q, P, B with copy-number function C, the weighted parameters satisfy:

1. **Query Size Transformation:**
   ```
   k_w = Σ_{g∈Q} C(g)
   ```

2. **Pathway Size Transformation:**
   ```
   m_w = Σ_{g∈P∩B} C(g)
   ```

3. **Non-pathway Size Transformation:**
   ```
   n_w = Σ_{g∈B\P} C(g) = Σ_{g∈B} C(g) - Σ_{g∈P∩B} C(g)
   ```

4. **Overlap Transformation:**
   ```
   q_w = Σ_{g∈Q∩P} C(g)
   ```

**Proof:** Direct from Definition 1.4 and set theory. □

### 2.2 Conservation Properties

**Theorem 2.2 (Total Instance Conservation):** 
The total instance count is preserved across different partitions:
```
k_w + (|I(B)| - k_w) = |I(B)|
m_w + n_w = |I(B)|
```

**Proof:** 
```
m_w + n_w = Σ_{g∈P∩B} C(g) + Σ_{g∈B\P} C(g)
          = Σ_{g∈B} C(g)    [since (P∩B) ∪ (B\P) = B and disjoint]
          = |I(B)|
```
□

## 3. Mathematical Equivalence with Instance Expansion

### 3.1 Instance Expansion Method

**Definition 3.1 (Instance Expansion):**
For gene sets with copy numbers, define expanded sets:
```
Q_exp = {(g,i) : g ∈ Q, i ∈ {1,...,C(g)}}
P_exp = {(g,i) : g ∈ P, i ∈ {1,...,C(g)}}
B_exp = {(g,i) : g ∈ B, i ∈ {1,...,C(g)}}
```

The expansion method parameters are:
```
k_exp = |Q_exp|
m_exp = |P_exp ∩ B_exp|
n_exp = |B_exp| - |P_exp ∩ B_exp|
q_exp = |Q_exp ∩ P_exp|
```

### 3.2 Equivalence Theorem

**Theorem 3.1 (Method Equivalence):**
For any gene sets Q, P, B with copy-number function C:
```
k_w = k_exp,  m_w = m_exp,  n_w = n_exp,  q_w = q_exp
```

**Proof:**

*Step 1:* Show k_w = k_exp
```
k_exp = |Q_exp| = |{(g,i) : g ∈ Q, i ∈ {1,...,C(g)}}|
      = Σ_{g∈Q} |{i : i ∈ {1,...,C(g)}}|
      = Σ_{g∈Q} C(g)
      = k_w
```

*Step 2:* Show m_w = m_exp
```
m_exp = |P_exp ∩ B_exp|
      = |{(g,i) : g ∈ P ∩ B, i ∈ {1,...,C(g)}}|
      = Σ_{g∈P∩B} C(g)
      = m_w
```

*Step 3:* Show n_w = n_exp
```
n_exp = |B_exp| - |P_exp ∩ B_exp|
      = Σ_{g∈B} C(g) - Σ_{g∈P∩B} C(g)
      = Σ_{g∈B\P} C(g)
      = n_w
```

*Step 4:* Show q_w = q_exp
```
q_exp = |Q_exp ∩ P_exp|
      = |{(g,i) : g ∈ Q ∩ P, i ∈ {1,...,C(g)}}|
      = Σ_{g∈Q∩P} C(g)
      = q_w
```
□

**Corollary 3.1 (Statistical Equivalence):**
The hypergeometric probability computed using weighted parameters equals that computed using instance expansion:
```
P(X ≥ q_w + 1) = phyper(q_w, m_w, n_w, k_w, lower.tail=FALSE)
                = phyper(q_exp, m_exp, n_exp, k_exp, lower.tail=FALSE)
```

## 4. Parameter Constraint Analysis

### 4.1 Basic Feasibility Constraints

**Theorem 4.1 (Parameter Validity):**
For valid hypergeometric parameters, the following constraints must hold:

1. **Non-negativity:** q_w, m_w, n_w, k_w ∈ ℕ ∪ {0}

2. **Logical bounds:**
   ```
   q_w ≤ k_w           (overlap cannot exceed sample size)
   q_w ≤ m_w           (overlap cannot exceed pathway size)
   k_w ≤ m_w + n_w     (sample cannot exceed population)
   ```

3. **Feasibility range:**
   ```
   max(0, k_w - n_w) ≤ q_w ≤ min(k_w, m_w)
   ```

**Proof:** These follow directly from the hypergeometric model constraints and the interpretation of parameters as counts. □

### 4.2 Copy Number Model Constraints

**Definition 4.1 (Copy Number Model Consistency):**
A copy number function C is consistent with background B if:
```
∀g ∈ G: C(g) > 0 ⟺ g ∈ B
```

**Theorem 4.2 (Model Consistency Requirement):**
For valid statistical inference, the copy number functions used for query (C_Q) and background (C_B) must satisfy:
```
∀g ∈ Q ∩ B: C_Q(g) = C_B(g)
```

**Proof:** If copy numbers differ for the same gene in query and background, the hypergeometric model assumes different populations, violating the fundamental assumption of sampling without replacement from a fixed population. □

### 4.3 Extreme Case Analysis

**Theorem 4.3 (Degenerate Case Handling):**

1. **Empty pathway (P = ∅):** m_w = 0, q_w = 0, P(X ≥ 1) = 0

2. **Empty query (Q = ∅):** k_w = 0, q_w = 0, P(X ≥ 1) = 1

3. **Complete pathway (P = B):** n_w = 0, q_w = k_w, test becomes trivial

4. **Zero copy genes:** If C(g) = 0 for some g ∈ B, then g effectively removed from analysis

**Proof:** Each case follows from parameter definitions and hypergeometric distribution properties. □

## 5. Statistical Properties and Validation

### 5.1 Null Distribution Preservation

**Theorem 5.1 (Null Hypothesis Invariance):**
Under the null hypothesis H₀ (no enrichment), both standard and weighted approaches test the same hypothesis:

H₀: Query instances are randomly sampled from background instances without regard to pathway membership.

**Proof:** The hypergeometric distribution models random sampling without replacement. Both methods model this sampling process - standard method with gene units, weighted method with gene instance units. The null hypothesis of random sampling is preserved under parameter transformation. □

### 5.2 Convergence Properties

**Theorem 5.2 (Uniform Copy Number Convergence):**
If C(g) = c for all g ∈ B (uniform copy numbers), then:
```
lim_{c→1} (k_w/k₀, m_w/m₀, n_w/n₀, q_w/q₀) = (1, 1, 1, 1)
```

**Proof:**
```
k_w/k₀ = (Σ_{g∈Q} c) / |Q| = c|Q| / |Q| = c → 1 as c → 1
```
Similar for other parameters. □

**Corollary 5.1:** Weighted and standard methods converge as copy number variation decreases.

## 6. Computational Complexity Analysis

### 6.1 Algorithm Complexity

**Theorem 6.1 (Computational Complexity):**

1. **Instance Expansion Method:**
   - Time: O(Σ_{g∈B} C(g))
   - Space: O(Σ_{g∈B} C(g))

2. **Parameter Weighting Method:**
   - Time: O(|B|)
   - Space: O(|B|)

**Proof:** Instance expansion requires creating Σ_{g∈B} C(g) instance objects. Parameter weighting requires one pass through |B| genes to compute sums. □

**Corollary 6.1 (Efficiency Gain):**
Let κ = (Σ_{g∈B} C(g)) / |B| be the average copy number. Then parameter weighting achieves:
- Time speedup: O(κ)
- Space reduction: O(κ)

For PHR datasets with κ ≈ 34, this represents ~34× improvement.

## 7. Practical Implementation Considerations

### 7.1 Numerical Stability

**Theorem 7.1 (Overflow Bounds):**
For parameters (q_w, m_w, n_w, k_w), the hypergeometric computation is numerically stable if:
```
max(m_w, n_w, k_w) < 2³¹ - 1    (for 32-bit integers)
```

**Mitigation:** Use appropriate data types (64-bit integers) for large copy numbers.

### 7.2 Precision Requirements

**Theorem 7.2 (Floating Point Precision):**
For copy numbers C(g) ≤ C_max and background size |B| ≤ B_max, the relative error in parameter computation is bounded by:
```
ε_rel ≤ (B_max × C_max × ε_machine) / min(parameters)
```

where ε_machine is machine epsilon.

## 8. Conclusion

This mathematical formulation establishes the theoretical foundation for copy-number weighted hypergeometric testing. Key results include:

1. **Formal parameter transformation equations** with rigorous mathematical notation
2. **Proof of mathematical equivalence** between parameter weighting and instance expansion
3. **Complete constraint analysis** ensuring statistical validity
4. **Complexity analysis** demonstrating computational advantages

The approach provides a statistically sound and computationally efficient method for incorporating copy number information into hypergeometric enrichment testing, with strong theoretical guarantees and practical advantages over naive instance expansion methods.

## Appendix: Algorithm Implementation

### A.1 Formal Algorithm

```
Algorithm: WeightedHypergeometricTest
Input: Q (query genes), P (pathway genes), B (background genes), C (copy function)
Output: p-value

1. Validate inputs:
   - Ensure Q ⊆ B, P ⊆ B
   - Verify C(g) > 0 for all g ∈ B
   
2. Compute weighted parameters:
   k_w ← Σ_{g∈Q} C(g)
   m_w ← Σ_{g∈P∩B} C(g)
   n_w ← Σ_{g∈B\P} C(g)
   q_w ← Σ_{g∈Q∩P} C(g)
   
3. Validate constraints:
   - Check 0 ≤ q_w ≤ min(k_w, m_w)
   - Check k_w ≤ m_w + n_w
   
4. Compute hypergeometric probability:
   return phyper(q_w - 1, m_w, n_w, k_w, lower.tail = FALSE)
```

### A.2 Complexity Analysis Summary

| Aspect | Instance Expansion | Parameter Weighting | Improvement Factor |
|--------|-------------------|--------------------|--------------------|
| Time | O(Σ C(g)) | O(\|B\|) | κ (avg copy number) |
| Space | O(Σ C(g)) | O(\|B\|) | κ (avg copy number) |
| Accuracy | Exact | Exact | 1 (equivalent) |