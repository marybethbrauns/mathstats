# Vector and Matrix Norms: Theoretical Analysis

## Introduction

This document presents formal proofs and counterexamples relating to vector and matrix norms in finite-dimensional spaces. The analysis covers verification of norm properties, counterexample construction, limit characterization of the infinity norm, and transitivity of norm equivalence. These results are fundamental in functional analysis, numerical linear algebra, and approximation theory.

---

## Problem 1: Verification of the 1-Norm

### Problem Statement
Verify that the 1-norm of a vector x in R^n, defined as ||x||₁ = Σ(i=1 to n) |xᵢ|, satisfies all the necessary properties to be a valid norm.

### Background
A function ||·||: R^n → R is a norm if it satisfies three fundamental properties for all vectors in the space and all scalars. These properties establish the geometric and algebraic behavior of the norm and ensure its usefulness in analysis of vector spaces.

### Solution

To establish that ||x||₁ = Σ(i=1 to n) |xᵢ| is a valid norm, we must verify the following three properties:

#### 1.1 Positive Definiteness
We must show that:
- ||x||₁ ≥ 0 for all x in R^n
- ||x||₁ = 0 if and only if x = 0

**Proof:**
- Since |xᵢ| ≥ 0 for all i in {1,2,...,n}, it follows that Σ(i=1 to n) |xᵢ| ≥ 0 for all x in R^n.
- If x = 0, then |xᵢ| = 0 for all i in {1,2,...,n}, thus ||x||₁ = Σ(i=1 to n) |xᵢ| = 0.
- Conversely, if ||x||₁ = 0, then Σ(i=1 to n) |xᵢ| = 0. Since each term |xᵢ| ≥ 0, we must have |xᵢ| = 0 for all i in {1,2,...,n}, which implies xᵢ = 0 for all i, thus x = 0.

#### 1.2 Absolute Homogeneity
We must prove that ||αx||₁ = |α|·||x||₁ for all α in R and x in R^n.

**Proof:**
||αx||₁ = Σ(i=1 to n) |αxᵢ|
       = Σ(i=1 to n) |α|·|xᵢ|
       = |α|·Σ(i=1 to n) |xᵢ|
       = |α|·||x||₁

#### 1.3 Triangle Inequality
We must demonstrate that ||x + y||₁ ≤ ||x||₁ + ||y||₁ for all x,y in R^n.

**Proof:**
||x + y||₁ = Σ(i=1 to n) |xᵢ + yᵢ|
          ≤ Σ(i=1 to n) (|xᵢ| + |yᵢ|)  (by the triangle inequality for real numbers)
          = Σ(i=1 to n) |xᵢ| + Σ(i=1 to n) |yᵢ|
          = ||x||₁ + ||y||₁

Since all three properties are satisfied, ||x||₁ = Σ(i=1 to n) |xᵢ| is indeed a valid norm on R^n.

---

## Problem 2: Counterexample for Matrix Infinity Norm

### Problem Statement
Show by counterexample that the infinity norm of a matrix A in R^(m×n), defined as ||A||∞ = max(i,j) |aᵢⱼ|, does not define a valid matrix norm.

### Background
A matrix norm ||·||: R^(m×n) → R must satisfy the same three properties as vector norms (positive definiteness, absolute homogeneity, and triangle inequality), but additionally must satisfy the submultiplicative property: ||AB|| ≤ ||A||·||B|| for compatible matrices A and B. This property is essential for many applications in numerical analysis and linear algebra.

### Solution

To disprove that ||A||∞ = max(i,j) |aᵢⱼ| is a valid matrix norm, it suffices to find a counterexample to the submultiplicative property.

**Counterexample:**
Consider the matrices:
A = [1 1]
    [1 1]

and

B = [1 1]
    [1 1]

Computing the infinity norm as defined:
- ||A||∞ = max(i,j) |aᵢⱼ| = 1
- ||B||∞ = max(i,j) |bᵢⱼ| = 1

The matrix product is:
AB = [1 1] [1 1] = [2 2]
     [1 1] [1 1]   [2 2]

Thus ||AB||∞ = max(i,j) |(AB)ᵢⱼ| = 2

We can now verify whether the submultiplicative property holds:
||AB||∞ = 2 > 1 = ||A||∞ · ||B||∞

This directly contradicts the submultiplicative property. Therefore, ||A||∞ = max(i,j) |aᵢⱼ| does not define a valid matrix norm. It's worth noting that while this "max norm" fails as a matrix norm, it does define a valid vector norm when applied to vectors.

---

## Problem 3: Characterization of Vector Infinity Norm

### Problem Statement
Prove that the infinity norm of a vector x in R^n, defined as ||x||∞ = max(1≤i≤n) |xᵢ|, can be characterized as the limit: ||x||∞ = lim(p→∞) (Σ(i=1 to n) |xᵢ|^p)^(1/p).

### Background
The family of p-norms for vectors provides a flexible framework for measuring vector magnitudes with different characteristics. The relationship between different p-norms, especially as p approaches infinity, reveals important theoretical connections in functional analysis and provides insight into the geometry of normed spaces.

### Solution

We need to prove that:
||x||∞ = max(1≤i≤n) |xᵢ| = lim(p→∞) (Σ(i=1 to n) |xᵢ|^p)^(1/p)

Let M = max(1≤i≤n) |xᵢ| = ||x||∞. We can rewrite the p-norm as:

(Σ(i=1 to n) |xᵢ|^p)^(1/p) = M · (Σ(i=1 to n) (|xᵢ|/M)^p)^(1/p)

Observe that for each i, |xᵢ|/M ≤ 1 with equality for at least one index i (specifically, for any i where |xᵢ| = M).

Let S = {i : |xᵢ| = M} be the set of indices where the maximum is attained, and let k = |S| be the cardinality of this set. Then:

(Σ(i=1 to n) |xᵢ|^p)^(1/p) = M · (Σ(i∈S) (|xᵢ|/M)^p + Σ(i∉S) (|xᵢ|/M)^p)^(1/p)
                            = M · (k + Σ(i∉S) (|xᵢ|/M)^p)^(1/p)

For any i∉S, we have |xᵢ|/M < 1, so (|xᵢ|/M)^p → 0 as p → ∞.

Thus:
lim(p→∞) (Σ(i=1 to n) |xᵢ|^p)^(1/p) = M · lim(p→∞) (k + Σ(i∉S) (|xᵢ|/M)^p)^(1/p)
                                      = M · lim(p→∞) (k)^(1/p) = M · 1 = M = ||x||∞

since lim(p→∞) k^(1/p) = 1 for any fixed positive integer k.

This establishes that the infinity norm of a vector can indeed be characterized as the limit of its p-norms as p approaches infinity.

---

## Problem 4: Transitivity of Norm Equivalence

### Problem Statement
Prove that if two vector norms ||·||ₐ and ||·||ᵦ are each equivalent to a third vector norm ||·||ᵧ, then ||·||ₐ and ||·||ᵦ are equivalent to each other.

### Background
Norm equivalence is a fundamental concept in functional analysis. Two norms ||·||ₐ and ||·||ᵦ on a vector space V are equivalent if there exist positive constants c and C such that c||x||ₐ ≤ ||x||ᵦ ≤ C||x||ₐ for all x in V. Equivalent norms induce the same topology on the vector space, meaning that convergence in one norm implies convergence in the other. The transitivity of this relation has important implications for analyzing different norms in the same space.

### Solution

Suppose we have three norms ||·||ₐ, ||·||ᵦ, and ||·||ᵧ on a vector space V where:

1. ||·||ₐ is equivalent to ||·||ᵧ, so there exist positive constants c₁, C₁ such that:
   c₁||x||ₐ ≤ ||x||ᵧ ≤ C₁||x||ₐ for all x in V

2. ||·||ᵦ is equivalent to ||·||ᵧ, so there exist positive constants c₂, C₂ such that:
   c₂||x||ᵦ ≤ ||x||ᵧ ≤ C₂||x||ᵦ for all x in V

We need to establish that ||·||ₐ and ||·||ᵦ are equivalent, i.e., find positive constants c₃, C₃ such that:
c₃||x||ₐ ≤ ||x||ᵦ ≤ C₃||x||ₐ for all x in V

**Proof:**

From inequality (1), we have ||x||ᵧ ≤ C₁||x||ₐ for all x in V.

From inequality (2), we have c₂||x||ᵦ ≤ ||x||ᵧ for all x in V.

Combining these inequalities: 
c₂||x||ᵦ ≤ ||x||ᵧ ≤ C₁||x||ₐ

This gives us: 
||x||ᵦ ≤ (C₁/c₂)||x||ₐ

So we have established the upper bound with C₃ = C₁/c₂.

Similarly, from inequality (1), we have c₁||x||ₐ ≤ ||x||ᵧ for all x in V.

From inequality (2), we have ||x||ᵧ ≤ C₂||x||ᵦ for all x in V.

Combining these inequalities:
c₁||x||ₐ ≤ ||x||ᵧ ≤ C₂||x||ᵦ

This gives us:
||x||ₐ ≤ (C₂/c₁)||x||ᵦ

Inverting this inequality:
(c₁/C₂)||x||ₐ ≤ ||x||ᵦ

So we have established the lower bound with c₃ = c₁/C₂.

Combining our results, we have:
(c₁/C₂)||x||ₐ ≤ ||x||ᵦ ≤ (C₁/c₂)||x||ₐ for all x in V

Therefore, ||·||ₐ and ||·||ᵦ are equivalent norms, demonstrating that norm equivalence is indeed a transitive relation. This result has significant implications for the study of normed spaces, as it allows us to establish equivalence classes of norms that induce the same topological and analytical properties on the vector space.
