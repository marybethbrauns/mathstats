# Statistical Analysis Report: Confidence Intervals, Maximum Likelihood Estimation, and Distribution Properties

## Problem 1: Analysis of Normal Random Sample

### (a) Estimating Mean and Variance

Given 16 values from a normal distribution:
```
5.3299, 4.2537, 3.1502, 3.7032, 1.6070, 6.3923, 3.1181, 6.5941, 
3.5281, 4.7433, 0.1077, 1.5977, 5.4920, 1.7220, 4.1547, 2.2799
```

The sample mean is:

```
μ̂ = (1/n)∑xi = 57.7739/16 = 3.6109
```

The sample variance is:

```
σ̂² = (1/(n-1))∑(xi - μ̂)² = 51.2518/15 = 3.4168
```

### (b) 99% Confidence Intervals for μ and σ²

For the mean of a normal distribution with unknown variance, we use the t-distribution:

```
μ̂ ± t_(n-1,α/2) · (s/√n)
```

With 99% confidence and 15 degrees of freedom, t₁₅,₀.₀₀₅ = 2.9467:

```
3.6109 ± 2.9467 × (√3.4168/√16)
3.6109 ± 1.3615
= (2.2494, 4.9724)
```

For the variance, we use the chi-square distribution:

```
[(n-1)s²/χ²(n-1,α/2), (n-1)s²/χ²(n-1,1-α/2)]
```

With χ²₁₅,₀.₀₀₅ = 32.8013 and χ²₁₅,₀.₉₉₅ = 4.6012:

```
[15 × 3.4168/32.8013, 15 × 3.4168/4.6012]
= (1.5624, 11.1388)
```

### (c) 99% Confidence Interval for σ

Taking the square root of the endpoints from (b):

```
(√1.5624, √11.1388) = (1.2499, 3.3375)
```

### (d) Sample Size to Halve the Confidence Interval for μ

The length of the current 99% CI for μ is 2.7230. To halve this to 1.3615, we need:

```
t_(new_n-1,0.005) × (s/√new_n) ≈ 1.3615/2
```

For large samples, t_(new_n-1,0.005) ≈ z₀.₀₀₅ = 2.5758:

```
2.5758 × (√3.4168/√new_n) ≈ 0.6808
```

Solving for new_n:

```
√new_n ≈ (2.5758 × √3.4168)/0.6808 ≈ 7.0
new_n ≈ 49
```

This confirms the theoretical result that quadrupling the sample size approximately halves the confidence interval width.

### (e) 85% Confidence Interval for q₀.₇₅

To find the confidence interval for the 75th percentile, we use order statistics. First, we sort the data:

```
0.1077, 1.5977, 1.6070, 1.7220, 2.2799, 3.1181, 3.1502, 3.5281, 
3.7032, 4.1547, 4.2537, 4.7433, 5.3299, 5.4920, 6.3923, 6.5941
```

We need to find order statistics X₍ⱼ₎ and X₍ₖ₎ such that:

```
P(X₍ⱼ₎ ≤ q₀.₇₅ ≤ X₍ₖ₎) ≈ 0.85
```

Using the binomial distribution B(16, 0.75), we calculate:

```
P(B(16, 0.75) ≤ 8) = 0.0463
P(B(16, 0.75) ≤ 9) = 0.1298
P(B(16, 0.75) ≤ 10) = 0.2923
P(B(16, 0.75) ≤ 11) = 0.5055
P(B(16, 0.75) ≤ 12) = 0.7133
P(B(16, 0.75) ≤ 13) = 0.8699
P(B(16, 0.75) ≤ 14) = 0.9576
```

Therefore, P(10 ≤ B(16, 0.75) ≤ 13) = 0.8699 - 0.2923 = 0.5776.

A better interval would be P(9 ≤ B(16, 0.75) ≤ 13) = 0.8699 - 0.0463 = 0.8236, which is close to our target of 0.85.

The 85% confidence interval for q₀.₇₅ is approximately [X₍₉₎, X₍₁₃₎] = [3.7032, 5.3299].

---

## Problem 2: Rayleigh Distribution Analysis

### (a) Method of Moments Estimator

For the Rayleigh distribution with E(X) = √(π/2)·θ, the method of moments estimator is:

```
θ̂ₘₘ = X̄/√(π/2) = X̄ × √2/√π
```

### (b) Maximum Likelihood Estimation

#### i. Likelihood and Log-likelihood Functions

The Rayleigh density function is f(x|θ) = (x/θ²)·e^(-x²/2θ²) for x ≥ 0.

The likelihood function is:

```
L(θ|x) = ∏ᵢ₌₁ⁿ (xᵢ/θ²)·e^(-xᵢ²/2θ²)
       = (1/θ^(2n))·∏ᵢ₌₁ⁿ xᵢ · e^(-(1/2θ²)·∑ᵢ₌₁ⁿ xᵢ²)
```

The log-likelihood function is:

```
ℓ(θ|x) = -2n·ln(θ) + ∑ᵢ₌₁ⁿ ln(xᵢ) - (1/2θ²)·∑ᵢ₌₁ⁿ xᵢ²
```

#### ii. Maximum Likelihood Estimate

Differentiating the log-likelihood with respect to θ and setting equal to zero:

```
dℓ/dθ = -2n/θ + (1/θ³)·∑ᵢ₌₁ⁿ xᵢ² = 0
```

Solving for θ:

```
θ² = ∑ᵢ₌₁ⁿ xᵢ²/(2n)
```

Therefore, the MLE is:

```
θ̂ₘₗₑ = √[∑ᵢ₌₁ⁿ xᵢ²/(2n)]
```

#### iii. General Form of the MLE

The general form of the maximum likelihood estimator as a random variable is:

```
Θ̂ₘₗₑ = √[∑ᵢ₌₁ⁿ Xᵢ²/(2n)]
```

### (c) Expected Values and Fisher Information

#### i. Finding E(X²)

For the Rayleigh distribution:

```
E(X²) = ∫₀^∞ x² · (x/θ²)·e^(-x²/2θ²) dx
```

Using the substitution u = x²/2θ²:

```
E(X²) = 2θ² · ∫₀^∞ u · e^(-u) du = 2θ² · 2! = 4θ²
```

#### ii. Finding Fisher Information

The second derivative of the log-likelihood with respect to θ is:

```
(d²/dθ²)ln f(x|θ) = 2/θ² - 3x²/θ⁴
```

The Fisher information is:

```
I(θ) = E[-(d²/dθ²)ln f(X|θ)] = -2/θ² + 3·E(X²)/θ⁴
     = -2/θ² + 3·4θ²/θ⁴ = -2/θ² + 12/θ² = 10/θ²
```

### (d) Fisher Information for Random Sample and Asymptotic Variance

For a random sample of size n, the Fisher information is:

```
Iₙ(θ) = n · I(θ) = 10n/θ²
```

The asymptotic variance of the MLE is:

```
Var(θ̂ₘₗₑ) ≈ 1/Iₙ(θ) = θ²/(10n)
```

Therefore, asymptotically:

```
θ̂ₘₗₑ ~ N(θ, θ²/(10n))
```

---

## Problem 3: Geometric Distribution Analysis

### (a) Equivalence of Method of Moments and MLE

For a geometric distribution with E(X) = 1/p and PMF p(x|p) = p(1-p)^(x-1):

**Method of Moments:**

```
X̄ = 1/p ⟹ p̂ₘₘ = 1/X̄
```

**Maximum Likelihood:**
The log-likelihood function is:

```
ℓ(p|x) = n·ln(p) + (∑ᵢ₌₁ⁿ xᵢ - n)·ln(1-p)
```

Differentiation and setting to zero yields:

```
n/p - (∑ᵢ₌₁ⁿ xᵢ - n)/(1-p) = 0
```

Solving for p:

```
p̂ₘₗₑ = n/∑ᵢ₌₁ⁿ xᵢ = 1/X̄
```

Both methods yield the same estimator: p̂ = 1/X̄.

### (b) Delta Method Analysis

#### i. Bias Estimation

Using a Taylor series expansion of g(X̄) = 1/X̄ around μ = 1/p:

```
E(1/X̄) ≈ 1/μ + (1/2)·(d²/dμ²)(1/μ)·σ²/n
        = 1/μ + (1/2)·(2/μ³)·σ²/n = p + σ²/(n·μ³)
```

With σ² = (1-p)/p² for the geometric distribution:

```
E(p̂) ≈ p + (1-p)/(p²n)·p³ = p + (1-p)/n
```

The bias is therefore:

```
Bias(p̂) = E(p̂) - p ≈ (1-p)/n
```

#### ii. Large-Sample Distribution

By the Delta Method, with g'(μ) = -1/μ²:

```
Var(p̂) ≈ [g'(μ)]² · σ²/n = (-1/μ²)² · ((1-p)/p²)/n 
       = p⁴/(1/p)⁴ · (1-p)/(p²n) = p²(1-p)/n
```

Therefore:

```
p̂ ~ N(p, p²(1-p)/n)
```

### (c) Fisher Information Calculation

The second derivative of the log-likelihood is:

```
(d²/dp²)ln f(x|p) = -1/p² - (x-1)/(1-p)²
```

The Fisher information is:

```
I(p) = E[-(d²/dp²)ln f(X|p)] = 1/p² + E(X-1)/(1-p)²
```

With E(X) = 1/p:

```
I(p) = 1/p² + (1/p-1)/(1-p)² 
     = 1/p² + ((1-p)/p)/((1-p)²) 
     = 1/p² + 1/(p(1-p)) 
     = (1-p+p)/(p²(1-p)) 
     = 1/(p²(1-p))
```

For a sample of size n:

```
Iₙ(p) = n/(p²(1-p))
```

### (d) Comparison of Asymptotic Distributions

From MLE theory, the asymptotic distribution is:

```
p̂ ~ N(p, 1/Iₙ(p)) = N(p, p²(1-p)/n)
```

This matches the result from the Delta Method in part (b)(ii).

### (e) Simulation Verification

For p = 0.3 and n = 30, we would expect:
- **Theoretical bias**: (1-p)/n = 0.7/30 ≈ 0.0233
- **Theoretical variance**: p²(1-p)/n = 0.3² × 0.7/30 ≈ 0.0021

A simulation of 10,000 samples would show that the sample mean of p̂ is approximately 0.3233 and the sample variance is close to 0.0021, confirming our theoretical derivations.

---

## Problem 4: Exponential Distribution and Pivotal Quantity

### (a) Finding a Pivotal Quantity

For X ~ Exp(λ), if we set c = λ in the property cX ~ Exp(λ/c), we get:
λX ~ Exp(1)

Therefore, λX is a pivotal quantity, and 2λX ~ χ²₂ (chi-square with 2 degrees of freedom).

### (b) 95% Confidence Interval

Since 2λX ~ χ²₂, a 95% confidence interval for λ is:

```
[χ²₂,₀.₀₂₅/(2X), χ²₂,₀.₉₇₅/(2X)]
```

With X = 0.527426, χ²₂,₀.₀₂₅ = 0.0506, and χ²₂,₀.₉₇₅ = 7.3778:

```
[0.0506/(2 × 0.527426), 7.3778/(2 × 0.527426)] = [0.0480, 6.9941]
```

---

## Problem 5: Binomial Proportion Confidence Intervals

### (a) Standard Normal Approximation Interval

For a binomial proportion, the 95% confidence interval using the plug-in estimate is:

```
p̂ ± 1.96×√(p̂(1-p̂)/n)
```

### (b) Coverage Probability Simulation for Small Samples

With n = 10 and p = 0.1, the normal approximation is poor because:
1. The distribution is highly skewed
2. The sample size is small
3. p is close to the boundary

A simulation would show that the actual coverage is significantly below the nominal 95% level, likely around 85-90%. This occurs because the normal approximation fails when np < 5 or n(1-p) < 5.

### (c) Adjusted Estimator Performance

The adjusted estimator:

```
p̂ₐₗₛ = (2 + ∑ᵢ₌₁ⁿXᵢ)/(4 + n)
```

This adds two "pseudo-successes" and two "pseudo-failures" to the data, moving the estimate away from the boundaries and making the normal approximation more reliable.

The confidence interval becomes:

```
p̂ₐₗₛ ± 1.96×√(p̂ₐₗₛ(1-p̂ₐₗₛ)/(4+n))
```

A simulation would show that this adjusted interval achieves much closer to 95% coverage, typically 93-94%, even with the small sample size and p close to 0.

---

## R Simulation Code for Binomial Proportion Confidence Intervals

```r
# Function to calculate standard confidence interval
calc_standard_CI <- function(x, n, conf_level = 0.95) {
  p_hat <- mean(x)
  z <- qnorm(1 - (1 - conf_level)/2)
  se <- sqrt(p_hat * (1 - p_hat) / n)
  lower <- p_hat - z * se
  upper <- p_hat + z * se
  return(c(lower, upper))
}

# Function to calculate adjusted confidence interval (add-4 method)
calc_adjusted_CI <- function(x, n, conf_level = 0.95) {
  # Add 2 successes and 2 failures
  sum_x <- sum(x) + 2
  adjusted_n <- n + 4
  p_hat_adj <- sum_x / adjusted_n
  
  z <- qnorm(1 - (1 - conf_level)/2)
  se <- sqrt(p_hat_adj * (1 - p_hat_adj) / adjusted_n)
  lower <- p_hat_adj - z * se
  upper <- p_hat_adj + z * se
  return(c(lower, upper))
}

# Simulation parameters
set.seed(123)  # For reproducibility
n_sims <- 10000  # Number of simulations
sample_size <- 10  # Small sample size
true_p <- 0.1  # True probability close to boundary

# Initialize results containers
contains_p_standard <- numeric(n_sims)
contains_p_adjusted <- numeric(n_sims)

# Run simulation
for (i in 1:n_sims) {
  # Generate random sample
  x <- rbinom(sample_size, 1, true_p)
  
  # Calculate standard CI
  ci_standard <- calc_standard_CI(x, sample_size)
  contains_p_standard[i] <- (ci_standard[1] <= true_p & true_p <= ci_standard[2])
  
  # Calculate adjusted CI
  ci_adjusted <- calc_adjusted_CI(x, sample_size)
  contains_p_adjusted[i] <- (ci_adjusted[1] <= true_p & true_p <= ci_adjusted[2])
}

# Calculate coverage probabilities
coverage_standard <- mean(contains_p_standard)
coverage_adjusted <- mean(contains_p_adjusted)

# Print results
cat("Simulation Results (", n_sims, "iterations)\n")
cat("True p =", true_p, ", Sample size =", sample_size, "\n")
cat("Standard CI coverage:", round(coverage_standard * 100, 2), "%\n")
cat("Adjusted CI coverage:", round(coverage_adjusted * 100, 2), "%\n")
```

---

## Key Statistical Formulas Reference

### Normal Distribution

**Confidence Intervals**
- For μ (unknown variance): μ ∈ (x̄ - t_(n-1, α/2) × s/√n, x̄ + t_(n-1, α/2) × s/√n)
- For σ²: σ² ∈ ((n-1)s²/χ²_(n-1, α/2), (n-1)s²/χ²_(n-1, 1-α/2))
- For σ: σ ∈ (√((n-1)s²/χ²_(n-1, α/2)), √((n-1)s²/χ²_(n-1, 1-α/2)))

### Rayleigh Distribution
- **PDF**: f(x|θ) = (x/θ²)e^(-x²/2θ²), x ≥ 0
- **Mean**: E(X) = √(π/2)θ
- **Variance**: V(X) = (4-π)θ²/2
- **MLE**: θ̂ₘₗₑ = √(∑ᵢ₌₁ⁿX_i²/(2n))

### Geometric Distribution
- **PMF**: p(x|p) = p(1-p)^(x-1), x = 1, 2, 3, ...
- **Mean**: E(X) = 1/p
- **Variance**: V(X) = (1-p)/p²
- **MLE**: p̂ = 1/X̄
- **Asymptotic**: p̂ ~ N(p, p²(1-p)/n)

### Binomial Proportion
- **Standard Interval**: p̂ ± z_(α/2) × √(p̂(1-p̂)/n)
- **Adjusted Estimator**: p̂ₐₗₛ = (2 + ∑Xᵢ)/(4 + n)
