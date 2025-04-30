# Comprehensive Statistical Analysis with Complete Solutions

## Problem 1: Normal Distribution Hypothesis Testing

**Problem Statement:** Consider a random sample X₁, X₂, ..., Xₙ from a N(μ, 4) distribution. We want to test H₀: μ = 1 versus H₁: μ = 0.5.

### Part (a): Distributions of Xₙ under H₀ and H₁

Since X₁, X₂, ..., Xₙ ~ N(μ, 4) independently, the sample mean Xₙ = (1/n)∑ᵢ₌₁ⁿ Xᵢ follows a normal distribution:

- Under H₀ (μ = 1): Xₙ ~ N(1, 4/n)
- Under H₁ (μ = 0.5): Xₙ ~ N(0.5, 4/n)

**Statistical Connection:** This exact result demonstrates a fundamental sampling distribution property: when sampling from a normal population, the sample mean is precisely normally distributed with the same mean and reduced variance (by a factor of n). This is a special case of the Central Limit Theorem that holds exactly for any sample size when the underlying population is normal.

### Part (b): Finding rejection region and power with α = 0.02

**Step-by-step calculation:**
1) Since H₁ specifies μ = 0.5 which is less than μ = 1 under H₀, we need a one-sided lower-tail test. We reject H₀ when Xₙ < c.

2) For significance level α = 0.02, we need P(Xₙ < c | H₀ true) = 0.02

3) Under H₀: Xₙ ~ N(1, 4/20) = N(1, 0.2) [using n = 20]

4) Standardizing to get a standard normal Z:
   - P(Xₙ < c) = P((Xₙ - 1)/√0.2 < (c - 1)/√0.2) = 0.02
   - For a standard normal, P(Z < -2.0537) = 0.02
   - Therefore: (c - 1)/√0.2 = -2.0537
   - c - 1 = -2.0537 × √0.2 = -2.0537 × 0.4472 = -0.9184
   - c = 1 - 0.9184 = 0.0816 ≈ 0.0812

5) Power calculation (probability of correctly rejecting H₀ when H₁ is true):
   - Under H₁: Xₙ ~ N(0.5, 0.2)
   - Power = P(Xₙ < 0.0812 | H₁ true)
   - = P((Xₙ - 0.5)/√0.2 < (0.0812 - 0.5)/√0.2)
   - = P(Z < (0.0812 - 0.5)/0.4472)
   - = P(Z < -0.9366)
   - = 0.1745 ≈ 0.175 or 17.5%

**Statistical Connection:** This calculation demonstrates the Neyman-Pearson framework for hypothesis testing, showing the trade-off between Type I error (controlled at α = 0.02) and Type II error (β = 0.825). The low power indicates a high probability of missing a true effect, illustrating why researchers must carefully consider effect size, variance, and sample size during experimental design.

### Part (c): Increasing power to 0.9

**Step-by-step derivation:**
1) We need to determine n and c such that:
   - P(Xₙ < c | H₀ true) = 0.02 (control Type I error)
   - P(Xₙ < c | H₁ true) = 0.9 (achieve desired power)

2) Under H₀: Xₙ ~ N(1, 4/n)
   From condition 1:
   - P((Xₙ - 1)/√(4/n) < (c - 1)/√(4/n)) = 0.02
   - (c - 1)/√(4/n) = z₀.₀₂ = -2.054
   - c - 1 = -2.054 × √(4/n)
   - c = 1 - 4.108/√n

3) Under H₁: Xₙ ~ N(0.5, 4/n)
   From condition 2:
   - P((Xₙ - 0.5)/√(4/n) < (c - 0.5)/√(4/n)) = 0.9
   - (c - 0.5)/√(4/n) = z₀.₉ = 1.282
   - c - 0.5 = 1.282 × √(4/n)
   - c = 0.5 + 2.564/√n

4) Setting these equal (the critical value c must be the same for both conditions):
   - 1 - 4.108/√n = 0.5 + 2.564/√n
   - 0.5 = 6.672/√n
   - √n = 13.344
   - n = 178.06 ≈ 179

5) With n = 179, the critical value is:
   - c = 1 - 4.108/√179 = 1 - 4.108/13.38 = 1 - 0.307 = 0.693

**Statistical Connection:** This demonstrates how sample size determines statistical power. The substantial increase from n = 20 to n = 179 illustrates the "cost" of achieving high power when effect sizes are modest relative to variability. This connects to the Law of Large Numbers - as n increases, test statistics become more precise and capable of detecting smaller effects.

## Problem 2: Multinomial Probability and Maximum Likelihood

**Problem Statement:** Consider a genetic model with three genotypes having probabilities (1-θ)², 2θ(1-θ), and θ² respectively.

### Part (a): Log-likelihood function

For a multinomial distribution with three categories having probabilities p₁ = (1-θ)², p₂ = 2θ(1-θ), and p₃ = θ²:

**Step-by-step derivation:**
1) The probability mass function is:
   - P(O₁, O₂, O₃) = n!/(O₁!O₂!O₃!) × p₁^O₁ × p₂^O₂ × p₃^O₃

2) Taking natural logarithm:
   - ℓ(θ) = log[n!/(O₁!O₂!O₃!)] + O₁log((1-θ)²) + O₂log(2θ(1-θ)) + O₃log(θ²)

3) Using logarithm properties:
   - log((1-θ)²) = 2log(1-θ)
   - log(2θ(1-θ)) = log(2) + log(θ) + log(1-θ)
   - log(θ²) = 2log(θ)

4) Simplifying:
   - ℓ(θ) = constant + 2O₁log(1-θ) + O₂[log(2) + log(θ) + log(1-θ)] + 2O₃log(θ)
   - = constant + (2O₁ + O₂)log(1-θ) + (O₂ + 2O₃)log(θ) + O₂log(2)

**Statistical Connection:** This exemplifies Fisher's scoring method, where log-transformation converts multiplicative probabilities to additive terms. The expression reveals sufficient statistics - the entire dataset's information for estimating θ is captured in just the count combinations (2O₁ + O₂) and (O₂ + 2O₃).

### Part (b): Maximum likelihood estimator

**Step-by-step derivation:**
1) To find the MLE, differentiate ℓ(θ) with respect to θ and set equal to zero:
   - dℓ/dθ = -(2O₁ + O₂)/(1-θ) + (O₂ + 2O₃)/θ = 0

2) Solving for θ:
   - (2O₁ + O₂)/(1-θ) = (O₂ + 2O₃)/θ
   - θ(2O₁ + O₂) = (1-θ)(O₂ + 2O₃)
   - θ(2O₁ + O₂) = O₂ + 2O₃ - θO₂ - 2θO₃
   - θ(2O₁ + O₂) + θO₂ + 2θO₃ = O₂ + 2O₃
   - θ(2O₁ + 2O₂ + 2O₃) = O₂ + 2O₃
   - 2θ(O₁ + O₂ + O₃) = O₂ + 2O₃
   - 2nθ = O₂ + 2O₃     (since O₁ + O₂ + O₃ = n)

3) Therefore:
   - θ̂ = (O₂ + 2O₃)/(2n)

**Statistical Connection:** The MLE has a beautiful genetic interpretation - it estimates allele frequency by counting alleles (heterozygotes contribute one copy, homozygotes two copies) divided by total possible alleles (2n). This demonstrates how maximum likelihood translates a probabilistic model into an optimal parameter estimator.

### Part (c): Fisher information and asymptotic distribution

**Step-by-step derivation:**
1) Calculate the second derivative of log-likelihood:
   - d²ℓ/dθ² = -(2O₁ + O₂)/(1-θ)² - (O₂ + 2O₃)/θ²

2) The Fisher information is I(θ) = -E[d²ℓ/dθ²]. Using:
   - E[2O₁ + O₂] = 2n(1-θ)
   - E[O₂ + 2O₃] = 2nθ

3) Therefore:
   - I(θ) = 2n/(1-θ) + 2n/θ = 2n/(θ(1-θ))

4) By asymptotic theory of maximum likelihood estimation:
   - θ̂ ~ N(θ, 1/I(θ)) = N(θ, θ(1-θ)/(2n))

**Statistical Connection:** This connects to the Cramér-Rao lower bound theorem, which establishes the theoretical minimum variance for any unbiased estimator. The Fisher information I(θ) quantifies how much "statistical information" the data provides about the parameter θ, with larger values indicating more precise estimation potential.

### Part (d): MLE and confidence interval for Plato data

**Step-by-step calculation:**
1) Given data: O₁ = 10, O₂ = 68, O₃ = 112, n = 190

2) Calculate MLE:
   - θ̂ = (O₂ + 2O₃)/(2n) = (68 + 2×112)/(2×190) = (68 + 224)/380 = 292/380 = 0.7684

3) For a 99% confidence interval, we use the asymptotic variance from part (c):
   - Var(θ̂) = θ̂(1-θ̂)/(2n) = 0.7684×(1-0.7684)/(2×190) = 0.7684×0.2316/380 = 0.0004682
   - Standard error = √Var(θ̂) = 0.02164

4) For 99% confidence:
   - z₀.₉₉₅ = 2.576 (the z-value cutting off 0.5% in each tail)
   - 99% CI = θ̂ ± z₀.₉₉₅ × SE = 0.7684 ± 2.576 × 0.02164
   - = 0.7684 ± 0.0557 = (0.713, 0.824)

**Statistical Connection:** This application demonstrates how asymptotic theory enables practical inference when exact distributions are intractable. The confidence interval embodies the sampling variability principle - if we repeated the experiment many times, about 99% of similarly constructed intervals would contain the true parameter.

### Part (e): Likelihood ratio test

**Step-by-step calculation:**
1) Parameter spaces:
   - Θ₀ = {1/2} (under H₀: θ = 1/2, Hardy-Weinberg equilibrium)
   - Θ = (0,1) (full parameter space)

2) Expected cell counts under H₀ (θ = 0.5):
   - p₁ = (1-0.5)² = 0.25
   - p₂ = 2×0.5×0.5 = 0.5
   - p₃ = 0.5² = 0.25
   - E₁ = n×p₁ = 190×0.25 = 47.5
   - E₂ = n×p₂ = 190×0.5 = 95
   - E₃ = n×p₃ = 190×0.25 = 47.5

3) Pearson chi-square statistic:
   - χ² = Σ(Oᵢ - Eᵢ)²/Eᵢ
   - = (10-47.5)²/47.5 + (68-95)²/95 + (112-47.5)²/47.5
   - = (-37.5)²/47.5 + (-27)²/95 + (64.5)²/47.5
   - = 1406.25/47.5 + 729/95 + 4160.25/47.5
   - = 29.61 + 7.67 + 87.58 = 124.86

4) For the likelihood ratio test:
   - Under θ̂ = 0.7684:
     * p̂₁ = (1-θ̂)² = 0.2316² = 0.05364
     * p̂₂ = 2θ̂(1-θ̂) = 2(0.7684)(0.2316) = 0.3560
     * p̂₃ = θ̂² = 0.7684² = 0.5904
   
   - Under θ₀ = 0.5:
     * p₀,₁ = 0.25, p₀,₂ = 0.5, p₀,₃ = 0.25

   - −2ln(Λ) = 2∑ᵢ₌₁³ Oᵢ ln(p̂ᵢ/p₀,ᵢ)
   - O₁ ln(p̂₁/p₀,₁) = 10 × ln(0.05364/0.25) = 10 × ln(0.21456) = -15.394
   - O₂ ln(p̂₂/p₀,₂) = 68 × ln(0.3560/0.5) = 68 × ln(0.712) = -23.106
   - O₃ ln(p̂₃/p₀,₃) = 112 × ln(0.5904/0.25) = 112 × ln(2.3616) = 96.242
   - −2ln(Λ) = 2[-15.394 - 23.106 + 96.242] = 2[57.742] = 115.5

5) Both test statistics follow a χ² distribution with df = 1
   - P-values are extremely small (< 10⁻¹⁵)
   - Therefore, we strongly reject H₀: θ = 0.5

6) This is consistent with our 99% confidence interval (0.713, 0.824), which excludes 0.5.

**Statistical Connection:** This applies Wilks' theorem, which states that -2log(Λ) asymptotically follows a chi-square distribution with degrees of freedom equal to the difference in dimensionality between parameter spaces. The test provides overwhelming evidence against Hardy-Weinberg equilibrium (θ = 0.5) in this population.

## Problem 3: Chi-Square Statistic for Binomial

**Problem Statement:** Show that for a binomial distribution with two categories, the Pearson chi-square statistic can be expressed as the square of a standardized normal random variable.

**Step-by-step derivation:**
1) For a binomial with two categories and probabilities p₁ and p₂ = 1-p₁:
   - The Pearson chi-square statistic is:
   - χ² = Σ(Xᵢ - npᵢ)²/(npᵢ) = (X₁ - np₁)²/(np₁) + (X₂ - np₂)²/(np₂)

2) Using the constraint X₁ + X₂ = n:
   - X₂ = n - X₁
   
3) Using the constraint p₁ + p₂ = 1:
   - p₂ = 1 - p₁

4) Substituting these constraints:
   - χ² = (X₁ - np₁)²/(np₁) + ((n - X₁) - n(1 - p₁))²/(n(1 - p₁))
   - = (X₁ - np₁)²/(np₁) + (n - X₁ - n + np₁)²/(n(1 - p₁))
   - = (X₁ - np₁)²/(np₁) + (np₁ - X₁)²/(n(1 - p₁))

5) Noting that (np₁ - X₁) = -(X₁ - np₁):
   - χ² = (X₁ - np₁)²/(np₁) + (X₁ - np₁)²/(n(1 - p₁))
   - = (X₁ - np₁)² × [1/(np₁) + 1/(n(1 - p₁))]

6) Simplifying the common factor:
   - 1/(np₁) + 1/(n(1 - p₁)) = [(1-p₁) + p₁]/(np₁(1-p₁)) = 1/(np₁(1-p₁))

7) Therefore:
   - χ² = (X₁ - np₁)²/(np₁(1-p₁)) = [(X₁ - np₁)/√(np₁(1-p₁))]²

**Statistical Connection:** This elegant result reveals a fundamental relationship between the chi-square and normal distributions: the chi-square statistic with 1 degree of freedom is precisely the square of a standard normal random variable. For large n, X₁ approximately follows a normal distribution with mean np₁ and variance np₁(1-p₁). The constraints reduced the degrees of freedom from 2 categories to 1 parameter, explaining why the chi-square distribution has 1 degree of freedom in this case.

## Problem 4: Comparing Means of Two Normal Distributions

**Problem Statement:** Compare the means of two normal distributions using sample data.

### Part (a): Estimating means, difference, and variance

**Step-by-step calculation:**
1) For first group (x₁): (1.1650, 0.6268, 0.0751, 0.3516)
   - Sample mean:
   - x̄₁ = (1.1650 + 0.6268 + 0.0751 + 0.3516)/4 = 2.2185/4 = 0.554625

2) For second group (x₂): (0.3035, 2.6961, 1.0591, 2.7971, 1.2641)
   - Sample mean:
   - x̄₂ = (0.3035 + 2.6961 + 1.0591 + 2.7971 + 1.2641)/5 = 8.1199/5 = 1.62398

3) Difference of means:
   - x̄₂ - x̄₁ = 1.62398 - 0.554625 = 1.069355

4) Sample variances:
   - For first group:
     * s₁² = Σ(x₁ᵢ - x̄₁)²/(n₁-1)
     * = [(1.1650-0.5546)² + (0.6268-0.5546)² + (0.0751-0.5546)² + (0.3516-0.5546)²]/3
     * = [(0.6104)² + (0.0722)² + (-0.4795)² + (-0.2030)²]/3
     * = [0.3726 + 0.0052 + 0.2299 + 0.0412]/3 = 0.6489/3 = 0.2163

   - For second group:
     * s₂² = Σ(x₂ᵢ - x̄₂)²/(n₂-1)
     * = [(0.3035-1.6240)² + (2.6961-1.6240)² + (1.0591-1.6240)² + (2.7971-1.6240)² + (1.2641-1.6240)²]/4
     * = [(-1.3205)² + (1.0721)² + (-0.5649)² + (1.1731)² + (-0.3599)²]/4
     * = [1.7437 + 1.1494 + 0.3191 + 1.3762 + 0.1295]/4 = 4.7179/4 = 1.1795

5) Pooled variance (assuming equal population variances):
   - s²ₚ = [(n₁-1)s₁² + (n₂-1)s₂²]/(n₁+n₂-2)
   - = [(4-1)×0.2163 + (5-1)×1.1795]/(4+5-2)
   - = [0.6489 + 4.718]/7 = 5.3669/7 = 0.7667

**Statistical Connection:** The pooled variance estimator exemplifies the principle of efficiency in statistical estimation. By assuming homogeneity of variances, we combine information from both samples to get a more precise estimate with more degrees of freedom (n₁+n₂-2 instead of separate estimates with n₁-1 and n₂-1 degrees of freedom).

### Part (b): Standard error of difference of means

**Step-by-step calculation:**
- The standard error quantifies uncertainty in the difference of means
- For independent samples with equal variances:
  * SE(x̄₂ - x̄₁) = √[s²ₚ × (1/n₁ + 1/n₂)]
  * = √[0.7667 × (1/4 + 1/5)]
  * = √[0.7667 × (0.25 + 0.20)]
  * = √[0.7667 × 0.45]
  * = √0.345 = 0.5874

**Statistical Connection:** This formula demonstrates how sampling errors propagate when taking the difference of two independently estimated quantities. The expression (1/n₁ + 1/n₂) shows how both sample sizes affect precision, highlighting the statistical principle that larger samples provide more accurate estimates.

### Part (c): 90% confidence interval

**Step-by-step calculation:**
1) For a 90% confidence interval with unknown variance:
   - We use the t-distribution with degrees of freedom = n₁ + n₂ - 2 = 7
   - For a 90% interval, we need t₀.₉₅₍₇₎ (for a two-sided interval)
   - From t-tables or software: t₀.₉₅₍₇₎ = 1.8946

2) The confidence interval formula:
   - (x̄₂ - x̄₁) ± t₀.₉₅₍₇₎ × SE
   - = 1.069355 ± 1.8946 × 0.5874
   - = 1.069355 ± 1.112863
   - = (-0.0435, 2.1822) ≈ (-0.03, 2.17)

**Statistical Connection:** This employs Student's t-distribution (developed by William Gosset) which accounts for the additional uncertainty introduced by estimating the variance. Unlike the normal distribution, the t-distribution has heavier tails, resulting in wider confidence intervals especially for small samples. The interval contains zero, suggesting we lack strong evidence of a difference between means.

### Part (d): P-value for test of equal means

**Step-by-step calculation:**
1) The null hypothesis is H₀: μ₁ = μ₂ (or μ₂ - μ₁ = 0)
   The alternative is H₁: μ₁ ≠ μ₂ (two-sided)

2) Test statistic:
   - t = (x̄₂ - x̄₁)/SE = 1.069355/0.5874 = 1.820

3) P-value calculation for a two-tailed test:
   - p-value = 2 × P(t₇ > 1.820)
   - Using t-distribution with df = 7
   - p-value ≈ 0.112

4) Since p-value > 0.05 (conventional significance level), we fail to reject H₀.

**Statistical Connection:** This exemplifies the duality between hypothesis tests and confidence intervals - the 90% confidence interval contains zero, and correspondingly, the p-value exceeds 0.10. The result illustrates the challenge of drawing conclusions with small sample sizes - despite an estimated mean difference of 1.07, we lack sufficient evidence to declare it statistically significant.

### Part (e): Analysis with known variance of 1

**Step-by-step calculation:**
1) If the variance is known to be σ² = 1 (rather than estimated):
   - Standard error changes to:
   - SE = √(σ² × (1/n₁ + 1/n₂)) = √(1 × (1/4 + 1/5)) = √0.45 = 0.6708

2) With known variance, we use the standard normal distribution (Z) instead of t:
   - Z = (x̄₂ - x̄₁)/SE = 1.069355/0.6708 = 1.593

3) P-value:
   - p-value = 2 × P(Z > 1.593) ≈ 0.111

4) 90% confidence interval:
   - z₀.₉₅ = 1.6449
   - 90% CI = 1.069355 ± 1.6449 × 0.6708 = 1.069355 ± 1.103419 = (-0.034, 2.173)

**Statistical Connection:** This demonstrates how prior information changes our statistical approach - from a t-test to a z-test. Interestingly, despite using a z-test (which is typically more powerful), the standard error increases from 0.5874 to 0.6708 because the assumed variance (1) is larger than the estimated pooled variance (0.767). This illustrates how assumptions affect statistical results, sometimes in counter-intuitive ways.
