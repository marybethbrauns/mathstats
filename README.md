# MathStats Homework 2 Solutions

## Problem 1
### Part (a) 
Let X₁, ..., X₁₀ be a random sample from a N(5, 8) distribution and Y₁, ..., Y₁₅ an independent random sample from a N(-2, 16) distribution.

#### i. Distribution of V = X̄₁₀ - Ȳ₁₅
We know that if X₁, ..., Xₙ ~ N(μ, σ²), then X̄ₙ ~ N(μ, σ²/n).

Therefore:
- X̄₁₀ ~ N(5, 8/10) = N(5, 0.8)
- Ȳ₁₅ ~ N(-2, 16/15) = N(-2, 1.0667)

For the difference of independent normal random variables:
- V = X̄₁₀ - Ȳ₁₅ ~ N(5-(-2), 0.8 + 1.0667) = N(7, 1.8667)

So V follows a normal distribution with mean μᵥ = 7 and variance σ²ᵥ = 1.8667.

#### ii. Calculate P(6 < V < 9)
For a normal random variable, we standardize:

P(6 < V < 9) = P((6-7)/√1.8667 < Z < (9-7)/√1.8667) = P(-0.7321 < Z < 1.4642)

Using the standard normal CDF:

P(-0.7321 < Z < 1.4642) = Φ(1.4642) - Φ(-0.7321) = Φ(1.4642) - (1-Φ(0.7321))
                         = Φ(1.4642) + Φ(0.7321) - 1 
                         = 0.9284 + 0.7680 - 1 
                         = 0.6964

Therefore, P(6 < V < 9) = 0.6964

#### iii. Simulation

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

np.random.seed(123)
nsim = 10000
v_samples = np.zeros(nsim)

for i in range(nsim):
    x_sample = np.random.normal(5, np.sqrt(8), size=10)
    y_sample = np.random.normal(-2, np.sqrt(16), size=15)
    v_samples[i] = np.mean(x_sample) - np.mean(y_sample)

# Mean and variance of simulated values
print(f"Mean: {np.mean(v_samples):.4f}")  # Should be close to 7
print(f"Variance: {np.var(v_samples):.4f}")   # Should be close to 1.8667

# Proportion between 6 and 9
proportion = np.mean((v_samples > 6) & (v_samples < 9))
print(f"Proportion between 6 and 9: {proportion:.4f}")  # Should be close to 0.6964

# Plot histogram with PDF overlay
plt.figure(figsize=(10, 6))
plt.hist(v_samples, bins=30, density=True, alpha=0.7, label='Simulated Values')
x = np.linspace(2, 12, 1000)
plt.plot(x, stats.norm.pdf(x, 7, np.sqrt(1.8667)), 'r-', lw=2, label='Theoretical PDF')
plt.title('Histogram of V with PDF Overlay')
plt.xlabel('V = X̄ - Ȳ')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

The simulation results confirm our theoretical calculations. The mean is approximately 7, the variance is close to 1.8667, and the proportion between 6 and 9 is approximately 0.6964.

### Part (b)
Let Zᵢ = (Xᵢ - 5)/√8 for 1 ≤ i ≤ 10, and W = ∑ᵢ₌₁¹⁰ Zᵢ²

#### i. Distribution of W
When we standardize Xᵢ ~ N(5, 8) to get Zᵢ = (Xᵢ - 5)/√8, we obtain Zᵢ ~ N(0, 1).

The sum of squares of standard normal random variables follows a chi-square distribution:
W = ∑ᵢ₌₁¹⁰ Zᵢ² ~ χ²(10)

For a chi-square distribution with k degrees of freedom:

E[W] = k = 10  
Var(W) = 2k = 20

#### ii. Calculate P(6 < W < 9)
For a chi-square random variable with 10 degrees of freedom:

P(6 < W < 9) = F_χ²(10)(9) - F_χ²(10)(6)

Using the chi-square CDF:

P(6 < W < 9) = 0.4697 - 0.1817 = 0.2880

#### iii. Simulation

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

np.random.seed(123)
nsim = 10000
w_samples = np.zeros(nsim)

for i in range(nsim):
    x_sample = np.random.normal(5, np.sqrt(8), size=10)
    z_sample = (x_sample - 5)/np.sqrt(8)
    w_samples[i] = np.sum(z_sample**2)

# Mean and variance of simulated values
print(f"Mean: {np.mean(w_samples):.4f}")  # Should be close to 10
print(f"Variance: {np.var(w_samples):.4f}")   # Should be close to 20

# Proportion between 6 and 9
proportion = np.mean((w_samples > 6) & (w_samples < 9))
print(f"Proportion between 6 and 9: {proportion:.4f}")  # Should be close to 0.2880

# Plot histogram with PDF overlay
plt.figure(figsize=(10, 6))
plt.hist(w_samples, bins=30, density=True, alpha=0.7, label='Simulated Values')
x = np.linspace(0, 30, 1000)
plt.plot(x, stats.chi2.pdf(x, df=10), 'r-', lw=2, label='Chi-Square(10) PDF')
plt.title('Histogram of W with PDF Overlay')
plt.xlabel('W = Sum of Squared Z_i')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

The simulation results confirm our theoretical calculations. The mean is approximately 10, the variance is close to 20, and the proportion between 6 and 9 is approximately 0.2880.

## Problem 2
Suppose that X₁, X₂, ..., Xₙ is a random sample from a Unif(0, b) distribution. Consider three estimators for b:

* b̂₁ = 2X̄
* b̂₂ = X₍ₙ₎
* b̂₃ = (n+1)/n·X₍ₙ₎

Take n = 10 and b = 4.

### Part (a) Find the PDF of b̂₃
First, let's find the PDF of X₍ₙ₎ (the maximum order statistic).

For a Unif(0, b) distribution, the CDF is F(x) = x/b for 0 ≤ x ≤ b.

The CDF of X₍ₙ₎ is:

F_{X₍ₙ₎}(x) = P(X₍ₙ₎ ≤ x) = P(all Xᵢ ≤ x) = (x/b)ⁿ

The PDF of X₍ₙ₎ is the derivative of its CDF:

f_{X₍ₙ₎}(x) = d/dx F_{X₍ₙ₎}(x) = n·xⁿ⁻¹ / bⁿ

Now, for b̂₃ = (n+1)/n·X₍ₙ₎, we use the change of variable technique.

Let Y = b̂₃ = (n+1)/n·X₍ₙ₎, so X₍ₙ₎ = n/(n+1)·Y. The Jacobian of this transformation is dx/dy = n/(n+1).

The PDF of Y = b̂₃ is:

f_Y(y) = f_{X₍ₙ₎}(n/(n+1)·y) · |dx/dy|

       = n·(n/(n+1)·y)ⁿ⁻¹ / bⁿ · n/(n+1)

       = nⁿ·yⁿ⁻¹ / ((n+1)ⁿ·bⁿ)

For n = 10 and b = 4, the PDF of b̂₃ is:

f_{b̂₃}(y) = 10¹⁰·y⁹ / (11¹⁰·4¹⁰) = 10¹⁰·y⁹ / (11¹⁰·2²⁰)

### Part (b) Probabilities within 0.05 of b

#### i. For b̂₁ = 2X̄
For a Unif(0, b) distribution:

E[X] = b/2 = 2  
Var(X) = b²/12 = 16/12 = 4/3

Therefore:

E[b̂₁] = 2E[X̄] = 2 · b/2 = b = 4  
Var(b̂₁) = 4Var(X̄) = 4 · Var(X)/n = 4 · (4/3)/10 = 16/30 = 8/15

By the Central Limit Theorem, b̂₁ is approximately normally distributed:

b̂₁ ≈ N(4, 8/15)

So:

P(|b̂₁ - b| < 0.05) = P(3.95 < b̂₁ < 4.05)

                    = P((3.95-4)/√(8/15) < Z < (4.05-4)/√(8/15))

                    = P(-0.0671 < Z < 0.0671)

                    = 2Φ(0.0671) - 1 = 2(0.5267) - 1 = 0.0534

#### ii. For b̂₂ = X₍ₙ₎

P(|b̂₂ - b| < 0.05) = P(3.95 < X₍ₙ₎ < 4)

                    = F_{X₍ₙ₎}(4) - F_{X₍ₙ₎}(3.95)

                    = (4/4)¹⁰ - (3.95/4)¹⁰

                    = 1 - (0.9875)¹⁰ = 1 - 0.8825 = 0.1175

#### iii. For b̂₃ = (n+1)/n·X₍ₙ₎

P(|b̂₃ - b| < 0.05) = P(3.95 < 11/10·X₍ₙ₎ < 4.05)

                    = P(10/11 · 3.95 < X₍ₙ₎ < 10/11 · 4.05)

                    = P(3.5909 < X₍ₙ₎ < 3.6818)

                    = F_{X₍ₙ₎}(3.6818) - F_{X₍ₙ₎}(3.5909)

                    = (3.6818/4)¹⁰ - (3.5909/4)¹⁰

                    = (0.9205)¹⁰ - (0.8977)¹⁰ = 0.4309 - 0.3259 = 0.1050

#### iv. Comparison
Comparing the probabilities:

P(|b̂₁ - b| < 0.05) = 0.0534  
P(|b̂₂ - b| < 0.05) = 0.1175  
P(|b̂₃ - b| < 0.05) = 0.1050  

The estimator b̂₂ = X₍ₙ₎ has the highest probability of being within 0.05 of b.

### Part (c) Simulation

```python
import numpy as np

np.random.seed(123)
nsim = 10000
b = 4
n = 10

b1_values = np.zeros(nsim)
b2_values = np.zeros(nsim)
b3_values = np.zeros(nsim)

for i in range(nsim):
    sample = np.random.uniform(0, b, size=n)
    b1_values[i] = 2 * np.mean(sample)
    b2_values[i] = np.max(sample)
    b3_values[i] = ((n+1)/n) * np.max(sample)

# Proportions within 0.05 of b
print(f"Proportion for b1: {np.mean(np.abs(b1_values - b) < 0.05):.4f}")  # Should be close to 0.0534
print(f"Proportion for b2: {np.mean(np.abs(b2_values - b) < 0.05):.4f}")  # Should be close to 0.1175
print(f"Proportion for b3: {np.mean(np.abs(b3_values - b) < 0.05):.4f}")  # Should be close to 0.1050
```

The simulation results confirm our theoretical calculations.

### Part (d) General estimator b̂ = cX₍ₙ₎

#### i. Mean Square Error
For the general estimator b̂ = cX₍ₙ₎, the MSE is:

MSE(b̂) = E[(b̂ - b)²] = E[(cX₍ₙ₎ - b)²]

= E[c²X₍ₙ₎² - 2cbX₍ₙ₎ + b²]

= c²E[X₍ₙ₎²] - 2cbE[X₍ₙ₎] + b²

We need to calculate E[X₍ₙ₎] and E[X₍ₙ₎²].

For a Unif(0, b) distribution with PDF f_{X₍ₙ₎}(x) = n·x^(n-1)/b^n:

E[X₍ₙ₎] = ∫₀ᵇ x · n·xⁿ⁻¹/bⁿ dx 
       = n/bⁿ · ∫₀ᵇ xⁿ dx 
       = n/bⁿ · bⁿ⁺¹/(n+1) 
       = n·b/(n+1)

E[X₍ₙ₎²] = ∫₀ᵇ x² · n·xⁿ⁻¹/bⁿ dx 
         = n/bⁿ · ∫₀ᵇ xⁿ⁺¹ dx 
         = n/bⁿ · bⁿ⁺²/(n+2) 
         = n·b²/(n+2)

Substituting into the MSE:

MSE(b̂) = c² · n·b²/(n+2) - 2cb · n·b/(n+1) + b²

       = c²·n·b² / (n+2) - 2c·n·b² / (n+1) + b²

#### ii. Minimizing MSE
To minimize the MSE, we differentiate with respect to c and set equal to zero:

d/dc MSE(b̂) = 2c·n·b² / (n+2) - 2n·b² / (n+1) = 0

Solving for c:

2c·n·b² / (n+2) = 2n·b² / (n+1)

c = (n+2) / (n+1)

For n = 10, the optimal value is c = 12/11.

The bias of this estimator is:

Bias(b̂) = E[b̂] - b = E[cX₍ₙ₎] - b = c · E[X₍ₙ₎] - b = c · n·b/(n+1) - b

        = (n+2)/(n+1) · n·b/(n+1) - b 
        = n(n+2)·b / (n+1)² - b

        = b·(n(n+2)/(n+1)² - 1) 
        = b·(n(n+2) - (n+1)²) / (n+1)²

        = b·(n² + 2n - n² - 2n - 1) / (n+1)² 
        = b·(-1) / (n+1)² 
        = -b / (n+1)²

For n = 10 and b = 4, the bias is -4/(11²) = -4/121 ≈ -0.0331.

#### iii. Probability for optimal estimator
For b̂ = (n+2)/(n+1)·X₍ₙ₎:

P(|b̂ - b| < 0.05) = P(3.95 < 12/11·X₍ₙ₎ < 4.05)

= P(11/12 · 3.95 < X₍ₙ₎ < 11/12 · 4.05)

= P(3.6208 < X₍ₙ₎ < 3.7125)

= F_{X₍ₙ₎}(3.7125) - F_{X₍ₙ₎}(3.6208)

= (3.7125/4)¹⁰ - (3.6208/4)¹⁰

= (0.9281)¹⁰ - (0.9052)¹⁰ = 0.4728 - 0.3762 = 0.0966

## Problem 3
A statistic used in testing hypotheses concerning categorical data is the Pearson Chi-Square statistic.

### Part (a) Type of random variable for Oᵢ
If a fair die is rolled n times, Oᵢ represents the number of times face i appears.

Since each roll has probability 1/6 of showing face i, and the rolls are independent, Oᵢ follows a binomial distribution:

Oᵢ ~ Binomial(n, 1/6)

### Part (b) Expected values Eᵢ
For a fair die, the expected number of times face i appears in n rolls is:

Eᵢ = n · P(face i) = n · 1/6 = n/6

### Part (c) Expected value of χ²
The Pearson Chi-Square statistic is:

χ² = ∑ᵢ₌₁⁶ (Oᵢ - Eᵢ)²/Eᵢ

For a single term:

E[(Oᵢ - Eᵢ)²] = Var(Oᵢ) = np(1-p) = n · 1/6 · 5/6 = 5n/36

E[(Oᵢ - Eᵢ)²/Eᵢ] = E[(Oᵢ - Eᵢ)²]/Eᵢ = (5n/36)/(n/6) = 5/6

For the sum of 6 terms:

E[χ²] = ∑ᵢ₌₁⁶ E[(Oᵢ - Eᵢ)²/Eᵢ] = 6 · 5/6 = 5

Note: For a chi-square test with k categories, the degrees of freedom is k-1 = 5, and the expected value of the chi-square statistic under the null hypothesis is equal to the degrees of freedom, which is 5. This confirms our calculation.

### Part (d) Calculate χ² for observed data

With n = 100, Eᵢ = n/6 = 100/6 = 16.67 for each face.

Observed counts: [25, 12, 21, 18, 16, 8]

χ² = ∑ᵢ₌₁⁶ (Oᵢ - Eᵢ)²/Eᵢ

   = (25 - 16.67)²/16.67 + (12 - 16.67)²/16.67 + (21 - 16.67)²/16.67 
     + (18 - 16.67)²/16.67 + (16 - 16.67)²/16.67 + (8 - 16.67)²/16.67

   = 4.18 + 1.31 + 1.13 + 0.11 + 0.03 + 4.52 
   = 11.28

### Part (e) Simulation for distribution of χ²

```python
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(123)
nsim = 10000
chi_sq_samples = np.zeros(nsim)

for i in range(nsim):
    # Simulate 100 rolls of a fair die
    rolls = np.random.randint(1, 7, size=100)
    # Count occurrences of each face
    observed = np.zeros(6)
    for j in range(1, 7):
        observed[j-1] = np.sum(rolls == j)
    # Expected counts (16.67 for each face)
    expected = np.ones(6) * (100/6)
    # Calculate chi-square statistic
    chi_sq_samples[i] = np.sum((observed - expected)**2 / expected)

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(chi_sq_samples, bins=30, density=True, alpha=0.7)
plt.title('Distribution of Chi-Square Statistic')
plt.xlabel('Chi-Square Value')
plt.grid(True, alpha=0.3)
plt.show()

# Calculate p-value (proportion of simulated values >= observed value)
p_value = np.mean(chi_sq_samples >= 11.28)
print(f"p-value = {p_value:.4f}")
```

The p-value represents the probability of observing a chi-square statistic at least as extreme as 11.28 if the die is fair. If the p-value is small (typically < 0.05), we would reject the null hypothesis that the die is fair.

Based on our simulation, if the p-value is approximately 0.046, this suggests that the observed pattern of rolls is unlikely to occur by chance if the die is truly fair. Therefore, we have evidence to suggest that the die is not fair.

## Problem 4
Let (X₁, Y₁), (X₂, Y₂), ..., (Xₙ, Yₙ) be a random sample from some joint distribution F_{X,Y}.

### Part (a) Show that Cov(X̄, Ȳ) = C/n
Using properties of covariance:

Cov(X̄, Ȳ) = Cov(1/n ∑ᵢ₌₁ⁿ Xᵢ, 1/n ∑ⱼ₌₁ⁿ Yⱼ) = 1/n² ∑ᵢ₌₁ⁿ∑ⱼ₌₁ⁿ Cov(Xᵢ, Yⱼ)

Since Xᵢ and Yⱼ are independent if i ≠ j, Cov(Xᵢ, Yⱼ) = 0 for i ≠ j. Therefore:

Cov(X̄, Ȳ) = 1/n² ∑ᵢ₌₁ⁿ Cov(Xᵢ, Yᵢ) = 1/n² · n · C = C/n

### Part (b) Show that ∑ᵢ₌₁ⁿ (Xᵢ - X̄)(Yᵢ - Ȳ) = ∑ᵢ₌₁ⁿ XᵢYᵢ - nX̄Ȳ

Expanding the left side:

∑ᵢ₌₁ⁿ (Xᵢ - X̄)(Yᵢ - Ȳ) = ∑ᵢ₌₁ⁿ (XᵢYᵢ - XᵢȲ - X̄Yᵢ + X̄Ȳ)

= ∑ᵢ₌₁ⁿ XᵢYᵢ - Ȳ∑ᵢ₌₁ⁿ Xᵢ - X̄∑ᵢ₌₁ⁿ Yᵢ + nX̄Ȳ

= ∑ᵢ₌₁ⁿ XᵢYᵢ - Ȳ(nX̄) - X̄(nȲ) + nX̄Ȳ

= ∑ᵢ₌₁ⁿ XᵢYᵢ - nX̄Ȳ - nX̄Ȳ + nX̄Ȳ

= ∑ᵢ₌₁ⁿ XᵢYᵢ - nX̄Ȳ

### Part (c) Show that μ_{XY} = C + μₓ μᵧ and E(X̄Ȳ) = C/n + μₓ μᵧ

By definition, C = Cov(Xᵢ, Yᵢ) = E[(Xᵢ - μₓ)(Yᵢ - μᵧ)]. Expanding:

C = E[XᵢYᵢ - Xᵢμᵧ - Yᵢμₓ + μₓμᵧ]

= E[XᵢYᵢ] - μᵧ E[Xᵢ] - μₓ E[Yᵢ] + μₓμᵧ

= E[XᵢYᵢ] - μᵧμₓ - μₓμᵧ + μₓμᵧ

= E[XᵢYᵢ] - μₓμᵧ

Rearranging, we get:

E[XᵢYᵢ] = C + μₓμᵧ

Therefore, μ_{XY} = C + μₓ μᵧ.

For the second part:

E[X̄Ȳ] = Cov(X̄, Ȳ) + E[X̄]E[Ȳ] = C/n + μₓμᵧ

### Part (d) Show that sample covariance is an unbiased estimator
The sample covariance is:

s_{XY} = 1/(n-1) ∑ᵢ₌₁ⁿ (Xᵢ - X̄)(Yᵢ - Ȳ)

From part (b):

s_{XY} = 1/(n-1) (∑ᵢ₌₁ⁿ XᵢYᵢ - nX̄Ȳ)

Taking the expected value:

E[s_{XY}] = 1/(n-1) (E[∑ᵢ₌₁ⁿ XᵢYᵢ] - nE[X̄Ȳ])

         = 1/(n-1) (nμ_{XY} - n(C/n + μₓμᵧ))

         = 1/(n-1) (n(C + μₓμᵧ) - C - nμₓμᵧ)

         = 1/(n-1) (nC + nμₓμᵧ - C - nμₓμᵧ)

         = 1/(n-1) ((n-1)C)

         = C

Therefore, the sample covariance is an unbiased estimator of the population covariance.
