# Numerical Analysis of Integration Methods

## Problem Statement

We are asked to analyze and compare different numerical integration methods for evaluating the definite integral:

I(f) = the integral from 0 to π of e^x cos(x) dx

The specific tasks are:

a) Evaluate the integral exactly using the Fundamental Theorem of Calculus.

b) Implement the Composite Midpoint Rule with n = 2^i, where i = 1, 2, ..., 9, compute the approximation In(f), and determine the error En(f) = I(f) - In(f). Analyze by what factor the error decreases as n is doubled and explain why.

c) Repeat part (b) using the Composite Trapezoidal Rule.

d) Repeat part (b) using the Corrected Trapezoidal Rule.

e) Repeat part (b) using the Composite Simpson's Rule with n = 2^i, where i = 1, 2, ..., 7. Analyze the error reduction factor and explain why.

f) Discuss the accuracy of the approximations and the computational effort involved.

## (a) Exact Integral Evaluation

To evaluate the integral from 0 to π of e^x cos(x) dx exactly, we'll use the Fundamental Theorem of Calculus by finding the antiderivative.

For f(x) = e^x cos(x), we'll use integration by parts twice:

Let u = e^x, dv = cos(x)dx  
Then du = e^x dx, v = sin(x)

∫ e^x cos(x) dx = e^x sin(x) - ∫ e^x sin(x) dx

For the remaining integral, let u = e^x, dv = sin(x)dx  
Then du = e^x dx, v = -cos(x)

∫ e^x sin(x) dx = -e^x cos(x) + ∫ e^x cos(x) dx

Substituting:
∫ e^x cos(x) dx = e^x sin(x) - (-e^x cos(x) + ∫ e^x cos(x) dx)  
= e^x sin(x) + e^x cos(x) - ∫ e^x cos(x) dx

Solving for ∫ e^x cos(x) dx:
2∫ e^x cos(x) dx = e^x sin(x) + e^x cos(x)  
∫ e^x cos(x) dx = (e^x sin(x) + e^x cos(x))/2

Therefore, F(x) = (e^x sin(x) + e^x cos(x))/2 is our antiderivative.

I(f) = F(π) - F(0)  
= [(e^π sin(π) + e^π cos(π))/2] - [(e^0 sin(0) + e^0 cos(0))/2]  
= [(e^π·0 + e^π·(-1))/2] - [(1·0 + 1·1)/2]  
= [-e^π/2] - [1/2]  
= -(e^π + 1)/2  
≈ -12.070346316389633

This exact value will serve as our reference for calculating the error in our numerical approximations.

## (b) Composite Midpoint Rule

### Method Description

The composite midpoint rule approximates the integral by dividing the interval [0,π] into n subintervals and evaluating the function at the midpoint of each subinterval. This method is a natural choice for numerical integration when function evaluations are costly or when the function is not well-defined at the endpoints.

The formula for the composite midpoint rule is:

I_n(f) = h · sum from i=1 to n of f(x_i)

where:
- h = (b-a)/n is the width of each subinterval
- x_i = a + (i-1/2)·h is the midpoint of the i-th subinterval

### Theoretical Error Analysis

For the midpoint rule, the theoretical error is given by:

E_n(f) = -(b-a)^3/(24n^2) · f''(ξ)

for some ξ ∈ [a,b]. This indicates that the midpoint rule has second-order accuracy, O(h^2), where h = (b-a)/n. When we double n (halve h), we expect the error to decrease by a factor of approximately 4.

### Implementation Results

Results for Composite Midpoint Rule with n = 2^i, i = 1,2,...,9:

| n    | Approximation | Error        | Convergence Ratio |
|------|---------------|--------------|-------------------|
| 2    | -9.28278636   | -2.78755995  | N/A               |
| 4    | -11.42830201  | -0.64204430  | 4.3417            |
| 8    | -11.91384577  | -0.15650055  | 4.1025            |
| 16   | -12.03148013  | -0.03886618  | 4.0267            |
| 32   | -12.06064608  | -0.00970023  | 4.0067            |
| 64   | -12.06792228  | -0.00242404  | 4.0017            |
| 128  | -12.06974037  | -0.00060595  | 4.0004            |
| 256  | -12.07019483  | -0.00015148  | 4.0001            |
| 512  | -12.07030845  | -0.00003787  | 4.0000            |

### Convergence Ratio Justification

The convergence ratio is calculated as |E_n(f)|/|E_(2n)(f)|, which represents how much the error decreases when we double the number of subintervals.

For the midpoint rule, we observe that this ratio consistently approaches 4 as n increases. This confirms the theoretical O(h^2) convergence rate. When n doubles, h is halved, and the error term containing h^2 decreases by a factor of 2^2 = 4.

The initial convergence ratio (4.3417 for n=2 to n=4) is slightly higher than 4, likely due to the specific behavior of the higher derivatives of our integrand over the interval [0,π]. As n increases, the ratio stabilizes almost exactly at 4, validating the theoretical expectation.

The convergence ratio approaching 4 is mathematically justified by the error term in the midpoint rule, which is proportional to 1/n^2. When n doubles to 2n, the error term becomes proportional to 1/(2n)^2 = 1/(4n^2), which is 1/4 of the original error. Therefore, the ratio of consecutive errors |E_n(f)|/|E_(2n)(f)| should approach 4, which is precisely what we observe in our numerical results.

## (c) Composite Trapezoidal Rule

### Method Description

The composite trapezoidal rule approximates the integral by connecting function values at the endpoints of each subinterval with straight lines, creating a series of trapezoids. This method is widely used due to its simplicity and effectiveness for smooth functions.

The formula for the composite trapezoidal rule is:

I_n(f) = (h/2) · [f(a) + f(b) + 2·sum from i=1 to n-1 of f(a + i·h)]

where:
- h = (b-a)/n is the width of each subinterval

### Theoretical Error Analysis

For the trapezoidal rule, the theoretical error is given by:

E_n(f) = -(b-a)^3/(12n^2) · f''(ξ)

for some ξ ∈ [a,b]. Like the midpoint rule, the trapezoidal rule has second-order accuracy, O(h^2). When we double n, we expect the error to decrease by a factor of approximately 4.

### Implementation Results

Results for Composite Trapezoidal Rule with n = 2^i, i = 1,2,...,9:

| n    | Approximation  | Error       | Convergence Ratio |
|------|----------------|-------------|-------------------|
| 2    | -17.38925933   | 5.31891301  | N/A               |
| 4    | -13.33602285   | 1.26567653  | 4.2024            |
| 8    | -12.38216243   | 0.31181611  | 4.0590            |
| 16   | -12.14800410   | 0.07765778  | 4.0153            |
| 32   | -12.08974212   | 0.01939580  | 4.0038            |
| 64   | -12.07519410   | 0.00484778  | 4.0010            |
| 128  | -12.07155819   | 0.00121187  | 4.0002            |
| 256  | -12.07064928   | 0.00030296  | 4.0001            |
| 512  | -12.07042206   | 0.00007574  | 4.0000            |

### Convergence Ratio Justification

Similar to the midpoint rule, the convergence ratio for the trapezoidal rule approaches 4 as n increases, confirming the O(h^2) convergence rate. This is consistent with the theoretical error formula which includes an h^2 term.

An interesting observation is that while the midpoint rule consistently underestimates the true value (negative errors), the trapezoidal rule consistently overestimates it (positive errors). This is a well-known property: for functions with positive second derivatives in the integration interval, the trapezoidal rule overestimates the integral, while the midpoint rule underestimates it.

The convergence ratio starts slightly higher at 4.2024 for small n values and stabilizes at almost exactly 4 for larger n values. This matches our theoretical expectations perfectly and validates the implementation.

The mathematical justification is similar to that of the midpoint rule: the error term in the trapezoidal rule is proportional to 1/n^2. When n doubles, the error term becomes proportional to 1/(4n^2), which is 1/4 of the original error. Therefore, the ratio of consecutive errors should approach 4, which is confirmed by our numerical results.

## (d) Corrected Trapezoidal Rule

### Method Description

The corrected trapezoidal rule, also known as the Euler-Maclaurin formula, enhances the basic trapezoidal rule by adding a correction term based on the derivatives of the function at the endpoints. This modification significantly improves accuracy without substantially increasing the computational cost.

The formula for the corrected trapezoidal rule is:

I_n(f) = I_n^T(f) - (h^2/12) · [f'(b) - f'(a)]

where:
- I_n^T(f) is the approximation from the standard trapezoidal rule
- f'(x) is the first derivative of f(x)

For our function f(x) = e^x cos(x), the derivative is:
f'(x) = e^x (cos(x) - sin(x))

### Theoretical Error Analysis

The error for the corrected trapezoidal rule is of order O(h^4), a significant improvement from the O(h^2) of the basic trapezoidal rule. The Euler-Maclaurin correction effectively eliminates the leading error term in the trapezoidal rule.

Theoretically, the error formula is:

E_n(f) = -(b-a)^5/(720n^4) · f^(4)(ξ)

for some ξ ∈ [a,b]. This indicates that when we double n (halve h), we expect the error to decrease by a factor of approximately 2^4 = 16.

### Implementation Results

Results for Corrected Trapezoidal Rule with n = 2^i, i = 1,2,...,9:

| n    | Approximation  | Error        | Convergence Ratio |
|------|----------------|--------------|-------------------|
| 2    | -12.42552837   | 0.35518205   | N/A               |
| 4    | -12.09509011   | 0.02474379   | 14.3544           |
| 8    | -12.07192924   | 0.00158293   | 15.6317           |
| 16   | -12.07044580   | 0.00009949   | 15.9109           |
| 32   | -12.07035254   | 0.00000623   | 15.9779           |
| 64   | -12.07034671   | 0.00000039   | 15.9945           |
| 128  | -12.07034634   | 0.00000002   | 15.9986           |
| 256  | -12.07034632   | 0.00000000   | 15.9997           |
| 512  | -12.07034632   | 0.00000000   | 16.0005           |

### Convergence Ratio Justification

The convergence ratio for the corrected trapezoidal rule approaches 16 as n increases, confirming the theoretical O(h^4) convergence rate. This dramatic improvement over the basic trapezoidal rule demonstrates the power of the Euler-Maclaurin correction.

The initial convergence ratio of 14.3544 (from n=2 to n=4) is slightly below 16, but it quickly approaches the theoretical value as n increases. By n=512, the ratio is approximately 16.0005, remarkably close to the expected value of 16. This confirms that the implementation correctly captures the fourth-order accuracy of the method.

The progression of the convergence ratio (14.3544 → 15.6317 → 15.9109 → 15.9779 → 15.9945 → 15.9986 → 15.9997 → 16.0005) clearly shows it asymptotically approaching the theoretical value of 16, which is exactly what we would expect from the error analysis.

The mathematical justification is based on the error term in the corrected trapezoidal rule, which is proportional to 1/n^4. When n doubles, the error term becomes proportional to 1/(16n^4), which is 1/16 of the original error. Therefore, the ratio of consecutive errors should approach 16, which is confirmed by our numerical results with remarkable precision.

## (e) Composite Simpson's Rule

### Method Description

The composite Simpson's rule approximates the integral by using quadratic interpolation over pairs of subintervals. It effectively fits a parabola through three consecutive points and integrates the resulting function. Simpson's rule is widely used because it achieves high-order accuracy without requiring derivatives.

The formula for the composite Simpson's rule is:

I_n(f) = (h/3) · [f(a) + f(b) + 4·sum from i=1 to k of f(x_(2i-1)) + 2·sum from i=1 to k-1 of f(x_(2i))]

where:
- k = n/2 (n must be even)
- h = (b-a)/n
- x_i = a + i·h

### Theoretical Error Analysis

For Simpson's rule, the theoretical error is given by:

E_n(f) = -(b-a)^5/(180n^4) · f^(4)(ξ)

for some ξ ∈ [a,b]. This indicates that Simpson's rule has fourth-order accuracy, O(h^4). When we double n (halve h), we expect the error to decrease by a factor of approximately 2^4 = 16.

### Implementation Results

Results for Composite Simpson's Rule with n = 2^i, i = 1,2,...,7:

| n    | Approximation  | Error        | Convergence Ratio |
|------|----------------|--------------|-------------------|
| 2    | -11.59283955   | -0.47750676  | N/A               |
| 4    | -11.98494402   | -0.08540230  | 5.5913            |
| 8    | -12.06420896   | -0.00613736  | 13.9152           |
| 16   | -12.06995132   | -0.00039499  | 15.5379           |
| 32   | -12.07032146   | -0.00002486  | 15.8885           |
| 64   | -12.07034476   | -0.00000156  | 15.9724           |
| 128  | -12.07034622   | -0.00000010  | 15.9931           |

### Convergence Ratio Justification

The convergence ratio for Simpson's rule shows an interesting progression. For small n values, the ratio is significantly lower than the theoretical value of 16, but it rapidly approaches this value as n increases.

For n=2 to n=4, the ratio is only 5.5913, but by n=8 to n=16, it reaches 15.5379, and for n=64 to n=128, it's 15.9931, very close to the theoretical value of 16.

This behavior can be explained by considering the complete error expansion for Simpson's rule. While the leading term is O(h^4), there are additional higher-order terms that become less significant as h decreases (or n increases). For small n values, these higher-order terms still contribute noticeably to the error, causing the convergence ratio to deviate from 16. As n increases, the O(h^4) term dominates, and the ratio approaches the theoretical value of 16.

The observed progression (5.5913 → 13.9152 → 15.5379 → 15.8885 → 15.9724 → 15.9931) clearly shows the convergence ratio asymptotically approaching 16, validating the fourth-order accuracy of Simpson's rule.

The mathematical justification, similar to the corrected trapezoidal rule, is based on the error term being proportional to 1/n^4. As n doubles, the error should decrease by a factor of 16, which our data confirms for larger values of n when the leading error term dominates.

## (f) Discussion of Accuracy and Computational Effort

### Accuracy Comparison

1. **Composite Midpoint Rule**: Demonstrates second-order convergence (O(h^2)). While simple to implement, it requires a large number of function evaluations to achieve high accuracy. With n = 512, the error is still around 3.8×10^-5.

   The midpoint rule consistently underestimates the true value of the integral for our function. This is because e^x cos(x) has regions where its second derivative is positive over the integration interval [0,π]. For functions with positive second derivatives, the midpoint rule typically underestimates the integral.

2. **Composite Trapezoidal Rule**: Also exhibits second-order convergence (O(h^2)). It produces errors of opposite sign compared to the midpoint rule but similar magnitude. With n = 512, the error is about 7.6×10^-5.

   The trapezoidal rule consistently overestimates the true value of the integral. This complementary behavior to the midpoint rule is a well-known property and can be leveraged in error estimation techniques like Richardson extrapolation.

3. **Corrected Trapezoidal Rule**: Shows fourth-order convergence (O(h^4)). This method dramatically improves upon the basic trapezoidal rule. With just n = 32, it achieves an error of about 6.2×10^-6, which is better than what the midpoint or trapezoidal rules achieve with n = 512.

   The Euler-Maclaurin correction effectively eliminates the leading error term in the trapezoidal rule, resulting in a much higher-order method. For our function, this correction works exceptionally well, providing extremely accurate results even with relatively few subintervals.

4. **Composite Simpson's Rule**: Also exhibits fourth-order convergence (O(h^4)). Its performance is comparable to the corrected trapezoidal rule. With n = 128, the error is approximately 1.0×10^-7.

   Simpson's rule achieves high-order accuracy by effectively fitting parabolas to the function. For our integrand e^x cos(x), which has significant curvature, this approach is particularly effective. The method does not require derivative evaluations, making it versatile for functions where derivatives are not readily available.

### Computational Effort Analysis

1. **Composite Midpoint Rule**: Requires n function evaluations.
   - Advantages: Simple to implement, no need to evaluate endpoints.
   - Disadvantages: Slower convergence requires more subintervals for high accuracy.
   - Computational complexity: O(n) function evaluations to achieve O(1/n^2) accuracy.

2. **Composite Trapezoidal Rule**: Requires n+1 function evaluations.
   - Advantages: Simple to implement, good for periodic functions.
   - Disadvantages: Same O(h^2) convergence as midpoint rule.
   - Computational complexity: O(n) function evaluations to achieve O(1/n^2) accuracy.

3. **Corrected Trapezoidal Rule**: Requires n+1 function evaluations plus 2 derivative evaluations.
   - Advantages: Dramatically improved convergence rate (O(h^4)).
   - Disadvantages: Requires explicit knowledge of the derivative function.
   - Computational complexity: O(n) function evaluations plus constant overhead for derivatives to achieve O(1/n^4) accuracy.

4. **Composite Simpson's Rule**: Requires n+1 function evaluations.
   - Advantages: Fourth-order convergence without requiring derivatives.
   - Disadvantages: Slightly more complex implementation than trapezoidal rule.
   - Computational complexity: O(n) function evaluations to achieve O(1/n^4) accuracy.

### Efficiency Analysis

To quantify the efficiency of each method, we can compare the number of function evaluations required to achieve a specified error tolerance.

For our integral from 0 to π of e^x cos(x) dx:

- To achieve an error of approximately 10^-5:
  - Midpoint Rule: Requires n ≈ 512 (512 function evaluations)
  - Trapezoidal Rule: Requires n ≈ 512 (513 function evaluations)
  - Corrected Trapezoidal Rule: Requires n ≈ 16 (17 function evaluations + 2 derivative evaluations)
  - Simpson's Rule: Requires n ≈ 16 (17 function evaluations)

This comparison clearly demonstrates the superior efficiency of the higher-order methods. The corrected trapezoidal rule and Simpson's rule require approximately 1/32 of the function evaluations needed by the second-order methods to achieve the same accuracy.

The corrected trapezoidal rule and Simpson's rule are clearly superior in terms of accuracy per function evaluation. With n = 64, both methods achieve an error on the order of 10^-7, while the basic midpoint and trapezoidal rules would require n > 4096 to achieve similar accuracy.

### Theoretical Validation

We can validate our observed convergence rates by comparing them with the theoretical error terms for each method:

1. Midpoint Rule: E_n(f) = O(h^2) ∝ 1/n^2
   - Observed: Convergence ratio ≈ 4 = 2^2 when n doubles (h halves)
   - Theoretical: Error should decrease by a factor of 2^2 = 4 when n doubles
   - Conclusion: Perfect agreement with theory

2. Trapezoidal Rule: E_n(f) = O(h^2) ∝ 1/n^2
   - Observed: Convergence ratio ≈ 4 = 2^2 when n doubles
   - Theoretical: Error should decrease by a factor of 2^2 = 4 when n doubles
   - Conclusion: Perfect agreement with theory

3. Corrected Trapezoidal Rule: E_n(f) = O(h^4) ∝ 1/n^4
   - Observed: Convergence ratio approaches 16 = 2^4 when n doubles
   - Theoretical: Error should decrease by a factor of 2^4 = 16 when n doubles
   - Conclusion: Excellent agreement with theory

4. Simpson's Rule: E_n(f) = O(h^4) ∝ 1/n^4
   - Observed: Convergence ratio approaches 16 = 2^4 as n increases
   - Theoretical: Error should decrease by a factor of 2^4 = 16 when n doubles
   - Conclusion: Good agreement with theory, especially for larger n values

The observed convergence rates match the theoretical predictions remarkably well, providing strong validation for both our implementation and the theoretical error analysis of these numerical integration methods.
