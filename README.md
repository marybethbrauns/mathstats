# Numerical Analaysis HW #6 | Marybeth Brauns

## Problem 1: Evaluating $\int_{0}^{\pi} e^x \cos x \, dx$

### Problem Statement:
This problem requires the evaluation of an integral involving exponential and trigonometric functions. The purpose is to compare different integration techniques and their accuracy. The task involves finding $\int_{0}^{\pi} e^x \cos x \, dx$ using analytical methods and numerical approximations, followed by analysis of the results. This integral serves as an excellent case study because it combines two important functions (exponential and trigonometric) whose interplay creates an interesting challenge for numerical methods. The objective is to demonstrate how quickly different methods converge to the exact solution and to understand the computational tradeoffs involved.

### (a) Evaluation via Integration by Parts

To evaluate $I(f) = \int_{0}^{\pi} e^x \cos x \, dx$, we must find an antiderivative of $f(x) = e^x \cos x$.

Using integration by parts with $u = e^x$ and $dv = \cos x \, dx$:
$$\int e^x \cos x \, dx = e^x \sin x - \int e^x \sin x \, dx$$

For the remaining integral, applying integration by parts again with $u = e^x$ and $dv = \sin x \, dx$:
$$\int e^x \sin x \, dx = -e^x \cos x + \int e^x \cos x \, dx$$

Substituting this result back into our original equation:
$$\int e^x \cos x \, dx = e^x \sin x - \left(-e^x \cos x + \int e^x \cos x \, dx\right)$$

Simplifying:
$$\int e^x \cos x \, dx = e^x \sin x + e^x \cos x - \int e^x \cos x \, dx$$
$$2\int e^x \cos x \, dx = e^x (\sin x + \cos x)$$
$$\int e^x \cos x \, dx = \frac{e^x (\sin x + \cos x)}{2}$$

Now evaluating the definite integral:
$$I(f) = \int_{0}^{\pi} e^x \cos x \, dx = \left[ \frac{e^x (\sin x + \cos x)}{2} \right]_{0}^{\pi}$$
$$= \frac{e^{\pi}(\sin \pi + \cos \pi)}{2} - \frac{e^0(\sin 0 + \cos 0)}{2}$$
$$= \frac{e^{\pi}(0 + (-1))}{2} - \frac{1(0 + 1)}{2}$$
$$= -\frac{e^{\pi}}{2} - \frac{1}{2} = -\frac{e^{\pi} + 1}{2}$$

Calculating this value numerically:
$$e^{\pi} \approx 23.1407$$
$$-\frac{e^{\pi} + 1}{2} \approx -\frac{24.1407}{2} \approx -12.070346316389633$$

### (b) Gauss-Legendre Quadrature Approximation

To apply Gauss-Legendre quadrature, we must first transform the integration interval from $[0,\pi]$ to $[-1,1]$ using the substitution:
$$x = \frac{\pi}{2}(t+1)$$

This gives:
$$dx = \frac{\pi}{2}dt$$

The transformed integrand becomes:
$$f(t) = \frac{\pi}{2}e^{\frac{\pi}{2}(t+1)}\cos\left(\frac{\pi}{2}(t+1)\right)$$

The $n$-point Gauss-Legendre approximation is given by:
$$I_n(f) = \sum_{i=1}^{n} w_i f(t_i)$$

where $t_i$ are the Gauss-Legendre nodes and $w_i$ are the corresponding weights.

The Gauss-Legendre quadrature approximations for various values of $n$ are presented in Table 1.

**Table 1.** Gauss-Legendre Quadrature Results for the integral from 0 to π of e^x·cos(x)dx

| $n$ | $I_n(f)$ | Error |
|-----|----------|-------|
| 2   | -12.33621047 | 2.6586e-1 |
| 3   | -12.12742045 | 5.7074e-2 |
| 4   | -12.07018949 | 1.5683e-4 |
| 5   | -12.07032854 | 1.7781e-5 |
| 6   | -12.07034633 | 1.4712e-8 |

### (c) Analysis of Quadrature Methods

For the function $f(x) = e^x \cos x$, we observe:

1. **Error Convergence**:
   - The Gauss-Legendre quadrature demonstrates spectral convergence
   - The error decreases by approximately two orders of magnitude when increasing from $n$ to $n+1$ points
   - At $n=6$, the error is already on the order of $10^{-8}$

2. **Comparative Efficiency**:
   - Gauss-Legendre quadrature achieves higher accuracy with fewer function evaluations compared to Newton-Cotes formulas
   - For example, 6-point Gauss-Legendre quadrature achieves an error of approximately $10^{-8}$, which would require significantly more points with Newton-Cotes methods
   - Newton-Cotes formulas with $n$ points can exactly integrate polynomials of degree at most $n-1$ (open) or $n$ (closed)
   - Gauss-Legendre quadrature with $n$ points can exactly integrate polynomials of degree up to $2n-1$

3. **Computational Considerations**:
   - Newton-Cotes formulas require less preprocessing since they use evenly spaced nodes
   - Gauss-Legendre requires computation of specific nodes and weights, which involves finding roots of Legendre polynomials
   - However, these nodes and weights can be pre-computed and tabulated
   - For smooth, analytic functions like $e^x \cos x$, Gauss-Legendre is significantly more efficient in terms of function evaluations required for a given accuracy

## Problem 2: Integrating $\int_{-1}^{1} \frac{1}{1+25x^2} \, dx$

### Problem Statement:
This problem focuses on evaluating the integral of a rational function with a sharp peak. The function $\frac{1}{1+25x^2}$ forms a tall, narrow spike at x=0, which makes it particularly difficult for polynomial-based approximation methods. The purpose of this problem is to explore how different numerical integration techniques handle functions with localized features. By comparing polynomial interpolation with Gauss-Legendre quadrature, the phenomenon of Runge oscillation becomes observable, illustrating why adaptive methods might be preferable for certain classes of functions. The goal is to determine which method provides the best accuracy with the least computational effort for this challenging integrand.

### (a) Evaluation via Substitution Method

For the integrand $f(x) = \frac{1}{1+25x^2}$, we use the substitution $u = 5x$, which gives $du = 5\,dx$:

$$\int \frac{1}{1+25x^2} \, dx = \int \frac{1}{1+(5x)^2} \, dx = \frac{1}{5} \int \frac{1}{1+u^2} \, du = \frac{1}{5} \arctan(u) + C = \frac{1}{5} \arctan(5x) + C$$

Evaluating the definite integral:
$$\int_{-1}^{1} \frac{1}{1+25x^2} \, dx = \left[ \frac{1}{5} \arctan(5x) \right]_{-1}^{1}$$
$$= \frac{1}{5} \arctan(5) - \frac{1}{5} \arctan(-5)$$

Since $\arctan$ is an odd function, $\arctan(-5) = -\arctan(5)$, so:
$$= \frac{1}{5} \arctan(5) - \frac{1}{5} \cdot (-\arctan(5))$$
$$= \frac{2}{5} \arctan(5)$$

Computing the value:
$$\arctan(5) \approx 1.3734$$
$$\frac{2}{5} \arctan(5) \approx \frac{2}{5} \times 1.3734 \approx 0.5493603067780064$$

### (b) Integration of Interpolating Polynomials

Let us examine how polynomial interpolation performs for this function. We'll construct interpolation polynomials using different sets of points and integrate them analytically.

#### 2-point Equally-Spaced Interpolation:

Using points $x_0 = -1$ and $x_1 = 1$:

1. Function values:
   - $f(-1) = \frac{1}{1+25(-1)^2} = \frac{1}{1+25} = \frac{1}{26} \approx 0.0384615385$
   - $f(1) = \frac{1}{1+25(1)^2} = \frac{1}{1+25} = \frac{1}{26} \approx 0.0384615385$

2. Lagrange basis polynomials:
   - $L_0(x) = \frac{x-1}{-1-1} = \frac{x-1}{-2} = \frac{1-x}{2}$
   - $L_1(x) = \frac{x-(-1)}{1-(-1)} = \frac{x+1}{2}$

3. Interpolation polynomial:
   - $P_1(x) = f(-1)L_0(x) + f(1)L_1(x)$
   - $P_1(x) = \frac{1}{26} \cdot \frac{1-x}{2} + \frac{1}{26} \cdot \frac{x+1}{2}$
   - $P_1(x) = \frac{1}{52}(1-x) + \frac{1}{52}(x+1)$
   - $P_1(x) = \frac{1}{52}(1-x+x+1) = \frac{2}{52} = \frac{1}{26}$

4. Exact integration:
   - $\int_{-1}^{1} P_1(x) \, dx = \int_{-1}^{1} \frac{1}{26} \, dx = \frac{1}{26} \cdot 2 = \frac{2}{26} = \frac{1}{13} \approx 0.0769230769$

5. Error:
   - Exact value: $\frac{2}{5}\arctan(5) \approx 0.5493603068$
   - Error: $|0.5493603068 - 0.0769230769| \approx 0.4724$

#### 3-point Equally-Spaced Interpolation:

Using points $x_0 = -1$, $x_1 = 0$, and $x_2 = 1$:

1. Function values:
   - $f(-1) = \frac{1}{26} \approx 0.0384615385$
   - $f(0) = \frac{1}{1+25(0)^2} = \frac{1}{1} = 1$
   - $f(1) = \frac{1}{26} \approx 0.0384615385$

2. Lagrange basis polynomials:
   - $L_0(x) = \frac{(x-0)(x-1)}{(-1-0)(-1-1)} = \frac{x(x-1)}{(-1)(-2)} = \frac{x(x-1)}{2}$
   - $L_1(x) = \frac{(x+1)(x-1)}{(0+1)(0-1)} = \frac{(x+1)(x-1)}{1(-1)} = -(x+1)(x-1) = -(x^2-1) = -x^2+1$
   - $L_2(x) = \frac{(x+1)(x-0)}{(1+1)(1-0)} = \frac{(x+1)x}{2}$

3. Interpolation polynomial:
   - $P_2(x) = f(-1)L_0(x) + f(0)L_1(x) + f(1)L_2(x)$
   - $P_2(x) = \frac{1}{26} \cdot \frac{x(x-1)}{2} + 1 \cdot (-x^2+1) + \frac{1}{26} \cdot \frac{x(x+1)}{2}$
   
   Expanding and simplifying:
   - $P_2(x) = \frac{1}{52}x(x-1) + (-x^2+1) + \frac{1}{52}x(x+1)$
   - $P_2(x) = \frac{1}{52}(x^2-x) + (-x^2+1) + \frac{1}{52}(x^2+x)$
   - $P_2(x) = \frac{1}{52}(2x^2) - x^2 + 1$
   - $P_2(x) = \frac{2x^2}{52} - x^2 + 1$
   - $P_2(x) = -\frac{50x^2}{52} + 1$
   - $P_2(x) = -\frac{25x^2}{26} + 1$

4. Exact integration:
   - $\int_{-1}^{1} P_2(x) \, dx = \int_{-1}^{1} (-\frac{25x^2}{26} + 1) \, dx$
   - $= -\frac{25}{26} \int_{-1}^{1} x^2 \, dx + \int_{-1}^{1} 1 \, dx$
   - $= -\frac{25}{26} \cdot [\frac{x^3}{3}]_{-1}^{1} + [x]_{-1}^{1}$
   - $= -\frac{25}{26} \cdot (\frac{1}{3} - \frac{-1}{3}) + (1 - (-1))$
   - $= -\frac{25}{26} \cdot \frac{2}{3} + 2$
   - $= -\frac{50}{78} + \frac{156}{78}$
   - $= \frac{106}{78} \approx 1.3590$

5. Error:
   - Error: $|0.5493603068 - 1.3590| \approx 0.8096$

For all interpolation degrees from 2 to 10, the results are summarized in Table 2.

**Table 2.** Polynomial Interpolation Results for the integral from -1 to 1 of 1/(1+25x²)dx

| $n$ | Equally Spaced | Error | Chebyshev | Error |
|-----|----------------|-------|-----------|-------|
| 2   | 0.0769230769   | 4.7244e-1 | 0.1481481481 | 4.0121e-1 |
| 3   | 1.3589743590   | 8.0961e-1 | 1.1561181435 | 6.0676e-1 |
| 4   | 0.4162895928   | 1.3307e-1 | 0.3393357343 | 2.1002e-1 |
| 5   | 0.4748010610   | 7.4559e-2 | 0.7366108212 | 1.8725e-1 |
| 6   | 0.4615384615   | 8.7822e-2 | 0.4422623071 | 1.0710e-1 |
| 7   | 0.7740897347   | 2.2473e-1 | 0.6363602552 | 8.7000e-2 |
| 8   | 0.5797988819   | 3.0439e-2 | 0.4995830749 | 4.9777e-2 |
| 9   | 0.3000977814   | 2.4926e-1 | 0.5839263513 | 3.4566e-2 |
| 10  | 0.4797235796   | 6.9637e-2 | 0.5259711610 | 2.3389e-2 |

#### Observations on Interpolating Polynomials:

The results highlight why polynomial interpolation struggles with this function:

1. **Low-degree approximation (2 points)**:
   - The linear interpolant forms a flat line at height $\frac{1}{26}$, completely missing the peak at $x=0$
   - This leads to significant underestimation of the integral (only about 14% of the true value)

2. **Adding the peak (3 points)**:
   - Including the point at $x=0$ creates a parabola that matches the function at $x = 0$ (value = 1)
   - However, the quadratic approximation decreases much more slowly than the original function
   - This leads to significant overestimation (about 247% of the true value)

3. **Runge's phenomenon challenge**:
   - As we increase the degree of interpolation and add more equally spaced points, oscillations appear
   - These oscillations can lead to over/underestimation in different regions
   - This explains the non-monotonic convergence pattern seen in the data table
   - For example, the error with 7 points is worse than with 6 points

4. **Chebyshev points advantage**:
   - Chebyshev points distribute the interpolation points with higher density near the endpoints
   - This helps reduce oscillations from Runge's phenomenon
   - At $n=10$, Chebyshev points give a more accurate result (error ≈ 2.3%) than equally spaced points (error ≈ 7.0%)
   - However, for this function with a central peak, Chebyshev points aren't optimally distributed either

### (c) Gauss-Legendre Quadrature

For this integral, we apply Gauss-Legendre quadrature with $n = 2, 4, 6$ nodes:

The complete Gauss-Legendre quadrature results for this integral are presented in Table 3.

**Table 3.** Gauss-Legendre Quadrature Results for the integral from -1 to 1 of 1/(1+25x²)dx

| $n$ | Approximation | Absolute Error | Relative Error |
|-----|---------------|----------------|----------------|
| 2   | 0.2142857143  | 3.3507e-1      | 6.0994e-1      |
| 4   | 0.3709273183  | 1.7843e-1      | 3.2480e-1      |
| 6   | 0.4617005584  | 8.7660e-2      | 1.5957e-1      |

#### Observations:

1. Convergence is much slower than for the first integral
2. Even with $n=6$ nodes, the relative error remains approximately 16%
3. The function's sharp peak presents a significant challenge for polynomial-based quadrature
4. The performance is still better than polynomial interpolation but requires more points for high accuracy
5. For this type of function with localized features, adaptive methods would be more appropriate

## Problem 3: Developing a Hermite-Type Quadrature Formula

### Problem Statement:
This problem requires the construction of a specialized numerical integration formula that leverages both function values and derivatives. The goal is to determine specific coefficients a, b, c, and d for the formula:
$$\int_{-1}^{1} f(x) \, dx = a \cdot f(-1) + b \cdot f(1) + c \cdot f'(-1) + d \cdot f'(1)$$

The purpose of this problem is to demonstrate how derivative information can enhance the accuracy of quadrature formulas. While standard methods like Newton-Cotes and Gauss-Legendre use only function values, this Hermite-type approach incorporates the function's rate of change at specific points. By requiring exactness for polynomials up to degree 3, the formula achieves higher-order accuracy with minimal evaluation points. This highlights an important principle in numerical analysis: utilizing additional information about a function (beyond simple values) can significantly improve approximation quality.

To find these constants, we need the formula to be exact for the basis $\{1, x, x^2, x^3\}$.

### Testing on $f(x) = 1$:

Function and derivative values:
- $f(-1) = 1$
- $f(1) = 1$
- $f'(-1) = 0$
- $f'(1) = 0$

Exact integral:
$$\int_{-1}^{1} 1 \, dx = [x]_{-1}^{1} = 1 - (-1) = 2$$

This gives our first equation:
$$a \cdot 1 + b \cdot 1 + c \cdot 0 + d \cdot 0 = 2$$
$$a + b = 2 \quad \text{(Equation 1)}$$

### Testing on $f(x) = x$:

Function and derivative values:
- $f(-1) = -1$
- $f(1) = 1$
- $f'(-1) = 1$ (derivative of $x$ is 1)
- $f'(1) = 1$

Exact integral:
$$\int_{-1}^{1} x \, dx = \left[\frac{x^2}{2}\right]_{-1}^{1} = \frac{1}{2} - \frac{1}{2} = 0$$

This gives our second equation:
$$a \cdot (-1) + b \cdot 1 + c \cdot 1 + d \cdot 1 = 0$$
$$-a + b + c + d = 0 \quad \text{(Equation 2)}$$

### Testing on $f(x) = x^2$:

Function and derivative values:
- $f(-1) = (-1)^2 = 1$
- $f(1) = 1^2 = 1$
- $f'(-1) = 2(-1) = -2$ (derivative of $x^2$ is $2x$)
- $f'(1) = 2(1) = 2$

Exact integral:
$$\int_{-1}^{1} x^2 \, dx = \left[\frac{x^3}{3}\right]_{-1}^{1} = \frac{1}{3} - \frac{-1}{3} = \frac{2}{3}$$

This gives our third equation:
$$a \cdot 1 + b \cdot 1 + c \cdot (-2) + d \cdot 2 = \frac{2}{3}$$
$$a + b - 2c + 2d = \frac{2}{3} \quad \text{(Equation 3)}$$

### Testing on $f(x) = x^3$:

Function and derivative values:
- $f(-1) = (-1)^3 = -1$
- $f(1) = 1^3 = 1$
- $f'(-1) = 3(-1)^2 = 3$ (derivative of $x^3$ is $3x^2$)
- $f'(1) = 3(1)^2 = 3$

Exact integral:
$$\int_{-1}^{1} x^3 \, dx = \left[\frac{x^4}{4}\right]_{-1}^{1} = \frac{1}{4} - \frac{1}{4} = 0$$

This gives our fourth equation:
$$a \cdot (-1) + b \cdot 1 + c \cdot 3 + d \cdot 3 = 0$$
$$-a + b + 3c + 3d = 0 \quad \text{(Equation 4)}$$

### Solving the System of Equations:

From Equation 1:
$$b = 2 - a \quad \text{(Equation 5)}$$

Substituting Equation 5 into Equation 2:
$$-a + (2 - a) + c + d = 0$$
$$-a + 2 - a + c + d = 0$$
$$-2a + 2 + c + d = 0$$
$$c + d = 2a - 2 \quad \text{(Equation 6)}$$

Substituting Equation 5 into Equation 3:
$$a + (2 - a) - 2c + 2d = \frac{2}{3}$$
$$a + 2 - a - 2c + 2d = \frac{2}{3}$$
$$2 - 2c + 2d = \frac{2}{3}$$
$$-c + d = \frac{2/3 - 2}{2} = \frac{2/3 - 6/3}{2} = \frac{-4/3}{2} = -\frac{2}{3} \quad \text{(Equation 7)}$$

From Equations 6 and 7:
$$c + d = 2a - 2 \quad \text{(Equation 6)}$$
$$-c + d = -\frac{2}{3} \quad \text{(Equation 7)}$$

Adding these equations:
$$2d = 2a - 2 - \frac{2}{3}$$
$$2d = 2a - \frac{6}{3} - \frac{2}{3}$$
$$2d = 2a - \frac{8}{3}$$
$$d = a - \frac{4}{3} \quad \text{(Equation 8)}$$

Substituting Equation 8 into Equation 6:
$$c + \left(a - \frac{4}{3}\right) = 2a - 2$$
$$c = 2a - 2 - a + \frac{4}{3}$$
$$c = a - 2 + \frac{4}{3}$$
$$c = a - \frac{6}{3} + \frac{4}{3}$$
$$c = a - \frac{2}{3} \quad \text{(Equation 9)}$$

Now, substituting Equations 5, 8, and 9 into Equation 4:
$$-a + (2-a) + 3\left(a - \frac{2}{3}\right) + 3\left(a - \frac{4}{3}\right) = 0$$
$$-a + 2 - a + 3a - 2 + 3a - 4 = 0$$
$$-a - a + 3a + 3a + 2 - 2 - 4 = 0$$
$$-2a + 6a - 4 = 0$$
$$4a = 4$$
$$a = 1$$

Substituting back:
$$b = 2 - a = 2 - 1 = 1 \quad \text{(from Equation 5)}$$
$$c = a - \frac{2}{3} = 1 - \frac{2}{3} = \frac{3}{3} - \frac{2}{3} = \frac{1}{3} \quad \text{(from Equation 9)}$$
$$d = a - \frac{4}{3} = 1 - \frac{4}{3} = \frac{3}{3} - \frac{4}{3} = -\frac{1}{3} \quad \text{(from Equation 8)}$$

Therefore, the quadrature formula is:
$$\int_{-1}^{1} f(x) \, dx \approx 1 \cdot f(-1) + 1 \cdot f(1) + \frac{1}{3} \cdot f'(-1) - \frac{1}{3} \cdot f'(1)$$

### Verification:

Let's verify this formula works for polynomials of degree $\leq 3$:

For $f(x) = 1$:
$$1 \cdot 1 + 1 \cdot 1 + \frac{1}{3} \cdot 0 - \frac{1}{3} \cdot 0 = 2 \quad \text{(Correct!)}$$

For $f(x) = x$:
$$1 \cdot (-1) + 1 \cdot 1 + \frac{1}{3} \cdot 1 - \frac{1}{3} \cdot 1 = -1 + 1 + \frac{1}{3} - \frac{1}{3} = 0 \quad \text{(Correct!)}$$

For $f(x) = x^2$:
$$1 \cdot 1 + 1 \cdot 1 + \frac{1}{3} \cdot (-2) - \frac{1}{3} \cdot 2 = 1 + 1 - \frac{2}{3} - \frac{2}{3} = 2 - \frac{4}{3} = \frac{6}{3} - \frac{4}{3} = \frac{2}{3} \quad \text{(Correct!)}$$

For $f(x) = x^3$:
$$1 \cdot (-1) + 1 \cdot 1 + \frac{1}{3} \cdot 3 - \frac{1}{3} \cdot 3 = -1 + 1 + 1 - 1 = 0 \quad \text{(Correct!)}$$

The formula is therefore verified to be exact for polynomials of degree $\leq 3$.

## Conclusion

This analysis demonstrates the strengths and limitations of different numerical integration methods:

1. **Gauss-Legendre Quadrature**:
   - Highly efficient for smooth functions (like $e^x \cos x$)
   - Requires more points for functions with sharp features (like $\frac{1}{1+25x^2}$)
   - Achieves accuracy with fewer function evaluations than Newton-Cotes methods

2. **Polynomial Interpolation**:
   - Struggles with Runge's phenomenon, especially for equally spaced points
   - Chebyshev points help mitigate oscillations but aren't always optimal
   - May exhibit non-monotonic convergence for challenging functions

3. **Hermite-Type Quadrature**:
   - Incorporates derivative information to enhance accuracy
   - Can be designed to be exact for specific polynomial degrees
   - Provides a balance between accuracy and computational effort

The choice of method depends on the characteristics of the integrand and the desired accuracy. For functions with sharp features or singularities, adaptive methods that focus computational effort where the function varies rapidly are generally more effective.
