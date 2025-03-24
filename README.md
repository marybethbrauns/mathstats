
///

///

///

///

///

///

///

///

///

//

///

///



///

///

///

///

///




# Numerical Analysis Homework #6

## Problem 1: Evaluating $\int_{0}^{\pi} e^x \cos x\,dx$

### Problem Statement

Evaluate the integral
$$I(f)=\int_{0}^{\pi} e^x \cos x\,dx$$

using the following approaches:

1. **(a) Analytical Evaluation:** Use the Fundamental Theorem of Calculus (via integration by parts) to obtain an exact result.
2. **(b) Gauss–Legendre Quadrature:** Approximate the integral using Gauss–Legendre quadrature for $n=2,3,4,5,6$ nodes. Compute the approximations $I_n(f)$ and the corresponding errors $E_n(f)$.
3. **(c) Discussion:** Discuss the accuracy of Gauss–Legendre quadrature in comparison to Newton–Cotes formulas and comment on the computational effort involved.

### (a) Analytical Evaluation

To evaluate the integral $I(f)=\int_{0}^{\pi} e^x \cos x\,dx$, we'll use integration by parts twice.

First integration by parts with $u = e^x$ and $dv = \cos x\,dx$:
- $du = e^x\,dx$
- $v = \sin x$

This gives us:
$$\int e^x\cos x\,dx = e^x \sin x - \int e^x \sin x\,dx$$

For the remaining integral, we apply integration by parts again with $u = e^x$ and $dv = \sin x\,dx$:
- $du = e^x\,dx$
- $v = -\cos x$

This gives us:
$$\int e^x \sin x\,dx = -e^x\cos x + \int e^x\cos x\,dx$$

Substituting this result back into our first equation:
$$\int e^x\cos x\,dx = e^x \sin x - \left(-e^x\cos x + \int e^x\cos x\,dx\right)$$

Simplifying:
$$\int e^x\cos x\,dx = e^x (\sin x + \cos x) - \int e^x\cos x\,dx$$

Adding $\int e^x\cos x\,dx$ to both sides:
$$2\int e^x\cos x\,dx = e^x (\sin x + \cos x)$$

Therefore:
$$\int e^x\cos x\,dx = \frac{e^x (\sin x + \cos x)}{2}$$

Now, evaluating the definite integral:
$$I(f)=\left[\frac{e^x (\sin x + \cos x)}{2}\right]_0^{\pi}$$

At the limits:
- At $x=\pi$: $\sin \pi=0$ and $\cos \pi=-1$, giving $\frac{e^\pi (0-1)}{2} = -\frac{e^\pi}{2}$
- At $x=0$: $\sin 0=0$ and $\cos 0=1$, giving $\frac{e^0(0+1)}{2} = \frac{1}{2}$

Therefore:
$$I(f) = -\frac{e^\pi}{2} - \frac{1}{2} = -\frac{e^\pi+1}{2}$$

Numerically, with $e^\pi \approx 23.1407$:
$$I(f) \approx -\frac{24.1407}{2} \approx -12.07035$$

### (b) Gauss–Legendre Quadrature

To apply Gauss–Legendre quadrature, we need to transform the integration interval from $[0,\pi]$ to $[-1,1]$ using the substitution:

$$x = \frac{\pi}{2}(t+1), \quad dx = \frac{\pi}{2}\,dt$$

The transformed integrand becomes:
$$f(t)=\frac{\pi}{2}\,e^{\frac{\pi}{2}(t+1)}\cos\left(\frac{\pi}{2}(t+1)\right)$$

The $n$-point Gauss-Legendre quadrature approximation is:
$$I_n(f)=\sum_{i=1}^{n} w_i\,f(t_i)$$

where $t_i$ are the nodes and $w_i$ are the corresponding weights.

**Results for different values of $n$:**

| $n$ | $I_n(f)$ | Error $E_n(f)$ |
|-----|-----------|----------------|
| 2   | -12.33621047 | $2.66\times10^{-1}$ |
| 3   | -12.12742045 | $5.71\times10^{-2}$ |
| 4   | -12.07018949 | $1.57\times10^{-4}$ |
| 5   | -12.07032854 | $1.78\times10^{-5}$ |
| 6   | -12.07034633 | $1.47\times10^{-8}$ |

### (c) Discussion

The Gauss–Legendre quadrature method demonstrates impressive efficiency for this integral:

- **Rapid convergence:** The error drops dramatically from $2.66\times10^{-1}$ with just 2 nodes to $1.47\times10^{-8}$ with 6 nodes. This exponential error reduction is characteristic of Gauss quadrature for smooth functions.

- **Comparison with Newton-Cotes:** Newton–Cotes formulas (which use evenly spaced nodes) would require significantly more function evaluations to achieve similar accuracy. This is because Gauss-Legendre positions its nodes optimally to maximize accuracy.

- **Computational efficiency:** While Gauss–Legendre quadrature requires special nodes and weights, these can be precomputed and tabulated. For well-behaved functions like $e^x\cos(x)$, the method achieves high precision with minimal computational effort.

- **Error behavior:** The errors follow a pattern consistent with the theoretical error bound for Gauss-Legendre quadrature, which is proportional to the $(2n)$th derivative of the function for $n$ nodes.

## Problem 2: Integration of $\int_{-1}^{1}\frac{1}{1+25x^2}\,dx$

### Problem Statement

(a) **Analytical Evaluation:** Evaluate the integral $\int_{-1}^{1}\frac{1}{1+25x^2}\,dx$ using the Fundamental Theorem of Calculus.

(b) **Global Polynomial Interpolation via Newton's Divided Differences:** Approximate the integral by exactly integrating two interpolating polynomials:
   - A linear polynomial (using 2 nodes)
   - A quadratic polynomial (using 3 nodes)

(c) **Gauss–Legendre Quadrature:** Approximate the same integral using Gauss–Legendre quadrature with $n=2$, $4$, and $6$ node-point formulas.

### (a) Analytical Evaluation

We use the substitution $u = 5x$ to simplify the integral:
- $u = 5x$
- $du = 5\,dx$
- $dx = \frac{du}{5}$

This transforms our integral:
$$\int_{-1}^{1} \frac{1}{1+25x^2}\,dx = \int_{-5}^{5} \frac{1}{5} \cdot \frac{1}{1+u^2}\,du = \frac{1}{5}\int_{-5}^{5} \frac{1}{1+u^2}\,du$$

The standard form $\int \frac{1}{1+u^2}\,du = \arctan(u) + C$ gives us:
$$\int_{-1}^{1}\frac{1}{1+25x^2}\,dx = \frac{1}{5}\left[\arctan(u)\right]_{-5}^{5} = \frac{1}{5}[\arctan(5) - \arctan(-5)]$$

Since $\arctan(-x) = -\arctan(x)$, we get:
$$\int_{-1}^{1}\frac{1}{1+25x^2}\,dx = \frac{1}{5}[\arctan(5) - (-\arctan(5))] = \frac{2}{5}\arctan(5)$$

Numerically, with $\arctan(5) \approx 1.3734$:
$$\frac{2}{5}\arctan(5) \approx 0.54936$$

### (b) Integration via Newton's Divided Difference Interpolating Polynomials

This approach uses Newton's divided difference method to construct interpolating polynomials that are then integrated exactly.

#### Global 2-Point (Linear) Interpolation

Using nodes $x_0 = -1$ and $x_1 = 1$, with function values:
- $f(-1) = \frac{1}{1+25(-1)^2} = \frac{1}{26}$
- $f(1) = \frac{1}{1+25(1)^2} = \frac{1}{26}$

First, we construct the Newton's divided difference polynomial:

1. Zeroth divided difference: $f[x_0] = f(-1) = \frac{1}{26}$

2. First divided difference:
   $$f[x_0,x_1] = \frac{f(1) - f(-1)}{1-(-1)} = \frac{\frac{1}{26} - \frac{1}{26}}{2} = 0$$

3. The linear interpolation polynomial in Newton form:
P₁(x) = f[x₀] + f[x₀,x₁](x + 1) = 1/26 + 0 · (x + 1) = 1/26


Integrating this constant polynomial:
$$\int_{-1}^{1} P_1(x)\,dx = \frac{1}{26} \cdot (1-(-1)) = \frac{2}{26} = \frac{1}{13} \approx 0.07692$$

This value significantly underestimates the true integral (0.54936) because the constant interpolant fails to capture the sharp peak at x=0.

#### Global 3-Point (Quadratic) Interpolation

Using nodes $x_0 = -1$, $x_1 = 0$, and $x_2 = 1$, with function values:
- $f(-1) = \frac{1}{26}$
- $f(0) = \frac{1}{1+25(0)^2} = 1$
- $f(1) = \frac{1}{26}$

Computing the divided differences:

1. First divided differences:
   $$f[x_0,x_1] = \frac{f(0) - f(-1)}{0-(-1)} = \frac{1 - \frac{1}{26}}{1} = \frac{25}{26}$$

   $$f[x_1,x_2] = \frac{f(1) - f(0)}{1-0} = \frac{\frac{1}{26} - 1}{1} = -\frac{25}{26}$$

2. Second divided difference:
   $$f[x_0,x_1,x_2] = \frac{f[x_1,x_2] - f[x_0,x_1]}{x_2-x_0} = \frac{-\frac{25}{26} - \frac{25}{26}}{2} = -\frac{25}{26}$$

3. The quadratic interpolation polynomial in Newton form:
   P₂(x) = f[x₀] + f[x₀,x₁](x - x₀) + f[x₀,x₁,x₂](x - x₀)(x - x₁)

   Substituting values:
   $$P_2(x) = \frac{1}{26} + \frac{25}{26}(x+1) - \frac{25}{26}(x+1)(x)$$

   Simplifying and noting that $(x+1)(x) = x^2 + x$:
   $$P_2(x) = \frac{1}{26} + \frac{25}{26}(x+1) - \frac{25}{26}(x^2+x) = \frac{1}{26} + \frac{25}{26} + \frac{25x}{26} - \frac{25x^2}{26} - \frac{25x}{26}$$

   This further simplifies to:
   $$P_2(x) = \frac{26}{26} - \frac{25x^2}{26} = 1 - \frac{25}{26}x^2$$

Integrating the quadratic polynomial:
$$\int_{-1}^{1} P_2(x)\,dx = \int_{-1}^{1} \left(1 - \frac{25}{26}x^2\right) dx = \left[x - \frac{25}{26}\frac{x^3}{3}\right]_{-1}^{1}$$

$$= \left(1 - \frac{25}{26}\frac{1}{3}\right) - \left(-1 - \frac{25}{26}\frac{-1}{3}\right) = 2 - 2\frac{25}{26}\frac{1}{3} = 2 - \frac{50}{78} = \frac{78-50}{39} = \frac{28}{39} \approx 0.7179$$

This value overestimates the true integral (0.54936), but is closer than the linear approximation. The quadratic interpolant captures the central peak but doesn't properly model the function's decay away from the origin.

### (c) Gauss–Legendre Quadrature

For the integral $\int_{-1}^{1}\frac{1}{1+25x^2}\,dx$, Gauss-Legendre quadrature yields the following results:

| $n$ | Approximation | Absolute Error | Relative Error |
|-----|---------------|----------------|----------------|
| 2   | 0.2142857143  | 0.33507        | 60.99%         |
| 4   | 0.3709273183  | 0.17843        | 32.48%         |
| 6   | 0.4617005584  | 0.08766        | 15.96%         |

**Observations:**
- Increasing the number of nodes steadily improves the approximation
- Even with 6 nodes, the quadrature result (0.4617) underestimates the true value (0.54936)
- This is due to the function's sharp peak at x=0, which is challenging to capture with global methods
- The function $f(x)=\frac{1}{1+25x^2}$ has a Runge phenomenon-like behavior in the sense that it has a very steep gradient near x=0

## Problem 3: Determining a Hermite-Type Quadrature Formula

### Problem Statement

Determine constants $a$, $b$, $c$, and $d$ so that the quadrature formula
$$\int_{-1}^{1} f(x)\,dx = a\,f(-1) + b\,f(1) + c\,f'(-1) + d\,f'(1)$$
is exact for all polynomials $f(x)$ of degree 3 or less.

### Derivation

For the formula to be exact for all polynomials up to degree 3, it must produce exact results for $f(x)=1$, $f(x)=x$, $f(x)=x^2$, and $f(x)=x^3$. This gives us four equations to determine our four unknowns.

#### Equation 1: $f(x)=1$
- Function values: $f(-1)=1$, $f(1)=1$
- Derivatives: $f'(-1)=0$, $f'(1)=0$
- Exact integral: $\int_{-1}^{1}1\,dx = 2$
- Resulting equation: $a + b = 2$

#### Equation 2: $f(x)=x$
- Function values: $f(-1)=-1$, $f(1)=1$
- Derivatives: $f'(-1)=1$, $f'(1)=1$
- Exact integral: $\int_{-1}^{1}x\,dx = 0$
- Resulting equation: $-a + b + c + d = 0$

#### Equation 3: $f(x)=x^2$
- Function values: $f(-1)=1$, $f(1)=1$
- Derivatives: $f'(-1)=-2$, $f'(1)=2$
- Exact integral: $\int_{-1}^{1}x^2\,dx = \frac{2}{3}$
- Resulting equation: $a + b - 2c + 2d = \frac{2}{3}$

#### Equation 4: $f(x)=x^3$
- Function values: $f(-1)=-1$, $f(1)=1$
- Derivatives: $f'(-1)=3$, $f'(1)=3$
- Exact integral: $\int_{-1}^{1}x^3\,dx = 0$
- Resulting equation: $-a + b + 3c + 3d = 0$

#### Solving the System of Equations

From equation 1: $b = 2 - a$

Substituting into equation 2:
$-a + (2-a) + c + d = 0$
$-2a + 2 + c + d = 0$
$c + d = 2a - 2$

Substituting $b = 2-a$ into equation 3:
$a + (2-a) - 2c + 2d = \frac{2}{3}$
$2 - 2c + 2d = \frac{2}{3}$
$-c + d = \frac{1}{3} - 1 = -\frac{2}{3}$

From these two equations:
$c + d = 2a - 2$
$-c + d = -\frac{2}{3}$

Adding them:
$(c + d) + (-c + d) = (2a - 2) + (-\frac{2}{3})$
$2d = 2a - 2 - \frac{2}{3}$
$2d = 2a - \frac{8}{3}$
$d = a - \frac{4}{3}$

Substituting back:
$c + \left(a - \frac{4}{3}\right) = 2a - 2$
$c = 2a - 2 - a + \frac{4}{3} = a - \frac{2}{3}$

Now we substitute these expressions for $c$ and $d$ into equation 4:
$-a + (2-a) + 3\left(a - \frac{2}{3}\right) + 3\left(a - \frac{4}{3}\right) = 0$

Simplifying:
$-a + 2 - a + 3a - 2 + 3a - 4 = 0$
$-2a + 6a - 4 = 0$
$4a = 4$
$a = 1$

Substituting back:
$b = 2 - 1 = 1$
$c = 1 - \frac{2}{3} = \frac{1}{3}$
$d = 1 - \frac{4}{3} = -\frac{1}{3}$

Therefore, the Hermite-type quadrature formula is:
$$\int_{-1}^{1} f(x)\,dx \approx f(-1) + f(1) + \frac{1}{3}f'(-1) - \frac{1}{3}f'(1)$$

#### Verification

Let's verify this formula with our test functions:

For $f(x) = 1$:
- Exact: $\int_{-1}^{1} 1\,dx = 2$
- Our formula: $1 \cdot 1 + 1 \cdot 1 + \frac{1}{3} \cdot 0 - \frac{1}{3} \cdot 0 = 2$ ✓

For $f(x) = x$:
- Exact: $\int_{-1}^{1} x\,dx = 0$
- Our formula: $1 \cdot (-1) + 1 \cdot 1 + \frac{1}{3} \cdot 1 - \frac{1}{3} \cdot 1 = 0$ ✓

For $f(x) = x^2$:
- Exact: $\int_{-1}^{1} x^2\,dx = \frac{2}{3}$
- Our formula: $1 \cdot 1 + 1 \cdot 1 + \frac{1}{3} \cdot (-2) - \frac{1}{3} \cdot 2 = 2 - \frac{2}{3} - \frac{2}{3} = \frac{2}{3}$ ✓

For $f(x) = x^3$:
- Exact: $\int_{-1}^{1} x^3\,dx = 0$
- Our formula: $1 \cdot (-1) + 1 \cdot 1 + \frac{1}{3} \cdot 3 - \frac{1}{3} \cdot 3 = 0$ ✓

The formula correctly integrates all polynomials of degree 3 or less.

## Conclusion

This homework explored three different numerical integration problems, highlighting the strengths and limitations of various integration techniques:

1. **Problem 1** demonstrated that Gauss-Legendre quadrature converges extremely rapidly for smooth functions ($e^x\cos x$), achieving high precision with relatively few function evaluations.

2. **Problem 2** revealed the challenges of integrating functions with sharp peaks ($\frac{1}{1+25x^2}$):
   - Global polynomial interpolation using Newton's divided differences can lead to significant under- or overestimation
   - Even Gauss-Legendre quadrature requires more nodes to achieve acceptable accuracy
   - The function's shape creates difficulties for all global approximation methods

3. **Problem 3** showed how incorporating derivative information in Hermite-type quadrature can improve the degree of precision, achieving exactness for cubic polynomials with just two function evaluations and two derivative evaluations.

The primary insights from this homework are:

- **Function smoothness** greatly impacts the convergence rate of numerical integration methods
- **Gauss-Legendre quadrature** generally outperforms Newton-Cotes formulas for the same number of function evaluations
- **Global polynomial interpolation** may not be suitable for functions with sharp features
- **Derivative information** can enhance the precision of quadrature formulas

These concepts are fundamental to choosing appropriate numerical integration strategies for different types of functions and required accuracy levels.
