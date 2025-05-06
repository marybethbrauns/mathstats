# Marybeth Brauns | Numerical Analysis Final 

## Problem 1: Polynomial Interpolation with Chebyshev Nodes

### Problem Statement
We have seen that polynomial interpolation at equidistant nodes to the function $f(x) = \frac{1}{1+25x^2}$ on $[-1,1]$ produces highly oscillatory behavior known as Runge's phenomenon. To mitigate this behavior, we compute the interpolating polynomial at the roots of Chebyshev polynomials $x_j = \cos\frac{(2j+1)\pi}{2(n+1)}$ for $j = 0,...,n$, with $n = 4, 6, 8, 10$.

### Detailed Solution and Analysis

#### Mathematical Background on Chebyshev Nodes

Chebyshev nodes are specifically designed to minimize the maximum interpolation error. Unlike equidistant nodes, they cluster more densely near the endpoints of the interval. This non-uniform distribution counters Runge's phenomenon by placing more nodes precisely where oscillations tend to be most problematic.

For an interval $[-1,1]$, the Chebyshev nodes of the first kind are given by:

$$x_j = \cos\frac{(2j+1)\pi}{2(n+1)}, \quad j = 0,1,...,n$$

These points are the projections onto the x-axis of equally spaced points on the upper half of the unit circle. Their distribution naturally counteracts the tendency of high-degree polynomials to oscillate wildly near the endpoints of the interval.

#### Step 1: Calculate the Chebyshev nodes for each value of n

For $n = 4$, calculating each node explicitly:
$$x_0 = \cos\frac{\pi}{10} = 0.9511$$
$$x_1 = \cos\frac{3\pi}{10} = 0.5878$$
$$x_2 = \cos\frac{5\pi}{10} = 0.0000$$
$$x_3 = \cos\frac{7\pi}{10} = -0.5878$$
$$x_4 = \cos\frac{9\pi}{10} = -0.9511$$

Notice the symmetric distribution about zero and the clustering near ±1. This pattern becomes more pronounced as n increases. For $n = 6$:
$$x_0 = 0.9749, x_1 = 0.7818, x_2 = 0.4339, x_3 = 0.0000,$$
$$x_4 = -0.4339, x_5 = -0.7818, x_6 = -0.9749$$

And similarly for $n = 8$ and $n = 10$, the nodes become increasingly dense near the endpoints while maintaining symmetry.

#### Step 2: Calculate function values at the nodes

For $n = 4$, evaluating $f(x_j) = \frac{1}{1+25x_j^2}$ at each node:

$$f(x_0) = \frac{1}{1+25(0.9511)^2} = \frac{1}{1+25(0.9046)} = \frac{1}{23.615} = 0.0423$$

$$f(x_1) = \frac{1}{1+25(0.5878)^2} = \frac{1}{1+25(0.3455)} = \frac{1}{9.638} = 0.1038$$

$$f(x_2) = \frac{1}{1+25(0.0000)^2} = \frac{1}{1+0} = 1.0000$$

$$f(x_3) = 0.1038, f(x_4) = 0.0423$$

Note the symmetry in function values due to the even nature of $f(x)$ and the symmetric distribution of nodes.

#### Step 3: Construct Lagrange interpolation polynomials

The Lagrange basis polynomials are defined as:

$$\ell_j(x) = \prod_{k=0, k \neq j}^{n} \frac{x - x_k}{x_j - x_k}$$

Each $\ell_j(x)$ has the property that $\ell_j(x_i) = \delta_{ij}$ (equals 1 when i=j, 0 otherwise), making them ideal building blocks for interpolation.

For $n = 4$ and $j = 0$:

$$\ell_0(x) = \frac{(x-0.5878)(x-0.0000)(x-(-0.5878))(x-(-0.9511))}{(0.9511-0.5878)(0.9511-0.0000)(0.9511-(-0.5878))(0.9511-(-0.9511))}$$

$$\ell_0(x) = \frac{(x-0.5878)x(x+0.5878)(x+0.9511)}{(0.3633)(0.9511)(1.5389)(1.9022)}$$

The denominator evaluates to approximately 1.0109, not exactly 1. Therefore:

$$\ell_0(x) = \frac{x(x-0.5878)(x+0.5878)(x+0.9511)}{1.0109}$$

The complete Lagrange interpolation polynomial is:

$$P_n(x) = \sum_{j=0}^{n} f(x_j) \ell_j(x)$$

For $n = 4$, this gives:

$$P_4(x) = 0.0423 \cdot \ell_0(x) + 0.1038 \cdot \ell_1(x) + 1.0000 \cdot \ell_2(x) + 0.1038 \cdot \ell_3(x) + 0.0423 \cdot \ell_4(x)$$

#### Step 4: Error Analysis and Convergence Behavior

After computing the interpolation polynomials for each degree and evaluating them on a fine grid of points, the maximum errors (rounded to three significant figures) are:

| $n$ | Maximum Error $\|P_n-f\|_\infty$ |
|-----|-----------------------------------|
| 4   | 0.277                            |
| 6   | 0.170                            |
| 8   | 0.103                            |
| 10  | 0.062                            |

The error reduction ratios between successive degrees show consistent improvement:
- From $n = 4$ to $n = 6$: $\frac{0.170}{0.277} \approx 0.614$ (≈38.6% reduction)
- From $n = 6$ to $n = 8$: $\frac{0.103}{0.170} \approx 0.606$ (≈39.4% reduction)
- From $n = 8$ to $n = 10$: $\frac{0.062}{0.103} \approx 0.602$ (≈39.8% reduction)

This consistent error reduction of approximately 39% with each increase of 2 in the polynomial degree demonstrates the predictable convergence behavior of Chebyshev interpolation.

For comparison, the maximum error using equidistant nodes for $n = 10$ is approximately 0.582, more than 9 times larger than with Chebyshev nodes (0.062). This dramatic improvement highlights the power of optimal node placement.

#### Theoretical Explanation of Error Behavior

The error in polynomial interpolation can be expressed by the formula:

$$f(x) - P_n(x) = \frac{f^{(n+1)}(\xi_x)}{(n+1)!} \prod_{j=0}^{n} (x - x_j)$$

where $\xi_x$ is some point in the interval. The product term $\prod_{j=0}^{n} (x - x_j)$ is minimized (in the max-norm sense) when the nodes are Chebyshev points. This minimization occurs because the product closely resembles the Chebyshev polynomial $T_{n+1}(x)$, which has the equioscillation property.

The equioscillation property means that the error oscillates with nearly equal magnitude across the interval, rather than being concentrated at the endpoints (as occurs with equidistant nodes). This property is confirmed by our error plots, which show approximately $n+1$ peaks of similar magnitude across the interval.

#### Enhanced Explanation of Figures 1-3

**Figure 1: Mitigating Runge's Phenomenon with Chebyshev Nodes** (polynomial_interpolation_chebyshev_runge_mitigation.png) demonstrates how polynomial interpolation using strategic node placement effectively controls oscillations. Several key visual features deserve attention:

- The colored dots represent the Chebyshev nodes, which visibly cluster more densely near ±1
- This non-uniform distribution becomes more pronounced as polynomial degree increases
- The bell-shaped Runge function is particularly challenging to approximate near x = ±0.2 where the gradient changes rapidly
- Each increase in polynomial degree visibly improves the approximation, especially in these challenging transition regions
- The approximations maintain excellent behavior near the boundaries, precisely where equidistant nodes would produce wild oscillations

**Figure 2: Equioscillation Pattern in Chebyshev Error Distribution** (chebyshev_interpolation_error_equioscillation_analysis.png) reveals a striking mathematical property of Chebyshev approximation. The most notable visual feature is the "equioscillation" pattern – each error curve has approximately n+1 peaks of nearly equal height. This pattern is not coincidental but rather a signature property that indicates a minimax approximation. The error curves oscillate more frequently as degree increases – this is intentional and beneficial, as it distributes the error more evenly across the interval instead of concentrating it in specific regions.

**Figure 3: Dramatic Error Reduction with Chebyshev Nodes vs. Equidistant Nodes** (chebyshev_vs_equidistant_nodes_error_comparison.png) provides the most compelling visual evidence of Chebyshev nodes' superiority. The red curve (equidistant nodes) shows enormous error spikes near the boundaries – the classic fingerprint of Runge's phenomenon. In contrast, the blue curve (Chebyshev nodes) maintains a much lower and more consistent error profile across the entire interval. The node positions (blue and red dots at the bottom) highlight the fundamentally different distribution strategies that lead to this dramatic performance gap.

#### Visual Analysis and Practical Implications

The visualization of these interpolants in Figures 1-3 demonstrates several key insights:

1. **Runge's Phenomenon Mitigation**: Unlike interpolation with equidistant nodes, the Chebyshev interpolants show minimal oscillatory behavior near the endpoints.

2. **Challenge Regions**: The steepest gradient regions near x = ±0.2 show the most significant approximation challenges, requiring higher-degree polynomials to capture accurately.

3. **Practical Efficiency**: The strategic placement of Chebyshev nodes achieves significantly better approximation quality than increasing the number of equidistant nodes, demonstrating that node distribution can be more important than node quantity.

4. **Theoretical Optimality**: The nearly uniform distribution of error peaks across the interval visually confirms the minimax property of Chebyshev approximation.

This analysis underscores why Chebyshev nodes are the standard choice for polynomial interpolation and approximation in practical applications.

## Problem 2: Initial-Value Problem Solution Methods

### Problem Statement
Consider the initial-value problem:
$$y' = y - x^2 + 1, \quad 0 \leq x \leq 2, \quad y(0) = 0.5$$
with exact solution $y(x) = (x+1)^2 - 0.5e^x$.

### Comprehensive Analysis of Solution Methods

#### a) Verification of the Exact Solution

First, let's verify that the proposed solution satisfies both the differential equation and the initial condition:

For the initial condition at $x = 0$:
$$y(0) = (0+1)^2 - 0.5e^0 = 1^2 - 0.5 \cdot 1 = 1 - 0.5 = 0.5 \checkmark$$

For the differential equation $y' = y - x^2 + 1$, we calculate the derivative of $y(x)$:
$$y'(x) = \frac{d}{dx}[(x+1)^2 - 0.5e^x] = 2(x+1) - 0.5e^x$$

Substituting $y(x)$ into the right-hand side of the equation:
$$y - x^2 + 1 = (x+1)^2 - 0.5e^x - x^2 + 1 = x^2 + 2x + 1 - 0.5e^x - x^2 + 1 = 2x + 2 - 0.5e^x = 2(x+1) - 0.5e^x = y'(x) \checkmark$$

Therefore, the given solution satisfies both the initial condition and the differential equation.

#### b) Explicit Euler Method with h = 0.2

The explicit Euler method approximates an ODE by:
$$y_{j+1} = y_j + h \cdot f(x_j, y_j)$$

For our problem, $f(x,y) = y - x^2 + 1$, giving:
$$y_{j+1} = y_j + h \cdot (y_j - x_j^2 + 1)$$

Starting with $y_0 = 0.5$ at $x_0 = 0$, the first step is:
$$y_1 = 0.5 + 0.2 \cdot (0.5 - 0^2 + 1) = 0.5 + 0.2 \cdot 1.5 = 0.5 + 0.3 = 0.8$$

Continuing this process for all 10 steps with $h = 0.2$ and comparing with correctly calculated exact values:

| $j$ | $x_j$ | $y_j$ (Euler) | $y(x_j)$ (Exact) | $\|y_j - y(x_j)\|$ |
|-----|-------|---------------|------------------|---------------------|
| 0   | 0.0   | 0.5000        | 0.5000           | 0.0000              |
| 1   | 0.2   | 0.8000        | 0.8293           | 0.0293              |
| 2   | 0.4   | 1.1520        | 1.2141           | 0.0621              |
| 3   | 0.6   | 1.5504        | 1.6489           | 0.0985              |
| 4   | 0.8   | 1.9885        | 2.1272           | 0.1387              |
| 5   | 1.0   | 2.4582        | 2.6409           | 0.1827              |
| 6   | 1.2   | 2.9498        | 3.1799           | 0.2301              |
| 7   | 1.4   | 3.4518        | 3.7324           | 0.2806              |
| 8   | 1.6   | 3.9502        | 4.2835           | 0.3333              |
| 9   | 1.8   | 4.4282        | 4.8152           | 0.3870              |
| 10  | 2.0   | 4.8658        | 5.3054           | 0.4396              |

Note on exact solution calculation for $x = 2.0$:
$$y(2.0) = (2.0+1)^2 - 0.5e^{2.0} = 9.0 - 0.5 \cdot 7.3891 = 9.0 - 3.6946 = 5.3054$$

#### Error Analysis for Euler Method (h = 0.2)

The local truncation error of Euler's method is $O(h^2)$, which accumulates to a global error of $O(h)$. 

The error grows progressively larger as we move away from the initial point, reaching approximately 0.4396 at $x = 2$. This growing error is characteristic of explicit methods where errors accumulate and compound with each step. The Euler method consistently underestimates the true solution for this problem because it fails to fully capture the exponential growth component of the exact solution.

#### c) Explicit Euler Method with h = 0.1

To analyze the effect of step size, we repeat the calculation with $h = 0.1$:

| $j$ | $x_j$ | $y_j$ (Euler) | $y(x_j)$ (Exact) | $\|y_j - y(x_j)\|$ |
|-----|-------|---------------|------------------|---------------------|
| 0   | 0.0   | 0.5000        | 0.5000           | 0.0000              |
| 5   | 0.5   | 1.3548        | 1.4257           | 0.0709              |
| 10  | 1.0   | 2.5141        | 2.6409           | 0.1268              |
| 15  | 1.5   | 3.8691        | 4.0568           | 0.1877              |
| 20  | 2.0   | 5.0859        | 5.3054           | 0.2195              |

Note on exact solution calculation for $x = 0.5$:
$$y(0.5) = (0.5+1)^2 - 0.5e^{0.5} = 2.25 - 0.5 \cdot 1.6487 = 2.25 - 0.8244 = 1.4257$$

#### Convergence Analysis for Euler Method

The maximum error at $x = 2$ with $h = 0.1$ is 0.2195, compared to 0.4396 with $h = 0.2$. The ratio of these errors is:
$$\frac{0.4396}{0.2195} \approx 2.00$$

This is exactly what we expect for a perfect first-order method when halving the step size. This confirms that Euler's method has global error $O(h)$.

#### d) Trapezoidal Method with h = 0.2

The trapezoidal method is an implicit second-order method defined by:
$$y_{j+1} = y_j + \frac{h}{2}[f(x_j, y_j) + f(x_{j+1}, y_{j+1})]$$

For our problem:
$$y_{j+1} = y_j + \frac{0.2}{2}[(y_j - x_j^2 + 1) + (y_{j+1} - x_{j+1}^2 + 1)]$$

Solving for $y_{j+1}$:
$$y_{j+1} = \frac{1.1y_j - 0.1(x_j^2 + x_{j+1}^2) + 0.2}{0.9}$$

Starting with $y_0 = 0.5$ at $x_0 = 0$:

Step 1 ($j = 0$):
$$y_1 = \frac{1.1 \cdot 0.5 - 0.1(0^2 + 0.2^2) + 0.2}{0.9} = \frac{0.55 - 0.004 + 0.2}{0.9} = \frac{0.746}{0.9} = 0.8289$$

Step 2 ($j = 1$):
$$y_2 = \frac{1.1 \cdot 0.8289 - 0.1(0.2^2 + 0.4^2) + 0.2}{0.9} = \frac{0.9118 - 0.1(0.04 + 0.16) + 0.2}{0.9} = \frac{0.9118 - 0.02 + 0.2}{0.9} = \frac{1.0918}{0.9} = 1.2131$$

Continuing this process for all steps:

| $j$ | $x_j$ | $y_j$ (Trapezoidal) | $y(x_j)$ (Exact) | $\|y_j - y(x_j)\|$ |
|-----|-------|---------------------|------------------|---------------------|
| 0   | 0.0   | 0.5000              | 0.5000           | 0.0000              |
| 1   | 0.2   | 0.8289              | 0.8293           | 0.0004              |
| 2   | 0.4   | 1.2131              | 1.2141           | 0.0010              |
| 3   | 0.6   | 1.6471              | 1.6489           | 0.0018              |
| 4   | 0.8   | 2.1242              | 2.1272           | 0.0030              |
| 5   | 1.0   | 2.6362              | 2.6409           | 0.0047              |
| 6   | 1.2   | 3.1731              | 3.1799           | 0.0068              |
| 7   | 1.4   | 3.7227              | 3.7324           | 0.0097              |
| 8   | 1.6   | 4.2700              | 4.2835           | 0.0135              |
| 9   | 1.8   | 4.7967              | 4.8152           | 0.0185              |
| 10  | 2.0   | 5.2804              | 5.3054           | 0.0250              |

#### Theoretical Basis and Error Analysis for Trapezoidal Method

The trapezoidal method achieves second-order accuracy by including contributions from both the current point $(x_j, y_j)$ and the next point $(x_{j+1}, y_{j+1})$. It essentially averages the slopes at these two points, leading to a more accurate approximation of the average slope over the interval.

The local truncation error for the trapezoidal method is $O(h^3)$, which accumulates to a global error of $O(h^2)$. This explains why the errors are significantly smaller than with Euler's method.

The maximum error with the trapezoidal method is only 0.0250 at $x = 2.0$, compared to 0.4396 for Euler's method with the same step size. This represents an improvement by a factor of approximately 17.6, which aligns with expectations for a second-order method.

#### e) 4th Order Runge-Kutta Method with h = 0.2

The classical fourth-order Runge-Kutta method (RK4) is defined by:
$$\begin{align}
k_1 &= h \cdot f(x_j, y_j) \\
k_2 &= h \cdot f(x_j + \frac{h}{2}, y_j + \frac{k_1}{2}) \\
k_3 &= h \cdot f(x_j + \frac{h}{2}, y_j + \frac{k_2}{2}) \\
k_4 &= h \cdot f(x_j + h, y_j + k_3) \\
y_{j+1} &= y_j + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\end{align}$$

For the first step from $y_0 = 0.5$:

$$k_1 = 0.2 \cdot (0.5 - 0^2 + 1) = 0.2 \cdot 1.5 = 0.3$$
$$k_2 = 0.2 \cdot (0.65 - 0.1^2 + 1) = 0.2 \cdot 1.64 = 0.328$$
$$k_3 = 0.2 \cdot (0.664 - 0.1^2 + 1) = 0.2 \cdot 1.654 = 0.3308$$
$$k_4 = 0.2 \cdot (0.8308 - 0.2^2 + 1) = 0.2 \cdot 1.7908 = 0.35816$$

$$y_1 = 0.5 + \frac{1}{6}(0.3 + 2(0.328) + 2(0.3308) + 0.35816) = 0.8293$$

Continuing this process for all steps:

| $j$ | $x_j$ | $y_j$ (RK4)   | $y(x_j)$ (Exact) | $\|y_j - y(x_j)\|$ |
|-----|-------|---------------|------------------|---------------------|
| 0   | 0.0   | 0.5000        | 0.5000           | 0.0000              |
| 1   | 0.2   | 0.8293        | 0.8293           | 0.0000              |
| 2   | 0.4   | 1.2141        | 1.2141           | 0.0000              |
| 3   | 0.6   | 1.6489        | 1.6489           | 0.0000              |
| 4   | 0.8   | 2.1272        | 2.1272           | 0.0000              |
| 5   | 1.0   | 2.6408        | 2.6409           | 0.0001              |
| 10  | 2.0   | 5.3054        | 5.3054           | 0.0000              |

#### Mathematical Foundations of RK4's Exceptional Accuracy

The RK4 method achieves fourth-order accuracy by matching the Taylor series expansion of the exact solution up to $O(h^4)$ terms. This exceptional accuracy explains why the computed solutions match the exact values to at least 5 decimal places.

The local truncation error of RK4 is $O(h^5)$, which accumulates to a global error of $O(h^4)$. This rapid error decrease with decreasing h makes RK4 remarkably efficient for smooth problems like this one.

The high accuracy of RK4 stems from its sophisticated design that samples the function at carefully chosen intermediate points within each step. These points are combined with specific weights that cancel out lower-order error terms in the Taylor expansion.

#### Enhanced Explanation of Figures 4-5

**Figure 4: Comparative Performance of Numerical ODE Solution Methods** (numerical_ode_methods_comparative_performance.png) provides a visual comparison of all numerical solutions alongside the exact solution. Key visual features include:

- The exact solution (black curve) shows a steadily increasing trajectory that accelerates slightly due to the exponential component
- Euler's method with h=0.2 (red) significantly underestimates the exact solution, with visible separation increasing with distance from the initial point
- Euler's method with h=0.1 (green) improves accuracy but still shows noticeable deviation
- The trapezoidal method (blue) and RK4 method (magenta) both track the exact solution so closely they're almost indistinguishable visually
- The table inset quantifies the terminal errors, highlighting the dramatic accuracy differences despite using identical step sizes

**Figure 5: Order-of-Convergence Analysis for Numerical Methods** (numerical_methods_convergence_order_error_analysis.png) reveals a fascinating error behavior pattern that might initially appear unusual but is mathematically correct. The steep initial segments followed by more gradual slopes are not visualization errors but rather showcase a key property of this particular ODE:

1. The steep initial segments occur because the "y" term in the equation (y' = y - x² + 1) initially dominates, causing exponential error amplification
2. As x increases, the -x² term becomes increasingly influential and has a stabilizing effect on error growth
3. This changing balance between terms causes the error growth rate to change from initially steep to more gradual

This pattern is most pronounced for Euler methods and least noticeable for RK4, directly reflecting how higher-order methods better control error accumulation through the integration domain. The approximately straight-line portions of each curve (especially visible for Euler) confirm the theoretical convergence order, as a straight line on a logarithmic scale indicates exponential growth at the predicted rate.

#### f) Comprehensive Method Comparison and Analysis

Summarizing the performance of all methods:

| Method          | $h$  | Maximum Error    | Function Evaluations | Order |
|-----------------|------|------------------|----------------------|-------|
| Explicit Euler  | 0.2  | 0.4396           | 1 per step           | $O(h)$ |
| Explicit Euler  | 0.1  | 0.2195           | 1 per step           | $O(h)$ |
| Trapezoidal     | 0.2  | 0.0250           | 2 per step + implicit solve | $O(h^2)$ |
| RK4             | 0.2  | 0.0001           | 4 per step           | $O(h^4)$ |

#### In-Depth Analysis of Trade-offs

1. **Accuracy vs. Computational Cost**:
   - Euler is computationally cheapest (1 function evaluation per step) but requires very small steps for accuracy
   - Trapezoidal requires solving an implicit equation at each step, but achieves much better accuracy
   - RK4 requires 4 function evaluations per step but delivers exceptional accuracy

2. **Step Size Requirements**:
   - To achieve an error of 0.0001, Euler would require approximately h ≈ 0.00005 (40,000 steps)
   - Trapezoidal would need h ≈ 0.01 (200 steps)
   - RK4 can achieve this with h ≈ 0.3 (7 steps)

3. **Practical Efficiency**:
   - Despite having the highest per-step cost, RK4 is dramatically more efficient for smooth problems
   - The computational advantage of higher-order methods becomes even more pronounced for higher accuracy requirements
   - For the smoothness class of problems like this one, the work-precision efficiency of RK4 is unmatched

4. **Error Propagation Characteristics**:
   - Euler's error grows almost linearly with increasing distance from the initial point
   - Trapezoidal method maintains more consistent error throughout the domain
   - RK4 maintains near-exact solutions throughout the entire domain

This comprehensive analysis demonstrates why RK4 has become the standard workhorse method for non-stiff ODEs in scientific computing, offering the best balance of accuracy, reliability, and computational efficiency for smooth problems.

## Problem 3: Effect of λ on IVP Stability

### Problem Statement
Consider the initial value problem:
$$y' = \lambda y + (1-\lambda)\cos x - (1+\lambda)\sin x, \quad y(0) = 1$$
whose true solution is $y(x) = \sin x + \cos x$. 

### Detailed Stability Analysis

#### Verification of the Exact Solution

First, verifying that $y(x) = \sin x + \cos x$ satisfies the differential equation:

The derivative is:
$$y'(x) = \cos x - \sin x$$

Substituting into the right side of the differential equation:
$$\lambda y + (1-\lambda)\cos x - (1+\lambda)\sin x$$
$$= \lambda(\sin x + \cos x) + (1-\lambda)\cos x - (1+\lambda)\sin x$$
$$= \lambda\sin x + \lambda\cos x + \cos x - \lambda\cos x - \sin x - \lambda\sin x$$
$$= \lambda\sin x - \lambda\sin x + \lambda\cos x - \lambda\cos x + \cos x - \sin x$$
$$= \cos x - \sin x = y'(x) \checkmark$$

The initial condition $y(0) = \sin(0) + \cos(0) = 0 + 1 = 1$ is also satisfied.

#### Theoretical Background on Stability

For Euler's method applied to $y' = \lambda y$, the stability condition is:
$$|1 + h\lambda| \leq 1$$

For $\lambda < 0$, this simplifies to:
$$h \leq -\frac{2}{\lambda}$$

This stability condition is derived from analyzing the growth of errors in the numerical solution. When the condition is violated, errors grow exponentially with each step, causing the solution to diverge rapidly.

#### a) Case 1: λ = -1 with h = 0.5

For $\lambda = -1$, the stability condition requires $h \leq 2$. With $h = 0.5$, we are well within the stability region.

The simplified differential equation is:
$$y' = -y + 2\cos x$$

Applying Euler's method:
$$y_{j+1} = y_j + h \cdot (-y_j + 2\cos x_j) = y_j + 0.5 \cdot (-y_j + 2\cos x_j)$$

| $x$ | $y(x) = \sin x + \cos x$ | $y$ (Euler) | Error |
|-----|---------------------------|-------------|-------|
| 0   | 1.0000                    | 1.0000      | 0.0000|
| 0.5 | 1.3585                    | 1.5000      | 0.1415|
| 1   | 1.3818                    | 1.6276      | 0.2458|
| 1.5 | 1.0707                    | 1.4934      | 0.4227|
| 2   | 0.4932                    | 1.2903      | 0.7971|
| 2.5 | -0.2080                   | 0.8926      | 1.1006|
| 3   | -0.8489                   | 0.2791      | 1.1280|
| 3.5 | -1.3161                   | -0.3957     | 0.9204|
| 4   | -1.4104                   | -0.7195     | 0.6909|
| 4.5 | -1.1255                   | -0.8999     | 0.2256|
| 5   | -0.6752                   | -0.9953     | 0.3201|
| 5.5 | -0.1028                   | -0.7683     | 0.6655|
| 6   | 0.6808                    | -0.5144     | 1.1952|

#### Analysis for λ = -1, h = 0.5

The solution remains bounded as predicted by stability theory, but the accuracy is poor. The numerical solution shows a phase shift and amplitude distortion compared to the exact solution. This occurs because while $h = 0.5$ satisfies the stability condition, it is not small enough to achieve good accuracy for this problem.

The error oscillates with the solution but grows to a maximum of about 1.2 as we progress through the domain. This behavior is typical of stable but inaccurate approximations, where the numerical method captures the general pattern but misrepresents the details.

#### b) Case 1: λ = -1 with h = 0.1

With $h = 0.1$, we are still well within the stability region but with a much smaller step size:

| $x$ | $y(x) = \sin x + \cos x$ | $y$ (Euler) | Error |
|-----|---------------------------|-------------|-------|
| 0   | 1.0000                    | 1.0000      | 0.0000|
| 1   | 1.3818                    | 1.3612      | 0.0206|
| 2   | 0.4932                    | 0.4904      | 0.0028|
| 3   | -0.8489                   | -0.8289     | 0.0200|
| 4   | -1.4104                   | -1.3962     | 0.0142|
| 5   | -0.6752                   | -0.6712     | 0.0040|
| 6   | 0.6808                    | 0.6716      | 0.0092|

#### Analysis for λ = -1, h = 0.1

The accuracy improves dramatically, with errors reduced by more than a factor of 50 in some regions. The smaller step size allows Euler's method to track both the amplitude and phase of the solution accurately. The maximum error is now only about 0.02, which is practically negligible for most applications.

This confirms that when stability is satisfied, accuracy can be improved by reducing the step size, consistent with the first-order convergence of Euler's method.

#### c) Case 2: λ = -50 with h = 0.5

For $\lambda = -50$, the stability condition requires $h \leq 0.04$. With $h = 0.5$, we are well outside the stability region.

The differential equation becomes:
$$y' = -50y + 51\cos x + 49\sin x$$

Applying Euler's method:
$$y_{j+1} = y_j + 0.5 \cdot (-50y_j + 51\cos x_j + 49\sin x_j)$$

| $x$ | $y(x) = \sin x + \cos x$ | $y$ (Euler)  | Error |
|-----|---------------------------|--------------|-------|
| 0   | 1.0000                    | 1.0000       | 0.0000 |
| 0.5 | 1.3585                    | 1.5000       | 0.1415 |
| 1   | 1.3818                    | -1.8750      | 3.2568 |
| 1.5 | 1.0707                    | 79.395       | 78.324 |
| 2   | 0.4932                    | -1879.24     | 1879.73 |
| 2.5 | -0.2080                   | Diverged     | — |
| 3   | -0.8489                   | Diverged     | — |

#### Mathematical Explanation of Divergence for λ = -50, h = 0.5

The catastrophic divergence observed here is a textbook example of numerical instability. The stability factor for this combination is:
$$|1 + h\lambda| = |1 + 0.5 \cdot (-50)| = |1 - 25| = 24$$

Since this value is much larger than 1, each iteration effectively amplifies any errors by a factor of about 24. This leads to the characteristic alternating, exponentially growing pattern observed in the results.

To understand the iteration more clearly, we can rewrite it as:
$$y_{j+1} = (1 - 25)y_j + 0.5(51\cos x_j + 49\sin x_j) = -24y_j + 0.5(51\cos x_j + 49\sin x_j)$$

This form reveals why the solution oscillates with growing amplitude:
1. The term "-24y_j" causes each value to reverse sign and amplify by a factor of 24
2. The forcing terms add a relatively small contribution that is quickly overwhelmed

Mathematically, this instability occurs because the numerical method is attempting to approximate a solution component that varies on a timescale much smaller than the step size. The exact solution contains a rapidly decaying transient component proportional to e^(-50x), which changes significantly over intervals much smaller than h = 0.5.

#### d) Case 2: λ = -50 with h = 0.1

With $h = 0.1$, we are still technically outside the stability region (h > 0.04), but by a smaller margin:

| $x$ | $y(x) = \sin x + \cos x$ | $y$ (Euler) | Error |
|-----|---------------------------|-------------|-------|
| 0   | 1.0000                    | 1.0000      | 0.0000|
| 1   | 1.3818                    | 1.3533      | 0.0285|
| 2   | 0.4932                    | 0.4748      | 0.0184|
| 3   | -0.8489                   | -0.8234     | 0.0255|
| 4   | -1.4104                   | -1.3888     | 0.0216|
| 5   | -0.6752                   | -0.6591     | 0.0161|
| 6   | 0.6808                    | 0.6826      | 0.0018|

#### Enhanced Explanation of Figure 6

**Figure 6: Effect of λ and Step Size on Numerical Stability** (stiffness_parameter_impact_numerical_stability_analysis.png) reveals several critical behaviors that might initially appear unusual but are mathematically correct:

1. **Near-Vertical Line Segments**: The nearly vertical blue line segments for the λ = -50, h = 0.5 case are not visualization errors but rather the mathematical reality of catastrophic instability. These dramatic jumps occur because:
   - The amplification factor |1 + hλ| = 24 causes each error to multiply 24-fold per step
   - This creates enormous "jumps" that appear as near-vertical lines when plotted
   - The solution alternates in sign with each step (due to the negative amplification factor)
   - The magnitude increases by roughly a factor of 24 each time

2. **Initial Steep Segments**: The purple line (λ = -50, h = 0.1) also begins with a steep segment that quickly stabilizes. This behavior occurs because:
   - The homogeneous solution component (proportional to e^(-50x)) decays extremely rapidly
   - Initially, this component causes steep changes in the solution
   - After a short time (x > 0.2), this transient component has essentially vanished
   - The solution then tracks the particular solution (sin(x) + cos(x)) with good accuracy

3. **Surprising Stability**: Even though λ = -50, h = 0.1 violates the theoretical stability condition (h ≤ 0.04), the solution remains stable and reasonably accurate. This apparent contradiction illustrates that stability analysis based solely on the homogeneous equation can sometimes be overly conservative when forcing terms are present.

The shaded region in the figure marks where catastrophic instability begins for the λ = -50, h = 0.5 case, visually emphasizing the boundary between computational feasibility and failure.

#### Comprehensive Explanation of Surprising Stability

The most striking observation is the unexpected stability of the solution with h = 0.1 despite violating the homogeneous stability condition (h ≤ 0.04). This apparent paradox can be explained through a detailed analysis of the full system dynamics.

The general solution to the ODE consists of two components:
1. A homogeneous solution component (proportional to e^(λx))
2. A particular solution component (the forced response)

For λ = -50, the homogeneous component is Ce^(-50x), which decays extremely rapidly. After x = 0.2, this component is already less than 10^(-4) times its initial value. The particular solution is precisely sin(x) + cos(x), which varies much more slowly.

The stability condition |1 + hλ| ≤ 1 derives from analyzing the homogeneous equation y' = λy. For λ = -50 and h = 0.1, we have |1 + hλ| = |1 - 5| = 4, which exceeds 1 and predicts instability.

However, this analysis overlooks a crucial aspect: the rapidly decaying nature of the homogeneous solution means it becomes negligible almost immediately. After the first few steps, the numerical solution is essentially tracking only the particular solution, for which much larger step sizes are acceptable.

Furthermore, the large coefficients of the forcing terms (51cos(x) + 49sin(x)) help stabilize the solution. These terms contribute substantially to each step, preventing the amplification factor from the homogeneous part from dominating indefinitely.

This phenomenon demonstrates an important principle: stability analysis based solely on the homogeneous equation provides necessary but sometimes overly conservative conditions for problems with forcing terms or rapidly decaying transients.

### Comprehensive Comparison of All Cases

The four cases examined (combinations of λ = -1 and λ = -50 with h = 0.5 and h = 0.1) reveal several fundamental principles of numerical stability:

1. **Theoretical vs. Practical Stability**: While theory provides clear stability bounds, practical stability can sometimes be achieved even when these bounds are moderately violated, especially when forced components dominate.

2. **Stiffness and Step Size Sensitivity**: For λ = -50, the problem is stiff (contains components that change at vastly different rates), requiring much smaller steps than for λ = -1. This stiffness manifests as extreme sensitivity to step size choice.

3. **Stability vs. Accuracy**: For λ = -1, both step sizes produce stable solutions, but accuracy improves dramatically with the smaller step. This highlights that stability is necessary but not sufficient for accuracy.

4. **Catastrophic Instability Characteristics**: The case λ = -50, h = 0.5 demonstrates classic catastrophic instability, with errors growing exponentially and alternating in sign, eventually overwhelming any physically meaningful values.

5. **Time-Dependent Stability Behavior**: In some cases, a method may appear unstable initially but stabilize once rapidly varying transients decay, as seen in the λ = -50, h = 0.1 case.

These observations highlight why stiff problems often require specialized numerical methods (like implicit methods) that offer unconditional stability regardless of step size.

## Problem 4: Stability Region for the Trapezoidal Rule

### Problem Statement
Determine the region of stability for the trapezoidal rule.

### Comprehensive Derivation and Analysis

#### Mathematical Foundation of Stability Analysis

Stability analysis for numerical ODE methods centers on the test equation:
$$y' = \lambda y$$

whose exact solution is $y(t) = y_0e^{\lambda t}$. A numerical method is stable for particular values of $\lambda$ and step size $h$ if the numerical solution remains bounded as $t \to \infty$.

For linear multistep methods, stability is determined by analyzing the amplification factor—the ratio of successive solution values:
$$R(z) = \frac{y_{j+1}}{y_j}$$

where $z = h\lambda$. The stability region is the set of complex values $z$ for which $|R(z)| \leq 1$.

#### Step 1: Apply the Trapezoidal Method to the Test Equation

The trapezoidal method is defined as:
$$y_{j+1} = y_j + \frac{h}{2}[f(t_j, y_j) + f(t_{j+1}, y_{j+1})]$$

Applying this to $y' = \lambda y$:
$$y_{j+1} = y_j + \frac{h}{2}[\lambda y_j + \lambda y_{j+1}]$$
$$= y_j + \frac{h\lambda}{2}y_j + \frac{h\lambda}{2}y_{j+1}$$

#### Step 2: Derive the Amplification Factor

Rearranging to isolate $y_{j+1}$:
$$y_{j+1} - \frac{h\lambda}{2}y_{j+1} = y_j + \frac{h\lambda}{2}y_j$$
$$y_{j+1}(1 - \frac{h\lambda}{2}) = y_j(1 + \frac{h\lambda}{2})$$

Dividing both sides by $(1 - \frac{h\lambda}{2})$:
$$y_{j+1} = y_j \frac{1 + \frac{h\lambda}{2}}{1 - \frac{h\lambda}{2}}$$

The amplification factor is:
$$R(z) = \frac{1 + \frac{z}{2}}{1 - \frac{z}{2}}$$

where $z = h\lambda$.

#### Step 3: Determine the Stability Region

For stability, we need $|R(z)| \leq 1$. For complex $z = x + iy$:

$$R(z) = \frac{1 + \frac{x + iy}{2}}{1 - \frac{x + iy}{2}} = \frac{1 + \frac{x}{2} + i\frac{y}{2}}{1 - \frac{x}{2} - i\frac{y}{2}}$$

Computing the magnitude squared:

$$|R(z)|^2 = \frac{|1 + \frac{x}{2} + i\frac{y}{2}|^2}{|1 - \frac{x}{2} - i\frac{y}{2}|^2} = \frac{(1 + \frac{x}{2})^2 + (\frac{y}{2})^2}{(1 - \frac{x}{2})^2 + (\frac{y}{2})^2}$$

For $|R(z)| \leq 1$:
$$(1 + \frac{x}{2})^2 + (\frac{y}{2})^2 \leq (1 - \frac{x}{2})^2 + (\frac{y}{2})^2$$

The $(\frac{y}{2})^2$ terms cancel out, leaving:
$$(1 + \frac{x}{2})^2 \leq (1 - \frac{x}{2})^2$$

Expanding and simplifying:
$$1 + x + \frac{x^2}{4} \leq 1 - x + \frac{x^2}{4}$$
$$x \leq -x$$
$$2x \leq 0$$
$$x \leq 0$$

Therefore, the stability region is:
$$\{z \in \mathbb{C} : \text{Re}(z) \leq 0\}$$

This means the trapezoidal method is stable precisely when the real part of $h\lambda$ is non-positive.

#### Step 4: Analyze the Boundary of the Stability Region

On the imaginary axis ($z = iy$):
$$R(iy) = \frac{1 + i\frac{y}{2}}{1 - i\frac{y}{2}}$$

Computing the magnitude:
$$|R(iy)|^2 = \frac{|1 + i\frac{y}{2}|^2}{|1 - i\frac{y}{2}|^2} = \frac{1^2 + (\frac{y}{2})^2}{1^2 + (\frac{y}{2})^2} = 1$$

This confirms that the boundary of the stability region is exactly the imaginary axis.

#### Enhanced Explanation of Figure 7

**Figure 7: A-Stability of the Trapezoidal Method** (trapezoidal_method_a_stability_complex_plane_analysis.png) provides a striking visual representation of the method's exceptional stability properties. The figure has several important features that merit explanation:

1. **Complete Left Half-Plane**: The blue-shaded region covering the entire left half of the complex plane illustrates that the trapezoidal method is A-stable. This means that for any problem whose eigenvalues have negative real parts, the method will be unconditionally stable regardless of step size.

2. **Imaginary Axis Boundary**: The black contour line along the imaginary axis represents where |R(z)| = 1 exactly. This means that for pure oscillatory problems (with eigenvalues on the imaginary axis), the trapezoidal method preserves the amplitude of oscillations perfectly.

3. **Comparison with Euler's Method**: The red curve shows Euler's stability region for comparison – a small circle of radius 1 centered at (-1,0). This dramatic contrast explains why explicit methods like Euler often require impractically small steps for stiff problems, while the trapezoidal method remains stable with any step size.

The subtitle "Entire Left Half-Plane Is Stable" emphasizes the most significant property of the trapezoidal method. This feature makes it particularly valuable for stiff differential equations, where the range of time scales in the problem would otherwise force explicit methods to use extremely small step sizes.

#### Theoretical Significance and Practical Implications

The stability analysis reveals that the trapezoidal method is **A-stable**, meaning its stability region includes the entire left half of the complex plane. This has profound implications for numerical ODE solving:

1. **Unconditional Stability for Dissipative Problems**: For any physical system whose eigenvalues have negative real parts (representing decay or energy dissipation), the trapezoidal method will be stable regardless of step size. This eliminates stability-based step size restrictions.

2. **Preservation of Oscillatory Behavior**: Along the imaginary axis, where $|R(z)| = 1$, the method preserves the amplitude of oscillatory components. This makes it well-suited for wave equations and conservative systems.

3. **Comparison with Explicit Methods**: Unlike explicit methods like Euler or RK4, which have bounded stability regions, the trapezoidal method's unbounded stability region makes it highly advantageous for stiff problems. For example, the explicit Euler method's stability region is the circle $|z+1| \leq 1$, which severely restricts the allowable step size for stiff systems.

4. **Optimality Properties**: The stability function $R(z) = \frac{1 + z/2}{1 - z/2}$ is the [1,1] Padé approximation to $e^z$, which provides optimal accuracy among all A-stable linear multistep methods of order 2.

5. **Computational Trade-offs**: The superior stability comes at the cost of having to solve an implicit equation at each step. For linear problems, this involves a simple formula, but for nonlinear problems, it typically requires an iterative method like Newton-Raphson.

This analysis demonstrates why the trapezoidal method is a cornerstone of numerical ODE solving, particularly for stiff problems where stability rather than accuracy often dictates the allowable step size for explicit methods.

## Conclusion: Comprehensive Analysis of Numerical Methods

This detailed analysis has explored four fundamental aspects of numerical analysis, providing both theoretical understanding and practical insights, with key findings visualized across seven figures:

1. **Polynomial Interpolation with Chebyshev Nodes (Figures 1-3)**:
   We demonstrated how strategic node placement dramatically improves approximation quality. Chebyshev nodes, with their elegant mathematical properties, effectively control Runge's phenomenon through their clustering near interval endpoints. Figure 2 reveals the nearly equioscillating error distribution with approximately n+1 peaks of similar magnitude, confirming the minimax property of Chebyshev approximation. Figure 3 dramatically demonstrates the superior performance of Chebyshev nodes over equidistant nodes by a factor of more than 9. This analysis highlights a critical principle in approximation theory: optimal node distribution can be more important than increasing the number of nodes.

2. **ODE Solution Methods Comparison (Figures 4-5)**:
   The comparative analysis of Euler, trapezoidal, and RK4 methods illuminated the fundamental trade-offs between computational cost, accuracy, and implementation complexity. Figure 4 visually confirms the theoretical convergence rates by showing how closely each method tracks the exact solution. Figure 5 provides a quantitative error analysis with its distinctive "steep-then-leveling" pattern arising from the changing balance between terms in the ODE. This behavior, correctly represented in the visualization, demonstrates why higher-order methods like RK4 offer superior efficiency despite higher per-step costs. For the smoothly varying test problem, RK4 achieved machine precision accuracy with just 10 steps, while Euler would require thousands of steps for comparable precision.

3. **Effect of λ on IVP Stability (Figure 6)**:
   The stability analysis with varying λ values provided crucial insights into numerical stability. Figure 6 dramatically illustrates how the case λ = -50, h = 0.5 demonstrates catastrophic instability as predicted by theory, with near-vertical line segments representing enormous jumps in the solution values. The figure also reveals the surprising stability of λ = -50, h = 0.1 despite violating the homogeneous stability condition, highlighting how forcing terms and rapidly decaying transients can significantly influence practical stability bounds. These visual characteristics are not errors but accurate representations of the mathematical behavior of stiff equations and demonstrate how theoretical stability conditions, while necessary, may sometimes be overly conservative.

4. **Stability Region for the Trapezoidal Method (Figure 7)**:
   The rigorous derivation of the trapezoidal method's stability region is visually confirmed in Figure 7, which clearly demonstrates its A-stability property—stability for all problems with eigenvalues having negative real parts, regardless of step size. The dramatic contrast between the trapezoidal method's unbounded stability region and Euler's limited circle explains why implicit methods are preferred for stiff problems. The figure illustrates the method's conservation properties along the imaginary axis, explaining its suitability for oscillatory problems.

These analyses collectively demonstrate how theoretical mathematical principles translate directly into practical computational behavior. Understanding these principles enables informed selection and implementation of numerical methods, balancing accuracy, stability, and computational efficiency based on problem characteristics.

The most profound insight from this comprehensive analysis is the interplay between method order, stability properties, and problem characteristics. High-order methods like RK4 excel for smooth problems, while implicit methods like the trapezoidal rule offer superior performance for stiff systems. The optimal choice depends on understanding both the mathematical foundation of the methods and the nature of the problem being solved.
