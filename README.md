# Numerical Analysis Homework | Marybeth Brauns

## Problem 1: Evaluating \(\displaystyle \int_{0}^{\pi} e^x \cos x\,dx\)

### Problem Statement

Evaluate the integral

\[
I(f)=\int_{0}^{\pi} e^x \cos x\,dx
\]

using the following approaches:

1. **(a) Analytical Evaluation:** Use the Fundamental Theorem of Calculus (via integration by parts) to obtain an exact result.
2. **(b) Gauss–Legendre Quadrature:** Approximate the integral using Gauss–Legendre quadrature for \(n=2,3,4,5,6\) nodes. Compute the approximations \(I_n(f)\) and the corresponding errors \(E_n(f)\).
3. **(c) Discussion:** Discuss the accuracy of Gauss–Legendre quadrature in comparison to Newton–Cotes formulas and comment on the computational effort involved.

### (a) Analytical Evaluation

To evaluate
\[
I(f)=\int_{0}^{\pi} e^x \cos x\,dx,
\]
we use integration by parts twice.

1. **First integration by parts:**  
   Let
   - \(u = e^x\)  so \(du = e^x\,dx\),
   - \(dv = \cos x\,dx\)  so \(v = \sin x\).

   Then,
   \[
   \int e^x\cos x\,dx = e^x \sin x - \int e^x \sin x\,dx.
   \]

2. **Second integration by parts:**  
   For the remaining integral, let
   - \(u = e^x\)  so \(du = e^x\,dx\),
   - \(dv = \sin x\,dx\)  so \(v = -\cos x\).

   Then,
   \[
   \int e^x \sin x\,dx = -e^x\cos x + \int e^x\cos x\,dx.
   \]

Substitute the second result back into the first:
\[
\begin{aligned}
\int e^x\cos x\,dx &= e^x \sin x - \Bigl(-e^x\cos x + \int e^x\cos x\,dx\Bigr)\\[1mm]
&= e^x (\sin x + \cos x) - \int e^x\cos x\,dx.
\end{aligned}
\]

Adding \(\int e^x\cos x\,dx\) to both sides:
\[
2\int e^x\cos x\,dx = e^x (\sin x + \cos x),
\]
so that
\[
\int e^x\cos x\,dx = \frac{e^x (\sin x + \cos x)}{2}.
\]

Now evaluate the definite integral:
\[
I(f)=\left[\frac{e^x (\sin x + \cos x)}{2}\right]_0^{\pi}.
\]

At the limits:

- **At \(x=\pi\):**  
  \(\sin \pi=0\) and \(\cos \pi=-1\) so
  \[
  \frac{e^\pi (0-1)}{2} = -\frac{e^\pi}{2}.
  \]
- **At \(x=0\):**  
  \(\sin 0=0\) and \(\cos 0=1\) so
  \[
  \frac{e^0 (0+1)}{2} = \frac{1}{2}.
  \]

Thus,
\[
I(f) = \left(-\frac{e^\pi}{2}\right) - \frac{1}{2} = -\frac{e^\pi+1}{2}.
\]

Numerically, with \(e^\pi \approx 23.1407\),
\[
I(f) \approx -\frac{24.1407}{2} \approx -12.07035.
\]

**Observations on the Analytical Solution:**
- Integration by parts neatly handles the combination of exponential and trigonometric functions.
- The large magnitude arises due to the factor \(e^\pi\).
- This exact value serves as our benchmark for assessing numerical approximations.

---

### (b) Gauss–Legendre Quadrature

To apply Gauss–Legendre quadrature on the interval \([0,\pi]\), we transform the interval to \([-1,1]\) by letting

\[
x = \frac{\pi}{2}(t+1), \quad dx = \frac{\pi}{2}\,dt.
\]

The transformed integrand becomes

\[
f(t)=\frac{\pi}{2}\,e^{\frac{\pi}{2}(t+1)}\cos\!\Bigl(\frac{\pi}{2}(t+1)\Bigr).
\]

The \(n\)-point Gauss–Legendre quadrature approximation is given by

\[
I_n(f)=\sum_{i=1}^{n} w_i\,f(t_i),
\]

where the \(t_i\) are the nodes and \(w_i\) the corresponding weights.

**Results for different values of \(n\):**

| \(n\) | \(I_n(f)\)      | Absolute Error       | Relative Error |
|:-----:|:---------------:|:--------------------:|:--------------:|
| 2     | \(-12.33621047\)| \(2.66\times10^{-1}\)| 2.20%          |
| 3     | \(-12.12742045\)| \(5.71\times10^{-2}\)| 0.47%          |
| 4     | \(-12.07018949\)| \(1.57\times10^{-4}\)| 0.0013%        |
| 5     | \(-12.07032854\)| \(1.78\times10^{-5}\)| 0.00015%       |
| 6     | \(-12.07034633\)| \(1.47\times10^{-8}\)| 0.00000012%    |

**Observations on Gauss–Legendre Convergence:**
- The error decreases by roughly one to two orders of magnitude with each additional node.
- With \(n=4\) nodes, about 4–5 significant digits of accuracy are achieved.
- The rapid (exponential) convergence is typical for analytic integrands such as \(e^x\cos x\).

---

### (c) Discussion

The Gauss–Legendre quadrature method demonstrates impressive efficiency:

- **Rapid Convergence:** The error drops from roughly \(2.66\times10^{-1}\) (with 2 nodes) to \(1.47\times10^{-8}\) (with 6 nodes).
- **Comparison with Newton–Cotes:** For this problem, methods using equally spaced nodes (e.g., Simpson's or Boole's rule) would require many more subintervals to reach the same accuracy.
- **Computational Efficiency:** Although Gauss–Legendre requires precomputed nodes and weights, its optimal node placement leads to a high degree of precision with few function evaluations.
- **Error Behavior:** The error is tied to the \((2n)\)th derivative of the integrand; for our smooth function, this produces an exponential error decrease.
- **Handling Challenging Integrands:** Even though \(e^x\cos x\) combines exponential growth with oscillations, Gauss–Legendre quadrature captures the behavior extremely well.

---

## Problem 2: Integration of \(\displaystyle \int_{-1}^{1}\frac{1}{1+25x^2}\,dx\)

### Problem Statement

(a) **Analytical Evaluation:** Evaluate

\[
\int_{-1}^{1}\frac{1}{1+25x^2}\,dx
\]

using the Fundamental Theorem of Calculus.

(b) **Global Polynomial Interpolation via Newton's Divided Differences:** Approximate the integral by exactly integrating the interpolating polynomial using:
   - A linear polynomial (2 nodes),
   - A quadratic polynomial (3 nodes),
   - A quartic polynomial (5 nodes).

(c) **Gauss–Legendre Quadrature:** Approximate the same integral using Gauss–Legendre quadrature with 2, 4, and 6 nodes.

### (a) Analytical Evaluation

Let \(u = 5x\); then \(du = 5\,dx\) or \(dx = \frac{du}{5}\). Changing limits:

- When \(x=-1\), \(u=-5\); when \(x=1\), \(u=5\).

Thus, the integral becomes

\[
\begin{aligned}
\int_{-1}^{1} \frac{1}{1+25x^2}\,dx 
&= \int_{-5}^{5} \frac{1}{1+u^2}\cdot\frac{du}{5} \\
&= \frac{1}{5}\int_{-5}^{5} \frac{du}{1+u^2}.
\end{aligned}
\]

Since

\[
\int \frac{du}{1+u^2} = \arctan u + C,
\]

we have

\[
\int_{-1}^{1}\frac{1}{1+25x^2}\,dx = \frac{1}{5}\Bigl[\arctan(5)-\arctan(-5)\Bigr] = \frac{2}{5}\arctan(5).
\]

Numerically, with \(\arctan(5)\approx 1.3734\),
\[
\frac{2}{5}\arctan(5) \approx 0.54936.
\]

**Observations:**
- The substitution reduces the integral to a standard arctangent form.
- The denominator’s large coefficient (25) causes a sharp peak at \(x=0\) and very small values at \(x=\pm1\).

---

### (b) Integration via Newton's Divided Difference Interpolating Polynomials

This method constructs an interpolating polynomial from the function values at selected nodes and then integrates that polynomial exactly.

#### Global 2-Point (Linear) Interpolation

**Nodes:** \(x_0 = -1\) and \(x_1 = 1\)

**Function values:**
\[
f(-1)=\frac{1}{1+25}=\frac{1}{26},\quad f(1)=\frac{1}{26}.
\]

Since both values are equal, the linear (first–degree) interpolant is the constant function

\[
P_1(x) = \frac{1}{26}.
\]

**Integration:**
\[
\int_{-1}^{1} P_1(x)\,dx = \frac{1}{26}\cdot(1-(-1))=\frac{2}{26}=\frac{1}{13}\approx 0.07692.
\]

*Observations:*  
This approximation captures only about 14% of the true area because the two endpoints miss the central peak.

---

#### Global 3-Point (Quadratic) Interpolation

**Nodes:** \(x_0=-1\), \(x_1=0\), and \(x_2=1\)

**Function values:**
\[
f(-1)=\frac{1}{26},\quad f(0)=1,\quad f(1)=\frac{1}{26}.
\]

**Divided differences:**

- Zeroth order:
  \[
  f[x_0]=\frac{1}{26},\quad f[x_1]=1,\quad f[x_2]=\frac{1}{26}.
  \]
- First order:
  \[
  f[x_0,x_1]=\frac{1-1/26}{0-(-1)}=\frac{25}{26},\quad
  f[x_1,x_2]=\frac{1/26-1}{1-0}=-\frac{25}{26}.
  \]
- Second order:
  \[
  f[x_0,x_1,x_2]=\frac{-\frac{25}{26}-\frac{25}{26}}{1-(-1)}=-\frac{25}{26}.
  \]

**Newton form:**
\[
\begin{aligned}
P_2(x)&=f(-1)+f[x_0,x_1](x+1)+f[x_0,x_1,x_2](x+1)(x-0)\\[1mm]
&=\frac{1}{26}+\frac{25}{26}(x+1)-\frac{25}{26}(x+1)x.
\end{aligned}
\]

Notice that
\[
(x+1)-x(x+1)= (x+1)(1-x)=1-x^2.
\]
Thus,
\[
P_2(x)=\frac{1}{26}+\frac{25}{26}(1-x^2)=1-\frac{25}{26}x^2.
\]

**Integration:**
\[
\begin{aligned}
\int_{-1}^{1}P_2(x)\,dx &= \int_{-1}^{1}\Bigl(1-\frac{25}{26}x^2\Bigr)dx\\[1mm]
&= \Bigl[x\Bigr]_{-1}^{1} - \frac{25}{26}\Bigl[\frac{x^3}{3}\Bigr]_{-1}^{1}\\[1mm]
&= \Bigl[(1-(-1))\Bigr] - \frac{25}{26}\Bigl(\frac{1^3-(-1)^3}{3}\Bigr)\\[1mm]
&= 2 - \frac{25}{26}\cdot\frac{2}{3}\\[1mm]
&= 2 - \frac{50}{78} = 2 - \frac{25}{39}\\[1mm]
&= \frac{78-25}{39}=\frac{53}{39}\approx 1.359.
\end{aligned}
\]

*Observations:*  
Even though the quadratic interpolant exactly matches the values at \(-1\), 0, and 1, it grossly overestimates the integral (by about 147%) because the interpolating parabola does not capture the rapid decay away from the central peak.

---

#### Global 5-Point (Quartic) Interpolation

**Nodes:** Equally spaced at
\[
x_0=-1,\quad x_1=-0.5,\quad x_2=0,\quad x_3=0.5,\quad x_4=1.
\]

**Function values:**
\[
\begin{aligned}
f(-1)&=\frac{1}{1+25}= \frac{1}{26}\approx 0.03846,\\[1mm]
f(-0.5)&=\frac{1}{1+25(0.25)}=\frac{1}{7.25}\approx 0.13793,\\[1mm]
f(0)&=1,\\[1mm]
f(0.5)&\approx 0.13793,\\[1mm]
f(1)&\approx 0.03846.
\end{aligned}
\]

**Due to the even symmetry** of the function, the interpolating polynomial must be even; thus, it takes the form
\[
P_4(x)=Ax^4+Bx^2+C,
\]
with \(C=P_4(0)=1\).

Matching the endpoints:
\[
P_4(1)=A+B+1=f(1)=\frac{1}{26}\quad\Longrightarrow\quad A+B=\frac{1}{26}-1=-\frac{25}{26}.
\]

Also, matching at \(x=0.5\):
\[
P_4(0.5)=A(0.5)^4+B(0.5)^2+1=\frac{A}{16}+\frac{B}{4}+1\approx 0.13793.
\]
Thus,
\[
\frac{A}{16}+\frac{B}{4} \approx 0.13793-1=-0.86207.
\]
Multiplying by 16:
\[
A+4B\approx -13.79312.
\]

Now, subtract the equation \(A+B=-\frac{25}{26}\approx -0.96154\) from this:
\[
(A+4B) - (A+B)=3B\approx -13.79312+0.96154=-12.83158,\quad B\approx -4.27719.
\]
Then,
\[
A\approx -0.96154 - B\approx 3.31565.
\]

Thus, the quartic polynomial (written in simplified form) is

\[
P_4(x)=\frac{1250}{377}x^4-\frac{3225}{754}x^2+1,
\]

since
\[
\frac{1250}{377}\approx 3.316\quad\text{and}\quad\frac{3225}{754}\approx 4.278.
\]

**Integration:**

Because \(P_4(x)\) is even, we have

\[
\int_{-1}^{1}P_4(x)\,dx=2\int_{0}^{1}\left(\frac{1250}{377}x^4-\frac{3225}{754}x^2+1\right)dx.
\]

Carrying out the integration yields

\[
\int_{-1}^{1}P_4(x)\,dx=\frac{179}{377}\approx 0.4748.
\]

*Observations:*  
- The quartic interpolation (using 5 nodes) produces a much better approximation than the linear or quadratic cases.
- With an approximation of \(0.4748\) (relative error of about 13.57%), the quartic polynomial better captures the sharp central peak as well as the rapid decay.
- Still, it underestimates the exact value \(0.54936\), highlighting the challenge posed by the function’s sharp features.

---

### (c) Gauss–Legendre Quadrature

For the integral

\[
\int_{-1}^{1}\frac{1}{1+25x^2}\,dx,
\]

the following results were obtained:

| \(n\) | Approximation  | Absolute Error | Relative Error |
|:-----:|:--------------:|:--------------:|:--------------:|
| 2     | 0.2142857143   | 0.33507        | 60.99%         |
| 4     | 0.3709273183   | 0.17843        | 32.48%         |
| 6     | 0.4617005584   | 0.08766        | 15.96%         |

**Observations on Gauss–Legendre Performance:**
- Convergence is slower than in Problem 1 because the sharp peak near \(x=0\) challenges the quadrature.
- Even with 6 nodes, the relative error remains around 16%.
- This example underscores how the smoothness of the integrand affects the efficiency of Gauss–Legendre quadrature.

---

### Summary of Polynomial Approximations

Below is a table comparing the approximations against the exact value \(0.54936\):

| Method                   | Polynomial                                  | Approximation | Absolute Error | Relative Error |
|--------------------------|---------------------------------------------|---------------|----------------|----------------|
| **Exact Value**          | -                                           | 0.54936       | -              | -              |
| **Linear (2-point)**     | \(P_1(x)=\tfrac{1}{26}\)                      | 0.07692       | 0.47244        | 86.00%         |
| **Quadratic (3-point)**  | \(P_2(x)=1-\tfrac{25}{26}x^2\)                | 1.35900       | 0.80964        | 147.38%        |
| **Quartic (5-point)**    | \(P_4(x)=\tfrac{1250}{377}x^4-\tfrac{3225}{754}x^2+1\) | 0.47480       | 0.07456        | 13.57%         |
| **Gauss–Legendre (2-node)** | -                                        | 0.21429       | 0.33507        | 60.99%         |
| **Gauss–Legendre (4-node)** | -                                        | 0.37093       | 0.17843        | 32.48%         |
| **Gauss–Legendre (6-node)** | -                                        | 0.46170       | 0.08766        | 15.96%         |

**Comparative Analysis:**
- The linear interpolant grossly underestimates the integral.
- The quadratic interpolant overestimates by a large margin.
- The quartic interpolant (with nodes strategically including \(x=0\)) achieves the best accuracy among the polynomial methods.
- Gauss–Legendre quadrature, although optimal for many smooth functions, also underestimates the integral here due to the sharp peak at \(x=0\).

---

## Problem 3: Determining a Hermite-Type Quadrature Formula

### Problem Statement

Determine constants \(a\), \(b\), \(c\), and \(d\) so that the quadrature formula

\[
\int_{-1}^{1} f(x)\,dx = a\,f(-1) + b\,f(1) + c\,f'(-1) + d\,f'(1)
\]

is exact for all polynomials \(f(x)\) of degree 3 or less.

### Derivation

For exactness on \(f(x)=1\), \(x\), \(x^2\), and \(x^3\), we set up:

1. **For \(f(x)=1\):**  
   - \(f(-1)=1\), \(f(1)=1\); \(f'(-1)=0\), \(f'(1)=0\).  
   - \(\int_{-1}^{1}1\,dx = 2\).  
   - **Equation:** \(a+b=2\).

2. **For \(f(x)=x\):**  
   - \(f(-1)=-1\), \(f(1)=1\); \(f'(-1)=1\), \(f'(1)=1\).  
   - \(\int_{-1}^{1}x\,dx = 0\).  
   - **Equation:** \(-a+b+c+d=0\).

3. **For \(f(x)=x^2\):**  
   - \(f(-1)=1\), \(f(1)=1\); \(f'(-1)=-2\), \(f'(1)=2\).  
   - \(\int_{-1}^{1}x^2\,dx = \frac{2}{3}\).  
   - **Equation:** \(a+b-2c+2d=\frac{2}{3}\).

4. **For \(f(x)=x^3\):**  
   - \(f(-1)=-1\), \(f(1)=1\); \(f'(-1)=3\), \(f'(1)=3\).  
   - \(\int_{-1}^{1}x^3\,dx = 0\).  
   - **Equation:** \(-a+b+3c+3d=0\).

**Solving the system:**

- From Equation 1:  
  \[
  b=2-a.
  \]
- Substitute into Equation 2:
  \[
  -a+(2-a)+c+d=0\quad\Longrightarrow\quad c+d=2a-2.
  \]
- Substitute into Equation 3:
  \[
  a+(2-a)-2c+2d=\frac{2}{3}\quad\Longrightarrow\quad 2-2c+2d=\frac{2}{3}.
  \]
  Divide by 2:
  \[
  -c+d=-\frac{2}{3}\quad\Longrightarrow\quad d-c=-\frac{2}{3}.
  \]
  (Equivalently, \( -c+d = -\frac{2}{3}\).)
- Now, adding \(c+d=2a-2\) and \(-c+d=-\frac{2}{3}\) gives:
  \[
  2d=2a-2-\frac{2}{3}=2a-\frac{8}{3}\quad\Longrightarrow\quad d=a-\frac{4}{3}.
  \]
- Then,
  \[
  c=(2a-2)-d=2a-2-\left(a-\frac{4}{3}\right)=a-\frac{2}{3}.
  \]
- Substitute \(b=2-a\), \(c=a-\frac{2}{3}\), and \(d=a-\frac{4}{3}\) into Equation 4:
  \[
  -a+(2-a)+3\Bigl(a-\frac{2}{3}\Bigr)+3\Bigl(a-\frac{4}{3}\Bigr)=0.
  \]
  Simplify:
  \[
  -2a+2+3a-2+3a-4=0\quad\Longrightarrow\quad 4a-4=0,\quad\text{so}\quad a=1.
  \]
- Then,
  \[
  b=2-1=1,\quad c=1-\frac{2}{3}=\frac{1}{3},\quad d=1-\frac{4}{3}=-\frac{1}{3}.
  \]

Thus, the Hermite-type quadrature formula is

\[
\int_{-1}^{1} f(x)\,dx \approx f(-1)+f(1)+\frac{1}{3}\,f'(-1)-\frac{1}{3}\,f'(1).
\]

**Verification:**

- For \(f(x)=1\):  
  Exact: \(2\);  
  Formula: \(1+1+0-0=2\).

- For \(f(x)=x\):  
  Exact: \(0\);  
  Formula: \((-1)+1+\frac{1}{3}(1)-\frac{1}{3}(1)=0\).

- For \(f(x)=x^2\):  
  Exact: \(\frac{2}{3}\);  
  Formula: \(1+1+\frac{1}{3}(-2)-\frac{1}{3}(2)=2-\frac{4}{3}=\frac{2}{3}\).

- For \(f(x)=x^3\):  
  Exact: \(0\);  
  Formula: \((-1)+1+\frac{1}{3}(3)-\frac{1}{3}(3)=0\).

---

## Conclusion

This homework explored three different numerical integration problems:

1. **Problem 1:**  
   The analytical evaluation of \(\int_{0}^{\pi} e^x\cos x\,dx\) yielded the exact result \(-\frac{e^\pi+1}{2}\) (approximately \(-12.07035\)). Gauss–Legendre quadrature was shown to converge extremely rapidly for this smooth integrand.

2. **Problem 2:**  
   The integral \(\int_{-1}^{1}\frac{1}{1+25x^2}\,dx\) was evaluated analytically to be \(\frac{2}{5}\arctan(5)\approx0.54936\).  
   - **Global polynomial interpolation** using Newton's divided differences was performed with 2, 3, and 5 nodes. The linear interpolant grossly underestimated the area, the quadratic interpolant overestimated it, and the quartic interpolant (with nodes strategically including \(x=0\)) achieved a much better (though still not perfect) approximation.  
   - **Gauss–Legendre quadrature** applied directly on \([-1,1]\) showed moderate improvement with increased nodes; even with 6 nodes the relative error was around 16%.

3. **Problem 3:**  
   A Hermite-type quadrature formula incorporating derivative information was derived. The final formula,
   \[
   \int_{-1}^{1} f(x)\,dx \approx f(-1)+f(1)+\frac{1}{3}\,f'(-1)-\frac{1}{3}\,f'(1),
   \]
   is exact for all polynomials up to degree 3.

**Key Insights:**

- **Function Smoothness:** Smooth functions (like \(e^x\cos x\)) allow Gauss–Legendre quadrature to converge exponentially fast, whereas functions with sharp features (like \(\frac{1}{1+25x^2}\)) challenge both global interpolation and quadrature methods.
- **Node Placement:** The strategic placement of nodes (especially including the region where the function attains its peak) is critical. For instance, the quartic interpolant that includes \(x=0\) performs significantly better.
- **Use of Derivative Information:** Incorporating derivatives (as in the Hermite-type formula) can dramatically improve accuracy with fewer evaluation points.
- **No Universal Method:** The optimal numerical integration method depends on the integrand’s properties and the desired accuracy.

These principles are fundamental when selecting numerical integration strategies for various applications in scientific computing and engineering.
