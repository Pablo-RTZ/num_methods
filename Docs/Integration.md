# Integration

## Introduction

This section is dedicated to explaining the integration methods included in this repository. It is designed more as a companion to using them, rather than a full explanation or proof. This section also includes some strategies to extend those methods to iterative ones, or to adapt them for improper or infinite integrals.

As with all the code in this repository, it is meant to be academic, so possible efficiency improvements are traded for readability and clarity. Most integration methods consist of a vector of nodes, where the function is evaluated, and a vector of weights associated to those nodes. Then, the integral is approximated by $I=\int_a^b f(x)dx\approx \sum W\odot f(X)$, where $\odot$ denotes element by element multiplication. Therefore, these methods can be easily extended to multiple integrals, extending the vectors to matrices or tensors.

On the implementation of those methods in this MatLab repository, they are mostly useable with base MATLAB, but some Gaussian quadrature methods require the Symbolic Math Toolbox. Also, numerical stability is to be considered with these methods, as in many cases the methods will generate and work with high degree polynomials (for big $n$), so rounding errors will significantly affect the result.

## Trapezoidal

The trapezoidal method is one of the simplest methods to approximate definite integrals. It divides the interval into regularly spaced subintervals, separated by nodes. Taking functional evaluations on those nodes, it approximates the function by linear functions for each pair of nodes. It then sums them. It is given by $I\approx \sum W\odot f(x_i)$ with $x_i$ being equispaced nodes in the integration interval, $W=[1\ 2\ 1\ ...\ 2\ 1]$, and $h=\frac{(b-a)}{n}$.

## Simpson

The Simpson method is conceptually an extension of the trapezoidal one. It divides the interval into regularly spaced subintervals, separated by nodes. Taking functional evaluations on those nodes, it approximates the function by quadratic functions for each pair of nodes. It then sums them. It is given by $I\approx \frac{h}{3}\sum W\odot f(x_i)$, $x_i$ being $n$ equispaced nodes in the integration interval, and $W=[1\ 4\ 2\ 4\ ...\ 4\ 1]$.

## Monte-Carlo

If, on the contrary, the nodes are chosen to be random, the Monte-Carlo method is obtained. The weights are all equal for the nodes, and the average of the functional evaluations, multiplied by the length of the integration interval is taken as an approximation of the integral, $I\approx \frac{(b-a)}{n}\sum f(x_i)$.

## Romberg

Romberg integration combines the trapezoidal rule with Richardson extrapolation to achieve high accuracy with fewer function evaluations. The method constructs a triangular array $R_{i,j}$:

- $R_{i,0}$: trapezoidal rule with $2^i$ subintervals
- $R_{i,j} = \frac{4^j R_{i,j-1} - R_{i-1,j-1}}{4^j - 1}$ for $j \ge 1$

It computes trapezoidal approximations with doubling step sizes: $h, h/2, h/4, \dots$, then applies extrapolation successively to cancel error terms, obtaining the diagonal entries $R_{i,i}$ as improved estimates. This way, it converges faster than trapezoidal or Simpson for smooth functions. It stops when the difference between iterations is lower than the desired tolerance.

## Gaussian quadrature

### Legendre

Gauss–Legendre quadrature approximates integrals on $[-1, 1]$ (or transformed to $[a, b]$) using roots of Legendre polynomials as nodes. Weights are derived from the Lagrange polynomials. It exactly integrates polynomials of degree up to $2n-1$ with $n$ nodes. It is suitable for smooth functions on finite intervals.

### Chebyshev

Gauss–Chebyshev quadrature uses the roots of Chebyshev polynomials of the first kind $T_n(x)$ as nodes on $[-1, 1]$, with constant weights $w_i = \frac{\pi}{n}$. Designed for integrals of the form $\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} dx$. In this implementation, only $f(x)$ has to be inputed, so without the $\sqrt{1-x^2}$ term.

### Laguerre

Gauss–Laguerre quadrature handles integrals on $[0, \infty)$ of the form $\int_0^\infty e^{-x} f(x) dx$. Nodes are roots of Laguerre polynomials $L_n(x)$, with appropriate weights. It is useful for semi-infinite domains with exponential decay. In this implementation, only $f(x)$ has to be inputed.

### Hermite

Gauss–Hermite quadrature targets integrals on $(-\infty, \infty)$ of the form $\int_{-\infty}^\infty e^{-x^2} f(x) dx$. Nodes are roots of Hermite polynomials $H_n(x)$. It is suitable for problems with Gaussian weighting. In this implementation, only $f(x)$ has to be inputed

### Lobatto

Lobatto quadrature is a variant of Gaussian quadrature that includes the endpoints $a$ and $b$ as nodes. It integrates polynomials of degree up to $2n-3$ exactly with $n$ nodes. It is useful when boundary values are known or required.

## Iterative schemes

## Infinite and improper integrals
