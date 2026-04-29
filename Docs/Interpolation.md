# Interpolation

## Introduction

The programs in this section are meant for both function fitting as well as interpolation. They allow to find the best fitting line or polynomial for a set of points, to interpolate a set of points using polynomials or piecewise functions, and to aproximate a continuous function using various families of orthogonal functions. Best fit lines are limited to polynomials, as for other types of fitting (exponential, rational,...), a variable change can be normally used to linearize the problem.

Except for the linear and polynomical fitting programs, the rest return symbolic expresisons, as it makes it easier to display neatly. Nonetheless, the code could easily be rewritten to use handles or any other format. As with all, the code is only one of many possible implementations, this one with readability as a focus. Another design choice is to use the backlash operator to solve linear systems, for simplicity. In many cases, the systems that arise in this section are tridiagonal, where much faster algorithms exist (Crout).

Finally, one of the main limitations of these methods is numerical stability when high degree polynomials are used, as coefficients tend to get really small, and some oscillatory behavious might appear as a result of the high degree terms (GIbbs phenomenon).

## Linear interpolation

Linear interpolation connects data points with straight lines. For fitting a line to data, the LinearFit program performs least squares linear regression, returning coefficients `[a0, a1]` for the line `y = a0 + a1*x`, as well as the quadratic error on the interpolation nodes as an optional argument.

## Polynomical interpolation

Polynomial interpolation fits a polynomial of degree $n-1$ through $n$ points. It uses least squares to find coefficients for the best-fit polynomial, returning the coefficients and error.

## Lagrange and Newton polynomials

For exact interpolation through given points, `LagrangeFit(xi, fi)` returns a symbolic Lagrange polynomial that passes through all points. `NewtonFit(xi, fi)` uses Newton's divided difference method to construct the interpolating polynomial symbolically. The unicity for the $n$ degree polynomial that interpolates $n+1$ points can be proven, which implies that both Lagrange and Newton polynomials return the same polynomial. However, the grouping of terms is different, so the expressions are visually different, despite being the exact same polynomial.

## Splines

Splines provide piecewise polynomial interpolation for smoother curves. `LinearSpline(xi, fi)` creates linear segments between points, with the resulting function being continous (but not $\mathcal C^1$ generally). `QuadraticSpline(xi, fi)` uses quadratic pieces, while keeping the first derivative continuous at piecewise changes ($\mathcal C^1$). `CubicSpline(xi, fi, cond)` generates cubic splines with optional boundary conditions for derivatives (defaults to natural splines with zero derivatives at endpoints). Due to the extra degrees of freedom, continuity at the first and second derivative can be enforced ($\mathcal C^2$).

## Continuous aproximation

These methods approximate functions using series expansions over intervals. All these methods start with a family of orthogonal functions (in the Linear Algebra sense), and by proyecting the function to be fit onto the corresponding subspace, aproximate it as a linear combination of the family's members, using the corresponding Fourier coefficients.

### Legendre polynomials

`LegendreTSeries(f, N, a, b)` computes the truncated Legendre series expansion of a symbolic function `f` over `[a, b]` with `N` terms, returning the symbolic approximation and coefficients.

### Chebyshev polynomials

`ChebyshevTSeries(f, N, a, b)` approximates a function using Chebyshev polynomials, providing better convergence properties. Returns the symbolic series and coefficients.

### Trigonometric series

`TrigSeries(f, N, a, b)` computes the Fourier series expansion of a function over `[a, b]`. The method considers only one period to get the coefficients, but if the function is periodic in the $[a,b]$ interval, and all of $\mathbb R$ is considered, it will obtain its truncated Fourier series.
