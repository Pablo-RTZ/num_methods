# Root Finding methods

## Introduction

This section is dedicated to explaining the root finding methods included on this repository. It will contain the general expression of the iterative method, convergence criterion where applicable, and some edge cases. It is by no means exhaustive or academic, as full proofs for the convergence and order of convergence of these methods will not be included. However, it can serve as a basic companion when exploring this repository.

Additionally, the implementation of these methods is meant to be academic, prioritizing the code being readable over efficiency. Several programming choices, such as the general preference to use step size $|x_{n+1}-x_n|$ over error $|f(x_n)|$ as a stopping criterion, or the use of options for additional arguments, are not relevant for the method being valid.

Except for the bisection method, the methods described here are fixed-point iterative methods. Let $f(x)$ be a function whose root is sought, so that $f(x_{\text{root}}) = 0$. A fixed-point function is a function $g(x)$ such that $g(x_{\text{root}}) = x_{\text{root}}$. In that case, the sequence $\{ x_n \}_{n=0}^\infty$, defined by $x_{n+1} = g(x_n)$, converges to its fixed point, which is also a root of the original function. Consequently, both $|x_{n+1} - x_n|$ and $|f(x_n)|$ converge to 0. The methods described here are of the form $x_{n+1} = G(x_n, f(x_n), f'(x_n), \dots)$, and, for a sufficiently smooth function and an initial guess $x_0$ sufficiently close to the root, they will converge. Generally, convergence criterion for the methods will not get more specific than this.

## Bisection

This is usually the first root finding method taught. It is a bracketing methods, that starts in the interval $(a,b)$ whose endpoints must have oposite signs ($sign(f(a))\ne sign(f(b))$). The midpoint is calculated, and depending on the sign of the functional evaluation at the midpoint, a new bracket is chosen. It reduces the searchspace by half on each iteration, so the error is approximately $e_{n+1}=\frac{e_n}2\implies e_n\approx |b-a|\cdot\left(\frac 12\right)^n$, making it a linear method (order of convergence 1). If the signs of the functional evaluations of the inital interval are opposite, it will always converge to a root.

## Newton

The Newton-Raphson method, is one of the most widely used and studied methods for root finding. It is given by:

$$x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}$$

It can be geometrically interpreted as the intersection between the tangent to the point and the y axis. It converges cuadratically, $e_{n+1}=c_2 e_k^2+\mathcal O(e_k^3)$ for simple roots, given sufficiently smooth functions and sufficiently close initial guesses. It will diverge if $f'(x_n)=0$ for some $n$.

## Secant

The secant method is an implementation of the Newton-Raphson method, but using an approximation for the derivative. It is given by:

$$x_{n+1}=x_n-\frac{f(x_n)(x_n-x_{n-1})}{f(x_n)-f(x_{n-1})}$$

It requires two initial guesses, and has a convergence order of $\varphi\approx 1.6$.

## NS1

The Newton-Raphson method converges linealry (rather than cuadratically) for multiple roots. However, this method, which is really similar, converges cuadratically for multiple roots, given the multiplicity is known. It is given by:

$$x_{n+1}=x_n-m\frac{f(x_n)}{f'(x_n)}$$

## NS2

If the root multiplicity is not known, the NS2 method can be used, having a cuadratic convergence order, but requiring an evaluation of the second derivative. It is given by:

$$x_{n+1}=x_n-\frac{f(x_n)f'(x_n)}{(f'(x_n))^2-f(x_n)f''(x_n)}$$

As all the other fixed point methods, it converges for a sufficiently smooth function, and a sufficiently close initial guess.

## Chebyshev

Chebyshev's method is a third-order iterative method, derived from Taylor expansion or from the concept of tangent hyperbolas. It is given by:

$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)} - \frac{f''(x_n)}{2f'(x_n)} \left( \frac{f(x_n)}{f'(x_n)} \right)^2$$

The method converges cubically for simple roots, requiring the evaluation of $f$, $f'$, and $f''$ per iteration. It is more efficient than Newton when higher accuracy is needed and second derivatives are inexpensive to compute.

## ChebyshevHalley

The Chebyshev-Halley method is a family of iterative methods that generalizes both Chebyshev's method and Halley's method. The general expression is:

$$x_{n+1} = x_n - \left(1 + \frac{L_f(x_n)}{2(1 - \beta L_f(x_n))}\right) \frac{f(x_n)}{f'(x_n)}$$

where $L_f(x_n) = \frac{f''(x_n)f(x_n)}{(f'(x_n))^2}$ and $\beta$ is a parameter. Special cases include:

- $\beta = 0$: Chebyshev's method (cubic convergence)
- $\beta = \frac{1}{2}$: Halley's method (cubic convergence)
- $\beta = 1$: Newton's method (quadratic convergence)
beta
The Chebyshev-Halley family converges cubically for simple roots when $\alpha \ne 1$, requiring $f$, $f'$, and $f''$ per iteration. It is particularly useful when dealing with functions where the second derivative is available and higher-order convergence is desired.

## AX19

The AX19 method is a more recent root-finding scheme (proposed around 2019). The iterative scheme is given by:

$$x_{n+1} = x_n - \frac{(f(x_n) \cdot f'(x_n)) }{ (f'(x_n)^2 + 2 \cdot f(x_n)^2)}$$

As with other fixed-point methods, it requires a sufficiently close initial guess and $f'(x_n) \ne 0$.

## Steffensen

Steffensen's method is a derivative-free method that achieves quadratic convergence, similar to Newton's method, but without requiring the derivative. It replaces $f'(x_n)$ with a finite difference approximation using the same point:

$$x_{n+1} = x_n - \frac{f(x_n)^2}{f(x_n + f(x_n)) - f(x_n)}$$

It can be seen as Newton's method with $f'(x_n) \approx \frac{f(x_n + f(x_n)) - f(x_n)}{f(x_n)}$. The method converges quadratically for simple roots, provided $x_0$ is sufficiently close to the root and $f'(x_{\text{root}}) \ne 0$. It requires two function evaluations per iteration. One edge case is when $f(x_n)$ is very large, causing the step $f(x_n)$ in the argument to be extreme.

## Traub

Traub's method is a method obtained from iterating Newton.

$$y_n = x_n - \frac{f(x_n)}{f'(x_n)}$$
$$x_{n+1} = y_n - \frac{f(y_n)}{f'(x_n)}$$

This method uses two evaluations of $f$ ($f(x_n)$ and $f(y_n)$) and one evaluation of $f'$ per full iteration, achieving fourth-order convergence. It is a precursor to many modern optimal methods. For simple roots and sufficiently smooth functions, the order of convergence is 4. Edge cases include $f'(x_n) = 0$ (causing division by zero) and cases where $y_n$ is not a better approximation (e.g., near oscillatory behavior).
