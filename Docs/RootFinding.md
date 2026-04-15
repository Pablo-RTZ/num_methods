# Root Finding methods

## Introduction

This section is dedicated to explaining the root finding methods included on this repository. It will contain the general expression of the iterative method, convergence criterion where applicable, and some edge cases. It is by no means exhaustive or academic, as full proofs for the convergence and order of convergence of these methods will not be included. However, it can serve as a basic companion when exploring this repository.

Additionally, the implementation of these methods is meant to be academic, prioritizing the code being readable over efficiency. Several programming choices, such as the general preference to use step size $|x_{n+1}-x_n|$ over error $|f(x_n)|$ as a stopping criterion, or the use of options for additional arguments, are not relevant for the method being valid.

Except for the bisection method, the methods described here are fixed-point iterative methods. Let $f(x)$ be a function whose root is sought, so that $f(x^*) = 0$. A fixed-point function is a function $g(x)$ such that $g(x^*) = x^*$. In that case, the sequence $\{ x_n \}_{n=0}^\infty$, defined by $x_{n+1} = g(x_n)$, converges to its fixed point, which is also a root of the original function. Consequently, both $|x_{n+1} - x_n|$ and $|f(x_n)|$ converge to 0. The methods described here are of the form $x_{n+1} = G(x_n, f(x_n), f'(x_n), \dots)$, and, for a sufficiently smooth function and an initial guess $x_0$ sufficiently close to the root, they will converge. Generally, convergence criterion for the methods will not get more specific than this.

## Bisection

This is usually the first root finding method taught. It is a bracketing methods, that starts in the interval $(a,b)$ whose endpoints must have oposite signs ($sign(f(a))\ne sign(f(b))$). The midpoint is calculated, and depending on the sign of the functional evaluation at the midpoint, a new bracket is chosen. It reduces the searchspace by half on each iteration, so the error is approximately $e_{n+1}=\frac{e_n}2\implies e_n\approx |b-a|\cdot(\frac 12)^n$, making it a linear method (order of convergence 1). If the signs of the functional evaluations of the inital interval are opposite, it will always converge to a root.

## Newton

The Newton-Raphson method, is one of the most widely used and studied methods for root finding. It is given by $x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}$, and can be geometrically interpreted as the intersection between the tangent to the point and the y axis. It converges cuadratically, $e_{n+1}=c_2 e_k^2+\mathcal O(e_k^3)$ for simple roots, given sufficiently smooth functions and sufficiently close initial guesses. It will diverge if $f'(x_n)=0$ for some $n$.

## Secant

The secant method is an implementation of the Newton-Raphson method, but using an approximation for the derivative. It is given by $x_{n+1}=x_n-\frac{f(x_n)(x_n-x_{n-1})}{f(x_n)-f(x_{n-1})}$. It requires two initial guesses, and has a convergence order of $\varphi\approx 1.6$.

## NS1

The Newton-Raphson method converges linealry (rather than cuadratically) for multiple roots. However, this method, which is really similar, converges cuadratically for multiple roots, given the multiplicity is known. It is given by $x_{n+1}=x_n-m\frac{f(x_n)}{f'(x_n)}$.

## NS2

If the root multiplicity is not known, the NS2 method can be used, having a cuadratic convergence order, but requiring an evaluation of the second derivative. It is given by $x_{n+1}=x_n-\frac{f(x_n)f'(x_n)}{(f'(x_n))^2-f(x_n)f''(x_n)}$. As all the other fixed point methods, it converges for a sufficiently smooth function, and a sufficiently close initial guess.

## Chebyshev

## ChebyshevHalley

## AX19

## Steffensen

## Traub
