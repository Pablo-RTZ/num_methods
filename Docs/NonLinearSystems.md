# Non Linear Systems methods

## Introduction

This section is dedicated to explaining the nonlinear systems methods included in this repository. These methods extend root-finding techniques to solve systems of nonlinear equations of the form $\mathbf{F}(\mathbf{x}) = \mathbf{0}$, where $\mathbf{F}: \mathbb{R}^n \to \mathbb{R}^n$ and $\mathbf{x} = (x_1, x_2, \ldots, x_n)^T$.

The documentation will contain the general expression of each iterative method, convergence properties where applicable, and relevant edge cases. As with root-finding methods, implementations prioritize readability over efficiency and are intended for academic purposes.

Most methods described here are Newton-like iterative methods, which extend the classical Newton-Raphson approach for single equations to systems. They generally take the form:

$$\mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{s}_n$$

where $\mathbf{s}_n$ is a step computed using the Jacobian matrix $J_F(\mathbf{x}_n)$ and function evaluations. For sufficiently smooth functions and initial guesses sufficiently close to the solution, these methods will converge. The primary computational cost at each iteration involves solving linear systems using the Jacobian matrix (typically via LU decomposition or similar).

## Newton

The Newton-Raphson method extended to systems is the most fundamental approach for solving nonlinear systems. The iterative scheme is:

$$\mathbf{x}_{n+1} = \mathbf{x}_n - J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$

where $J_F(\mathbf{x}_n)$ is the Jacobian matrix of $\mathbf{F}$ evaluated at $\mathbf{x}_n$. In practice, this is solved as a linear system:

$$J_F(\mathbf{x}_n) \mathbf{d}_n = \mathbf{F}(\mathbf{x}_n), \quad \mathbf{x}_{n+1} = \mathbf{x}_n - \mathbf{d}_n$$

It converges quadratically for sufficiently smooth functions and initial guesses sufficiently close to the solution. However, it requires the Jacobian matrix (all $n^2$ partial derivatives) and solving an $n \times n$ linear system at each iteration. The method fails if the Jacobian is singular at any iteration.

## Traub

Traub's method for systems extends the single-variable Traub method by using two function evaluations and a single Jacobian evaluation per iteration. The scheme is:

$$\mathbf{y}_n = \mathbf{x}_n - J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$
$$\mathbf{x}_{n+1} = \mathbf{y}_n - J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{y}_n)$$

This method achieves fourth-order convergence for simple roots when the Jacobian is sufficiently accurate and the initial approximation is sufficiently close to the solution. Computationally, it requires two evaluations of $\mathbf{F}$ and one Jacobian evaluation per iteration, plus two linear system solves. Edge cases include singular Jacobians and cases where $\mathbf{y}_n$ is not a better approximation.

## Trapezoidal

The Trapezoidal method for systems uses an average of the Jacobians at two points to improve the iterative process. The iterative scheme combines an initial Newton step with a corrector step:

$$\mathbf{y}_n = \mathbf{x}_n - J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$
$$\mathbf{x}_{n+1} = \mathbf{x}_n - \frac{1}{2}\left(J_F(\mathbf{x}_n) + J_F(\mathbf{y}_n)\right)^{-1} \mathbf{F}(\mathbf{x}_n)$$

This method uses a trapezoidal approximation of the average Jacobian between two consecutive points. It requires evaluations of both $\mathbf{F}$ and $J_F$ at $\mathbf{x}_n$ and $\mathbf{y}_n$, providing improved convergence properties compared to Newton in some cases. The method achieves super-linear convergence and typically has better stability properties than Newton's method for certain problems. The main computational drawback is the need to evaluate the Jacobian at two points and solve a different linear system.

## Jarrat

Jarrat's method is a three-step method designed to achieve higher convergence rates with efficient computational cost. The iterative scheme is:

$$\mathbf{y}_n = \mathbf{x}_n - \frac{2}{3} J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$
$$\mathbf{z}_n = \mathbf{x}_n - \frac{1}{2} (3J_F(\mathbf{y}_n) - J_F(\mathbf{x}_n))^{-1} (3J_F(\mathbf{y}_n) + J_F(\mathbf{x}_n)) J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$
$$\mathbf{x}_{n+1} = \mathbf{z}_n$$

Jarrat's method achieves fourth-order convergence per iteration while requiring only one Jacobian evaluation at the initial point and one at an intermediate point. The method is particularly efficient for systems where the Jacobian is expensive to compute. However, it requires careful numerical handling of the matrix operations and is sensitive to the initial guess.

## HMT

The HMT method is a three-step method that uses a weighted combination of function evaluations to improve convergence. The iterative scheme is:

$$\mathbf{y}_n = \mathbf{x}_n - J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$
$$\mathbf{z}_n = \mathbf{x}_n - J_F(\mathbf{x}_n)^{-1} (\mathbf{F}(\mathbf{y}_n) + p \mathbf{F}(\mathbf{x}_n))$$
$$\mathbf{x}_{n+1} = \mathbf{x}_n - J_F(\mathbf{x}_n)^{-1} (\mathbf{F}(\mathbf{z}_n) + \mathbf{F}(\mathbf{y}_n) + p \mathbf{F}(\mathbf{x}_n))$$

where $p$ is a weighting parameter. This method achieves eighth-order convergence with a single Jacobian evaluation, making it computationally efficient despite the multiple function evaluations. The method is particularly useful when the Jacobian is expensive but function evaluations are relatively cheap. The convergence order is very high, but the method requires careful initialization and is sensitive to the parameter $p$ and the initial guess.

## M5

The M5 method is an optimal fifth-order method with carefully chosen parameters to achieve high efficiency. The iterative scheme involves three sub-steps:

$$\mathbf{y}_n = \mathbf{x}_n - \alpha J_F(\mathbf{x}_n)^{-1} \mathbf{F}(\mathbf{x}_n)$$
$$\mathbf{z}_n = \mathbf{y}_n - J_F(\mathbf{x}_n)^{-1} (\beta \mathbf{F}(\mathbf{x}_n) + \gamma \mathbf{F}(\mathbf{y}_n))$$
$$\mathbf{x}_{n+1} = \mathbf{z}_n - J_F(\mathbf{x}_n)^{-1} (\rho \mathbf{F}(\mathbf{x}_n) + \mu \mathbf{F}(\mathbf{y}_n) + \sigma \mathbf{F}(\mathbf{z}_n))$$

with optimal parameters $\alpha = 1$, $\beta = 0$, $\gamma = 5$, $\rho = 0$, $\mu = -16/5$, $\sigma = 1/5$. This method achieves fifth-order convergence with only one Jacobian evaluation per iteration, using three function evaluations. It represents an optimal balance between convergence rate and computational cost, making it highly efficient for systems where evaluating the Jacobian is expensive relative to function evaluations. Edge cases include singular Jacobians and cases where the optimal parameter values are not appropriate for the specific problem structure.