# Linear Systems

## Introduction

This section is dedicated to explaining the linear systems methods included in this repository. It will contain the general approach of each method, convergence properties where applicable, and some requirements. It is by no means exhaustive or academic, as full proofs for convergence will not be included. However, it can serve as a basic companion when exploring this repository.

The implementation of these methods is meant to be academic, prioritizing code readability over efficiency. For iterative methods, a different stopping criteria is used in comparison with other sections, using the norm of the residual vector $||Ax_n - b||$ added to the change in the solution vector $||x_{n+1} - x_n||$. The evolution of the error is also shown as vectors, allowing the tracking of the residual norm or solution change over iterations.

Direct methods factorize the matrix A, then solve the resulting simpler systems via substitution. They are typically $O(n^3)$ but accurate and non-iterative. Iterative methods gradually improve an initial guess, useful for large sparse systems where factorization is expensive.

## Direct Methods

### DS

Direct substitution solves a lower triangular system $Ly = b$ via forward substitution. Used internally by many factorization methods (LU, Cholesky, QR) to complete the solve after factorization.

### IS

Inverse substitution solves an upper triangular system $Ux = y$ via backward substitution. Completes the factorization-based solve after forward substitution.

### Gauss Elimination

Gaussian elimination with partial pivoting solves $Ax = b$ by computing the LU decomposition $A = LU$. The algorithm eliminates variables column by column through row operations, maintaining numerical stability via pivoting. Forward substitution solves $Ly = b$, then backward substitution solves $Ux = y$. Works for any non-singular square matrix, though ill-conditioned matrices may suffer numerical issues.

### Cholesky Decomposition

Cholesky factorization is specialized for symmetric positive definite (SPD) matrices. It decomposes $A = LL^T$ where L is lower triangular, using roughly half the operations of LU. Then solves the triangular systems $Ly = b$ and $L^Tx = y$. More stable and efficient than Gaussian elimination for SPD matrices.

### CroutFact

Crout's LU decomposition variant computes LU factorization where U has 1s on its diagonal. Suitable for general square matrices and can be more numerically stable in certain situation.

### Crout

Crout's method to solve tridiagonal systems. It is significantly faster than Gauss, being linear ($\mathcal O(n)$), so it is widely used when tridiagonal systems appear.

### CroutPenta

Crout method optimized specifically for pentadiagonal matrices (5 non-zero diagonals). Exploits sparsity for efficiency on banded systems.

### Householder

QR decomposition using Householder reflections. Decomposes $A = QR$ where Q is orthogonal and R is upper triangular. Generally more numerically stable than Gram-Schmidt, solving being then done by forward and backward substitution.

### QR

QR decomposition. After factorization, solving $Ax = b$ becomes $R(Q^Tx) = Q^Tb$, which is numerically stable. Can also be useful for least squares and eigenvalue problems.

## Iterative Methods

### Jacobi

The Jacobi method is a stationary iterative method given by:

$$x_{n+1} = D^{-1}(b - (L + U)x_n)$$

where D, L, U are the diagonal, lower, and upper triangular parts of A. Each component of $x_{n+1}$ is computed independently from all components of $x_n$, making it easily parallelizable. It converges for strictly diagonally dominant matrices or symmetric positive definite matrices, though convergence may be slow. Convergence is linear (first-order).

### Gauss-Seidel

Gauss-Seidel improves on Jacobi by using the most recent values:

$$x_{n+1} = (D + L)^{-1}(b - U x_n)$$

As each component is updated, new values are used immediately in subsequent component calculations. This typically converges faster than Jacobi for the same class of matrices, with roughly twice the convergence rate. Still first-order convergence.

Based on the code you provided, here's the corrected documentation for SOR1 and SOR2:

## SOR1 and SOR2

Successive Over-Relaxation accelerates convergence by extrapolating the Gauss-Seidel step. The method can be formulated in two different ways:

**SOR1** uses the formulation:
$$x_{n+1} = (D + \omega L)^{-1}[(1-\omega)D - \omega U]x_n + \omega(D + \omega L)^{-1}b$$

**SOR2** uses the alternative formulation:
$$x_{n+1} = (D + \omega U)^{-1}[(1-\omega)D - \omega L]x_n + \omega(D + \omega U)^{-1}b$$

For $\omega > 2$, the method will diverge. The optimal $\omega$ depends on the problem and can significantly accelerate convergence. SOR1 requires a positively defined or strictly diagonally dominant matrix. SOR2 can work for more general complex matrices, but convergence isn't ensured.

### JOR

Jacobi Over-Relaxation applies relaxation to the Jacobi method:

$$x_{n+1} = (1 - \omega)x_n + \omega D^{-1}(b - (L + U)x_n)$$

A parallel-friendly alternative to SOR that combines Jacobi's parallelizability with over-relaxation.

### Gradient

Steepest descent (or gradient descent) is a basic iterative method that moves in the direction of steepest decrease of the quadratic form. Slower than conjugate gradient but simpler and more general (works on non-SPD systems if viewed as optimization). Convergence is linear and depends on the condition number.

### Conjugate Gradient

For symmetric positive definite A, conjugate gradient solves $Ax = b$ by minimizing the quadratic form $\frac{1}{2}x^TAx - b^Tx$. It is a Krylov subspace method that converges in at most n iterations (theoretically), though practically converges much faster. Excellent for large sparse SPD systems and does not require storing explicit LU factors.
