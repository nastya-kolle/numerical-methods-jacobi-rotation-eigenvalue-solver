# Jacobi Rotation Method for Symmetric Eigenvalue Problems — C++

A C++ implementation of the **Jacobi rotation method** for computing all eigenvalues and eigenvectors of a real symmetric matrix, along with a computational experiment that validates accuracy and convergence speed against matrices with known, controlled eigenvalues.

## Background

For a real **symmetric** matrix `A`, the Jacobi eigenvalue algorithm iteratively applies a sequence of orthogonal rotation (Givens/Jacobi) transformations:

```
A_(k+1) = J(p, q, θ)ᵀ · A_k · J(p, q, θ)
```

Each rotation `J(p, q, θ)` is chosen to zero out the **largest off-diagonal element** `A[p][q]` at that step. As rotations accumulate, the matrix converges to a diagonal matrix `Λ` whose diagonal entries are the eigenvalues of `A`, while the product of all rotation matrices `T = J₁ · J₂ · ... · Jₖ` converges to the matrix of eigenvectors, satisfying:

```
A · T ≈ T · Λ
```

The algorithm stops either when the largest off-diagonal element drops below a chosen tolerance `eps`, or after a maximum number of rotations.

## Project Structure

| File | Purpose |
|---|---|
| `Vector.h` / `Vector.cpp` | Custom vector class with manual memory management, 1-based bounds-checked indexing, and arithmetic operators. |
| `SquareMatrix.h` / `SquareMatrix.cpp` | Dense square matrix stored as an array of `Vector` rows. Supports matrix arithmetic, transpose, symmetry check, identity/zero fill, and string/file output. |
| `JacobiSolver.h` / `JacobiSolver.cpp` | The core algorithm: `JacobiSolver::solve(A, eps, max_rotations)` runs the rotation sweeps and returns eigenvalues, eigenvectors, rotation count, and accuracy measures. |
| `LinearAlgebra.h` / `LinearAlgebra.cpp` | Vector utilities used to build test matrices: Gram–Schmidt-style orthogonalization (`orthogonalizeAgainst`) and normalization. |
| `ComputationalExperiment.h` / `.cpp` | Generates random symmetric matrices with **known eigenvalues** (by conjugating a diagonal matrix with a random orthogonal matrix), runs the solver, and reports convergence/accuracy statistics. |
| `Main.cpp` | Runs the experiment for a set of matrix sizes and eigenvalue ranges and prints a summary table. |

## How eigenvalue accuracy is validated

Since directly-generated random symmetric matrices have unknown eigenvalues, the experiment builds test matrices with **eigenvalues chosen in advance**:

1. Pick a set of target eigenvalues `λ₁ … λₙ` (randomly, within a chosen range).
2. Build a random **orthogonal** matrix `V` via modified Gram–Schmidt orthogonalization (with re-orthogonalization for numerical stability) applied to a random matrix.
3. Construct `A = V · Λ · Vᵀ`, where `Λ = diag(λ₁, …, λₙ)`. Since `V` is orthogonal, `A` is symmetric and has exactly the chosen eigenvalues.
4. Run the Jacobi solver on `A` and compare its computed eigenvalues (sorted) against the true ones (sorted) to get a **maximum eigenvalue error**.

This gives a ground-truth way to measure the solver's accuracy, rather than only checking internal consistency.

## What `JacobiResult` reports

| Field | Meaning |
|---|---|
| `a` | Computed eigenvalues (1-based `Vector`). |
| `T` | Matrix of eigenvectors (columns), accumulated from all rotations. |
| `K` | Number of rotations performed until convergence (or the cap). |
| `IER` | `0` if converged before `max_rotations`; otherwise the rotation count at which the cap was hit. |
| `r` | Accuracy measure: `max_ij |A·T − T·Λ|`, i.e. how well the eigen-decomposition reconstructs the original system. |
| `final_max_off_diag` | Largest remaining off-diagonal element when iteration stopped. |

## Build & Run

The project is set up as a Visual Studio solution.

1. Open the `.sln` file in **Visual Studio** (2019+ recommended, C++17 or later).
2. Select **Debug** or **Release**, platform **x64**.
3. Build and run (`Ctrl+F5`).

No external dependencies are required — everything uses the C++ standard library.

## Sample Results

Running the program with the default configuration (sizes `{10, 30}` crossed with eigenvalue ranges `{[-2, 2], [-50, 50]}`, convergence tolerance `1e-9`) produces:

| # | N | Range λ | Max \|A_ij\| | Avg. rotations | Avg. λ error | Avg. measure r |
|---|---|---|---|---|---|---|
| 1 | 10 | [-2; 2] | 1.00e-09 | 134 | 0.000000 | 1.67e-06 |
| 2 | 10 | [-50; 50] | 1.00e-09 | 136 | 0.000000 | 1.40e-06 |
| 3 | 30 | [-2; 2] | 1.00e-09 | 1398 | 0.000000 | 2.66e-06 |
| 4 | 30 | [-50; 50] | 1.00e-09 | 1467 | 0.000000 | 2.57e-06 |

### Key takeaways

- **Eigenvalue accuracy is essentially exact.** The average eigenvalue error rounds to `0.000000` in all four configurations — the Jacobi method recovers the true eigenvalues used to construct the test matrix to very high precision.
- **The accuracy measure `r`** (how well `A·T ≈ T·Λ` holds) stays in the `1e-6` range across all configurations, showing the eigen-decomposition is self-consistent regardless of matrix size or eigenvalue spread.
- **Rotation count scales strongly with matrix size, not with eigenvalue range.** Going from `N = 10` to `N = 30` increases the average number of rotations roughly **tenfold** (from ~135 to ~1400+), while widening the eigenvalue range from `[-2, 2]` to `[-50, 50]` barely changes the rotation count for a fixed `N`. This matches the theoretical expectation that Jacobi's convergence rate depends on the number of off-diagonal elements needing to be zeroed (which grows roughly as `N²`), not on the magnitude of the eigenvalues themselves.
