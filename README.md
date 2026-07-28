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

## Output

Running the program prints a summary table to the console, one row per `(N, [λ_min, λ_max])` configuration, with columns for:

- matrix size `N`,
- the eigenvalue range tested,
- the tolerance (`eps`) reached,
- average number of rotations to converge,
- average maximum eigenvalue error (computed vs. true),
- average accuracy measure `r`.

By default, `Main.cpp` runs the experiment for sizes `{10, 30}` crossed with eigenvalue ranges `{[-2, 2], [-50, 50]}`, testing convergence tolerances `{1e-5, 1e-7, 1e-9}` for each configuration.

## What the experiment measures

For each combination of **matrix size** and **eigenvalue range**, the experiment checks how the Jacobi method's:

- **convergence speed** (number of rotations), and
- **eigenvalue accuracy** (deviation from the known true eigenvalues),

scale as the problem gets larger or the eigenvalues span a wider range — giving an empirical picture of the method's practical behavior beyond its theoretical guarantees.
