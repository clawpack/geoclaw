# SGN Boussinesq — PETSc-free OpenMP linear solver (`bouss_solver = 1`)

This directory implements the implicit SGN (Serre–Green–Naghdi) Boussinesq update
for GeoClaw. Alongside the existing PETSc path it now has a **self-contained,
OpenMP-threaded linear solver** that needs neither PETSc nor MPI.

## Solver selection (`bouss.data` / `setrun.py`)

- `bouss_solver = 1` — OpenMP GMRES + directional block Gauss–Seidel preconditioner (this work)
- `bouss_solver = 2` — Pardiso (expired / unsupported)
- `bouss_solver = 3` — PETSc (retained; unchanged)

`isolver = bouss_solver` throughout the code.

## The linear system

The implicit update solves a sparse system stored in **CSR (compressed row)**
format with a **2×2 block** structure (`MatSetBlockSize 2` on the PETSc side).
Unknowns are interleaved per cell: `u0, v0, u1, v1, …` (even rows = x-momentum
`u`, odd rows = y-momentum `v`). Globally cells are numbered i-fastest.

```
A = [ I - D11    -D12   ]
    [  -D21    I - D22   ]
```

- **`I - D11`** (u-block) couples `U(i,j) ↔ U(i±1,j)` only ⇒ a direct sum of
  independent **1-D tridiagonals in x**.
- **`I - D22`** (v-block) couples `V(i,j) ↔ V(i,j±1)` only ⇒ independent
  **tridiagonals in y**.
- `-D12`, `-D21` are the u↔v cross-coupling blocks.

## How the solver works

1. `prepBuildSparseMatrixSGNcrs` assembles the CSR matrix for the *union* of
   Bouss grids at the level (threaded over grids), then `compressOut` compacts it.
2. `build_line_decomp` (called twice — u-lines parity 0, v-lines parity 1) walks
   the matrix graph and extracts the maximal tridiagonal **lines**. Each field
   dof lands in **exactly one line** (a true partition), which is what makes the
   threaded line solve race-free.
3. `gmres_solve` runs **right-preconditioned restarted GMRES(50)**, `rtol 1e-9`,
   up to 4 cycles.
4. Preconditioner = **directional block Gauss–Seidel**: solve the u-block
   (tridiagonal lines in x, *exact* Thomas), update, solve the v-block. Forward
   sweep by default; symmetric (both cross blocks) available via
   `bgs_symmetric = .true.`.

## New source files

| File | Role |
|------|------|
| `bouss_tridiag_module.f90` | line decomposition (`build_line_decomp`), threaded Thomas line solve (`apply_line_solve`/`solve_one_path`), block-GS preconditioner (`apply_block_gs`), CSR sub-matvec |
| `bouss_gmres_module.f90` | OpenMP GMRES (`gmres_solve`), CSR matvec, threaded dot/axpy/scale |
| `bouss_pcshell.f90` | optional PETSc `PCShell` wrapper for the tridiagonal PC (only compiled with `HAVE_PETSC`) |
| `test_tridiag.f90`, `test_gmres.f90`, `test_compress.f90` | standalone unit tests |

Modified: `implicit_update_bouss_2Calls.f90`, `prepBuildSparseMatrixSGNcrs.f90`,
`compressOut.f`, `amr2.f90`, `amr_module.f90`, `stst1.f`, `bouss_module.f90`,
`Makefile.bouss`.

## Building and running

- Builds **with or without `-fopenmp`** (serial fallback runs one thread).
- Builds **with or without `HAVE_PETSC`** — `bouss_solver = 1` needs neither PETSc
  nor MPI. Run single-process with `OMP_NUM_THREADS` set.
- Unit tests (not in the Makefile — compile by hand):
  ```
  gfortran -O2 -fopenmp bouss_tridiag_module.f90 bouss_gmres_module.f90 test_gmres.f90 -o t && ./t
  ```
  Tests are **bit-identical at 1 and N threads** — use them to check any change.

## Performance (crater_westport, 10 threads)

- `bouss_solver = 1`: **~219 s** linear solve — faster than PETSc redundant-LU
  (~230 s) and GAMG (~352 s), PETSc-free.
- Best with **`max1d = 60`**: more/smaller grids ⇒ better grid-parallelism in the
  (grid-threaded) matrix build. `compressOut` and `line-decomp` are per-level and
  grid-count independent.
- Diagnostics printed at end of run (`bouss_solver=1`): `line-decomp` and
  `compressOut` timers, and the longest tridiagonal line seen (`maxLineLen`).

## Hard-won lessons / gotchas (read before editing)

1. **OpenMP without a hard dependency.** Do **not** `use omp_lib`. Declare
   `omp_get_max_threads`/`omp_get_thread_num` as `integer` externals and call them
   **only on `!$` sentinel lines**, defaulting thread count/id to 1/0. (Matches
   `amr2`/`advanc`/`flagger`.) This is what lets the code build without `-fopenmp`.

2. **Persistent `save` work buffers.** The big arrays (GMRES Krylov basis `V`
   ≈ 80 MB, block-GS `Ay`/`rw`, Thomas `thom_cp`, `compressOut` scratch) are
   module-`save`, grown-as-needed, reused across solves. `gmres_solve` /
   `apply_block_gs` are entered **only on the master thread** (`implicit_update`
   is called serially per level), so `save` is thread-safe. Eliminating per-call
   allocation of `V` saved ~27 s on crater.

3. **Worker-thread stacks.** OpenMP *worker* threads get a small stack
   (`OMP_STACKSIZE`); only the *primary* thread has the big main stack. Automatic
   (stack) arrays sized by problem size overflow worker stacks on large cases even
   when a single-thread run is fine. That is why the Thomas scratch is on the
   **heap** (`thom_cp(maxline, nthreads)`), not an automatic array. Tridiagonal
   lines cross grid boundaries, so line length is **not** bounded by `max1d`.
   On **macOS**, `limit stacksize unlimited` is silently capped (~64 MB main
   thread), and `OMP_STACKSIZE` only affects workers — so bumping the stack may
   not "fix" an overflow.

4. **The line decomposition must stay a true partition** (each dof in exactly one
   line). The parallel line solve writes disjoint dofs and reads read-only shared
   data; correctness depends on the partition. `build_line_decomp` guarantees it
   via a `visited` flag (max degree 2 asserted).

5. **Debugging NaNs.** A *quiet* NaN from uninitialized memory propagates silently
   through `+ − × ÷` and only traps at an **ordered comparison** under
   `-ffpe-trap=invalid` — so the trap fires where it's *used*, not created. Use
   **`-finit-real=snan`** (with `-ffpe-trap=invalid`) to trap at first *touch* and
   find the source. (A global-case NaN was tracked this way — it was in the
   coarse ghost-cell fill, **not** the solver.)

6. **Debugging on macOS.** `gfortran -fbacktrace` prints only raw addresses, and
   ASLR makes them unsymbolizable after the fact. Use **lldb**
   (`lldb --batch -o run -k bt -- ./xgeoclaw`) for symbolized traces.

## Provenance

GMRES (Saad & Schultz 1986), the Thomas tridiagonal solve (1949), graph traversal
for line extraction, and stream compaction are all standard algorithms, written
from scratch here — not copied from PETSc, SPARSKIT, Numerical Recipes, etc.
