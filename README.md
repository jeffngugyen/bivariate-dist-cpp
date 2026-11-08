# Bivariate / CSS code distance solver

An exact, exhaustive solver for small CSS quantum codes. It uses packed binary
vectors, Gaussian elimination over GF(2), row-space membership tests, and optional
OpenMP parallel search. A Python interface is available through pybind11.

## Quick start

Requirements: a C++20 compiler and GNU Make. Python 3.10+ is needed for the tests;
the standalone solver has no Python or pybind11 dependency.

```sh
make
make example
make test
```

`make example` uses `H_X = H_Z = [1 1 1 1]` and prints:

```text
d_Z=2 d_X=2 distance=2
```

This is a four-qubit code with two logical qubits and distance 2. Each kernel
contains the even-weight vectors. The stabilizer row space contains only `0000`
and `1111`, so the smallest nontrivial logical vectors have weight 2.

The executable is `build/dist`; there is no precompiled platform-specific binary.
The old top-level `H_X.txt` and `H_Z.txt` research inputs are retained, but are not
used for the quick start because their exhaustive searches can be expensive.

## Run your own matrices

```sh
./build/dist path/to/H_X.txt path/to/H_Z.txt 1 2
```

The last two arguments are positive integers `l` and `m`; each matrix row must
contain exactly `n = 2*l*m` space-separated binary entries. Files can have different
row counts, redundant checks, or no rows. Blank lines are ignored. Rows of `H_X`
and `H_Z` denote X and Z stabilizer checks respectively. The solver verifies
`H_X H_Z^T = 0` over GF(2).

The CLI reports `d_Z`, `d_X`, and their minimum. It prints an error and exits
nonzero for malformed matrices, noncommuting checks, unsupported sizes, or a
code with no logical qubits (whose logical distance is undefined).

### Parallel and serial builds

On Linux, the default build enables OpenMP (`-fopenmp`). To use a compiler without
OpenMP, including the default Apple Clang build on macOS:

```sh
make OPENMP=0
make test OPENMP=0
```

macOS defaults to `OPENMP=0`. If you have an OpenMP-capable compiler, set `CXX`,
`OPENMP=1`, and, if needed, `OPENMP_FLAGS`. Set `OMP_NUM_THREADS` to control the
parallel search, for example `OMP_NUM_THREADS=4 make example`.

## Python interface (optional)

Use the repository's pinned pybind11 headers:

```sh
git submodule update --init extern/pybind11
make python
make python-test
```

Alternatively, install `pybind11` into your Python environment. Python development
headers are also required. Set `PYTHON=/path/to/python` to select an interpreter.
The extension is built in the repository root; run the example from there:

```python
import dist_py

checks = [[1, 1, 1, 1]]
d_z, d_x = dist_py.compute_dist(checks, checks, l=1, m=2)
assert (d_z, d_x) == (2, 2)
```

The function releases the Python GIL only while executing native computation.
Invalid matrices and unsupported sizes raise `ValueError`.

**Compatibility:** `compute_dist(X, Z, l=4, m=7)` preserves the original numerical
tuple order: first the minimum weight in `ker(X) \ row(Z)`, then in
`ker(Z) \ row(X)`. With X/Z stabilizer naming, that is `(d_Z, d_X)`.
The old source called these `dx` and `dz`; these labels were reversed. The CLI
labels and Python documentation now state the convention explicitly.

## Why the row-space test works

1. XOR is vector addition over GF(2).
2. Gaussian elimination produces a basis for each kernel. Enumerating all
   nonzero XOR combinations visits each nonzero kernel vector once.
3. The stabilizer basis has distinct highest set bits, ordered highest first.
   For each basis row, XOR it into the candidate when doing so clears the
   candidate's bit at that row's leading position.
4. Later rows cannot reintroduce an earlier leading bit. The remainder is zero
   exactly when the candidate is a linear combination of the stabilizer rows.
5. Keep kernel vectors with a nonzero remainder and minimize their **original**
   Hamming weight. The reduced remainder is only a membership test, not the
   logical operator whose weight should be measured.

For example, stabilizer rows `1100` and `0011` reduce `1111` to zero, while `1000`
leaves a nonzero remainder. The first belongs to their row space; the second does
not. This replaces materializing every element of the stabilizer space with a
compact basis and a sequence of XOR operations.

A useful explanation exercise: why would measuring the remainder's weight give
the wrong objective? Also, why does the order of the distinct leading bits matter?

## Correctness checks and limits

`make test` checks the executable, using both `OMP_NUM_THREADS=1` and `4`.
`make python-test` runs the same mathematical tests against the Python API.
The suite includes:

- the known four-qubit distance-2 example;
- asymmetric distances to catch swapped tuple order;
- empty, redundant, and zero checks;
- 200 deterministic random commuting CSS instances with 2, 4, 6, or 8 columns,
  compared against an independent Python oracle that enumerates bit strings
  directly (no Gaussian elimination or optimized membership routine);
- malformed inputs, noncommuting checks, zero logical qubits, and the 64-bit
  enumeration boundary, including a large second nullspace.

GitHub Actions runs serial and OpenMP builds on Linux and a serial build on macOS,
including the Python extension. Local sanitizer check:

```sh
make test BUILD_DIR=build-sanitize OPENMP=0 \
  CXXFLAGS='-O1 -g -Wall -Wextra -Wpedantic -fsanitize=undefined -fno-sanitize-recover=all'
```

The representation stores at most 128 columns. Each nullspace dimension must be
less than 64 because enumeration uses a 64-bit mask. This is a representation
limit, **not a practical performance promise**: the search visits `2^r - 1`
candidates for a nullspace of dimension `r`, and becomes impractical well before
64. OpenMP reduces wall time without changing that exponential growth. Passing
small-instance tests does not prove correctness or feasibility for every input.
