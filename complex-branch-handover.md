# `complex` branch — handover notes

Working notes for the complex-number/matrix work on the `complex` branch
(tracks `origin/complex` at `light-matters/fastmath`). Kept up to date as
the work progresses so a new session (human or AI) can pick this up
without replaying the whole history. See also `complex-matrix-notes.org`
for the original design exploration/FAQ this work grew out of.

## Context

Branch originates from Thomas Clark's Clojurists Together proposal to add
complex-number support to fastmath - started as "add complex matrices",
grew into a broader generic-algebra API (the `fastmath.protocol.algebra.*`
structure protocol stack: semigroup -> monoid -> group -> ring -> field,
plus the matrix protocols: `D2`, `MatrixExtra`, `SquareMatrix`,
`MatrixSolve`, `MatrixPredicate`, `ComplexMatrix`/`ComplexPredicate`).

Three real-matrix backends now exist side by side: EJML (the original),
ojAlgo, and `commons` (fastmath's own pre-existing real matrix
implementation at `fastmath.matrix.dense.real.original`, backed by Apache
Commons Math's `RealMatrix`, wired in as a third backend rather than
reimplemented). Complex matrices only have EJML and ojAlgo so far - the
old implementation has no complex-number support to wire in. None of the
three are wired into `fastmath.algebra.object.matrix.create`'s
backend-selection dispatcher - that dispatcher still only ever builds
EJML instances. Making it backend-selectable is a real open decision, not
done.

## Completed work (chronological)

1. **Fixed broken/dead code** blocking the test suite: stale require paths
   left over from a `protocol.algebra.*` -> `protocol.algebra.structure.*`
   rename, a dangling symbol in `fastmath/constant.clj`, a dash/underscore
   filename mismatch that silently aborted the *entire* `lein test` run.
2. **Rewrote the complex scalar** (`fastmath.algebra.object.number.complex.primitive`,
   `ComplexNumber`) as two plain doubles instead of wrapping EJML's
   `Complex_F64` - decouples the scalar from any one matrix backend so
   each backend converts to/from it at its own boundary. Also fixed a
   missing `hashCode` (the type overrode `.equals` but not `.hashCode`,
   violating the Java contract - broke sets/map keys for complex numbers).
3. **`RealDense` brought to parity with `ComplexDense`** (EJML backend):
   added `diagonal`, `fmap`, `->vectors`, `Seqable`, `SquareMatrix`,
   `MatrixSolve`, `MatrixPredicate`. Found and fixed real bugs while
   writing tests for the previously-untested parity surface:
   - `RealDense.kronecker`/`.inner` only ever worked for column vectors
     (wrong index arithmetic / read raw fields that didn't exist).
   - `ComplexDense.multiply--e` called the wrong EJML `elementMultiply`
     overload (the 4-arg scalar-multiply form, not the 3-arg
     matrix-multiply form) - crashed with a `ClassCastException`.
   - `symmetric?`/`normal?` on both backends compared raw EJML matrices
     with `=`, but `DMatrixRMaj`/`ZMatrixRMaj` don't override `.equals`
     (identity comparison only) - these predicates were silently always
     `false`. Fixed to use `MatrixFeatures_*DRM/isEquals` with tolerance.
   - `cholesky` on both backends passed `M` itself (not a copy) to
     `.decompose`, which EJML documents as potentially mutating its
     input - silently corrupted the receiver. Fixed to decompose a copy.
   - `fmap` on `ComplexDense` mutated `M` in place instead of operating on
     a copy like every other method - risky since `fmap` was being called
     on `def`s shared across test namespaces.
4. **Extracted backend-independent matrix ops** into
   `fastmath.algebra.function.matrix` (`square?`, `symmetric?`, `normal?`,
   `singular?`, `hermitian?`, `inner`, `kronecker-rows`) - written once
   against protocol methods (`d2/shape`, `d2/element`, `emat/transpose`,
   `emat/multiply`, `sq/determinant`, `cmat/adjoint`) rather than any
   backend's raw arrays, mirroring how `fastmath.algebra.function.field`
   already derives `subtract`/`divide`/`square` generically. Bulk/perf-critical
   ops (add, multiply, transpose, the decompositions) stay backend-native.
5. **Built the ojAlgo backend** (`RealDense`/`ComplexDense` wrapping
   `MatrixR064`/`MatrixC128`) on top of that smaller surface, using
   ojAlgo's high-level `BasicMatrix` API (`Access2D/wrap` +
   `FACTORY.makeWrapper`) rather than the low-level
   `GenericStore`/`PhysicalStore` layer the abandoned
   `fastmath.matrix.dense.scratch` experiments got stuck on (that layer
   needs reflection to work around a Java access restriction).
6. **Built `fastmath.benchmark.matrix-backends`**, comparing EJML vs
   ojAlgo through fastmath's own protocol API (not the raw libraries)
   across both domains and sizes 8-256. See results below.
7. **Wired fastmath's own pre-existing real matrix implementation** in as
   a third backend (`fastmath.algebra.object.matrix.rectangular.real.commons`,
   backed by Apache Commons Math `RealMatrix` via the existing
   `fastmath.matrix.dense.real.original` functions - `old/add`, `old/mulm`,
   `old/det`, `old/cholesky`, etc.) rather than reimplementing anything.
   Same protocol surface as EJML/ojAlgo's `RealDense`, same generic ops
   from point 4 above reused unmodified (third independent confirmation
   they're genuinely backend-independent). `row`/`column`/`diagonal`/
   `multiply--v`/`solve` needed small conversions since the old API
   returns Commons Math `RealVector`s for those, not matrices. One real
   inconsistency found and fixed: the old API's `Seqable` would have
   produced nested per-row lists (`((1.0 2.0) (3.0 4.0))`), not the flat
   sequence (`(1.0 2.0 3.0 4.0)`) EJML/ojAlgo's `RealDense` produce -
   flattened for consistency. No complex-number equivalent exists in the
   old implementation, so this backend is real-only.

## Bugs found in external libraries (ojAlgo 56.1.0)

- **`BasicMatrix.conjugate()` is actually the Hermitian adjoint** (conjugate
  *transpose*), not elementwise conjugation. Confirmed empirically: calling
  it on a 1x2 matrix returns a 2x1 result. Worked around: `ComplexCoordinate/conjugate`
  is now an explicit elementwise loop; `cmat/adjoint` just calls
  `.conjugate()` directly (which is what it actually computes).
- **`BasicMatrix.isHermitian()` doesn't check that diagonal entries are
  real.** A diagonal matrix with a non-zero imaginary diagonal (e.g.
  `5+10i` on the diagonal) passes `.isHermitian()` as `true`, which is
  mathematically wrong. Confirmed EJML's `MatrixFeatures_ZDRM.isHermitian`
  does *not* have this problem (correctly rejects the same input), so
  it's specific to ojAlgo. Worked around: added a verified-correct
  `fastmath.algebra.function.matrix/hermitian?` (`A == adjoint(A)`), used
  for the ojAlgo backend only; EJML backend keeps its own correct native
  call.
- **`BasicMatrix.transpose()` returns a lazy view**, not an eager copy -
  confirmed by timing (near-zero regardless of size). Not a bug, but a
  real behavioral difference the benchmark now shows explicitly (see
  `:transpose-view` vs `:transpose` below).
- **`BasicMatrix` caches `determinant`/`cholesky`/`solve` on the
  (immutable) matrix instance.** Confirmed empirically: a repeat call on
  an *unchanged* matrix is 10-400x faster than the first call; EJML shows
  no such effect. Not a bug - a genuine advantage of ojAlgo's immutable
  design for workloads that query the same matrix more than once. See
  benchmark design notes below for how this is now measured and reported
  rather than designed away.

## Benchmark design notes

- **The sandbox this runs in is heavily CPU-throttled** - confirmed: a
  trivial double-array fill with *zero* library code involved takes the
  same order-of-magnitude time (~850ms for a 256x256 fill) as the actual
  library calls. Absolute times/GFLOP-s from this environment are **not**
  representative of real hardware. Only the EJML-vs-ojAlgo *ratio* at a
  given size is meaningful.
- **Had to fix ~35 reflection warnings per ojAlgo backend file** before
  the comparison meant anything: calling `doubleValue(int, int)` with
  `long` arguments (that overload doesn't exist, only `(int, int)` and a
  single-arg `(long)` do) silently falls back to Java reflection instead
  of a compile error. Measured directly: 78.7 seconds vs 41 milliseconds
  for the same loop with vs. without the correct hint (~1900x). Left
  unfixed, this would have made ojAlgo look catastrophically slower for
  reasons having nothing to do with the library itself. All four backend
  files (`real`/`complex` x `ejml`/`ojalgo`) are now warning-clean or
  close to it (a few low-severity "boxed math" warnings remain in paths
  the benchmark doesn't exercise - see e.g. `real/ejml.clj:263`,
  `complex/ejml.clj:213-214,345`).
- **`determinant`/`cholesky`/`solve` each have two benchmark variants**:
  the plain one builds a fresh matrix per trial (`measure-fresh`, cost of
  computing on data the backend hasn't seen before); the `-repeated` one
  queries the *same* matrix repeatedly on purpose, to surface and report
  ojAlgo's caching advantage explicitly rather than have it either
  silently inflate the "fresh" numbers or get averaged away.
- Also fixed a real, independent bug while building this: the pre-existing
  `fastmath.benchmark.apache-ejml/measure` function had backwards `nth`
  arguments (`(->> arr seq sort (nth (quot trials 2)))` threads the sorted
  array into `nth`'s *index* position) - threw `ClassCastException` on
  every call, never caught before because nothing used it.
- **A third real backend, `real/commons.clj`, is wired in too** (fastmath's
  own pre-existing implementation, backed by Apache Commons Math). It's
  absent from every `:complex` row (no complex equivalent exists) -
  `domain-backend-pairs` derives the iteration from `backend-info`'s own
  keys rather than a blind domain x backend cartesian product, so an
  asymmetric backend doesn't need special-casing anywhere else.
- **Found and fixed a real inefficiency in EJML's `->arrays`** while
  comparing it to `commons`: at n=256 `commons`'s materialized `transpose`
  was ~1000x "faster" than EJML's, which was too dramatic to be a genuine
  algorithmic difference. Root cause: `Array2DRowRealMatrix.getData()`
  (which `commons` calls directly) bulk-copies row-by-row internally,
  while EJML's `D2/->arrays` (code *we* wrote) read one `.get(i, j)` cell
  at a time. Confirmed EJML's `DMatrixRMaj` stores data as a flat
  row-major `double[]` already, so rewrote `->arrays` to use
  `System/arraycopy` per row instead - dropped from ~850ms to ~0.15ms at
  n=256 in this throttled sandbox. ojAlgo's `BasicMatrix` has no public
  bulk-access equivalent, so its `->arrays` is stuck doing it cell-by-cell
  - that's now a genuine, reportable cost of its abstraction level, not
  an artifact.

## Key results so far (relative only - see environment caveat above)

- EJML is faster than ojAlgo across the board (add, multiply, determinant,
  solve, norm), by a margin that grows with size, and **the gap is
  consistently wider for complex matrices than for real ones** at nearly
  every operation.
- `cholesky` is the one interesting exception: for **real** matrices
  ojAlgo is competitive and sometimes *faster* (e.g. beats EJML at n=32
  and n=128 in one run). For **complex** matrices that inverts hard -
  EJML's complex Cholesky was ~28x faster than ojAlgo's at n=128 in one
  run. Comparing each backend against itself: EJML's complex Cholesky is
  only ~1.4x slower than its own real Cholesky, while ojAlgo's complex
  Cholesky was ~51x slower than its own real Cholesky in that same run -
  not just "complex costs more arithmetic", but ojAlgo's complex
  decomposition path specifically carrying some disproportionate
  inefficiency its real path doesn't have. This is circumstantial but
  consistent with the `conjugate()`/`isHermitian()` findings above -
  ojAlgo's complex-matrix support looks like the less mature part of the
  library. Worth re-running before relying on it for anything
  complex-and-performance-sensitive (e.g. the FFT/Fourier-analysis
  integration from the original proposal).
- The `cholesky-repeated` caching effect showed up much more clearly for
  `determinant` than for `cholesky` itself in testing so far (determinant:
  clean 10-400x; cholesky: sometimes only ~1.5x) - possibly because our
  `solve/cholesky` wraps the decomposition result in a fresh matrix object
  each call even when the underlying decomposition is cached internally.
  Not fully explained yet; noted here rather than chased further given
  the environment's noise floor.
- **Real-only, three-way comparison (EJML / ojAlgo / commons) at n=128**,
  after the `->arrays` fix above: EJML fastest almost everywhere. `commons`
  is consistently the clear second choice for add/transpose/norm/
  determinant/cholesky/solve (typically 1.1-5x slower than EJML) - except
  `multiply`, where `commons` was actually the *slowest* of the three
  (0.035s vs EJML's 0.002s and ojAlgo's 0.017s), consistent with Apache
  Commons Math's matrix multiply being a known-unoptimized naive loop (no
  blocking/vectorization) rather than something specific to this wiring.
  ojAlgo trails on everything except the `-repeated` (caching) variants,
  where it wins dramatically, and `cholesky`, where all three land within
  a similar range.

## Commits on this branch (pushed to `origin/complex`)

- `f1db339` "Basic bug fixes and RealDense/ComplexDense parity"
- `497df6b` "Extract backend-independent matrix ops, add an ojAlgo backend"
- `7121fc5` "Add EJML vs ojAlgo matrix backend benchmark, fix reflection
  overhead" - the benchmark namespace, the `measure` bugfix, the
  type-hint/reflection fixes, and this handover doc itself.
- `e4fc563` - small follow-up fixing this doc's own commit list after the
  above push (meta, but illustrates: update this file *after* pushing,
  not before, or it goes stale immediately).

## Uncommitted as of this note

- `src/fastmath/algebra/object/matrix/rectangular/real/commons.clj` (new)
  and `test/.../real/commons_test.clj` (new) - the third backend
  described above.
- `src/fastmath/benchmark/matrix_backends.clj` - wired `commons` into the
  benchmark (`domain-backend-pairs` change).
- `src/fastmath/algebra/object/matrix/rectangular/real/ejml.clj` - the
  `->arrays` bulk-copy fix described above.

## Not done / explicitly out of scope so far

- None of the three real backends wired into `create.clj`'s dispatcher
  (separate decision - see Context above).
- `test/fastmath/matrix_test.clj` (tests for the old implementation
  directly, not through the new API) still requires `fastmath.matrix`,
  which no longer exists under that name - it was renamed to
  `fastmath.matrix.dense.real.original` before this session (commit
  `c7b247d` was the last touch on the test file, predating the rename).
  Fails with `FileNotFoundException`. Pre-existing, unrelated to the work
  in this doc, left alone - found while confirming the new `commons`
  backend didn't break the old implementation's own test coverage.
- `LBFGSBJava` git submodule is uninitialized on this machine
  (`git submodule status` shows it unchecked-out) - blocks a full,
  unfiltered `lein test` run with `ClassNotFoundException:
  org.generateme.lbfgsb.Parameters`. Pre-existing, unrelated to this
  branch's work, left alone.
- A few low-severity boxed-math warnings left in non-hot-path code (see
  benchmark design notes above).
- Three self-labeled scratch/draft files identified early on were left
  alone rather than fixed, since "fixing" them would mean guessing at
  unfinished design: `src/fastmath/algebra/object/matrix/rectangular/complex/ejmlDraft.clj`
  (marked `^:deprecated`, docstring says "temporary and dev-only"),
  `src/fastmath/api/v1/algebra/matrix.clj` (signature sketches with no
  bodies), `src/fastmath/matrix/dense/scratch.clj` and
  `src/fastmath/matrix/dense/ojalg-constructors--bench.clj` (both
  docstring-labeled "a playground for working through ideas").
