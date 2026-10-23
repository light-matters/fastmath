(ns fastmath.api.v2.algebra.fast
  "Explicit functions for performance-sensitive code: no multimethod dispatch
  and no domain coercion (compare `general`). Currently complex dense
  matrices (`ComplexDense`, EJML backend).

  Entries are complex: [re im] pairs on the way in (real numbers are
  accepted too) and on the way out (`->vectors`). Coefficients for `scale`
  and `linear-combination` are [re im] pairs or real numbers.

  All operations return new matrices; inputs are never modified. `solve`
  reuses one LU solver per thread and matrix size, so it is safe to call
  from many threads at once (e.g. `pmap`)."
  (:refer-clojure :exclude [identity])
  (:require
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as cd]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.solve :as solve]
   [fastmath.protocol.algebra.structure.additive.semigroup :as asg]
   [fastmath.protocol.algebra.structure.coordinate.complex :as cc]
   [fastmath.protocol.algebra.structure.module :as module]
   [fastmath.protocol.representation.d2 :as d2]))

;; --- construction ---------------------------------------------------------------

(defn zeros
  "nrows x ncols zero matrix."
  [nrows ncols] (cd/complexdense (long nrows) (long ncols)))

(defn identity
  "n x n identity."
  [n]
  (let [d (double-array (* 2 (long n)))]
    (dotimes [i n] (aset d (* 2 i) 1.0))
    (cd/complexdense d)))

(defn diagonal
  "Diagonal matrix from entries ([re im] or real)."
  [entries]
  (cd/complexdense (double-array (mapcat #(if (number? %) [% 0.0] %) entries))))

(def <-rows
  "From rows of entries: [re im] pairs, real numbers or complex numbers."
  cd/<-rows)

(def <-interleaved
  "From nrows, ncols and a row-major interleaved double[] (copied)."
  cd/<-interleaved)

;; --- access ---------------------------------------------------------------------

(def shape d2/shape)
(def element
  "Entry (i, j) as a fastmath complex number."
  d2/element)
(def ->vectors
  "Rows of [re im] pairs."
  d2/->vectors)
(def ->interleaved
  "Row-major interleaved double[] (a copy)."
  cd/->interleaved)

;; --- arithmetic -----------------------------------------------------------------

(def add asg/add)
(def subtract emat/subtract)
(defn scale
  "c A for a complex coefficient c ([re im] or real)."
  [A c] (module/scale A (if (number? c) [c 0.0] c)))
(def linear-combination
  "sum_k c_k A_k in one pass (one allocation). Terms are dense matrices or
  `sparse` ones (only their nonzero entries are visited)."
  cd/linear-combination)
(def sparse
  "Sparse copy of a matrix (nonzero entries), as a `linear-combination` term.
  Use it for fixed operators that are combined many times."
  cd/sparse)
(def add-diagonal
  "A + c diag(v) for square A and a length-n vector v (n x 1 or 1 x n)."
  cd/add-diagonal)
(def multiply
  "Matrix product A B."
  emat/multiply)
(def kronecker
  "Kronecker product A (x) B."
  emat/kronecker)
(def transpose emat/transpose)
(def adjoint
  "Conjugate transpose."
  cmat/adjoint)
(def conjugate
  "Elementwise conjugate."
  cc/conjugate)

;; --- linear systems -------------------------------------------------------------

(def solve
  "X with A X = B. Square A: LU (per-thread reused solver; throws if
  singular). Otherwise least squares."
  solve/solve)

(def solve-shifted
  "X with (A + c diag(v)) X = B: like (solve (add-diagonal A c v) B) without
  allocating the shifted matrix (per-thread scratch), for many shifts of
  one A."
  cd/solve-shifted)
