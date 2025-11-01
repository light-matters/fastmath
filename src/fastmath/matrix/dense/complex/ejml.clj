(ns fastmath.matrix.dense.complex.ejml
  "Implementing real matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml].

   Legend:
   m - matrix type
   i - row index
   j - column index

   TODO:
   - Should we implement a real vector type as well (consistency)?
  "

  (:require [fastmath.vector :as v]
            [fastmath.core :as fm]
            [fastmath.protocols.matrix2 :as proto])
  (:import
   (org.ejml.data ZMatrixRMaj)
   (org.ejml.dense.row
    CommonOps_DDRM
    MatrixFeatures_DDRM
    NormOps_DDRM)
   (org.ejml.dense.row.factory
    DecompositionFactory_DDRM)

   (org.ejml.interfaces.decomposition CholeskyDecomposition_F64)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^:private ^double tolerance--default 1.0e-10)

;; ==================================================
;; Functions 
;; ==================================================
(defn- ^ZMatrixRMaj extract
;; TODO: Add to type?
  "Extracts a submatrix m[i0:i1, j0:j1) into a new ZMatrixRMaj.
   Indices are zero-based and end-exclusive."
  [^ZMatrixRMaj m
   ^longs [i0 i1]
   ^longs [j0 j1]]
  (let [^ZMatrixRMaj out (ZMatrixRMaj. (- ^long i1 ^long i0) (- ^long j1 ^long j0))]
    (CommonOps_DDRM/extract m i0 i1 j0 j1 out 0 0)
    out))

(defn ^ZMatrixRMaj extract-multiple
  ;; TODO: efficiency implementaton
  "Reuses the destination memory for multiple access calls."
  [])

(deftype ComplexDense [^ZMatrixRMaj M]
  proto/MatrixReal

  ;; -------- Transformations --------
  (->seq [_]
    (seq (.data M)))

  (->array [_]
    (let [d (.data M) out (double-array (alength d))]
      (System/arraycopy d 0 out 0 (alength d))
      out))

  (->array--2d [_]
    (let [r (.numRows M) c (.numCols M)
          out (make-array Double/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-double out i j (.get M i j))))
      out))

  (->array--float [_]
    ;; TODO: faster
    (float-array (map float (seq (.data M)))))

  (->array--2dfloat [_]
    (let [r (.numRows M) c (.numCols M)
          out (make-array Float/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-float out i j (float (.get M i j)))))
      out))

  ;; -------- Retrieval --------
  (columns [_]
    (mapv (fn [^long ic]
            (extract M [0 (.numRows M)] [ic (+ 1 ic)]))
          (range (.numCols M))))

  (rows [_]
    (mapv (fn [^long ir]
            (extract M [ir (+ ir 1)] [0 (.numCols M)]))
          (range (.numRows M))))

  (diagonal [_]
    (let [n (min (.numRows M) (.numCols M))
          out (ZMatrixRMaj. n 1)]
      (CommonOps_DDRM/extractDiag M out)
      (ComplexDense. out)))

  (element [_ i j]
    (.get M (long i) (long j)))

  (column [_ j]
    (extract M [0 (.numRows M)] [j (+ 1 (long j))]))

  (row [_ i]
    (extract M [i (+ 1 (long i))] [0 (.numCols M)]))

  (num-rows [_] (.numRows M))
  (num-cols [_] (.numCols M))
  (shape    [_] [(.numRows M) (.numCols M)])

  ;; -------- Operations --------
  (add [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/add 1.0 M 1.0 B out)
      (ComplexDense. out)))

  (add--s [_ s]
    (let [A (.copy M) d (.data A) s (double s)]
      (dotimes [k (alength d)]
        (aset-double d k (+ (aget d k) s)))
      (ComplexDense. A)))

  (sub
    [this]
    (let [A (.copy M)]
      (CommonOps_DDRM/scale -1.0 A)
      (ComplexDense. A)))

  (sub [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/add 1.0 M -1.0 B out)
      (ComplexDense. out)))

  (map--m [_ f]
    (let [A (.copy M)]
      (dotimes [i (.numRows A)]
        (dotimes [j (.numCols A)]
          (.set A i j (double (f (.get A i j))))))
      (ComplexDense. A)))

  (mul
    [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols B))]
      (CommonOps_DDRM/mult M B out)
      (ComplexDense. out)))

  (mul [this tA? other tB?]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          rA (if tA? (.numCols M) (.numRows M))
          cA (if tA? (.numRows M) (.numCols M))
          rB (if tB? (.numCols B) (.numRows B))
          cB (if tB? (.numRows B) (.numCols B))]
      (when (not= cA rB)
        (throw (ex-info "Incompatible shapes for matmul"
                        {:A [rA cA] :B [rB cB] :tA? tA? :tB? tB?})))
      (let [out (ZMatrixRMaj. rA cB)]
        (cond
          (and tA? tB?) (CommonOps_DDRM/multTransAB M B out)
          tA?           (CommonOps_DDRM/multTransA  M B out)
          tB?           (CommonOps_DDRM/multTransB  M B out)
          :else         (CommonOps_DDRM/mult        M B out))
        (ComplexDense. out))))

  (mul--e [this other]
    (let [A (.copy M)
          ^ZMatrixRMaj B (.M ^ComplexDense other)]
      (CommonOps_DDRM/elementMult A B)
      (ComplexDense. A)))

  (mul--v [this v]
    ;; v is expected to be (n×1) ComplexDense
    (let [^ZMatrixRMaj x (.M ^ComplexDense v)
          out (ZMatrixRMaj. (.numRows M) 1)]
      (CommonOps_DDRM/mult M x out)
      (ComplexDense. out)))

  (mul--vt [this v]
    ;; v treated as row (1×n); caller ensures shape
    (let [^ZMatrixRMaj vt (.M ^ComplexDense v)
          out (ZMatrixRMaj. (.numRows M) (.numCols vt))]
      (CommonOps_DDRM/mult M vt out)
      (ComplexDense. out)))

  (mul--s [_ s]
    (let [A (.copy M)]
      (CommonOps_DDRM/scale (double s) A)
      (ComplexDense. A)))

  (cholesky [_]
    (let [n (.numRows ^ZMatrixRMaj M)
          _ (when (not= n (.numCols ^ZMatrixRMaj M))
              (throw (ex-info "Cholesky requires square matrix" {:shape [n (.numCols M)]})))
          ^CholeskyDecomposition_F64 chol (DecompositionFactory_DDRM/chol n true)]
      (when-not (.decompose chol M)
        (throw (ex-info "Cholesky failed (matrix not SPD)"
                        {:shape [n n]})))
      (let [L (.getT chol (ZMatrixRMaj. n n))]   ; lower if 'true' above
        {:L (ComplexDense. L) :lower? true :spd? true})))

  (determinant [_]
    (CommonOps_DDRM/det M))

  (inverse [_]
    (let [out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/invert M out)
      (ComplexDense. out)))

  (norm [_ t]
    (case t
      :fro (double (NormOps_DDRM/normF M))
      :two (double (NormOps_DDRM/normP2 M))
      :one (double (NormOps_DDRM/normP1 M))
      :inf (double (NormOps_DDRM/normPInf M))
      (throw (ex-info "Unknown norm kind" {:t t}))))

  (solve [this v]
    (let [^ZMatrixRMaj B (.M ^ComplexDense v)
          X (ZMatrixRMaj. (.numRows B) (.numCols B))]
      (CommonOps_DDRM/solve M B X)
      (ComplexDense. X)))

  (trace [_]
    (CommonOps_DDRM/trace M))

  (transpose [_]
    (let [T (ZMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_DDRM/transpose M T)
      (ComplexDense. T)))

;; ;; -------- Predicates --------
  (normal? [_]
    ;; AᵀA ≈ AAᵀ
    (let [AtA (ZMatrixRMaj. (.numCols M) (.numCols M))
          AAt (ZMatrixRMaj. (.numRows M) (.numRows M))]
      (CommonOps_DDRM/multTransA M M AtA)  ;; Aᵀ A
      (CommonOps_DDRM/multTransB M M AAt)  ;; A Aᵀ
      (MatrixFeatures_DDRM/isIdentical AtA AAt tolerance--default)))

  (singular?
    ;; "Approximate, but robust numerical method."
    [_]
    (let [cond (NormOps_DDRM/conditionP2 M)]
      (> cond 1.0e12)))

  (square? [_]
    (= (.numRows M) (.numCols M)))

  (symmetric? [_] (MatrixFeatures_DDRM/isSymmetric M))
  (symmetric? [_ tol] (MatrixFeatures_DDRM/isSymmetric M (double tol)))

  (unitary? [_]
    ;; For real matrices, “unitary” == orthogonal.
    (MatrixFeatures_DDRM/isOrthogonal M tolerance--default)))

;; -------------------------------------------------------------------
;; Constructors 
;; -------------------------------------------------------------------

;; (defn realdense
;;   "Construct:
;;    (realdense m n)                -> zeros m×n
;;    (realdense ^ZMatrixRMaj M)     -> wrap EJML matrix
;;    (realdense m n ^doubles data)  -> from row-major data"
;;   (^ComplexDense [^long m ^long n] (ComplexDense. (ZMatrixRMaj. m n)))
;;   (^ComplexDense [^ZMatrixRMaj M]  (ComplexDense. M))
;;   (^ComplexDense [^long m ^long n ^doubles data]
;;    (ComplexDense. (ZMatrixRMaj. m n false data))))

(comment (println "test"))
