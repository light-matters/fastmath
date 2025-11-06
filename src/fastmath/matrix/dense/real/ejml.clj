(ns fastmath.matrix.dense.real.ejml
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
            [fastmath.matrix.dense.constant :as constant]
            [fastmath.protocols.matrix2 :as mat])
  (:import
   (org.ejml.data DMatrixRMaj)
   (org.ejml.dense.row
    CommonOps_DDRM
    MatrixFeatures_DDRM
    NormOps_DDRM)
   (org.ejml.dense.row.factory
    DecompositionFactory_DDRM)

   (org.ejml.interfaces.decomposition CholeskyDecomposition_F64)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; ==================================================
;; Functions 
;; ==================================================
(defn- ^DMatrixRMaj extract
;; TODO: Add to type?
  "Extracts a submatrix m[i0:i1, j0:j1) into a new DMatrixRMaj.
   Indices are zero-based and end-exclusive."
  [^DMatrixRMaj m
   ^longs [i0 i1]
   ^longs [j0 j1]]
  (let [^DMatrixRMaj out (DMatrixRMaj. (- ^long i1 ^long i0) (- ^long j1 ^long j0))]
    (CommonOps_DDRM/extract m i0 i1 j0 j1 out 0 0)
    out))

(defn ^DMatrixRMaj extract-multiple
  ;; TODO: efficiency implementaton
  "Reuses the destination memory for multiple access calls."
  [])

(deftype RealDense [^DMatrixRMaj M]
  mat/MatrixReal

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
          out (DMatrixRMaj. n 1)]
      (CommonOps_DDRM/extractDiag M out)
      (RealDense. out)))

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
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/add 1.0 M 1.0 B out)
      (RealDense. out)))

  (add--s [_ s]
    (let [A (.copy M) d (.data A) s (double s)]
      (dotimes [k (alength d)]
        (aset-double d k (+ (aget d k) s)))
      (RealDense. A)))

  (sub
    [this]
    (let [A (.copy M)]
      (CommonOps_DDRM/scale -1.0 A)
      (RealDense. A)))

  (sub [this other]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/add 1.0 M -1.0 B out)
      (RealDense. out)))

  (map--m [_ f]
    (let [A (.copy M)]
      (dotimes [i (.numRows A)]
        (dotimes [j (.numCols A)]
          (.set A i j (double (f (.get A i j))))))
      (RealDense. A)))

  (mul
    [this other]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols B))]
      (CommonOps_DDRM/mult M B out)
      (RealDense. out)))

  (mul [this tA? other tB?]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          rA (if tA? (.numCols M) (.numRows M))
          cA (if tA? (.numRows M) (.numCols M))
          rB (if tB? (.numCols B) (.numRows B))
          cB (if tB? (.numRows B) (.numCols B))]
      (when (not= cA rB)
        (throw (ex-info "Incompatible shapes for matmul"
                        {:A [rA cA] :B [rB cB] :tA? tA? :tB? tB?})))
      (let [out (DMatrixRMaj. rA cB)]
        (cond
          (and tA? tB?) (CommonOps_DDRM/multTransAB M B out)
          tA?           (CommonOps_DDRM/multTransA  M B out)
          tB?           (CommonOps_DDRM/multTransB  M B out)
          :else         (CommonOps_DDRM/mult        M B out))
        (RealDense. out))))

  (mul--e [this other]
    (let [A (.copy M)
          ^DMatrixRMaj B (.M ^RealDense other)]
      (CommonOps_DDRM/elementMult A B)
      (RealDense. A)))

  (mul--v [this v]
    ;; v is expected to be (n×1) RealDense
    (let [^DMatrixRMaj x (.M ^RealDense v)
          out (DMatrixRMaj. (.numRows M) 1)]
      (CommonOps_DDRM/mult M x out)
      (RealDense. out)))

  (mul--vt [this v]
    ;; v treated as row (1×n); caller ensures shape
    (let [^DMatrixRMaj vt (.M ^RealDense v)
          out (DMatrixRMaj. (.numRows M) (.numCols vt))]
      (CommonOps_DDRM/mult M vt out)
      (RealDense. out)))

  (mul--s [_ s]
    (let [A (.copy M)]
      (CommonOps_DDRM/scale (double s) A)
      (RealDense. A)))

  (cholesky [_]
    (let [n (.numRows ^DMatrixRMaj M)
          _ (when (not= n (.numCols ^DMatrixRMaj M))
              (throw (ex-info "Cholesky requires square matrix" {:shape [n (.numCols M)]})))
          ^CholeskyDecomposition_F64 chol (DecompositionFactory_DDRM/chol n true)]
      (when-not (.decompose chol M)
        (throw (ex-info "Cholesky failed (matrix not SPD)"
                        {:shape [n n]})))
      (let [L (.getT chol (DMatrixRMaj. n n))]   ; lower if 'true' above
        {:L (RealDense. L) :lower? true :spd? true})))

  (determinant [_]
    (CommonOps_DDRM/det M))

  (inverse [_]
    (let [out (DMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/invert M out)
      (RealDense. out)))

  (norm [_ t]
    (case t
      :fro (double (NormOps_DDRM/normF M))
      :two (double (NormOps_DDRM/normP2 M))
      :one (double (NormOps_DDRM/normP1 M))
      :inf (double (NormOps_DDRM/normPInf M))
      (throw (ex-info "Unknown norm kind" {:t t}))))

  (solve [this v]
    (let [^DMatrixRMaj B (.M ^RealDense v)
          X (DMatrixRMaj. (.numRows B) (.numCols B))]
      (CommonOps_DDRM/solve M B X)
      (RealDense. X)))

  (trace [_]
    (CommonOps_DDRM/trace M))

  (transpose [_]
    (let [T (DMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_DDRM/transpose M T)
      (RealDense. T)))

;; ;; -------- Predicates --------
  (normal? [_]
    ;; AᵀA ≈ AAᵀ
    (let [AtA (DMatrixRMaj. (.numCols M) (.numCols M))
          AAt (DMatrixRMaj. (.numRows M) (.numRows M))]
      (CommonOps_DDRM/multTransA M M AtA)  ;; Aᵀ A
      (CommonOps_DDRM/multTransB M M AAt)  ;; A Aᵀ
      (MatrixFeatures_DDRM/isIdentical AtA AAt constant/tolerance--default)))

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
    (MatrixFeatures_DDRM/isOrthogonal M constant/tolerance--default))

  Object
  (toString [_]
    (.toString M)))

;; -------------------------------------------------------------------
;; Constructors 
;; -------------------------------------------------------------------

(defn realdense
  ;; I
  ;; (^RealDense [^long n]
  ;;  (RealDense. (CommonOps_DDRM/identity n)))

  ;; zero
  (^RealDense [^long n ^long o]
   (RealDense. (DMatrixRMaj. n o)))

  ;; diagonal
  (^RealDense [^doubles data]
   (-> (CommonOps_DDRM/diag (double-array data))
       RealDense.))

;; elements
  (^RealDense [^long n ^long o ^doubles data]
   (RealDense. (DMatrixRMaj. n o false data))))

(defn realdense<-rows [rows]
  (let [nrows (count rows)
        ncols (count (first rows))
        data  (double-array (apply concat rows))]
    (RealDense. (DMatrixRMaj. nrows ncols true data))))

(defn realdense<-cols [cols]
  (let [M (.M (realdense<-rows cols))
        out (DMatrixRMaj. (.numCols M) (.numRows M))]
    (CommonOps_DDRM/transpose M out)
    (RealDense. out)))

;; (defn realdense<-cols [cols]
;;   (let [nrows (count cols)
;;         ncols (count (first cols))
;;         data  (double-array (apply concat cols))
;;         M (DMatrixRMaj. nrows ncols true data)
;;         out (DMatrixRMaj. (.numCols M) (.numRows M))]
;;     (CommonOps_DDRM/transpose M out)
;;     (RealDense. out)))

(comment (println "test")
         (-> (->RealDense
              (DMatrixRMaj. 3 3 false (double-array [1 0 0 0 1 0 0 0 1])))
             println)

         (-> (RealDense. (DMatrixRMaj. 3 3))
             println))
