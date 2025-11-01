(ns fastmath.matrix.dense.real.ejml
  "Implementing real matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml]."

  (:require [fastmath.vector :as v]
            [fastmath.core :as fm]
            [fastmath.protocols.matrix2 :as proto])
  (:import
   (org.ejml.data DMatrixRMaj)
   (org.ejml.dense.row
    CommonOps_DDRM
    MatrixFeatures_DDRM
    NormOps_DDRM
    ;; DecompositionFactory_DDRM
    )
   (org.ejml.interfaces.decomposition CholeskyDecomposition_F64)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; ------------------------------------------------------------
;; Storage type
;; ------------------------------------------------------------

(deftype ComplexDense [^DMatrixRMaj M])

;; (defn realdense
;;   "Constructors:
;;    (realdense m n)                         ;; zeros m×n
;;    (realdense ^DMatrixRMaj M)              ;; wrap an EJML matrix
;;    (realdense m n ^doubles data)           ;; from row-major data"
;;   (^RealDense [^long m ^long n] (RealDense. (DMatrixRMaj. m n)))
;;   (^RealDense [^DMatrixRMaj M]  (RealDense. M))
;;   (^RealDense [^long m ^long n ^doubles data]
;;    (RealDense. (DMatrixRMaj. m n false data))))

(defn wrap
  "Wrap an existing EJML matrix."
  ^RealDense [^DMatrixRMaj M]
  (RealDense. M))

;; ------------------------------------------------------------
;; Utility 
;; ------------------------------------------------------------

(defn- ensure-R ^DMatrixRMaj [x]
  (cond
    (instance? DMatrixRMaj x) x
    (instance? RealDense x)      (.M ^RealDense x)
    :else (throw (ex-info "Expected RealDense or DMatrixRMaj" {:got (type x)}))))

(defn- copy ^DMatrixRMaj [^DMatrixRMaj A]
  (.copy A))

(defn- negate! ^DMatrixRMaj [^DMatrixRMaj A]
  (CommonOps_DDRM/scale -1.0 A) A)

(def ^:private tolerance--default 1.0e-10)

;; ------------------------------------------------------------
;; Implementation 
;; ------------------------------------------------------------

(extend-type ComplexDense
  proto/Matrix

  ;; -------- Transformations --------

  (->seq [m]
    (let [A (ensure-R m)
          r (.numRows A) c (.numCols A)]
      (for [i (range r), j (range c)] (.get A i j))))

  (->array [m]
    (let [A (ensure-R m)
          d (.data A)
          out (double-array (alength d))]
      (System/arraycopy d 0 out 0 (alength d))
      out))

  (->array--2d [m]
    (let [A (ensure-R m)
          r (.numRows A) c (.numCols A)
          out (make-array Double/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-double out i j (.get A i j))))
      out))

  (->array--float [m]
    (float-array (map float (seq (.data ^DMatrixRMaj (ensure-R m))))))

  (->array--2dfloat [m]
    (let [A (ensure-R m)
          r (.numRows A) c (.numCols A)
          out (make-array Float/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-float out i j (float (.get A i j)))))
      out))

  (->matrix--real [m] (ensure-R m))

  (->matrix--complex [m]
    (throw (ex-info "Real→complex promotion is not provided in the real-only backend" {})))

  ;; -------- Retrieval --------

  (cols [m] (.numCols (ensure-R m)))
  (rows [m] (.numRows (ensure-R m)))

  (diagonal [m]
    (let [A (ensure-R m)
          n (min (.numRows A) (.numCols A))
          out (DMatrixRMaj. n 1)]
      (dotimes [i n] (.set out i 0 (.get A i i)))
      (RealDense. out)))

  (element [m i j]
    (.get ^DMatrixRMaj (ensure-R m) (long i) (long j)))

  (column [m j]
    (let [A (ensure-R m)
          r (.numRows A)
          out (DMatrixRMaj. r 1)]
      (CommonOps_DDRM/extract A 0 r j (inc j) out)
      (RealDense. out)))

  (row [m i]
    (let [A (ensure-R m)
          c (.numCols A)
          out (DMatrixRMaj. 1 c)]
      (CommonOps_DDRM/extract A i (inc i) 0 c out)
      (RealDense. out)))

  (num-rows [m] (.numRows (ensure-R m)))
  (num-cols [m] (.numCols (ensure-R m)))
  (shape    [m] (let [A (ensure-R m)] [(.numRows A) (.numCols A)]))

  ;; -------- Operations --------

  (add [m1 m2]
    (let [A (ensure-R m1) B (ensure-R m2)
          out (DMatrixRMaj. (.numRows A) (.numCols A))]
      (CommonOps_DDRM/add 1.0 A 1.0 B out)
      (RealDense. out)))

  (add--s [m s]
    (let [A (copy (ensure-R m))
          d (.data A)
          s (double s)]
      (dotimes [k (alength d)]
        (aset-double d k (+ (aget d k) s)))
      (RealDense. A)))

  (sub
    ([m]
     (let [A (copy (ensure-R m))]
       (negate! A)
       (RealDense. A)))
    ([m1 m2]
     (let [A (ensure-R m1) B (ensure-R m2)
           out (DMatrixRMaj. (.numRows A) (.numCols A))]
       (CommonOps_DDRM/add 1.0 A -1.0 B out)
       (RealDense. out))))

  (map--m [f m]
    (let [A (copy (ensure-R m))
          r (.numRows A) c (.numCols A)]
      (dotimes [i r]
        (dotimes [j c]
          (.set A i j (double (f (.get A i j))))))
      (RealDense. A)))

  (mul
    ([m1 m2]
     (let [A (ensure-R m1) B (ensure-R m2)
           out (DMatrixRMaj. (.numRows A) (.numCols B))]
       (CommonOps_DDRM/mult A B out)
       (RealDense. out)))
    ([m1 tA? m2 tB?]
     (let [A (ensure-R m1) B (ensure-R m2)
           rA (if tA? (.numCols A) (.numRows A))
           cA (if tA? (.numRows A) (.numCols A))
           rB (if tB? (.numCols B) (.numRows B))
           cB (if tB? (.numRows B) (.numCols B))]
       (when (not= cA rB)
         (throw (ex-info "Incompatible shapes for matmul"
                         {:A [rA cA] :B [rB cB] :tA? tA? :tB? tB?})))
       (let [out (DMatrixRMaj. rA cB)]
         (cond
           (and tA? tB?) (CommonOps_DDRM/multTransAB A B out)
           tA?           (CommonOps_DDRM/multTransA  A B out)
           tB?           (CommonOps_DDRM/multTransB  A B out)
           :else         (CommonOps_DDRM/mult        A B out))
         (RealDense. out)))))

  (mul--e [m1 m2]
    (let [A (copy (ensure-R m1))
          B (ensure-R m2)]
      (CommonOps_DDRM/elementMult A B)
      (RealDense. A)))

  (mul--v [m v]
    ;; v is expected to be (n×1)
    (let [A (ensure-R m) x (ensure-R v)
          out (DMatrixRMaj. (.numRows A) 1)]
      (CommonOps_DDRM/mult A x out)
      (RealDense. out)))

  (mul--vt [m v]
    ;; v treated as row (1×n); caller ensures shape
    (let [A (ensure-R m) vt (ensure-R v)
          out (DMatrixRMaj. (.numRows A) (.numCols vt))]
      (CommonOps_DDRM/mult A vt out)
      (RealDense. out)))

  (mul--s [m s]
    (let [A (copy (ensure-R m))]
      (CommonOps_DDRM/scale (double s) A)
      (RealDense. A)))

  (cholesky [m]
    ;; Return a map to make intent explicit.
    (let [A (ensure-R m)
          ^CholeskyDecomposition_F64 decomp (DecompositionFactory_DDRM/cholesky true) ; lower
          ok (.decompose decomp A)]
      (when-not ok
        (throw (ex-info "Cholesky failed (matrix not SPD)" {:shape [(.numRows A) (.numCols A)]})))
      (let [L (.getT decomp (DMatrixRMaj. (.numRows A) (.numCols A)))]
        {:L (RealDense. L) :lower? true :spd? true})))

  (determinant [m]
    (CommonOps_DDRM/det (ensure-R m)))

  (inverse [m]
    (let [A (ensure-R m)
          out (DMatrixRMaj. (.numRows A) (.numCols A))]
      (CommonOps_DDRM/invert A out)
      (RealDense. out)))

  (norm [m t]
    (let [A (ensure-R m)]
      (case t
        :fro (double (NormOps_DDRM/normF A))
        :two (double (NormOps_DDRM/normP2 A))
        :one (double (NormOps_DDRM/normP1 A))
        :inf (double (NormOps_DDRM/normPInf A))
        (throw (ex-info "Unknown norm kind" {:t t})))))

  (solve [m v]
    (let [A (ensure-R m)
          B (ensure-R v)
          X (DMatrixRMaj. (.numRows B) (.numCols B))]
      (CommonOps_DDRM/solve A B X)
      (RealDense. X)))

  (trace [m]
    (CommonOps_DDRM/trace (ensure-R m)))

  (transpose [m]
    (let [A (ensure-R m)
          T (DMatrixRMaj. (.numCols A) (.numRows A))]
      (CommonOps_DDRM/transpose A T)
      (RealDense. T)))

  ;; -------- Predicates --------

  (normal? [m]
    ;; AᵀA ≈ AAᵀ
    (let [A (ensure-R m)
          n (.numRows A) p (.numCols A)
          AtA (DMatrixRMaj. p p)
          AAt (DMatrixRMaj. n n)]
      (CommonOps_DDRM/multTransA A A AtA) ; Aᵀ A
      (CommonOps_DDRM/mult A (DMatrixRMaj. n n) AAt) ; fill AAt
      (MatrixFeatures_DDRM/isIdentical
       AtA
       (doto AAt identity)
       tolerance--default)))

  (singular? [m]
    (MatrixFeatures_DDRM/isSingular (ensure-R m)))

  (square? [m]
    (let [A (ensure-R m)] (= (.numRows A) (.numCols A))))

  (symmetric?
    ([m]
     (MatrixFeatures_DDRM/isSymmetric (ensure-R m)))
    ([m tol]
     (MatrixFeatures_DDRM/isSymmetric (ensure-R m) (double tol))))

  (unitary? [m]
    ;; For real matrices, unitary == orthogonal.
    (MatrixFeatures_DDRM/isOrthogonal (ensure-R m) tolerance--default)))
