
(ns fastmath.matrix.dense.complex.ejml
  "Implementing complex matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml].

   Legend:
   m - matrix type
   i - row index
   j - column index
  "

  (:require
   [clojure.string :as str]
   [fastmath.vector :as v]
   [fastmath.core :as fm]
   [fastmath.protocols.matrix2 :as mat]
   [fastmath.matrix.dense.real.ejml :as realdense])

  (:import
   (java.lang Math)
   (org.ejml.data
    ;; DMatrixRMaj
    ZMatrixRMaj
    Complex_F64)
   (org.ejml.interfaces.decomposition LUDecomposition_F64)
   (org.ejml.dense.row
    CommonOps_ZDRM
    MatrixFeatures_ZDRM
    NormOps_ZDRM)
   (org.ejml.dense.row.factory
    DecompositionFactory_ZDRM)

   (org.ejml.interfaces.decomposition CholeskyDecomposition_F64)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^:private tolerance--default 1.0e-10)
;; (def ^:private ^double tolerance--default 1.0e-10)

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
    (CommonOps_ZDRM/extract m i0 i1 j0 j1 out 0 0)
    out))

(defn ^ZMatrixRMaj extract-multiple
  ;; TODO: efficiency implementaton
  "Reuses the destination memory for multiple access calls."
  [])

(defn- getc
  "Utility for returning complex numbers from a matrix."
  [^ZMatrixRMaj m ^long i ^long j]
  (let [^Complex_F64 cnum   (Complex_F64. 0.0 0.0)]
    (.get m i j cnum)
    cnum))

(defn- setc
  ;; TODO: Should have a multiple element version 
  "Utility for assigning complex numbers to a matrix."
  [^ZMatrixRMaj m  i j ^Complex_F64 cnum]
  (let [A (.copy m)]
    (.set A i j (.-real cnum) (.-imaginary cnum)) A))

(defn ^String format--real [^double x ^long p ^double eps]
  (format (str "%." p "f") (if (< (Math/abs x) eps) 0.0 x)))
(defn ^String format--complex [^double re ^double im ^long p ^double eps]
  (let [rz (if (< (Math/abs re) eps) 0.0 re)
        iz (if (< (Math/abs im) eps) 0.0 im)]
    (str (format (str "%." p "f") rz)
         (format (str "%+" p "f") iz)
         "i")))

(defn- realZM?
  ([^ZMatrixRMaj M]
   (realZM? M tolerance--default))
  ([^ZMatrixRMaj M tolerance]
   (let [^"[D" d (.-data M)
         n (alength d)]
     (loop [i 1]
       (cond
         (>= i n) true
         (> (Math/abs (aget d i)) (double tolerance)) false
         :else (recur (+ i 2)))))))

(defn ->str
  "Pretty string for ZMatrixRMaj with aligned columns.
   Options: :precision (3), :eps (1e-12), :max-rows (12), :max-cols (12)"

  ([^ZMatrixRMaj A] (->str A {}))
  ([^ZMatrixRMaj A {:keys [precision eps max-rows max-cols]
                    :or   {precision 3 eps tolerance--default max-rows 12 max-cols 12}}]
   (let [num-rows (.numRows A) num-cols (.numCols A)
         rlim (min num-rows ^long max-rows) clim (min num-cols ^long max-cols)
         real? (realZM? A eps)
         cell (fn [i j]
                (if real?
                  (format--real (.getReal A i j) precision eps)
                  (format--complex (.getReal A i j) (.getImag A i j) precision eps)))
         mat (vec (for [i (range rlim)]
                    (vec (concat (for [j (range clim)] (cell i j))
                                 (when (< ^long clim num-cols) ["…"])))))
         w   (if (seq mat) (apply max (mapcat #(map count %) mat)) 0)
         pad (fn [s] (format (str "%-" w "s") s))
         row->s (fn [row] (str "[" (str/join " " (map pad row)) "]"))
         head (format ":shape [%d %d] :type %s" num-rows num-cols (if real? "float64" "complex128"))]
     (str head "\n"
          "["
          (str/join "\n" (map-indexed #(str (when (> ^long %1 0) " ") (row->s %2)) mat))
          (when (< ^long rlim num-rows) "\n …")
          "]"))))

(defn print!
  ([^ZMatrixRMaj A] (print! A {}))
  ([^ZMatrixRMaj A opts]
   (println (->str A opts))))

(deftype ComplexDense [^ZMatrixRMaj M]
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
          (aset-double out i j (getc M i j))))
      out))

  (->array--float [_]
    ;; TODO: faster
    (float-array (map float (seq (.data M)))))

  (->array--2dfloat [_]
    (let [r (.numRows M) c (.numCols M)
          out (make-array Float/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-float out i j (float (getc M i j)))))
      out))

  ;; -------- Retrieval --------
  (columns [_]
    (mapv (fn [^long ic]
            (ComplexDense. (extract M [0 (.numRows M)] [ic (+ 1 ic)])))
          (range (.numCols M))))

  (rows [_]
    (mapv (fn [^long ir]
            (ComplexDense. (extract M [ir (+ ir 1)] [0 (.numCols M)])))
          (range (.numRows M))))

  (diagonal [_]
    (let [n (min (.numRows M) (.numCols M))
          out (ZMatrixRMaj. n 1)]
      (CommonOps_ZDRM/extractDiag M out)
      (ComplexDense. out)))

  (element [_ i j]
    (getc M (long i) (long j)))

  (column [_ j]
    (ComplexDense. (extract M [0 (.numRows M)] [j (+ 1 (long j))])))

  (row [_ i]
    (ComplexDense. (extract M [i (+ 1 (long i))] [0 (.numCols M)])))

  (num-rows [_] (.numRows M))
  (num-cols [_] (.numCols M))
  (shape    [_] [(.numRows M) (.numCols M)])

  ;; -------- Operations --------
  (add [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_ZDRM/add M B out)
      (ComplexDense. out)))

  (add--s [_ s]
    (let [A (.copy M) d (.data A) s (double s)]
      (dotimes [k (alength d)]
        (aset-double d k (+ (aget d k) s)))
      (ComplexDense. A)))

  (sub
    [this]
    (let [A (.copy M)]
      (CommonOps_ZDRM/scale -1.0 0.0 A)
      (ComplexDense. A)))

  (sub [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_ZDRM/subtract M B out)
      (ComplexDense. out)))

  (map--m [_ f]
    ;; TODO: Prevent function from multiple calls to set
    ;; - Work on (.data A) directly
    (let [^ZMatrixRMaj A (.copy M)]
      (dotimes [i (.numRows A)]
        (dotimes [j (.numCols A)]
          (setc A i j (f (getc A i j)))))
      (ComplexDense. A)))

  (mul
    [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols B))]
      (CommonOps_ZDRM/mult M B out)
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
          (and tA? tB?) (CommonOps_ZDRM/multTransAB M B out)
          tA?           (CommonOps_ZDRM/multTransA  M B out)
          tB?           (CommonOps_ZDRM/multTransB  M B out)
          :else         (CommonOps_ZDRM/mult        M B out))
        (ComplexDense. out))))

  (mul--e [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)]
      (ComplexDense. (CommonOps_ZDRM/elementMultiply M B nil nil))))

  (mul--v [this v]
    ;; v is expected to be (n×1) ComplexDense
    (let [^ZMatrixRMaj x (.M ^ComplexDense v)
          out (ZMatrixRMaj. (.numRows M) 1)]
      (CommonOps_ZDRM/mult M x out)
      (ComplexDense. out)))

  (mul--vt [this v]
    ;; v treated as row (1×n); caller ensures shape
    (let [^ZMatrixRMaj vt (.M ^ComplexDense v)
          out (ZMatrixRMaj. (.numRows M) (.numCols vt))]
      (CommonOps_ZDRM/mult M vt out)
      (ComplexDense. out)))

  (mul--s [_ [r i]]
    (let [A (.copy M)]
      (CommonOps_ZDRM/scale (double r) (double i) A)
      (ComplexDense. A)))

  (cholesky [_]
    (let [n (.numRows ^ZMatrixRMaj M)
          _ (when (not= n (.numCols ^ZMatrixRMaj M))
              (throw (ex-info "Cholesky requires square matrix" {:shape [n (.numCols M)]})))
          ^CholeskyDecomposition_F64 chol (DecompositionFactory_ZDRM/chol n true)]
      (when-not (.decompose chol M)
        (throw (ex-info "Cholesky failed (matrix not SPD)"
                        {:shape [n n]})))
      (let [L (.getT chol (ZMatrixRMaj. n n))]   ; lower if 'true' above
        {:L (ComplexDense. L) :lower? true :spd? true})))

  (determinant [_]
    (CommonOps_ZDRM/det M))

  (inverse [_]
    (let [out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_ZDRM/invert M out)
      (ComplexDense. out)))

  (norm [_ t]
    (case t
      :fro (double (NormOps_ZDRM/normF M))
      (throw (ex-info "Unknown norm kind" {:t t}))))

  (solve [this v]
    (let [^ZMatrixRMaj B (.M ^ComplexDense v)
          X (ZMatrixRMaj. (.numRows B) (.numCols B))]
      (CommonOps_ZDRM/solve M B X)
      (ComplexDense. X)))

  (trace [_]
    (CommonOps_ZDRM/trace M nil))

  (transpose [_]
    (let [T (ZMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_ZDRM/transpose M T)
      (ComplexDense. T)))

;; ;; -------- Predicates --------
  (normal? [_]
    ;; AᵀA ≈ AAᵀ
    (let [AtA (ZMatrixRMaj. (.numCols M) (.numCols M))
          AAt (ZMatrixRMaj. (.numRows M) (.numRows M))]
      (CommonOps_ZDRM/multTransA M M AtA)  ;; Aᵀ A
      (CommonOps_ZDRM/multTransB M M AAt)  ;; A Aᵀ
      (MatrixFeatures_ZDRM/isIdentical AtA AAt tolerance--default)))

  (singular?
    ;; "Approximate, but robust numerical method."
    [_]
    (let [^LUDecomposition_F64 lu (DecompositionFactory_ZDRM/lu
                                   (.numRows M) (.numCols M))]
      (or (not (.decompose lu M))
          (.isSingular lu))))

  (square? [_]
    (= (.numRows M) (.numCols M)))

  (symmetric? [_]
    (if (= (.numRows M) (.numCols M))
      (MatrixFeatures_ZDRM/isEquals M
                                    (CommonOps_ZDRM/transpose M nil)
                                    tolerance--default)
      false))
  (symmetric? [_ tol]
    (if (= (.numRows M) (.numCols M))
      (MatrixFeatures_ZDRM/isEquals M
                                    (CommonOps_ZDRM/transpose M nil)
                                    tol)
      false))

  (unitary? [_]
    ;; For real matrices, “unitary” == orthogonal.
    (MatrixFeatures_ZDRM/isUnitary M tolerance--default))

  mat/MatrixComplex

  (adjoint [_]
    (let [^ZMatrixRMaj out (ZMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_ZDRM/transposeConjugate M out)
      (ComplexDense. out)))
  (conjugate [_] (ComplexDense. (CommonOps_ZDRM/conjugate M nil)))
  (imag [_] (realdense/->RealDense (CommonOps_ZDRM/imaginary M nil)))
  (real [_] (realdense/->RealDense (CommonOps_ZDRM/real M nil)))

  (hermitian? [_]
    (MatrixFeatures_ZDRM/isHermitian M tolerance--default))

  (real? [_]
    (realZM? M))

  Object
  (toString [_]
    (print! M)))

;; -------------------------------------------------------------------
;; Constructors 
;; -------------------------------------------------------------------

(defn complexdense
  ;; I
  ;; (^ComplexDense [^long n]
  ;;  (ComplexDense. (CommonOps_ZDRM/identity n)))

  ;; zero
  (^ComplexDense [^long n ^long o]
   (ComplexDense. (ZMatrixRMaj. n o)))

  ;; diagonal
  (^ComplexDense [^doubles data]
   (-> (CommonOps_ZDRM/diag (double-array data))
       ComplexDense.))

;; elements
  (^ComplexDense [^long n ^long o ^doubles data]
   (ComplexDense. (ZMatrixRMaj. n o false data))))

(defn ^ComplexDense complexdense<-rows [rows]
  (let [nrows (count rows)
        ncols (count (first rows))
        data  (double-array (apply concat rows))]
    (ComplexDense. (ZMatrixRMaj. nrows ncols true data))))

(defn ^ComplexDense complexdense<-cols [cols]
  (let [^ZMatrixRMaj A (.M (complexdense<-rows cols))
        ^ZMatrixRMaj out (ZMatrixRMaj. (.numCols A) (.numRows A))]
    (CommonOps_ZDRM/transpose A out)

    (ComplexDense. out)))

(comment (println "test")
         (def A--test (complexdense 2 2 (double-array [1 1 0 0 1 0 0 0])))
         (def B--test (complexdense 2 2 (double-array [1 0 1 0 1 0 0 1])))
         (-> (mat/add A--test B--test)
             println)

         (-> (ComplexDense. (ZMatrixRMaj. 3 3))
             println))

