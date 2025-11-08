(ns fastmath.matrix.dense.complex.ojalgo
  "Implementing complex matrices using ojAlgo as a backend [https://github.com/optimatika/ojAlgo/tree/bd85bf0fe485016817fad945a5975eba8164b369?tab=readme-ov-file].

[https://javadoc.io/doc/org.ojalgo/ojalgo/latest/ojalgo/module-summary.html]
   NOTE:
   Ojalgo uses column-major ordering.
   TODO:
   - Understand the API! [...]
   - Get something working [x]
   - 'efficient' conversion from Java arrays [x] 
   - bench against ejml [ ]

   Legend:
   m - matrix type
   i - row index
   j - column index
  "
  (:import
   (org.ojalgo.matrix.store GenericStore PhysicalStore)
   (org.ojalgo.scalar ComplexNumber)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                               Helper Functions                              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Type hint the return value
(defn complex ^ComplexNumber
  ([real] (ComplexNumber/of (double real) 0.0))
  ([real imag] (ComplexNumber/of (double real) (double imag))))

;; Print function - use double cast instead of .doubleValue
(defn print-complex-matrix
  "TODO: Not really necessary if we use the MatrixC128 class."
  [^PhysicalStore matrix]
  (println)
  (dotimes [i (.countRows matrix)]
    (print "|")
    (dotimes [j (.countColumns matrix)]
      (let [^ComplexNumber c (.get matrix (long i) (long j))
            real (double (.getReal c))
            imag (double (.getImaginary c))]
        (print (format "%7.2f%+.2fi" real imag))))
    (print "  |")
    (println)))

;; Convert matrix to Clojure data
(defn matrix->vectors [^PhysicalStore matrix]
  (let [rows (.countRows matrix)
        cols (.countColumns matrix)]
    (vec (for [i (range rows)]
           (vec (for [j (range cols)]
                  (let [^ComplexNumber c (.get matrix (long i) (long j))]
                    [(double (.getReal c))
                     (double (.getImaginary c))])))))))

;; Test

(defn ->carray1d
  "A 1D Java array of `ComplexNumber`s from a 3D Clojure vector system of double-valued pairs, e.g.

[[[1 0] [2 1]]
 [[3 -1] [4 0]]].
  "
  [columns--complex]
  (->> columns--complex
       (into [] (comp cat))
       (mapv #(apply complex %))
       (into-array ComplexNumber)))

(defn ^GenericStore ->GenericStore
  "Use reflection to call GenericStore.wrap, bypassing access restrictions.

   TODO: Remove the accessibility methods when the bug is fixed upstream. 
  "
  [array1d num-cols]
  (let [arr (into-array ComplexNumber array1d)
        factory (GenericStore/C128)

        ;; Find the wrap method with 3 parameters
        wrap-method (->> (.getMethods GenericStore)
                         (filter #(= "wrap" (.getName ^java.lang.reflect.Method %)))
                         (filter #(= 3 (alength (.getParameterTypes ^java.lang.reflect.Method %))))
                         first)]
    (when wrap-method
      (.setAccessible ^java.lang.reflect.Method wrap-method true)
      ;; invoke takes: (method, object-instance, arg-array)
      ;; For static methods, object-instance is nil
      (.invoke ^java.lang.reflect.Method wrap-method nil (into-array Object [factory arr (int num-cols)])))))

(defn ->ComplexMatrix
  "A 2D matrix of `ComplexNumber`s from a 3D Clojure vector system of double-valued pairs, e.g.

[[[1 0] [2 1]]
 [[3 -1] [4 0]]].
  "
  [columns]

  (-> (->carray1d columns)
      (->GenericStore  (count columns))
      (->> (.makeWrapper MatrixC128/FACTORY))))

(comment
  (def vcs [[[1 0] [2 1]]
            [[3 -1] [4 0]]])
  (def result (->ComplexMatrix vcs))

  (println result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                  Protocols                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; (deftype ComplexDense [^ZMatrixRMaj M]
;;   mat/MatrixReal
;;   ;; -------- Transformations --------
;;   (->seq [_]
;;     (seq (.data M)))

;;   (->array [_]
;;     (let [d (.data M) out (double-array (alength d))]
;;       (System/arraycopy d 0 out 0 (alength d))
;;       out))

;;   (->array--2d [_]
;;     (let [r (.numRows M) c (.numCols M)
;;           out (make-array Double/TYPE r c)]
;;       (dotimes [i r]
;;         (dotimes [j c]
;;           (aset-double out i j (getc M i j))))
;;       out))

;;   (->array--float [_]
;;     ;; TODO: faster
;;     (float-array (map float (seq (.data M)))))

;;   (->array--2dfloat [_]
;;     (let [r (.numRows M) c (.numCols M)
;;           out (make-array Float/TYPE r c)]
;;       (dotimes [i r]
;;         (dotimes [j c]
;;           (aset-float out i j (float (getc M i j)))))
;;       out))

;;   ;; -------- Retrieval --------
;;   (columns [_]
;;     (mapv (fn [^long ic]
;;             (->ComplexDense (extract M [0 (.numRows M)] [ic (+ 1 ic)])))
;;           (range (.numCols M))))

;;   (rows [_]
;;     (mapv (fn [^long ir]
;;             (->ComplexDense (extract M [ir (+ ir 1)] [0 (.numCols M)])))
;;           (range (.numRows M))))

;;   (diagonal [_]
;;     (let [n (min (.numRows M) (.numCols M))
;;           out (ZMatrixRMaj. n 1)]
;;       (CommonOps_ZDRM/extractDiag M out)
;;       (->ComplexDense out)))

;;   (element [_ i j]
;;     (getc M (long i) (long j)))

;;   (column [_ j]
;;     (->ComplexDense (extract M [0 (.numRows M)] [j (+ 1 (long j))])))

;;   (row [_ i]
;;     (->ComplexDense (extract M [i (+ 1 (long i))] [0 (.numCols M)])))

;;   (num-rows [_] (.numRows M))
;;   (num-cols [_] (.numCols M))
;;   (shape    [_] [(.numRows M) (.numCols M)])

;;   ;; -------- Operations --------
;;   (add [this other]
;;     (let [^ZMatrixRMaj B (.M ^ComplexDense other)
;;           out (ZMatrixRMaj. (.numRows M) (.numCols M))]
;;       (CommonOps_ZDRM/add M B out)
;;       (->ComplexDense out)))

;;   (add--s [_ s]
;;     (let [A (.copy M) d (.data A) s (double s)]
;;       (dotimes [k (alength d)]
;;         (aset-double d k (+ (aget d k) s)))
;;       (->ComplexDense A)))

;;   (sub
;;     [this]
;;     (let [A (.copy M)]
;;       (CommonOps_ZDRM/scale -1.0 0.0 A)
;;       (->ComplexDense A)))

;;   (sub [this other]
;;     (let [^ZMatrixRMaj B (.M ^ComplexDense other)
;;           out (ZMatrixRMaj. (.numRows M) (.numCols M))]
;;       (CommonOps_ZDRM/subtract M B out)
;;       (->ComplexDense out)))

;;   (map--m [_ f]
;;     ;; TODO: Prevent function from multiple calls to set
;;     ;; - Work on (.data A) directly
;;     (let [^ZMatrixRMaj A (.copy M)]
;;       (dotimes [i (.numRows A)]
;;         (dotimes [j (.numCols A)]
;;           (setc A i j (f (getc A i j)))))
;;       (->ComplexDense A)))

;;   (mul
;;     [this other]
;;     (let [^ZMatrixRMaj B (.M ^ComplexDense other)
;;           out (ZMatrixRMaj. (.numRows M) (.numCols B))]
;;       (CommonOps_ZDRM/mult M B out)
;;       (->ComplexDense out)))

;;   (mul [this tA? other tB?]
;;     (let [^ZMatrixRMaj B (.M ^ComplexDense other)
;;           rA (if tA? (.numCols M) (.numRows M))
;;           cA (if tA? (.numRows M) (.numCols M))
;;           rB (if tB? (.numCols B) (.numRows B))
;;           cB (if tB? (.numRows B) (.numCols B))]
;;       (when (not= cA rB)
;;         (throw (ex-info "Incompatible shapes for matmul"
;;                         {:A [rA cA] :B [rB cB] :tA? tA? :tB? tB?})))
;;       (let [out (ZMatrixRMaj. rA cB)]
;;         (cond
;;           (and tA? tB?) (CommonOps_ZDRM/multTransAB M B out)
;;           tA?           (CommonOps_ZDRM/multTransA  M B out)
;;           tB?           (CommonOps_ZDRM/multTransB  M B out)
;;           :else         (CommonOps_ZDRM/mult        M B out))
;;         (->ComplexDense out))))

;;   (mul--e [this other]
;;     (let [^ZMatrixRMaj B (.M ^ComplexDense other)]
;;       (->ComplexDense (CommonOps_ZDRM/elementMultiply M B nil nil))))

;;   (mul--v [this v]
;;     ;; v is expected to be (n×1) ComplexDense
;;     (let [^ZMatrixRMaj x (.M ^ComplexDense v)
;;           out (ZMatrixRMaj. (.numRows M) 1)]
;;       (CommonOps_ZDRM/mult M x out)
;;       (->ComplexDense out)))

;;   (mul--vt [this v]
;;     ;; v treated as row (1×n); caller ensures shape
;;     (let [^ZMatrixRMaj vt (.M ^ComplexDense v)
;;           out (ZMatrixRMaj. (.numRows M) (.numCols vt))]
;;       (CommonOps_ZDRM/mult M vt out)
;;       (->ComplexDense out)))

;;   (mul--s [_ [r i]]
;;     (let [A (.copy M)]
;;       (CommonOps_ZDRM/scale (double r) (double i) A)
;;       (->ComplexDense A)))

;;   (cholesky [_]
;;     (let [n (.numRows ^ZMatrixRMaj M)
;;           _ (when (not= n (.numCols ^ZMatrixRMaj M))
;;               (throw (ex-info "Cholesky requires square matrix" {:shape [n (.numCols M)]})))
;;           ^CholeskyDecomposition_F64 chol (DecompositionFactory_ZDRM/chol n true)]
;;       (when-not (.decompose chol M)
;;         (throw (ex-info "Cholesky failed (matrix not SPD)"
;;                         {:shape [n n]})))
;;       (let [L (.getT chol (ZMatrixRMaj. n n))]   ; lower if 'true' above
;;         {:L (->ComplexDense L) :lower? true :spd? true})))

;;   (determinant [_]
;;     (CommonOps_ZDRM/det M))

;;   (inverse [_]
;;     (let [out (ZMatrixRMaj. (.numRows M) (.numCols M))]
;;       (CommonOps_ZDRM/invert M out)
;;       (->ComplexDense out)))

;;   (norm [_ t]
;;     (case t
;;       :fro (double (NormOps_ZDRM/normF M))
;;       (throw (ex-info "Unknown norm kind" {:t t}))))

;;   (solve [this v]
;;     (let [^ZMatrixRMaj B (.M ^ComplexDense v)
;;           X (ZMatrixRMaj. (.numRows B) (.numCols B))]
;;       (CommonOps_ZDRM/solve M B X)
;;       (->ComplexDense X)))

;;   (trace [_]
;;     (CommonOps_ZDRM/trace M nil))

;;   (transpose [_]
;;     (let [T (ZMatrixRMaj. (.numCols M) (.numRows M))]
;;       (CommonOps_ZDRM/transpose M T)
;;       (->ComplexDense T)))

;; ;; ;; -------- Predicates --------
;;   (normal? [_]
;;     ;; AᵀA ≈ AAᵀ
;;     (let [AtA (ZMatrixRMaj. (.numCols M) (.numCols M))
;;           AAt (ZMatrixRMaj. (.numRows M) (.numRows M))]
;;       (CommonOps_ZDRM/multTransA M M AtA)  ;; Aᵀ A
;;       (CommonOps_ZDRM/multTransB M M AAt)  ;; A Aᵀ
;;       (MatrixFeatures_ZDRM/isIdentical AtA AAt constant/tolerance--default)))

;;   (singular?
;;     ;; "Approximate, but robust numerical method."
;;     [_]
;;     (let [^LUDecomposition_F64 lu (DecompositionFactory_ZDRM/lu
;;                                    (.numRows M) (.numCols M))]
;;       (or (not (.decompose lu M))
;;           (.isSingular lu))))

;;   (square? [_]
;;     (= (.numRows M) (.numCols M)))

;;   (symmetric? [_]
;;     (if (= (.numRows M) (.numCols M))
;;       (MatrixFeatures_ZDRM/isEquals M
;;                                     (CommonOps_ZDRM/transpose M nil)
;;                                     constant/tolerance--default)
;;       false))
;;   (symmetric? [_ tol]
;;     (if (= (.numRows M) (.numCols M))
;;       (MatrixFeatures_ZDRM/isEquals M
;;                                     (CommonOps_ZDRM/transpose M nil)
;;                                     tol)
;;       false))

;;   (unitary? [_]
;;     ;; For real matrices, “unitary” == orthogonal.
;;     (MatrixFeatures_ZDRM/isUnitary M constant/tolerance--default))

;;   mat/MatrixComplex

;;   (adjoint [_]
;;     (let [^ZMatrixRMaj out (ZMatrixRMaj. (.numCols M) (.numRows M))]
;;       (CommonOps_ZDRM/transposeConjugate M out)
;;       (->ComplexDense out)))
;;   (conjugate [_] (->ComplexDense (CommonOps_ZDRM/conjugate M nil)))
;;   (imag [_] (realdense/->RealDense (CommonOps_ZDRM/imaginary M nil)))
;;   (real [_] (realdense/->RealDense (CommonOps_ZDRM/real M nil)))

;;   (hermitian? [_]
;;     (MatrixFeatures_ZDRM/isHermitian M constant/tolerance--default))

;;   (real? [_]
;;     (realZM? M))

;;   Object
;;   (toString [_]
;;     (print! M)))
