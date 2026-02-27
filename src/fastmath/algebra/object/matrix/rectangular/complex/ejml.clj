(ns fastmath.algebra.object.matrix.rectangular.complex.ejml
  "Implementing complex matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml].
  "
  ;; TODO: implement inc and dec etc.?
  (:require
   [clojure.string :as str]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as realdense]
   [fastmath.default :as default]
   [fastmath.algebra.object.number.complex.ejml :as Craw]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.solve :as solve]
   [fastmath.protocol.algebra.object.matrix.square :as p-sqmat]
   [fastmath.protocol.algebra.object.matrix.predicate :as pred]
   [fastmath.protocol.algebra.object.matrix.complex-predicate :as cpred]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.structure.additive.group :as ag]
   [fastmath.protocol.algebra.structure.additive.monoid :as am]
   [fastmath.protocol.algebra.structure.additive.semigroup :as asg]
   [fastmath.protocol.algebra.structure.coordinate.complex :as cc]
   [fastmath.protocol.algebra.structure.module :as module]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]
   [fastmath.protocol.representation.d2 :as d2])
  (:import
   (java.lang Math)
   (org.ejml.interfaces.decomposition CholeskyDecomposition_F64)
   (org.ejml.dense.row.factory DecompositionFactory_ZDRM)
   (org.ejml.simple.ops SimpleOperations_ZDRM)
   (org.ejml.data Complex_F64 ZMatrixRMaj DMatrixRMaj)
   (org.ejml.dense.row CommonOps_ZDRM NormOps_ZDRM MatrixFeatures_ZDRM)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
;; (println "=== start ===")
;; ^^^ Used for debugging

;; ==================================================
;; Utilities 
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

(defn- ->C [CF64]
  (C/i (.getReal CF64) (.getImaginary CF64)))

(defn- getc
  "Utility for returning complex numbers from a matrix."
  [^ZMatrixRMaj m ^long i ^long j]
  (let [^Complex_F64 cnum   (Complex_F64. 0.0 0.0)]
    (.get m i j cnum)
    cnum))

(defn- setc
  ;; TODO: Should have a multiple element version 
  "Utility for assigning complex numbers to a matrix."
  ([^ZMatrixRMaj m  i j cr ci]
   (let [A (.copy m)]
     (.set A i j cr ci) A))
  ([^ZMatrixRMaj m  i j ^Complex_F64 cnum]
   (let [A (.copy m)]
     (.set A i j (.-real cnum) (.-imaginary cnum))
     A)))

(defn- multiply--m [m1 m2]
  (let [out (ZMatrixRMaj. (.numRows m1) (.numCols m2))]
    (CommonOps_ZDRM/mult m1 m2 out)
    out))

(defn- transpose [m]
  (let [T (ZMatrixRMaj. (.numCols m) (.numRows m))]
    (CommonOps_ZDRM/transpose m T)
    T))

(defn- conjugate-transpose [^ZMatrixRMaj m]
  (let [^ZMatrixRMaj out (ZMatrixRMaj. (.numCols m) (.numRows m))]
    (CommonOps_ZDRM/transposeConjugate m out)
    out))

(defn  format--real ^String [^double x ^long p ^double eps]
  (format (str "%." p "f") (if (< (Math/abs x) eps) 0.0 x)))
(defn  format--complex ^String [^double re ^double im ^long p ^double eps]
  (let [rz (if (< (Math/abs re) eps) 0.0 re)
        iz (if (< (Math/abs im) eps) 0.0 im)]
    (str (format (str "%." p "f") rz)
         (format (str "%+" p "f") iz)
         "i")))

(defn- realZM?
  ([^ZMatrixRMaj M]
   (realZM? M default/tolerance))
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
                    :or   {precision 3 eps default/tolerance max-rows 12 max-cols 12}}]
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 ComplexDense                                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; (deftype ComplexDense [^ZMatrixRMaj M]
;;   rmat/RectangularMatrix
;;   asg/AdditiveSemigroup
;;   am/AdditiveMonoid
;;   ag/AdditiveGroup
;;   nspace/NormedSpace

;;   module/Module
;;   emat/MatrixExtra
;;   d2/D2
;;   cc/ComplexCoordinate
;;   cmat/ComplexMatrix
;;   Object)
;;   
(deftype ComplexDense [^ZMatrixRMaj M]
;; ==================================================
  rmat/RectangularMatrix
;; ==================================================
  asg/AdditiveSemigroup
  (add [_ other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_ZDRM/add M B out)
      (ComplexDense. out)))

  am/AdditiveMonoid
  (zero [_] (ComplexDense. (ZMatrixRMaj. (.numRows M) (.numCols M))))

  ag/AdditiveGroup
  (negate [_]
    (let [A (.copy M)]
      (CommonOps_ZDRM/scale -1.0 0.0 A)
      (ComplexDense. A)))

  nspace/NormedSpace
  (norm [_]
    ;; frobenius matrix norm
    (double (NormOps_ZDRM/normF M)))

  module/Module
  (scale [_ z]
    (let [A (.copy M)
          [r i] (if (seqable? z) z [z 0.0])]
      (CommonOps_ZDRM/scale (double r) (double i) A)
      (ComplexDense. A)))

;; ==================================================
  emat/MatrixExtra
;; ==================================================
  (diagonal [_]
    (let [n (min (.numRows M) (.numCols M))
          out (ZMatrixRMaj. n 1)]
      (CommonOps_ZDRM/extractDiag M out)
      (ComplexDense. out)))

  (add--s [_ ^doubles [re im]]
    (let [A (.copy M)
          d (.data A)
          n (alength d)]
      (loop [i 0]
        (when (< i n)
          (aset-double d i       (+ (aget d i) re))
          (aset-double d (inc i) (+ (aget d (inc i)) im))
          (recur (+ i 2))))
      (ComplexDense. A)))

  (inner [_ other]
    (let [B (.M ^ComplexDense other)
          ^doubles da (.getData M)
          ^doubles db (.getData B)
          n (alength da)]
      (when-not (== n (alength db))
        (throw (ex-info "Shape mismatch for inner product"
                        {:this-len n :other-len (alength db)})))
      (loop [i 0
             rsum 0.0
             isum 0.0]
        (if (>= i n)
          (C/i rsum isum)
          (let [a-re (aget da i)
                a-im (aget da (inc i))
                b-re (aget db i)
                b-im (aget db (inc i))]
          ;; conj(a) * b = (a_re - i a_im) * (b_re + i b_im)
            (recur (+ i 2)
                   (+ rsum (+ (* a-re b-re) (* a-im b-im)))
                   (+ isum (- (* a-re b-im) (* a-im b-re)))))))))

  (kronecker [_ other]
    (let [B (.M ^ComplexDense other)
          a-rows (.getNumRows M)
          a-cols (.getNumCols M)
          b-rows (.getNumRows B)
          b-cols (.getNumCols B)
          out    (ZMatrixRMaj. (* a-rows b-rows) (* a-cols b-cols))]
      (dotimes [i a-rows]
        (dotimes [j a-cols]
          (let [cnum1 (Complex_F64.)
                ^org.ejml.data.Complex_F64 aij (.get M i j cnum1)
                ar (.-real cnum1)
                ai (.-imaginary cnum1)]
            (dotimes [p b-rows]
              (dotimes [q b-cols]
                (let [cnum2 (Complex_F64.)
                      ^org.ejml.data.Complex_F64 bij (.get B p q cnum2)
                      br (.-real cnum2)
                      bi (.-imaginary cnum2)]
                ;; (ar + i ai) * (br + i bi)
                  (.set out
                        (+ (* i b-rows) p)
                        (+ (* j b-cols) q)
                        (- (* ar br) (* ai bi))
                        (+ (* ar bi) (* ai br)))))))))
      (ComplexDense. out)))

  (fmap [_ f]
    ;; TODO: VERY INEFFICIENT!
    ;; - Prevent function from multiple calls to set
    ;; - Work on (.data A) directly
    ;; - check number of arguments and branch?
    (let [^doubles d (.getData M)
          n (alength d)]
      (loop [i 0]
        (when (< i n)
          (let [re (aget d i)
                im (aget d (unchecked-inc-int i))
                [re' im'] (f re im)]
            (aset-double d i (double re'))
            (aset-double d (unchecked-inc-int i) (double im'))
            (recur (unchecked-add-int i 2)))))
      (ComplexDense. M)))

  (multiply [_ other]
    (ComplexDense. (multiply--m M (.M ^ComplexDense other))))

  (multiply--s [_ [r i]]
    (let [A (.copy M)]
      (CommonOps_ZDRM/scale (double r) (double i) A)
      (ComplexDense. A)))
  (multiply--e [this other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)]
      (ComplexDense. (CommonOps_ZDRM/elementMultiply M B nil nil))))
  (multiply--v [_ v]
    ;; v is expected to be (n×1) ComplexDense
    (let [^ZMatrixRMaj x (.M ^ComplexDense v)
          out (ZMatrixRMaj. (.numRows M) 1)]
      (CommonOps_ZDRM/mult M x out)
      (ComplexDense. out)))
  (subtract [_ other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_ZDRM/subtract M B out)
      (ComplexDense. out)))

  (transpose [_]
    (ComplexDense. (transpose M)))

;; ==================================================
  p-sqmat/SquareMatrix
;; ==================================================
  (determinant [_]
    (CommonOps_ZDRM/det M))
  (trace [_]
    (CommonOps_ZDRM/trace M nil))

;; ==================================================
  solve/MatrixSolve
;; ==================================================
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

  (solve [this v]
    (let [^ZMatrixRMaj B (.M ^ComplexDense v)
          X (ZMatrixRMaj. (.numRows B) (.numCols B))]
      (CommonOps_ZDRM/solve M B X)
      (ComplexDense. X)))
;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape    [_] [(.numRows M) (.numCols M)])
  (nrows [_] (.numRows M))
  (ncols [_] (.numCols M))

  ;; -------- Retrieval --------
  (element [_ i j]
    (Craw/->ComplexNumber (getc M (long i) (long j))))

  (column [_ j]
    (ComplexDense. (extract M [0 (.numRows M)] [j (+ 1 (long j))])))

  (row [_ i]
    (ComplexDense. (extract M [i (+ 1 (long i))] [0 (.numCols M)])))
  (columns [_]
    (mapv (fn [^long ic]
            (ComplexDense. (extract M [0 (.numRows M)] [ic (+ 1 ic)])))
          (range (.numCols M))))

  (rows [_]
    (mapv (fn [^long ir]
            (ComplexDense. (extract M [ir (+ ir 1)] [0 (.numCols M)])))
          (range (.numRows M))))

  (->arrays ^doubles [_]
    (let [rows (.numRows M)
          cols (.numCols M)
          data (.getData M)
          out  (make-array Complex_F64 rows cols)]
      (dotimes [i rows]
        (dotimes [j cols]
          (let [idx (+ (* 2 (+ j (* i cols))))]
            (aset out i j
                  (Complex_F64.
                   (aget data idx)
                   (aget data (inc idx)))))))
      out))

  (->vectors [_]
    (let [rows (.numRows M)
          cols (.numCols M)
          data (.getData M)]
      (mapv
       (fn [i]
         (mapv
          (fn [j]
            (let [idx (+ (* 2 (+ j (* i cols))))]
              [(aget data idx)
               (aget data (inc idx))]))
          (range cols)))
       (range rows))))

;; ==================================================
  cc/ComplexCoordinate
;; ==================================================
  (re [_] (realdense/->RealDense (CommonOps_ZDRM/real M nil)))
  (im [_] (realdense/->RealDense (CommonOps_ZDRM/imaginary M nil)))
  (conjugate [_] (ComplexDense. (CommonOps_ZDRM/conjugate M nil)))

;; ==================================================
  cmat/ComplexMatrix
;; ==================================================
  (adjoint [_]
    (ComplexDense. (conjugate-transpose M)))

  (real? [_]
    (realZM? M))
;; ==================================================
  pred/MatrixPredicate
;; ==================================================
  (square? [_] (= (.numRows M) (.numCols M)))
  (normal? [_] (and (= (.numRows M) (.numCols M))
                    (= (multiply--m M (conjugate-transpose M))
                       (multiply--m (conjugate-transpose M) M))))
  (symmetric? [_] (and (= (.numRows M) (.numCols M))
                       (= M (transpose M))))
  (singular? [_] (and (= (.numRows M) (.numCols M))
                      (= (.determinantComplex SimpleOperations_ZDRM M) (C/i 0.0))))
  (unitary? [_] (.isUnitary M default/tolerance))

;; ==================================================
  cpred/ComplexPredicate
;; ==================================================
  (hermitian? [_]
    (MatrixFeatures_ZDRM/isHermitian M default/tolerance))
;; ==================================================
  clojure.lang.Seqable
;; ==================================================
  (seq [_] (seq (.getData M)))
;; ==================================================
  clojure.lang.Sequential
;; ==================================================
;; ==================================================
  Object
;; ==================================================

  (equals [_ B]
    (MatrixFeatures_ZDRM/isEquals M (.-M B) default/tolerance))

  (toString [_]
    (print! M)))

;; -------------------------------------------------------------------
;; Constructors 
;; -------------------------------------------------------------------

(defn complexdense
  ;; I
  ;; (^ComplexDense [^long n]
  ;;  (->ComplexDense (CommonOps_ZDRM/identity n)))

  ;; zero
  (^ComplexDense [^long n ^long o]
   (->ComplexDense (ZMatrixRMaj. n o)))

  ;; diagonal
  (^ComplexDense [^doubles data]
   (-> (CommonOps_ZDRM/diag (double-array (flatten data)))
       ->ComplexDense))

;; elements
  (^ComplexDense [^long n ^long o ^doubles data]
   (->ComplexDense (ZMatrixRMaj. n o true (double-array (flatten data))))))

(defn <-real ^ComplexDense [^fastmath.algebra.object.matrix.rectangular.real.ejml.RealDense M]
  (let [A  (.-M  M)
        Z (ZMatrixRMaj. ^long (.numRows A) ^long (.numCols A))]
    (CommonOps_ZDRM/convert A Z)
    (->ComplexDense Z)))

(defn <-rows ^ComplexDense [rows]
  (let [nrows (count rows)
        ncols (count (first rows))

        data  (double-array (flatten rows))]
    (->ComplexDense (ZMatrixRMaj. nrows ncols true data))))

(defn <-cols ^ComplexDense [cols]
  (let [ncols (count cols)
        nrows (count (first cols))
        data  (double-array (flatten cols))]
    (->ComplexDense (ZMatrixRMaj. nrows ncols false data))))

(comment (println "test")
         (def A--test (complexdense 2 2 (double-array [1 1 0 0 1 0 0 0])))
         (def B--test (complexdense 2 2 (double-array [1 0 1 0 1 0 0 1])))

         (d2/element (complexdense 3 3) 0 0)

         (println A--test)
         (println B--test)
         (-> (cmat/add A--test B--test)
             println)

         (->> (ComplexDense. (ZMatrixRMaj. 3 3))
              .-M
              .getData
              seq))

