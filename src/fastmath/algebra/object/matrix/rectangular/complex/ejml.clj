(ns fastmath.algebra.object.matrix.rectangular.complex.ejml
  "Implementing complex matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml].
  "
  (:require
   [clojure.string :as str]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as realdense]
   [fastmath.default :as default]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
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
   (org.ejml.data Complex_F64 ZMatrixRMaj DMatrixRMaj)
   (org.ejml.dense.row CommonOps_ZDRM NormOps_ZDRM MatrixFeatures_ZDRM)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(println "=== start ===")
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
          [r i] (if (coll? z) z [z 0.0])]
      (CommonOps_ZDRM/scale (double r) (double i) A)
      (ComplexDense. A)))

;; ==================================================
  emat/MatrixExtra
;; ==================================================
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

  (outer [_ other]
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

  (map [_ f]
    ;; TODO: Prevent function from multiple calls to set
    ;; - Work on (.data A) directly
    (let [^ZMatrixRMaj A (.copy M)]
      (dotimes [i (.numRows A)]
        (dotimes [j (.numCols A)]
          (setc A i j (f (getc A i j)))))
      (ComplexDense. A)))
  (multiply [_ other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols B))]
      (CommonOps_ZDRM/mult M B out)
      (ComplexDense. out)))
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
    (let [T (ZMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_ZDRM/transpose M T)
      (ComplexDense. T)))

  (square? [_]
    (= (.numRows M) (.numCols M)))

;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape    [_] [(.numRows M) (.numCols M)])
  (nrows [_] (.numRows M))
  (ncols [_] (.numCols M))

  ;; -------- Retrieval --------
  (element [_ i j]
    (getc M (long i) (long j)))

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

  ;; (diagonal [_]
  ;;   (let [n (min (.numRows M) (.numCols M))
  ;;         out (ZMatrixRMaj. n 1)]
  ;;     (CommonOps_ZDRM/extractDiag M out)
  ;;     (ComplexDense. out)))

  (array<- [_]
    (let [r (.numRows M) c (.numCols M)
          out (make-array Double/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-double out i j (getc M i j))))
      out))
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
    (let [^ZMatrixRMaj out (ZMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_ZDRM/transposeConjugate M out)
      (ComplexDense. out)))

  ;; (hermitian? [_]
  ;;   (MatrixFeatures_ZDRM/isHermitian M default/tolerance))

  (real? [_]
    (realZM? M))

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

         (complexdense 3 3)
         (println A--test)
         (println B--test)
         (-> (cmat/add A--test B--test)
             println)

         (-> (->ComplexDense (ZMatrixRMaj. 3 3))
             println))

