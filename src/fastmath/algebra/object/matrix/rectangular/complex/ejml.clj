(ns fastmath.algebra.object.matrix.rectangular.complex.ejml
  "Implementing complex matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml].
  "
  ;; TODO: implement inc and dec etc.?
  (:require
   [clojure.string :as str]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as realdense]
   [fastmath.algebra.function.matrix :as gmat]
   [fastmath.default :as default]
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
   (org.ejml.dense.row.factory DecompositionFactory_ZDRM LinearSolverFactory_ZDRM)
   (org.ejml.interfaces.linsol LinearSolverDense)
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

(defn- ->C [^Complex_F64 CF64]
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

(defn- multiply--m [^ZMatrixRMaj m1 ^ZMatrixRMaj m2]
  (let [out (ZMatrixRMaj. (.numRows m1) (.numCols m2))]
    (CommonOps_ZDRM/mult m1 m2 out)
    out))

(defn- transpose [^ZMatrixRMaj m]
  (let [T (ZMatrixRMaj. (.numCols m) (.numRows m))]
    (CommonOps_ZDRM/transpose m T)
    T))

(defn- conjugate-transpose [^ZMatrixRMaj m]
  (let [^ZMatrixRMaj out (ZMatrixRMaj. (.numCols m) (.numRows m))]
    (CommonOps_ZDRM/transposeConjugate m out)
    out))

(declare <-rows)

(defn- kron
  "Kronecker product A (x) B on the interleaved data arrays."
  ^ZMatrixRMaj [^ZMatrixRMaj a ^ZMatrixRMaj b]
  (let [ra (.numRows a) ca (.numCols a) rb (.numRows b) cb (.numCols b)
        cols (* ca cb)
        out (ZMatrixRMaj. (* ra rb) cols)
        ^doubles ad (.data a) ^doubles bd (.data b) ^doubles od (.data out)]
    (dotimes [i ra]
      (dotimes [j ca]
        (let [ia (* 2 (+ (* i ca) j))
              xr (aget ad ia) xi (aget ad (inc ia))]
          (when-not (and (zero? xr) (zero? xi))
            (dotimes [k rb]
              (let [row (+ (* i rb) k)]
                (dotimes [l cb]
                  (let [ib (* 2 (+ (* k cb) l))
                        yr (aget bd ib) yi (aget bd (inc ib))
                        io (* 2 (+ (* row cols) (+ (* j cb) l)))]
                    (aset od io (- (* xr yr) (* xi yi)))
                    (aset od (inc io) (+ (* xr yi) (* xi yr)))))))))))
    out))

(def ^:private ^ThreadLocal lu-solvers
  ;; One LU solver per thread and size: EJML solvers keep internal workspace,
  ;; so they are reused within a thread and never shared between threads.
  (ThreadLocal/withInitial
   (reify java.util.function.Supplier
     (get [_] (java.util.HashMap.)))))

(defn- lu-solver ^LinearSolverDense [^long n]
  (let [^java.util.HashMap cache (.get lu-solvers)
        k (Long/valueOf n)]
    (or (.get cache k)
        (let [s (LinearSolverFactory_ZDRM/lu (int n))]
          (.put cache k s)
          s))))

(defn- zsolve
  "X with A X = B. Square A: LU with a per-thread reused solver (throws if
  singular); otherwise EJML's least-squares `CommonOps_ZDRM/solve`."
  ^ZMatrixRMaj [^ZMatrixRMaj A ^ZMatrixRMaj B]
  (let [X (ZMatrixRMaj. (.numCols A) (.numCols B))]
    (if (= (.numRows A) (.numCols A))
      (let [^LinearSolverDense s (lu-solver (.numRows A))]
        (when-not (.setA s (if (.modifiesA s) (.copy A) A))
          (throw (ex-info "singular matrix" {:shape [(.numRows A) (.numCols A)]})))
        (.solve s (if (.modifiesB s) (.copy B) B) X)
        ;; setA can succeed on an exactly singular matrix; the solution then
        ;; contains Inf/NaN (cheaper to detect here than with quality())
        (let [^doubles xd (.data X) n (alength xd)]
          (loop [t 0]
            (when (< t n)
              (if (Double/isFinite (aget xd t))
                (recur (inc t))
                (throw (ex-info "singular matrix" {:shape [(.numRows A) (.numCols A)]})))))))
      (CommonOps_ZDRM/solve A B X))
    X))

(def ^:private ^ThreadLocal work-matrices
  ;; per-thread scratch n x n matrices for `solve-shifted` (never escape)
  (ThreadLocal/withInitial
   (reify java.util.function.Supplier
     (get [_] (java.util.HashMap.)))))

(defn- work-matrix ^ZMatrixRMaj [^long n]
  (let [^java.util.HashMap cache (.get work-matrices)
        k (Long/valueOf n)]
    (or (.get cache k)
        (let [w (ZMatrixRMaj. (int n) (int n))]
          (.put cache k w)
          w))))

(defn- check-finite! [^ZMatrixRMaj X ^ZMatrixRMaj A]
  (let [^doubles xd (.data X) n (alength xd)]
    (loop [t 0]
      (when (< t n)
        (if (Double/isFinite (aget xd t))
          (recur (inc t))
          (throw (ex-info "singular matrix" {:shape [(.numRows A) (.numCols A)]})))))))

(defn- zsolve-shifted
  "X with (A + c diag(v)) X = B, without allocating the shifted matrix."
  ^ZMatrixRMaj [^ZMatrixRMaj A c-re c-im ^ZMatrixRMaj V ^ZMatrixRMaj B]
  (let [cr (double c-re) ci (double c-im)
        n (.numRows A)
        _ (when-not (and (= n (.numCols A)) (= n (* (.numRows V) (.numCols V))))
            (throw (ex-info "solve-shifted needs square A and a length-n vector"
                            {:A [(.numRows A) (.numCols A)] :v [(.numRows V) (.numCols V)]})))
        ^ZMatrixRMaj W (work-matrix n)
        ^doubles wd (.data W) ^doubles vd (.data V)
        ^LinearSolverDense s (lu-solver n)
        X (ZMatrixRMaj. n (.numCols B))]
    (System/arraycopy (.data A) 0 wd 0 (* 2 n n))
    (dotimes [i n]
      (let [vr (aget vd (* 2 i)) vi (aget vd (inc (* 2 i)))
            idx (* 2 (+ (* i n) i))]
        (aset wd idx (+ (aget wd idx) (- (* cr vr) (* ci vi))))
        (aset wd (inc idx) (+ (aget wd (inc idx)) (+ (* cr vi) (* ci vr))))))
    (when-not (.setA s W)          ; W is scratch, so modifiesA does not matter
      (throw (ex-info "singular matrix" {:shape [n n]})))
    (.solve s (if (.modifiesB s) (.copy B) B) X)
    (check-finite! X A)
    X))

(defn- zvectors
  "Rows of [re im] pairs."
  [^ZMatrixRMaj M]
  (let [rows (.numRows M) cols (.numCols M) ^doubles d (.data M)
        out (object-array rows)]
    (dotimes [i rows]
      (let [row (object-array cols)]
        (dotimes [j cols]
          (let [idx (* 2 (+ (* i cols) j))]
            (aset row j [(aget d idx) (aget d (inc idx))])))
        (aset out i (clojure.lang.LazilyPersistentVector/createOwning row))))
    (clojure.lang.LazilyPersistentVector/createOwning out)))

(defn- fill-rows!
  "Fills Z from rows whose entries are [re im] pairs, real numbers, or
  objects with re/im (fastmath complex numbers)."
  ^ZMatrixRMaj [^ZMatrixRMaj Z rows]
  (let [cols (.numCols Z)]
    (loop [i 0 rs (seq rows)]
      (when rs
        (loop [j 0 es (seq (first rs))]
          (when es
            (let [e (first es)]
              (cond
                (number? e) (.set Z (int i) (int j) (double e) 0.0)
                (sequential? e) (.set Z (int i) (int j) (double (nth e 0)) (double (nth e 1)))
                :else (.set Z (int i) (int j) (double (cc/re e)) (double (cc/im e)))))
            (recur (inc j) (next es))))
        (recur (inc i) (next rs))))
    Z))

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
    (let [^ZMatrixRMaj A (.copy M)
          ^doubles d (.data A)
          n (alength d)]
      (loop [i 0]
        (when (< i n)
          (aset-double d i       (+ (aget d i) re))
          (aset-double d (inc i) (+ (aget d (inc i)) im))
          (recur (+ i 2))))
      (ComplexDense. A)))

  (inner [this other] (gmat/inner this other))

  (kronecker [_ other] (ComplexDense. (kron M (.M ^ComplexDense other))))

  (fmap [_ f]
    ;; Operates on a copy, like every other operation here - `M` itself is
    ;; never mutated, since it may be shared with other ComplexDense values.
    (let [^ZMatrixRMaj A (.copy M)
          ^doubles d (.getData A)
          n (alength d)]
      (loop [i 0]
        (when (< i n)
          (let [re (aget d i)
                im (aget d (unchecked-inc-int i))
                [re' im'] (f re im)]
            (aset-double d i (double re'))
            (aset-double d (unchecked-inc-int i) (double im'))
            (recur (unchecked-add-int i 2)))))
      (ComplexDense. A)))

  (multiply [_ other]
    (ComplexDense. (multiply--m M (.M ^ComplexDense other))))

  (multiply--s [_ [r i]]
    (let [A (.copy M)]
      (CommonOps_ZDRM/scale (double r) (double i) A)
      (ComplexDense. A)))
  (multiply--e [_ other]
    (let [^ZMatrixRMaj B (.M ^ComplexDense other)
          out (ZMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_ZDRM/elementMultiply M B out)
      (ComplexDense. out)))
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
    (->C (CommonOps_ZDRM/det M)))
  (trace [_]
    (->C (CommonOps_ZDRM/trace M nil)))

;; ==================================================
  solve/MatrixSolve
;; ==================================================
  (cholesky [_]
    (let [n (.numRows ^ZMatrixRMaj M)
          _ (when (not= n (.numCols ^ZMatrixRMaj M))
              (throw (ex-info "Cholesky requires square matrix" {:shape [n (.numCols M)]})))
          ^CholeskyDecomposition_F64 chol (DecompositionFactory_ZDRM/chol n true)]
      ;; `decompose` may modify its input, so hand it a copy - `M` itself
      ;; must stay untouched, like every other operation here.
      (when-not (.decompose chol (.copy M))
        (throw (ex-info "Cholesky failed (matrix not SPD)"
                        {:shape [n n]})))
      (let [L (.getT chol (ZMatrixRMaj. n n))]   ; lower if 'true' above
        {:L (ComplexDense. L) :lower? true :spd? true})))

  (solve [_ v]
    (ComplexDense. (zsolve M (.M ^ComplexDense v))))
;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape    [_] [(.numRows M) (.numCols M)])
  (nrows [_] (.numRows M))
  (ncols [_] (.numCols M))

  ;; -------- Retrieval --------
  (element [_ i j]
    (->C (getc M (long i) (long j))))

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
          ^doubles data (.getData M)
          out  (make-array Complex_F64 rows cols)]
      (dotimes [i rows]
        (dotimes [j cols]
          (let [idx (* 2 (+ j (* i cols)))]
            (aset out i j
                  (Complex_F64.
                   (aget data idx)
                   (aget data (inc idx)))))))
      out))

  (->vectors [_] (zvectors M))

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
  (square? [this] (gmat/square? this))
  (normal? [this] (gmat/normal? this))
  (symmetric? [this] (gmat/symmetric? this))
  (singular? [this] (gmat/singular? this))
  (unitary? [_] (MatrixFeatures_ZDRM/isUnitary M default/tolerance))

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
    (MatrixFeatures_ZDRM/isEquals M ^ZMatrixRMaj (.-M ^ComplexDense B) default/tolerance))

  (toString [_]
    (print! M)))

;; -------------------------------------------------------------------
;; Constructors 
;; -------------------------------------------------------------------

(defn complexdense
  ;; zero
  (^ComplexDense [^long n ^long o]
   (->ComplexDense (ZMatrixRMaj. n o)))

  ;; diagonal, from interleaved [re0 im0 re1 im1 ...]
  (^ComplexDense [data]
   (-> (CommonOps_ZDRM/diag (if (instance? (Class/forName "[D") data) ^doubles data (double-array (flatten data))))
       ->ComplexDense))

  ;; elements, row-major interleaved: a double[] (used as is) or any nested collection
  (^ComplexDense [^long n ^long o data]
   (->ComplexDense (ZMatrixRMaj. n o true (if (instance? (Class/forName "[D") data)
                                            ^doubles data
                                            (double-array (flatten data)))))))

(defn <-real ^ComplexDense [^fastmath.algebra.object.matrix.rectangular.real.ejml.RealDense M]
  (let [^DMatrixRMaj A (.-M M)
        Z (ZMatrixRMaj. (.numRows A) (.numCols A))]
    (CommonOps_ZDRM/convert A Z)
    (->ComplexDense Z)))

(defn <-rows
  "From rows whose entries are [re im] pairs, real numbers or complex numbers."
  ^ComplexDense [rows]
  (->ComplexDense (fill-rows! (ZMatrixRMaj. (count rows) (count (first rows))) rows)))

(defn <-interleaved
  "nrows x ncols matrix from a row-major interleaved double[] (copied)."
  ^ComplexDense [^long nrows ^long ncols ^doubles data]
  (->ComplexDense (ZMatrixRMaj. nrows ncols true (aclone data))))

(defn ->interleaved
  "Copy of the row-major interleaved data."
  ^doubles [^ComplexDense A]
  (aclone ^doubles (.data ^ZMatrixRMaj (.M A))))

(defn solve-shifted
  "X with (A + c diag(v)) X = B for square A, a length-n vector v and a
  complex coefficient c ([re im] or real). Same as
  (solve (add-diagonal A c v) B) but works in a per-thread scratch matrix,
  so repeated shifted solves (resolvents, Doppler classes) allocate only X."
  ^ComplexDense [^ComplexDense A c ^ComplexDense v ^ComplexDense B]
  (->ComplexDense (zsolve-shifted (.M A)
                                  (double (if (number? c) c (nth c 0)))
                                  (double (if (number? c) 0.0 (nth c 1)))
                                  (.M v) (.M B))))

(defn add-diagonal
  "A + c diag(v): a copy of square A with c v_i added to entry (i, i). v is an
  n x 1 or 1 x n matrix, c a complex coefficient ([re im] or real)."
  ^ComplexDense [^ComplexDense A c ^ComplexDense v]
  (let [^ZMatrixRMaj M (.copy ^ZMatrixRMaj (.M A))
        ^ZMatrixRMaj V (.M v)
        n (.numRows M)
        cr (double (if (number? c) c (nth c 0)))
        ci (double (if (number? c) 0.0 (nth c 1)))
        ^doubles md (.data M) ^doubles vd (.data V)]
    (when-not (and (= n (.numCols M)) (= n (* (.numRows V) (.numCols V))))
      (throw (ex-info "add-diagonal needs square A and a length-n vector"
                      {:A [(.numRows M) (.numCols M)] :v [(.numRows V) (.numCols V)]})))
    (dotimes [i n]
      (let [vr (aget vd (* 2 i)) vi (aget vd (inc (* 2 i)))
            idx (* 2 (+ (* i n) i))]
        (aset md idx (+ (aget md idx) (- (* cr vr) (* ci vi))))
        (aset md (inc idx) (+ (aget md (inc idx)) (+ (* cr vi) (* ci vr))))))
    (->ComplexDense M)))

(deftype ComplexSparse [^long nrows ^long ncols ^ints idx ^doubles vals]
  ;; Nonzero entries of a complex matrix: row-major linear indices `idx` and
  ;; interleaved values `vals`. Immutable; made by `sparse`, consumed by
  ;; `linear-combination` (which then visits only these entries).
  Object
  (toString [_] (format "#ComplexSparse[%d x %d, %d nonzeros]" nrows ncols (alength idx))))

(defn sparse
  "Sparse copy (nonzero entries only) of a ComplexDense, for use as a term in
  `linear-combination`."
  ^ComplexSparse [^ComplexDense A]
  (let [^ZMatrixRMaj M (.M A)
        ^doubles d (.data M)
        n (* (.numRows M) (.numCols M))
        nz (vec (filter #(not (and (zero? (aget d (* 2 (long %)))) (zero? (aget d (inc (* 2 (long %))))))) (range n)))]
    (ComplexSparse. (.numRows M) (.numCols M)
                    (int-array nz)
                    (double-array (mapcat (fn [k] [(aget d (* 2 (long k))) (aget d (inc (* 2 (long k))))]) nz)))))

(defn- term-rows ^long [t]
  (if (instance? ComplexSparse t) (.-nrows ^ComplexSparse t) (.numRows ^ZMatrixRMaj (.M ^ComplexDense t))))
(defn- term-cols ^long [t]
  (if (instance? ComplexSparse t) (.-ncols ^ComplexSparse t) (.numCols ^ZMatrixRMaj (.M ^ComplexDense t))))

(defn linear-combination
  "sum_k c_k A_k in one pass, for complex coefficients c_k ([re im] or a
  real number) and same-shape terms A_k: ComplexDense, or ComplexSparse (from
  `sparse`), for which only the nonzero entries are visited."
  ^ComplexDense [coeffs mats]
  (let [t0 (first mats)
        r (term-rows t0) c (term-cols t0)
        out (ZMatrixRMaj. (int r) (int c))
        ^doubles od (.data out)
        n (alength od)]
    (loop [cs (seq coeffs) ms (seq mats)]
      (when (and cs ms)
        (let [k (first cs)
              cr (double (if (number? k) k (nth k 0)))
              ci (double (if (number? k) 0.0 (nth k 1)))
              t (first ms)]
          (when-not (and (== r (term-rows t)) (== c (term-cols t)))
            (throw (ex-info "shape mismatch" {:expected [r c] :got [(term-rows t) (term-cols t)]})))
          (if (instance? ComplexSparse t)
            (let [^ints ix (.-idx ^ComplexSparse t) ^doubles vs (.-vals ^ComplexSparse t)]
              (dotimes [q (alength ix)]
                (let [o (* 2 (aget ix q)) a (aget vs (* 2 q)) b (aget vs (inc (* 2 q)))]
                  (aset od o (+ (aget od o) (- (* cr a) (* ci b))))
                  (aset od (inc o) (+ (aget od (inc o)) (+ (* cr b) (* ci a)))))))
            (let [^doubles ad (.data ^ZMatrixRMaj (.M ^ComplexDense t))]
              (if (zero? ci)
                (loop [q 0]
                  (when (< q n)
                    (aset od q (+ (aget od q) (* cr (aget ad q))))
                    (recur (inc q))))
                (loop [q 0]
                  (when (< q n)
                    (let [a (aget ad q) b (aget ad (inc q))]
                      (aset od q (+ (aget od q) (- (* cr a) (* ci b))))
                      (aset od (inc q) (+ (aget od (inc q)) (+ (* cr b) (* ci a)))))
                    (recur (+ q 2)))))))
          (recur (next cs) (next ms)))))
    (->ComplexDense out)))

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

