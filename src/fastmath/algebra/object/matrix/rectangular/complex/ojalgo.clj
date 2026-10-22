(ns fastmath.algebra.object.matrix.rectangular.complex.ojalgo
  "Implementing complex matrices using ojAlgo as a backend
  [https://github.com/optimatika/ojAlgo].

  Uses ojAlgo's high-level `MatrixC128`/`BasicMatrix` API (immutable,
  functional) rather than the low-level `PhysicalStore`/`GenericStore`
  layer - the latter needs reflection to work around access restrictions
  (see the abandoned `fastmath.matrix.dense.scratch` experiments) and
  buys nothing here.
  "
  (:require
   [clojure.string :as str]
   [fastmath.algebra.object.matrix.rectangular.real.ojalgo :as realdense]
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
   (org.ojalgo.matrix MatrixC128)
   (org.ojalgo.matrix.decomposition Cholesky)
   (org.ojalgo.scalar ComplexNumber)
   (org.ojalgo.structure Access2D)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; ==================================================
;; Utilities
;; ==================================================

(defn- ->C
  "ojAlgo's own ComplexNumber -> fastmath's canonical ComplexNumber."
  [^ComplexNumber c]
  (C/i (.getReal c) (.getImaginary c)))

(defn- ->oj
  "A [re im] pair -> ojAlgo's own ComplexNumber."
  ^ComplexNumber [[re im]]
  (ComplexNumber/of (double re) (double im)))

(defn- ^MatrixC128 ->matrix
  "Wrap nested rows of [re im] pairs as a MatrixC128."
  [rows]
  (->> rows
       (map (fn [row] (into-array ComplexNumber (map ->oj row))))
       into-array
       Access2D/wrap
       (.makeWrapper MatrixC128/FACTORY)))

(declare <-rows)

(defn- format--real ^String [^double x ^long p ^double eps]
  (format (str "%." p "f") (if (< (Math/abs x) eps) 0.0 x)))

(defn- format--complex ^String [^double re ^double im ^long p ^double eps]
  (let [rz (if (< (Math/abs re) eps) 0.0 re)
        iz (if (< (Math/abs im) eps) 0.0 im)]
    (str (format (str "%." p "f") rz)
         (format (str "%+" p "f") iz)
         "i")))

(defn- real?--matrix
  "Whether every imaginary part is within tolerance of zero."
  ([^MatrixC128 M] (real?--matrix M default/tolerance))
  ([^MatrixC128 M tolerance]
   (every? #(< (Math/abs ^double %) ^double tolerance)
           (seq (realdense/->RealDense (.getImaginary M))))))

(defn ->str
  "Pretty string for MatrixC128 with aligned columns.
   Options: :precision (3), :eps (1e-12), :max-rows (12), :max-cols (12)"
  ([^MatrixC128 A] (->str A {}))
  ([^MatrixC128 A {:keys [precision eps max-rows max-cols]
                   :or   {precision 3 eps default/tolerance max-rows 12 max-cols 12}}]
   (let [num-rows (.countRows A) num-cols (.countColumns A)
         rlim (min num-rows ^long max-rows) clim (min num-cols ^long max-cols)
         real? (real?--matrix A eps)
         cell (fn [^long i ^long j]
                (let [^ComplexNumber c (.get A i j)]
                  (if real?
                    (format--real (.getReal c) precision eps)
                    (format--complex (.getReal c) (.getImaginary c) precision eps))))
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
  ([^MatrixC128 A] (print! A {}))
  ([^MatrixC128 A opts]
   (println (->str A opts))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                ComplexDense                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype ComplexDense [^MatrixC128 M]
;; ==================================================
  rmat/RectangularMatrix
;; ==================================================
  asg/AdditiveSemigroup
  (add [_ other] (ComplexDense. (.add M (.-M ^ComplexDense other))))

  am/AdditiveMonoid
  (zero [_] (ComplexDense. (.make MatrixC128/FACTORY (.countRows M) (.countColumns M))))

  ag/AdditiveGroup
  (negate [_] (ComplexDense. (.negate M)))

  module/Module
  (scale [_ z]
    (let [[r i] (if (seqable? z) z [z 0.0])]
      (ComplexDense. (.multiply M (ComplexNumber/of (double r) (double i))))))

  nspace/NormedSpace
  (norm [_] (.norm M))

;; ==================================================
  emat/MatrixExtra
;; ==================================================
  (diagonal [_]
    (let [n (min (.countRows M) (.countColumns M))]
      (<-rows (for [i (range n)]
                [[(.getReal ^ComplexNumber (.get M i i))
                  (.getImaginary ^ComplexNumber (.get M i i))]]))))

  (add--s [_ [re im]]
    (ComplexDense. (.add M (ComplexNumber/of (double re) (double im)))))

  (inner [this other] (gmat/inner this other))

  (kronecker [this other] (<-rows (gmat/kronecker-rows this other)))

  (fmap [_ f]
    (let [rows (.countRows M) cols (.countColumns M)]
      (<-rows (for [i (range rows)]
                (for [j (range cols)]
                  (let [^ComplexNumber c (.get M i j)]
                    (f (.getReal c) (.getImaginary c))))))))

  (multiply [_ other] (ComplexDense. (.multiply M (.-M ^ComplexDense other))))
  (multiply--s [this z] (module/scale this z))
  (multiply--e [_ other]
    (let [^MatrixC128 B (.-M ^ComplexDense other)
          rows (.countRows M) cols (.countColumns M)]
      (<-rows (for [i (range rows)]
                (for [j (range cols)]
                  (let [^ComplexNumber a (.get M i j)
                        ^ComplexNumber b (.get B i j)
                        ^ComplexNumber p (.multiply a b)]
                    [(.getReal p) (.getImaginary p)]))))))
  (multiply--v [_ v] (ComplexDense. (.multiply M (.-M ^ComplexDense v))))
  (subtract [_ other] (ComplexDense. (.subtract M (.-M ^ComplexDense other))))
  (transpose [_] (ComplexDense. (.transpose M)))

;; ==================================================
  p-sqmat/SquareMatrix
;; ==================================================
  (determinant [_] (->C (.getDeterminant M)))
  (trace [_] (->C (.getTrace M)))

;; ==================================================
  solve/MatrixSolve
;; ==================================================
  (cholesky [_]
    (let [n (.countRows M)
          _ (when (not= n (.countColumns M))
              (throw (ex-info "Cholesky requires square matrix" {:shape [n (.countColumns M)]})))
          chol (.decompose Cholesky/C128 M)]
      (when-not (.isSPD chol)
        (throw (ex-info "Cholesky failed (matrix not SPD)" {:shape [n n]})))
      {:L (ComplexDense. (.makeWrapper MatrixC128/FACTORY (.getL chol)))
       :lower? true :spd? true}))

  (solve [_ v] (ComplexDense. (.solve M (.-M ^ComplexDense v))))

;; ==================================================
  cc/ComplexCoordinate
;; ==================================================
  (re [_] (realdense/->RealDense (.getReal M)))
  (im [_] (realdense/->RealDense (.getImaginary M)))
  (conjugate [_]
    ;; NOTE: BasicMatrix's own `.conjugate()` is actually the Hermitian
    ;; adjoint (conjugate TRANSPOSE) - confirmed empirically (it changes
    ;; a 1x2 matrix's shape to 2x1). `ComplexCoordinate/conjugate` here
    ;; means plain elementwise conjugation with the shape unchanged, so
    ;; it can't be reused for that - see `adjoint` below instead.
    (let [rows (.countRows M) cols (.countColumns M)]
      (<-rows (for [i (range rows)]
                (for [j (range cols)]
                  (let [^ComplexNumber c (.get M i j)]
                    [(.getReal c) (- (.getImaginary c))]))))))

;; ==================================================
  cmat/ComplexMatrix
;; ==================================================
  (adjoint [_] (ComplexDense. (.conjugate M)))
  (real? [_] (real?--matrix M))

;; ==================================================
  pred/MatrixPredicate
;; ==================================================
  (square? [this] (gmat/square? this))
  (normal? [this] (gmat/normal? this))
  (symmetric? [this] (gmat/symmetric? this))
  (singular? [this] (gmat/singular? this))
  (unitary? [this]
    ;; No dedicated primitive on BasicMatrix - compare M* M to the
    ;; identity, built locally (no cross-namespace dependency needed).
    (and (gmat/square? this)
         (let [n (.countRows M)
               I (.makeEye MatrixC128/FACTORY n n)]
           (= (ComplexDense. (.multiply (.conjugate M) M)) (ComplexDense. I)))))

;; ==================================================
  cpred/ComplexPredicate
;; ==================================================
  (hermitian? [this]
    ;; Not `.isHermitian()`: confirmed empirically to not check that
    ;; diagonal entries are real (see `gmat/hermitian?`'s docstring).
    (gmat/hermitian? this))

;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape [_] [(.countRows M) (.countColumns M)])
  (nrows [_] (.countRows M))
  (ncols [_] (.countColumns M))

  ;; -------- Retrieval --------
  (element [_ i j] (->C (.get M (long i) (long j))))

  (column [_ j] (ComplexDense. (.column M (long j))))
  (row [_ i] (ComplexDense. (.row M (long i))))
  (columns [_] (mapv #(ComplexDense. (.column M (long %))) (range (.countColumns M))))
  (rows [_] (mapv #(ComplexDense. (.row M (long %))) (range (.countRows M))))

  (->arrays [_]
    (let [rows (.countRows M) cols (.countColumns M)
          out  (make-array ComplexNumber rows cols)]
      (dotimes [i rows]
        (dotimes [j cols]
          (aset out i j (.get M i j))))
      out))

  (->vectors [_]
    (let [rows (.countRows M) cols (.countColumns M)]
      (mapv (fn [i]
              (mapv (fn [j]
                      (let [^ComplexNumber c (.get M i j)]
                        [(.getReal c) (.getImaginary c)]))
                    (range cols)))
            (range rows))))

;; ==================================================
  clojure.lang.Seqable
;; ==================================================
  (seq [_]
    (for [i (range (.countRows M)) j (range (.countColumns M))]
      (->C (.get M i j))))
;; ==================================================
  clojure.lang.Sequential
;; ==================================================
;; ==================================================
  Object
;; ==================================================
  (equals [_ B]
    (let [^MatrixC128 Bm (.-M ^ComplexDense B)]
      (and (= (.countRows M) (.countRows Bm))
           (= (.countColumns M) (.countColumns Bm))
           (every? (fn [[i j]]
                     (let [^ComplexNumber a (.get M i j)
                           ^ComplexNumber b (.get Bm i j)]
                       (and (< (Math/abs (- (.getReal a) (.getReal b))) ^double default/tolerance)
                            (< (Math/abs (- (.getImaginary a) (.getImaginary b))) ^double default/tolerance))))
                   (for [i (range (.countRows M)) j (range (.countColumns M))] [i j])))))
  (toString [_] (print! M)))

;; -------------------------------------------------------------------
;; Constructors
;; -------------------------------------------------------------------

(defn complexdense
  ;; zero
  (^ComplexDense [^long n ^long o] (ComplexDense. (.make MatrixC128/FACTORY n o)))
  ;; diagonal, from a flat collection of [re im] pairs
  (^ComplexDense [data]
   (let [n (count data)]
     (<-rows (for [i (range n)]
               (vec (for [j (range n)] (if (= i j) (nth data i) [0.0 0.0])))))))
  ;; elements, flat row-major collection of [re im] pairs
  (^ComplexDense [^long n ^long o data]
   (<-rows (partition o data))))

(defn <-rows ^ComplexDense [rows] (ComplexDense. (->matrix rows)))

(defn <-cols ^ComplexDense [cols]
  (ComplexDense. (.transpose (->matrix cols))))

(defn <-real ^ComplexDense [^fastmath.algebra.object.matrix.rectangular.real.ojalgo.RealDense M]
  (let [[rows cols] (d2/shape M)]
    (<-rows (for [i (range rows)]
              (for [j (range cols)]
                [(d2/element M i j) 0.0])))))
