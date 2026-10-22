(ns fastmath.algebra.object.matrix.rectangular.real.ojalgo
  "Implementing real matrices using ojAlgo as a backend
  [https://github.com/optimatika/ojAlgo].

  Uses ojAlgo's high-level `MatrixR064`/`BasicMatrix` API (immutable,
  functional) rather than the low-level `PhysicalStore`/`GenericStore`
  layer - the latter needs reflection to work around access restrictions
  (see the abandoned `fastmath.matrix.dense.scratch` experiments) and
  buys nothing here.
  "
  (:require
   [clojure.string :as str]
   [fastmath.default :as default]
   [fastmath.algebra.function.matrix :as gmat]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.object.matrix.solve :as solve]
   [fastmath.protocol.algebra.object.matrix.square :as p-sqmat]
   [fastmath.protocol.algebra.object.matrix.predicate :as pred]
   [fastmath.protocol.algebra.structure.additive.group :as ag]
   [fastmath.protocol.algebra.structure.additive.monoid :as am]
   [fastmath.protocol.algebra.structure.additive.semigroup :as asg]
   [fastmath.protocol.algebra.structure.module :as module]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]
   [fastmath.protocol.representation.d2 :as d2])
  (:import
   (java.lang Math)
   (org.ojalgo.matrix MatrixR064)
   (org.ojalgo.matrix.decomposition Cholesky)
   (org.ojalgo.structure Access2D)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; ==================================================
;; Utilities
;; ==================================================

(defn- ^MatrixR064 ->matrix
  "Wrap plain 2D data (a `double[][]`, or anything else `Access2D/wrap`
  accepts) as a `MatrixR064`."
  [data]
  (.makeWrapper MatrixR064/FACTORY (Access2D/wrap ^"[[D" data)))

(defn- rows->array ^"[[D" [rows]
  (into-array (map double-array rows)))

(declare <-rows)

(defn- format--real ^String [^double x ^long p ^double eps]
  (format (str "%." p "f") (if (< (Math/abs x) eps) 0.0 x)))

(defn ->str
  "Pretty string for MatrixR064 with aligned columns.
   Options: :precision (3), :eps (1e-12), :max-rows (12), :max-cols (12)"
  ([^MatrixR064 A] (->str A {}))
  ([^MatrixR064 A {:keys [precision eps max-rows max-cols]
                   :or   {precision 3 eps default/tolerance max-rows 12 max-cols 12}}]
   (let [num-rows (.countRows A) num-cols (.countColumns A)
         rlim (min num-rows ^long max-rows) clim (min num-cols ^long max-cols)
         cell (fn [^long i ^long j] (format--real (.doubleValue A i j) precision eps))
         mat (vec (for [i (range rlim)]
                    (vec (concat (for [j (range clim)] (cell i j))
                                 (when (< ^long clim num-cols) ["…"])))))
         w   (if (seq mat) (apply max (mapcat #(map count %) mat)) 0)
         pad (fn [s] (format (str "%-" w "s") s))
         row->s (fn [row] (str "[" (str/join " " (map pad row)) "]"))
         head (format ":shape [%d %d] :type %s" num-rows num-cols "float64")]
     (str head "\n"
          "["
          (str/join "\n" (map-indexed #(str (when (> ^long %1 0) " ") (row->s %2)) mat))
          (when (< ^long rlim num-rows) "\n …")
          "]"))))

(defn print!
  ([^MatrixR064 A] (print! A {}))
  ([^MatrixR064 A opts]
   (println (->str A opts))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 RealDense                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype RealDense [^MatrixR064 M]
;; ==================================================
  rmat/RectangularMatrix
;; ==================================================
  asg/AdditiveSemigroup
  (add [_ other] (RealDense. (.add M ^MatrixR064 (.-M ^RealDense other))))

  am/AdditiveMonoid
  (zero [_] (RealDense. (.make MatrixR064/FACTORY (.countRows M) (.countColumns M))))

  ag/AdditiveGroup
  (negate [_] (RealDense. (.negate M)))

  module/Module
  (scale [_ scalar] (RealDense. (.multiply M (double scalar))))

  nspace/NormedSpace
  (norm [_] (.norm M))

;; ==================================================
  emat/MatrixExtra
;; ==================================================
  (diagonal [_]
    (let [n (min (.countRows M) (.countColumns M))]
      (<-rows (for [^long i (range n)] [(.doubleValue M (int i) (int i))]))))

  (add--s [_ s] (RealDense. (.add M (double s))))

  (inner [this other] (gmat/inner this other))

  (kronecker [this other] (<-rows (gmat/kronecker-rows this other)))

  (multiply [_ other] (RealDense. (.multiply M ^MatrixR064 (.-M ^RealDense other))))
  (multiply--s [this s] (module/scale this s))
  (multiply--e [_ other]
    (let [^MatrixR064 B (.-M ^RealDense other)
          rows (.countRows M) cols (.countColumns M)]
      (<-rows (for [^long i (range rows)]
                (for [^long j (range cols)]
                  (* (.doubleValue M (int i) (int j)) (.doubleValue B (int i) (int j))))))))
  (multiply--v [_ v] (RealDense. (.multiply M ^MatrixR064 (.-M ^RealDense v))))
  (subtract [_ other] (RealDense. (.subtract M ^MatrixR064 (.-M ^RealDense other))))
  (transpose [_] (RealDense. (.transpose M)))

;; ==================================================
  p-sqmat/SquareMatrix
;; ==================================================
  (determinant [_] (.getDeterminant M))
  (trace [_] (.getTrace M))

;; ==================================================
  solve/MatrixSolve
;; ==================================================
  (cholesky [_]
    (let [n (.countRows M)
          _ (when (not= n (.countColumns M))
              (throw (ex-info "Cholesky requires square matrix" {:shape [n (.countColumns M)]})))
          ^Cholesky chol (.decompose Cholesky/R064 M)]
      (when-not (.isSPD chol)
        (throw (ex-info "Cholesky failed (matrix not SPD)" {:shape [n n]})))
      {:L (RealDense. (.makeWrapper MatrixR064/FACTORY (.getL chol)))
       :lower? true :spd? true}))

  (solve [_ v] (RealDense. (.solve M ^MatrixR064 (.-M ^RealDense v))))

;; ==================================================
  pred/MatrixPredicate
;; ==================================================
  (square? [this] (gmat/square? this))
  (normal? [this] (gmat/normal? this))
  (symmetric? [this] (gmat/symmetric? this))
  (singular? [this] (gmat/singular? this))
  (unitary? [this]
    ;; No dedicated primitive on BasicMatrix - compare Qᵀ Q to the
    ;; identity, built locally (no cross-namespace dependency needed).
    (and (gmat/square? this)
         (let [n (.countRows M)
               ^MatrixR064 I (.makeEye MatrixR064/FACTORY n n)
               ^MatrixR064 Mt (.transpose M)]
           (= (RealDense. (.multiply Mt M)) (RealDense. I)))))

;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape [_] [(.countRows M) (.countColumns M)])
  (nrows [_] (.countRows M))
  (ncols [_] (.countColumns M))

  ;; -------- Retrieval --------
  (element [_ i j] (.doubleValue M (int i) (int j)))

  (column [_ j] (RealDense. (.column M (long j))))
  (row [_ i] (RealDense. (.row M (long i))))
  (columns [_] (mapv #(RealDense. (.column M (long %))) (range (.countColumns M))))
  (rows [_] (mapv #(RealDense. (.row M (long %))) (range (.countRows M))))

  (->arrays [_]
    (let [r (.countRows M) c (.countColumns M)
          out (make-array Double/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-double out i j (.doubleValue M (int i) (int j)))))
      out))

  (->vectors [_]
    (let [r (.countRows M) c (.countColumns M)]
      (mapv (fn [^long i] (mapv (fn [^long j] (.doubleValue M (int i) (int j))) (range c))) (range r))))

  (fmap [_ f]
    (let [r (.countRows M) c (.countColumns M)]
      (<-rows (for [^long i (range r)]
                (for [^long j (range c)] (f (.doubleValue M (int i) (int j))))))))

;; ==================================================
  clojure.lang.Seqable
;; ==================================================
  (seq [_]
    (let [r (.countRows M) c (.countColumns M)]
      (for [^long i (range r) ^long j (range c)] (.doubleValue M (int i) (int j)))))
;; ==================================================
  clojure.lang.Sequential
;; ==================================================
;; ==================================================
  Object
;; ==================================================
  (equals [_ B]
    (let [^MatrixR064 Bm (.-M ^RealDense B)
          rows (.countRows M) cols (.countColumns M)]
      (and (= rows (.countRows Bm))
           (= cols (.countColumns Bm))
           (every? (fn [[^long i ^long j]]
                     (< (Math/abs (- (.doubleValue M (int i) (int j)) (.doubleValue Bm (int i) (int j))))
                        ^double default/tolerance))
                   (for [i (range rows) j (range cols)] [i j])))))
  (toString [_] (print! M)))

;; -------------------------------------------------------------------
;; Constructors
;; -------------------------------------------------------------------

(defn realdense
  ;; zero
  (^RealDense [^long n ^long o] (RealDense. (.make MatrixR064/FACTORY n o)))
  ;; diagonal
  (^RealDense [^doubles data]
   (RealDense. (.makeDiagonal MatrixR064/FACTORY (double-array data))))
  ;; elements, flat row-major
  (^RealDense [^long n ^long o ^doubles data]
   (<-rows (partition o data))))

(defn <-rows ^RealDense [rows]
  (->> rows (map vec) rows->array ->matrix RealDense.))

(defn <-cols ^RealDense [cols]
  (RealDense. (.transpose (->matrix (rows->array (map vec cols))))))
