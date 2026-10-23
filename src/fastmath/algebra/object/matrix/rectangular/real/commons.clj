(ns fastmath.algebra.object.matrix.rectangular.real.commons
  "Wires fastmath's own pre-existing real dense matrix implementation
  (`fastmath.matrix`, backed by Apache Commons Math's
  `RealMatrix`) into the new protocol-based API as a third backend,
  alongside EJML and ojAlgo - same protocol surface, same method shapes,
  just delegating to the already-existing `old/*` functions instead of
  reimplementing the arithmetic.
  "
  (:require
   [fastmath.algebra.function.matrix :as gmat]
   [fastmath.matrix :as old]
   [fastmath.default :as default]
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
   (org.apache.commons.math3.linear RealMatrix ArrayRealVector RealVector)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; ==================================================
;; Utilities
;; ==================================================

(declare <-rows)

(defn- format--real ^String [^double x ^long p ^double eps]
  (format (str "%." p "f") (if (< (Math/abs x) eps) 0.0 x)))

(defn ->str
  "Pretty string for RealMatrix with aligned columns.
   Options: :precision (3), :eps (1e-12), :max-rows (12), :max-cols (12)"
  ([^RealMatrix A] (->str A {}))
  ([^RealMatrix A {:keys [precision eps max-rows max-cols]
                   :or   {precision 3 eps default/tolerance max-rows 12 max-cols 12}}]
   (let [num-rows (old/nrow A) num-cols (old/ncol A)
         rlim (min num-rows ^long max-rows) clim (min num-cols ^long max-cols)
         cell (fn [^long i ^long j] (format--real (old/entry A i j) precision eps))
         mat (vec (for [i (range rlim)]
                    (vec (concat (for [j (range clim)] (cell i j))
                                 (when (< ^long clim num-cols) ["…"])))))
         w   (if (seq mat) (apply max (mapcat #(map count %) mat)) 0)
         pad (fn [s] (format (str "%-" w "s") s))
         row->s (fn [row] (str "[" (clojure.string/join " " (map pad row)) "]"))
         head (format ":shape [%d %d] :type %s" num-rows num-cols "float64")]
     (str head "\n"
          "["
          (clojure.string/join "\n" (map-indexed #(str (when (> ^long %1 0) " ") (row->s %2)) mat))
          (when (< ^long rlim num-rows) "\n …")
          "]"))))

(defn print!
  ([^RealMatrix A] (print! A {}))
  ([^RealMatrix A opts]
   (println (->str A opts))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 RealDense                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype RealDense [^RealMatrix M]
;; ==================================================
  rmat/RectangularMatrix
;; ==================================================
  asg/AdditiveSemigroup
  (add [_ other] (RealDense. (old/add M (.-M ^RealDense other))))

  am/AdditiveMonoid
  (zero [_] (RealDense. (old/zero (old/nrow M) (old/ncol M) true)))

  ag/AdditiveGroup
  (negate [_] (RealDense. (old/sub M)))

  module/Module
  (scale [_ scalar] (RealDense. (old/muls M (double scalar))))

  nspace/NormedSpace
  (norm [_] (old/norm M :frobenius))

;; ==================================================
  emat/MatrixExtra
;; ==================================================
  (diagonal [_]
    (<-rows (mapv vector (.toArray ^RealVector (old/diag M)))))

  (add--s [_ s] (RealDense. (old/adds M (double s))))

  (inner [this other] (gmat/inner this other))

  (kronecker [this other] (<-rows (gmat/kronecker-rows this other)))

  (multiply [_ other] (RealDense. (old/mulm M (.-M ^RealDense other))))
  (multiply--s [this s] (module/scale this s))
  (multiply--e [_ other] (RealDense. (old/emulm M (.-M ^RealDense other))))
  (multiply--v [_ v]
    (let [n (old/nrow M)
          ^doubles arr (double-array (for [i (range n)] (old/entry ^RealMatrix (.-M ^RealDense v) i 0)))
          ^RealVector result (old/mulv M (ArrayRealVector. arr))]
      (<-rows (mapv vector (.toArray result)))))
  (subtract [_ other] (RealDense. (old/sub M (.-M ^RealDense other))))
  (transpose [_] (RealDense. (old/transpose M)))

;; ==================================================
  p-sqmat/SquareMatrix
;; ==================================================
  (determinant [_] (old/det M))
  (trace [_] (old/trace M))

;; ==================================================
  solve/MatrixSolve
;; ==================================================
  (cholesky [_]
    (let [n (old/nrow M)]
      (when (not= n (old/ncol M))
        (throw (ex-info "Cholesky requires square matrix" {:shape [n (old/ncol M)]})))
      (let [L (try (old/cholesky M)
                   (catch Exception _
                     (throw (ex-info "Cholesky failed (matrix not SPD)" {:shape [n n]}))))]
        {:L (RealDense. L) :lower? true :spd? true})))

  (solve [_ v]
    (let [n (old/nrow M)
          ^doubles arr (double-array (for [i (range n)] (old/entry ^RealMatrix (.-M ^RealDense v) i 0)))
          ^RealVector result (old/solve M (ArrayRealVector. arr))]
      (<-rows (mapv vector (.toArray result)))))

;; ==================================================
  pred/MatrixPredicate
;; ==================================================
  (square? [this] (gmat/square? this))
  (normal? [this] (gmat/normal? this))
  (symmetric? [this] (gmat/symmetric? this))
  (singular? [this] (gmat/singular? this))
  (unitary? [this]
    ;; No dedicated orthogonality check in the old API - compare Qᵀ Q to
    ;; the identity, built locally (no cross-namespace dependency needed).
    (and (gmat/square? this)
         (let [n (old/nrow M)
               I (old/eye n true)]
           (= (RealDense. (old/mulm (old/transpose M) M)) (RealDense. I)))))

;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape [_] (old/shape M))
  (nrows [_] (old/nrow M))
  (ncols [_] (old/ncol M))

  ;; -------- Retrieval --------
  (element [_ i j] (old/entry M (long i) (long j)))

  (column [_ j] (<-rows (mapv vector (.toArray ^RealVector (old/col M (long j))))))
  (row [_ i] (<-rows [(vec (.toArray ^RealVector (old/row M (long i))))]))
  (columns [this] (mapv #(d2/column this %) (range (old/ncol M))))
  (rows [this] (mapv #(d2/row this %) (range (old/nrow M))))

  (->arrays [_] (old/mat->array2d M))

  (->vectors [_] (mapv vec (old/mat->array2d M)))

  (fmap [_ f] (RealDense. (old/fmap M f)))

;; ==================================================
  clojure.lang.Seqable
;; ==================================================
  (seq [_] (apply concat (old/mat->seq M)))
;; ==================================================
  clojure.lang.Sequential
;; ==================================================
;; ==================================================
  Object
;; ==================================================
  (equals [_ B]
    (let [^RealMatrix Bm (.-M ^RealDense B)
          rows (old/nrow M) cols (old/ncol M)]
      (and (= rows (old/nrow Bm))
           (= cols (old/ncol Bm))
           (every? (fn [[^long i ^long j]]
                     (< (Math/abs (- (old/entry M i j) (old/entry Bm i j)))
                        ^double default/tolerance))
                   (for [i (range rows) j (range cols)] [i j])))))
  (toString [_] (print! M)))

;; -------------------------------------------------------------------
;; Constructors
;; -------------------------------------------------------------------

(defn realdense
  ;; zero
  (^RealDense [^long n ^long o] (RealDense. (old/zero n o true)))
  ;; diagonal
  (^RealDense [^doubles data] (RealDense. (old/diagonal (vec data))))
  ;; elements, flat row-major
  (^RealDense [^long n ^long o ^doubles data]
   (<-rows (partition o data))))

(defn <-rows ^RealDense [rows]
  (RealDense. (old/rows->RealMatrix (mapv vec rows))))

(defn <-cols ^RealDense [cols]
  (RealDense. (old/cols->RealMatrix (mapv vec cols))))
