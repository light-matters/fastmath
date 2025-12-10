(ns fastmath.algebra.object.matrix.rectangular.real.ejml
  "Implementing real matrices using EJML as a backend [https://github.com/lessthanoptimal/ejml].
  "
  (:require
   [clojure.string :as str]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as realdense]
   [fastmath.default :as default]
   [fastmath.protocol.algebra.object.matrix.general :as gmat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.structure.additive.group :as ag]
   [fastmath.protocol.algebra.structure.additive.monoid :as am]
   [fastmath.protocol.algebra.structure.additive.semigroup :as asg]
   [fastmath.protocol.algebra.structure.module :as module]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]
   [fastmath.protocol.representation.d2 :as d2])
  (:import
   (java.lang Math)
   (org.ejml.data DMatrixRMaj)
   (org.ejml.dense.row CommonOps_DDRM NormOps_DDRM)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(println "=== start ===")
;; ^^^ Used for debugging

;; ==================================================
;; Utilities 
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

(defn  format--real ^String [^double x ^long p ^double eps]
  (format (str "%." p "f") (if (< (Math/abs x) eps) 0.0 x)))

(defn- realZM?
  ([^DMatrixRMaj M]
   (realZM? M default/tolerance))
  ([^DMatrixRMaj M tolerance]
   (let [^"[D" d (.-data M)
         n (alength d)]
     (loop [i 1]
       (cond
         (>= i n) true
         (> (Math/abs (aget d i)) (double tolerance)) false
         :else (recur (+ i 2)))))))

(defn ->str
  "Pretty string for DMatrixRMaj with aligned columns.
   Options: :precision (3), :eps (1e-12), :max-rows (12), :max-cols (12)"

  ([^DMatrixRMaj A] (->str A {}))
  ([^DMatrixRMaj A {:keys [precision eps max-rows max-cols]
                    :or   {precision 3 eps default/tolerance max-rows 12 max-cols 12}}]
   (let [num-rows (.numRows A) num-cols (.numCols A)
         rlim (min num-rows ^long max-rows) clim (min num-cols ^long max-cols)
         real? (realZM? A eps)
         cell (fn [i j]
                (format--real (.get A i j) precision eps))
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
  ([^DMatrixRMaj A] (print! A {}))
  ([^DMatrixRMaj A opts]
   (println (->str A opts))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 RealDense                                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype RealDense [^DMatrixRMaj M]
;; ==================================================
  rmat/RectangularMatrix
;; ==================================================
  asg/AdditiveSemigroup
  (add [_ other]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/add M B out)
      (RealDense. out)))

  am/AdditiveMonoid
  (zero [_] (.zero (DMatrixRMaj. (.numRows M) (.numCols M))))
  ;; TODO: Check if `.zero` is necessary

  ag/AdditiveGroup
  (negate [_]
    (let [A (.copy M)]
      (CommonOps_DDRM/scale -1.0 0.0 A)
      (RealDense. A)))

  nspace/NormedSpace
  (norm [_]
    ;; frobenius matrix norm
    (double (NormOps_DDRM/normF M)))

  module/Module
  (scale [_ [r i]]
    (let [A (.copy M)]
      (CommonOps_DDRM/scale (double r) (double i) A)
      (RealDense. A)))

;; ==================================================
  gmat/GeneralMatrix
;; ==================================================
  (add--s [_ s]
    (let [A (.copy M) d (.data A) s (double s)]
      (dotimes [k (alength d)]
        (aset-double d k (+ (aget d k) s)))
      (RealDense. A)))

  (inner [_ other]
    (let [^doubles da (.copy M)
          ^doubles db (.M ^RealDense other)
          n (.getNumElements M)]
      (loop [i 0
             acc 0.0]
        (if (< i n)
          (recur (inc i)
                 (+ acc (* (aget da i) (aget db i))))
          acc))))

  (outer [_ other]
    (let [m (.getNumRows ^DMatrixRMaj M)
          n (.getNumRows ^DMatrixRMaj (.M ^RealDense other))
          ^doubles ad (.getData ^DMatrixRMaj M)
          ^doubles bd (.getData ^DMatrixRMaj (.M ^RealDense other))
          out (DMatrixRMaj. m n)]
      (dotimes [i m]
        (let [ai (aget ad i)]
          (dotimes [j n]
            (.set out i j (* ai (aget bd j))))))
      (RealDense. out)))

  (map--m [_ f]
    ;; TODO: Prevent function from multiple calls to set
    ;; - Work on (.data A) directly
    (let [^DMatrixRMaj A (.copy M)]
      (dotimes [i (.numRows A)]
        (dotimes [j (.numCols A)]
          (.set A i j (f (.get A i j)))))
      (RealDense. A)))
  (multiply [_ other]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols B))]
      (CommonOps_DDRM/mult M B out)
      (RealDense. out)))
  (multiply--e [_ other]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols B))]
      (CommonOps_DDRM/elementMult M B out)
      (RealDense. out)))
  (multiply--v [_ v]
    ;; v is expected to be (n×1) RealDense
    (let [^DMatrixRMaj x (.M ^RealDense v)
          out (DMatrixRMaj. (.numRows M) 1)]
      (CommonOps_DDRM/mult M x out)
      (RealDense. out)))
  (subtract [_ other]
    (let [^DMatrixRMaj B (.M ^RealDense other)
          out (DMatrixRMaj. (.numRows M) (.numCols M))]
      (CommonOps_DDRM/subtract M B out)
      (RealDense. out)))
  (transpose [_]
    (let [T (DMatrixRMaj. (.numCols M) (.numRows M))]
      (CommonOps_DDRM/transpose M T)
      (RealDense. T)))

  (square? [_]
    (= (.numRows M) (.numCols M)))

;; ==================================================
  d2/D2
;; ==================================================
  ;; -------- Info --------
  (shape    [_] [^long (.numRows M) ^long (.numCols M)])
  (num-rows [_] ^long (.numRows M))
  (num-cols [_] ^long (.numCols M))

  ;; -------- Retrieval --------
  (element [_ i j]
    (.get M (long i) (long j)))

  (column [_ j]
    (RealDense. (extract M [0 (.numRows M)] [j (+ 1 (long j))])))

  (row [_ i]
    (RealDense. (extract M [i (+ 1 (long i))] [0 (.numCols M)])))
  (columns [_]
    (mapv (fn [^long ic]
            (RealDense. (extract M [0 (.numRows M)] [ic (+ 1 ic)])))
          (range (.numCols M))))

  (rows [_]
    (mapv (fn [^long ir]
            (RealDense. (extract M [ir (+ ir 1)] [0 (.numCols M)])))
          (range (.numRows M))))

  (diagonal [_]
    (let [n (min (.numRows M) (.numCols M))
          out (DMatrixRMaj. n 1)]
      (CommonOps_DDRM/extractDiag M out)
      (RealDense. out)))

  (array<- [_]
    (let [r (.numRows M) c (.numCols M)
          out (make-array Double/TYPE r c)]
      (dotimes [i r]
        (dotimes [j c]
          (aset-double out i j (.get M i j))))
      out))

;; ==================================================
  Object
;; ==================================================
  (toString [_]
    (print! M)))

;; -------------------------------------------------------------------
;; Constructors 
;; -------------------------------------------------------------------

(defn realdense
  ;; I
  ;; (^RealDense [^long n]
  ;;  (->RealDense (CommonOps_DDRM/identity n)))

  ;; zero
  (^RealDense [^long n ^long o]
   (->RealDense (DMatrixRMaj. n o)))

  ;; diagonal
  (^RealDense [^doubles data]
   (-> (CommonOps_DDRM/diag (double-array data))
       ->RealDense))

;; elements
  (^RealDense [^long n ^long o ^doubles data]
   (->RealDense (DMatrixRMaj. n o false data))))

(defn <-rows ^RealDense [rows]
  (let [nrows (count rows)
        ncols (count (first rows))
        data  (double-array (apply concat rows))]
    (->RealDense (DMatrixRMaj. nrows ncols true data))))

(defn <-cols ^RealDense [cols]
  (let [^DMatrixRMaj A (.M (<-rows cols))
        ^DMatrixRMaj out (DMatrixRMaj. (.numCols A) (.numRows A))]
    (CommonOps_DDRM/transpose A out)

    (->RealDense out)))
