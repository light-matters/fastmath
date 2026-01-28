(ns fastmath.api.v2.algebra.general
  "WARNING: experimental & WIP!

  Towards a unified, flexible, entrypoint for linear algebra across different mathematical objects, i.e. for those who don't want to have to think about types.

  Such is the case, this namespace prioritizes convenience over 'the law of minimum surprise'. Therefore, multiplying matrices of differing shapes will try to perform something meaningful rather than failing, e.g. perform Kronecker products when normal matrix multiplication wouldn't otherwise work.
  "
  ;; Currently implemented using EJML. 
  ;; TODO: Check for fixed, square etc. types as well.
  ;; - is it necessary to distinguish vectors from matrices?
  ;; - Should I go all the way and make a protocol for real numbers? (avoid special cases)
  ;; - look into type hierarchies (`extend`...) etc. w.r.t. type hints. Currently, square and rectangular matrices count as completely different types and so hinting is limited.
  (:refer-clojure :exclude [type + - / * vector?])
  (:require
   [fastmath.api.v2.algebra.complex.number]
   [fastmath.core :as fm]
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.algebra.object.matrix.create :as mat]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;              Structure relationships
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn linear-shape?
  "Does this collection have the shape of a mathematical vector?"
  [coll]
  (if-not (some #{1} (d2/shape coll))
    false
    true))

(defn scalar? [x]
  (or (number? x)
      (C/? x)))
(defn matrix? [x]
  (rmat/? x))
(defn vector? [x]
  ((every-pred matrix?
               linear-shape?)
   x))
(defn real? [x]
  (if (scalar? x)
    (not (C/? x))
    (if (rmat/? x)
      (if (cmat/? x)
        (cmat/real? x)
        true)
      (ex-info "Not a number or matrix!" {}))))

(derive ::matrix ::type)
(derive ::vector ::type)
(derive ::scalar ::type)

(derive ::matrix--real ::matrix)
(derive ::matrix--real ::real)
(derive ::matrix--complex ::matrix)
(derive ::matrix--complex ::complex)

(derive ::vector ::matrix)
(derive ::vector--real ::vector)
(derive ::vector--real ::real)
(derive ::vector--complex ::vector)
(derive ::vector--complex ::complex)

(derive ::scalar--real ::scalar)
(derive ::scalar--real ::real)
(derive ::scalar--complex ::scalar)
(derive ::scalar--complex ::complex)

(defn type
  "Classify an argument so the arithmetic multimethods can dispatch on it."
  [x]
  (if (rmat/? x)
    (if (linear-shape? x)
      (if (cmat/? x)
        ::vector--complex
        ::vector--real)
      (if (cmat/? x)
        ::matrix--complex
        ::matrix--real))
    (cond
      (number? x) ::scalar--real
      (C/? x) ::scalar--complex)))

(defn- rank--domain [x]
  (cond
    (isa? x ::real) 0
    (isa? x ::complex) 1
    :else (throw (ex-info "No recognised domain for " {:value x}))))

(defn- rank--shape [x]
  (cond
    (isa? x ::scalar) 0
    (isa? x ::vector) 1
    (isa? x ::matrix) 2
    :else (throw (ex-info "No recognised algebraic type for " {:value x}))))

(defn- promote-domain
  "Increases the set type of the element (if possible), e.g. a real number becomes a complex number or a real matrix becomes a complex matrix."
  [x]
  (case (type x)
    ::scalar--real (cc/<-real (double x))
    ::scalar--complex x
    ::matrix--real (mat/<-real x)
    ::matrix--complex x

    (throw (ex-info "No promotion rule for " {:value x}))))

(defn- ensure-domain-match
  "Takes a pair of arguments and promotes arguments where necessary to ensure compatible domains.

   e.g. when adding a real number to a complex number, the `double` becomes a complex number with 0 for the 'imaginary' part."
  [a b]
  (let [ta (type a)
        tb (type b)
        ra (rank--domain ta)
        rb (rank--domain tb)]
    (cond
      (= ta tb) [a b]
      (= ra rb) [a b]
      (< ra rb) [(promote-domain a) b]
      :else [a (promote-domain b)])))

(defn- same-shape? [m1 m2]
  (->> [m1 m2]
       (map d2/shape)
       (apply =)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ; Core Aliases  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def norm rmat/norm)
(def zero rmat/zero)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                   Addition                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defmulti add*
  "Addition that understands matrices, vectors, and scalars."
  ;; TODO: Check for shape in multiplication dispatch?
  (fn [a b]
    [(type a) (type b)]))

(defmethod add* [::matrix ::matrix] [m1 m2]
  (when (not (same-shape? m1 m2))
    (throw (ex-info "Shape mismatch!" {:m1 m1 :m2 m2})))
  (apply cmat/add (ensure-domain-match m1 m2)))

(defmethod add* [::matrix ::scalar] [a s]
  (apply emat/add--s (ensure-domain-match a s)))
(defmethod add* [::scalar ::matrix] [s a]
  (apply emat/add--s (ensure-domain-match a s)))

(defmethod add* [::scalar--real ::scalar--real] [a b]
  (fm/+ a b))
(defmethod add* [::scalar--complex ::scalar--complex] [a b]
  (apply C/add (ensure-domain-match a b)))
(defmethod add* [::scalar--complex ::scalar--real] [a b]
  (apply C/add (ensure-domain-match a b)))
(defmethod add* [::scalar--real ::scalar--complex] [a b]
  (apply C/add (ensure-domain-match a b)))

(defn +
  "Variadic entry point that reduces via the multimethod."
  ([x] x)
  ([x y] (add* x y))
  ([x y & more]
   (reduce add* (add* x y) more)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                   Subtraction                               ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmulti subtract*
  "Subtraction that understands matrices, vectors, and scalars."
  (fn
    ([a] [(type a)])
    ([a b]
     [(type a) (type b)])))

(defmethod subtract* [::matrix] [m]
  (cmat/negate m))
(defmethod subtract* [::scalar--complex] [s]
  (cmat/negate s))
(defmethod subtract* [::scalar--real] [s]
  (fm/- s))

(defmethod subtract* [::matrix ::matrix] [m1 m2]
  (when (not (same-shape? m1 m2))
    (throw (ex-info "Shape mismatch!" {:m1 m1 :m2 m2})))
  (apply emat/subtract (ensure-domain-match m1 m2)))

(defmethod subtract* [::matrix ::scalar] [m s]
  (apply emat/add--s (ensure-domain-match m (subtract* s))))
(defmethod subtract* [::scalar ::matrix] [s m]
  (apply emat/add--s (ensure-domain-match (cmat/negate m) s)))

(defmethod subtract* [::scalar--real ::scalar--real] [a b]
  (fm/- a b))
(defmethod subtract* [::scalar--complex ::scalar--complex] [a b]
  (apply C/subtract (ensure-domain-match a b)))
(defmethod subtract* [::scalar--complex ::scalar--real] [a b]
  (apply C/subtract (ensure-domain-match a b)))
(defmethod subtract* [::scalar--real ::scalar--complex] [a b]
  (apply C/subtract (ensure-domain-match a b)))

(defn -
  "Variadic entry point that reduces via the multimethod."
  ([x] (subtract* x))
  ([x y] (subtract* x y))
  ([x y & more]
   (reduce subtract* (subtract* x y) more)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;            Multiplication           ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- compatible-shapes?
  "Checks that the given matrices have compatible shapes for 'matrix multiplication'."
  ;; TODO: Make this check early in the variadic version.
  [m1 m2]
  (= (second (d2/shape m1)) (first (d2/shape m2))))

(defmulti multiply*
  "Multiplication that understands matrices, vectors, and scalars."
  (fn
    ([a b]
     [(type a) (type b)])))

(defmethod multiply* [::matrix ::matrix] [m1 m2]
  (when (not (compatible-shapes? m1 m2))
    (throw (ex-info "Shape mismatch!" {:m1 m1 :m2 m2})))
  (apply emat/multiply (ensure-domain-match m1 m2)))

(defmethod multiply* [::matrix ::scalar] [m s]
  (apply rmat/scale (ensure-domain-match m s)))
(defmethod multiply* [::scalar ::matrix] [s m]
  (apply rmat/scale (ensure-domain-match m s)))

(defmethod multiply* [::scalar--real ::scalar--real] [a b]
  (fm/* a b))
(defmethod multiply* [::scalar--complex ::scalar--complex] [a b]
  (apply C/multiply (ensure-domain-match a b)))
(defmethod multiply* [::scalar--complex ::scalar--real] [a b]
  (apply C/multiply (ensure-domain-match a b)))
(defmethod multiply* [::scalar--real ::scalar--complex] [a b]
  (apply C/multiply (ensure-domain-match a b)))

(defn *
  "Variadic entry point that reduces via the multimethod."
  ([x] x)
  ([x y] (multiply* x y))
  ([x y & more]
   (reduce multiply* (multiply* x y) more)))

