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
  ;; - replace linear shape with `vector `protocol?
  ;; - better word than 'general'?
  ;; -- all, algebra 
  (:refer-clojure :exclude [type + - / * vector?])
  (:require
   [fastmath.api.v2.algebra.complex :as C]
   [fastmath.api.v2.algebra.predicate :as pred]
   [fastmath.algebra.plumb :as plumb]
   [fastmath.core :as fm]
   [fastmath.protocol.algebra.function.field :as pf]
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.algebra.object.matrix.create :as mat]))

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
  (when (not (pred/same-shape? m1 m2))
    (throw (ex-info "Shape mismatch!" {:m1 m1 :m2 m2})))
  (apply cmat/add (plumb/ensure-domain-match m1 m2)))

(defmethod add* [::matrix ::scalar] [a s]
  (apply emat/add--s (plumb/ensure-domain-match a s)))
(defmethod add* [::scalar ::matrix] [s a]
  (apply emat/add--s (plumb/ensure-domain-match a s)))

(defmethod add* [::scalar--real ::scalar--real] [a b]
  (fm/+ a b))
(defmethod add* [::scalar--complex ::scalar--complex] [a b]
  (apply C/add (plumb/ensure-domain-match a b)))
(defmethod add* [::scalar--complex ::scalar--real] [a b]
  (apply C/add (plumb/ensure-domain-match a b)))
(defmethod add* [::scalar--real ::scalar--complex] [a b]
  (apply C/add (plumb/ensure-domain-match a b)))

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
  (when (not (pred/same-shape? m1 m2))
    (throw (ex-info "Shape mismatch!" {:m1 m1 :m2 m2})))
  (apply emat/subtract (plumb/ensure-domain-match m1 m2)))

(defmethod subtract* [::matrix ::scalar] [m s]
  (apply emat/add--s (plumb/ensure-domain-match m (subtract* s))))
(defmethod subtract* [::scalar ::matrix] [s m]
  (apply emat/add--s (plumb/ensure-domain-match (cmat/negate m) s)))

(defmethod subtract* [::scalar--real ::scalar--real] [a b]
  (fm/- a b))
(defmethod subtract* [::scalar--complex ::scalar--complex] [a b]
  (apply pf/subtract (plumb/ensure-domain-match a b)))
(defmethod subtract* [::scalar--complex ::scalar--real] [a b]
  (apply pf/subtract (plumb/ensure-domain-match a b)))
(defmethod subtract* [::scalar--real ::scalar--complex] [a b]
  (apply pf/subtract (plumb/ensure-domain-match a b)))

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
  (apply emat/multiply (plumb/ensure-domain-match m1 m2)))

(defmethod multiply* [::matrix ::scalar] [m s]
  (apply rmat/scale (plumb/ensure-domain-match m s)))
(defmethod multiply* [::scalar ::matrix] [s m]
  (apply rmat/scale (plumb/ensure-domain-match m s)))

(defmethod multiply* [::scalar--real ::scalar--real] [a b]
  (fm/* a b))
(defmethod multiply* [::scalar--complex ::scalar--complex] [a b]
  (apply C/multiply (plumb/ensure-domain-match a b)))
(defmethod multiply* [::scalar--complex ::scalar--real] [a b]
  (apply C/multiply (plumb/ensure-domain-match a b)))
(defmethod multiply* [::scalar--real ::scalar--complex] [a b]
  (apply C/multiply (plumb/ensure-domain-match a b)))

(defn *
  "Variadic entry point that reduces via the multimethod."
  ([x] x)
  ([x y] (multiply* x y))
  ([x y & more]
   (reduce multiply* (multiply* x y) more)))

