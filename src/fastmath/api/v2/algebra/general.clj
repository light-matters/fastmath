(ns fastmath.api.v2.algebra.general
  "WARNING: experimental & WIP!

  Towards a unified, flexible, entrypoint for linear algebra across different types, i.e. for those who don't want to have to think about types.

  "
  ;; Currently implemented using EJML. 
  ;; TODO: Check for fixed, square etc. types as well.
  ;; - is it necessary to distinguish vectors from matrices?
  (:refer-clojure :exclude [type])
  (:require
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.protocol.algebra.object.number.complex :as C]
   [fastmath.protocol.algebra.object.matrix.rectangular.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.algebra.object.matrix.create :as createm]
   [fastmath.interpolation.linear :as linear]))

(defn linear-shape
  "Boolean or nil."
  [coll]
  (some #{1} (d2/shape coll)))

;; Type predicates
(defn scalar? [x]
  (or number?
      C/?))
(defn matrix? [x]
  ((some-fn cmat/? rmat/?) x))
(defn vector? [x]
  ((every-pred matrix?
               linear-shape)
   x))

(defn real? [x])

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

(defn- type
  "Classify an argument so the arithmetic multimethods can dispatch on it."
  [x]
  (if (rmat/? x)
    (if (linear-shape x)
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
    :scalar--real (createcn/<-real (double x))
    :scalar--complex x
    :matrix--real (createm/<-real x)
    :matrix--complex x

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
      (< ra rb) [(promote-domain a) b]
      :else [a (promote-domain b)])))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                   Addition                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defmulti add*
  "Addition that understands matrices, vectors, and scalars."
  (fn [a b]
    (let [[pa pb] (ensure-domain-match a b)]
      [(type pa) (type pb)])))

(defmethod add* [::matrix ::matrix] [a b]
  (cmat/add a b))
(defmethod add* [::matrix ::scalar] [a s]
  (cmat/add--s a s))
(defmethod add* [::scalar ::matrix] [s a]
  (cmat/add--s a s))
(defmethod add* [::scalar ::scalar] [a b]
  (clojure.core/+ a b))

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
;; WIP:
(defmulti subtract*
  "Subtraction that understands matrices, vectors, and scalars."
  (fn [a b]
    (let [[pa pb] (ensure-domain-match a b)]
      [(type pa) (type pb)])))

(defmethod subtract* [::matrix ::matrix] [a b]
  (cmat/sub a b))

(defmethod subtract* [::matrix ::scalar] [a b]
  (cmat/add--s a b))
;; (defmethod subtract* [::scalar ::matrix] [a b]
;;   (cmat/add--s a b))
;;   TODO: Does this ^ make sense?
;;   
(defmethod subtract* [::scalar ::scalar] [a b]
  (clojure.core/- a b))

(defn -
  "Variadic entry point that reduces via the multimethod."
  ([x] x)
  ([x y] (subtract* x y))
  ([x y & more]
   (reduce subtract* (subtract* x y) more)))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;








