(ns fastmath.api.v1.algebra.linear
  "WARNING: experimental & WIP!

  Towards a unified, flexible, entrypoint for linear algebra across different types, i.e. for those who don't want to have to think about types.
  "
  (:refer-clojure :exclude [type])
  (:require

   [fastmath.protocols.complex.number :as cn]
   [fastmath.protocols.linear-algebra.complex.matrix :as cm]
   [fastmath.matrix.dense.complex.ejml :as cmat]
   ;; TODO:
   ;; - API change to make complex more natural
   ;; - put `i` in a more central, reasonable, place
   ))

;; Type predicates
(defn vector? [x])
(defn matrix? [x])
(defn complex-number? [x])

(derive ::matrix ::type)
(derive ::vector ::type)
(derive ::scalar ::type)
(derive ::vector ::matrix)
(derive ::matrix--real ::matrix)
(derive ::matrix--real ::real)
(derive ::matrix--complex ::matrix)
(derive ::matrix--complex ::complex)
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
  (cond
    (matrix? x) (if (cm/real? x) ;; TODO: Change this to check for type!
                  ::matrix--real
                  ::matrix--complex)
    (vector? x) (if (cm/real? x);; TODO: Change this to check for type!
                  ::vector--real
                  ::vector--complex)
    (number? x) ::scalar--real
    (complex-number? x) ::scalar--complex))

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
    :scalar--real (cmat/i (double x))
    :scalar--complex x
    :matrix--real (cmat/<-real x)
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
  (cm/add a b))
(defmethod add* [::matrix ::scalar] [a s]
  (cm/add--s a s))
(defmethod add* [::scalar ::matrix] [s a]
  (cm/add--s a s))
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
  (cm/sub a b))

(defmethod subtract* [::matrix ::scalar] [a b]
  (cm/add--s a b))
;; (defmethod subtract* [::scalar ::matrix] [a b]
;;   (cm/add--s a b))
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








