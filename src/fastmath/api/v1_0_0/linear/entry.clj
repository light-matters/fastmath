(ns fastmath.api.v1-0-0.linear.entry
  "WARNING: experimental & WIP!

  Towards a unified, flexible, entrypoint for linear algebra across different types, i.e. for those who don't want to have to think about types.
  "
  (:require [fastmath.protocols.linear :as mat]))

;; Type predicates
(defn vector? [x])
(defn matrix? [x])
(defn complex-number? [x])

(derive :matrix :linear-algebra)
(derive :vector :linear-algebra)
(derive :scalar :linear-algebra)
(derive :vector :matrix)
(derive :matrix--real :matrix)
(derive :matrix--complex :matrix)
(derive :vector--real :vector)
(derive :vector--complex :vector)
(derive :scalar--real :scalar)
(derive :scalar--complex :scalar)

(def rank
  ;; "For ordering operations according to complexity. Amongst other things, helps to reduce code duplicaton."
  ;; TODO: Only sort scalar-involved operations. vector * matrix != matrix*vector!
  (into {}
        (map-indexed (fn [i k] [k i])
                     [:scalar--real
                      :scalar--complex
                      :vector--real
                      :vector--complex
                      :matrix--real
                      :matrix--complex])))

(defn tag
  "Classify an argument so the arithmetic multimethods can dispatch on it."
  [x]
  (cond
    (matrix? x) (if (mat/real? x)
                  :matrix--real
                  :matrix--complex)
    (vector? x) (if (mat/real? x)
                  :vector--real
                  :vector--complex)
    (number? x) :scalar--real
    (complex-number? x) :scalar--complex))

(defmulti add*
  "Addition that understands matrices, vectors, and scalars."
  (fn [a b]
    (->> [(tag a) (tag b)]
         (sort-by rank >)
         vec)))

(defmethod add* [:matrix :matrix] [a b]
  (m/add a b))

(defmethod add* [:matrix :scalar] [a s]
  (m/add--s a s))

(defmethod add* [:scalar :matrix] [s a]
  (m/add--s a s))

(defmethod add* [:vector :vector] [a b]
  (prot/add a b))

(defmethod add* [:scalar--complex :scalar--complex] [a b]
  (clojure.core/+ a b))
(defmethod add* [:scalar--real :scalar--real] [a b]
  (clojure.core/+ a b))

(defmethod add* [:scalar--real :scalar--complex] [a b]
  (clojure.core/+ (->complex a) b))
(defmethod add* [:scalar--complex :scalar--real] [a b]
  (clojure.core/+ a (->complex b)))

(defn +
  "Variadic entry point that reduces via the multimethod."
  ([x] x)
  ([x y] (add* x y))
  ([x y & more]
   (reduce add* (add* x y) more)))



