(ns fastmath.api.v2.algebra.predicate
  "For all of your (yes-no) algebraic questions."
  (:refer-clojure :exclude [vector?])
  (:require
   [fastmath.algebra.object.type :as type]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.object.number.complex :as pC]
   [fastmath.protocol.representation.d2 :as d2]))

(defn linear-shape?
  "Does this collection have the shape of a mathematical vector?"
  [coll]
  (if-not (some #{1} (d2/shape coll))
    false
    true))

(defn scalar? [x]
  (or (number? x)
      (pC/? x)))

(defn complex-number? [x]
  (pC/? x))

(defn matrix? [x]
  (rmat/? x))
(defn vector? [x]
  ((every-pred matrix?
               linear-shape?)
   x))
(defn real? [x]
;; TODO: Make this more elegant - using general protocols
  (if (scalar? x)
    (not (pC/? x))
    (if (rmat/? x)
      (if (cmat/? x)
        (cmat/real? x)
        true)
      (ex-info "Not a number or matrix!" {}))))

(defn same-shape? [m1 m2]
  (if (and (matrix? m1) (matrix? m2))
    (->> [m1 m2]

         (map d2/shape)
         (apply =))
    false))

(def type type/?)
