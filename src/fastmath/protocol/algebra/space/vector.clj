(ns fastmath.protocol.algebra.space.vector
  (:require
   ;; [fastmath.protocol.algebra.field :as field]
   ;;          [fastmath.protocol.algebra.module :as module]
   [fastmath.protocol.algebra.number.real :as real]))

(defprotocol VectorSpace)

(def add field/add)
(def multiply field/multiply)
(def negate field/negate)
(def one field/one)
(def zero field/zero)
(def inverse field/inverse)

(def scale module/scale)

(defn ? [x]
  (and
   (field/?  x)
   (module/?  x)
   (satisfies? VectorSpace x)))
