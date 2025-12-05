(ns fastmath.protocol.algebra.structure.space.vector
  (:require
   [fastmath.protocol.algebra.structure.additive.group :as agroup]
   [fastmath.protocol.algebra.structure.module :as module]))

(defprotocol VectorSpace)

(def add agroup/add)
(def negate agroup/negate)
(def zero agroup/zero)

(def scale module/scale)

(defn ? [x]
  (and
   (agroup/?  x)
   (module/?  x)
   (satisfies? VectorSpace x)))
