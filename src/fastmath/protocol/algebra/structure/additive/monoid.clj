(ns fastmath.protocol.algebra.structure.additive.monoid
  (:require [fastmath.protocol.algebra.structure.additive.semigroup :as semigroup]))

(defprotocol AdditiveMonoid
  (zero [x]))

(def add semigroup/add)

(defn ? [x]
  (and
   (semigroup/? x)
   (satisfies? Monoid x)))
