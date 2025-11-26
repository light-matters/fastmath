(ns fastmath.protocol.algebra.additive.monoid
  (:require [fastmath.protocol.algebra.additive.semigroup :as semigroup]))

(defprotocol Monoid
  (zero [x]))

(def add semigroup/add)

(defn ? [x]
  (and
   (semigroup/? x)
   (satisfies? Monoid x)))
