(ns fastmath.protocol.algebra.structure.multiplicative.monoid
  (:require [fastmath.protocol.algebra.structure.multiplicative.semigroup :as semigroup]))

(defprotocol Monoid
  (one [x]))

(def multiply semigroup/multiply)

(defn ? [x]
  (and
   (semigroup/? x)
   (satisfies? Monoid x)))
