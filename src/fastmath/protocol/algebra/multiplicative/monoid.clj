(ns fastmath.protocol.algebra.multiplicative.monoid
  (:require [fastmath.protocol.algebra.multiplicative.semigroup :refer [Semigroup]]))

(defprotocol Monoid
  (one [x]))
