(ns fastmath.protocol.algebra.additive.monoid
  (:require [fastmath.protocol.algebra.additive.semigroup :refer [Semigroup]]))

(defprotocol Monoid
  (zero [x]))

