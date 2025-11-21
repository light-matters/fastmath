(ns fastmath.protocol.algebra.additive.group
  (:require [fastmath.protocol.algebra.additive.monoid :refer [Monoid]]))

(defprotocol Group
  (negate [x]))
