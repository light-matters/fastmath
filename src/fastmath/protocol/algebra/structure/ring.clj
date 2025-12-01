(ns fastmath.protocol.algebra.structure.ring
  (:require [fastmath.protocol.algebra.structure.additive.group :as agroup]
            [fastmath.protocol.algebra.structure.multiplicative.monoid :as mmonoid]))

(defprotocol Ring)

(def add agroup/add)
(def negate agroup/negate)
(def zero agroup/zero)

(def multiply mmonoid/multiply)
(def one mmonoid/one)

(defn ? [x]
  (and
   (agroup/?  x)
   (mmonoid/?  x)
   (satisfies? Ring x)))
