(ns fastmath.protocol.algebra.ring
  (:require [fastmath.protocol.algebra.additive.group :as agroup]
            [fastmath.protocol.algebra.multiplicative.monoid :as mmonoid]))

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
   (satisfies? r/Ring x)))
