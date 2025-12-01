(ns fastmath.protocol.algebra.structure.field
  (:require [fastmath.protocol.algebra.structure.ring :as ring]
            [fastmath.protocol.algebra.structure.multiplicative.group :as mgroup]))

(defprotocol  Field)

(def add ring/add)
(def multiply ring/multiply)
(def negate ring/negate)
(def one ring/one)
(def zero ring/zero)

(def inverse mgroup/inverse)

(defn ? [x]
  (and
   (ring/?  x)
   (mgroup/?  x)
   (satisfies? Field x)))

