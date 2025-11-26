(ns fastmath.protocol.algebra.field
  (:require [fastmath.protocol.algebra.ring :as ring]
            [fastmath.protocol.algebra.multiplicative.group :as mgroup]))

(defprotocol  Field)

(def add (ring/add))
(def multiply (ring/multiply))
(def negate (ring/negate))
(def one (ring/one))
(def zero (ring/zero))

(def one (mgroup/one))

(defn ? [x]
  (and
   (ring/?  x)
   (mgroup/?  x)
   (satisfies? Field x)))

