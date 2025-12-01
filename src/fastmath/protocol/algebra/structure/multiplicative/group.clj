(ns fastmath.protocol.algebra.structure.multiplicative.group
  (:require [fastmath.protocol.algebra.structure.multiplicative.monoid :as mmonoid]))

(defprotocol Group
  (inverse [x]))

(def multiply mmonoid/multiply)
(def one mmonoid/one)

(defn ? [x]
  (and
   (mmonoid/? x)
   (satisfies? Group x)))

