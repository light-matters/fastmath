(ns fastmath.protocol.algebra.additive.group
  (:require [fastmath.protocol.algebra.additive.monoid :as monoid]))

(defprotocol Group
  (negate [x]))

(def add monoid/add)
(def zero monoid/zero)

(defn ? [x]
  (and
   (monoid/? x)
   (satisfies? Group x)))
