(ns fastmath.protocol.algebra.structure.additive.group
  (:require [fastmath.protocol.algebra.structure.additive.monoid :as monoid]))

(defprotocol AdditiveGroup
  (negate [x]))

(def add monoid/add)
(def zero monoid/zero)

(defn ? [x]
  (and
   (monoid/? x)
   (satisfies? AdditiveGroup x)))
