(ns fastmath.api.v2.algebra.complex.number
  "All of your favourite complex number functions in one place."
  (:require [fastmath.protocol.algebra.object.number.complex :as C :refer [inverse conjugate]]
            [fastmath.algebra.object.number.complex.create :as cc]))

(def i cc/create)

(def inverse C/inverse)
(def conjugate C/conjugate)
(def re C/re)
(def im C/im)
(def angle C/angle)
(def magnitude C/magnitude)
(def polar C/polar)
(def ? C/?)





