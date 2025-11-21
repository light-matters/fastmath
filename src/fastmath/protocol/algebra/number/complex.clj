(ns fastmath.protocol.algebra.number.complex
  (:require [fastmath.protocol.algebra.additive.group :as ag]
            [fastmath.protocol.algebra.multiplicative.group :as mg]
            [fastmath.protocol.algebra.field :as f]
            [fastmath.protocol.algebra.normed-space :as ns]
            [fastmath.protocol.algebra.coordinate.complex :as cc]
            [fastmath.protocol.algebra.coordinate.polar :as p]))

;; (defprotocol NumberComplex
;;   (->seq [z])
;;   (->array [z])

;;   )
(defprotocol Complex)

(defn ? [x]
  (and (satisfies? ag/Group x)
       (satisfies? mg/Group x)
       (satisfies? f/Field x)
       (satisfies? ns/NormedSpace x)
       (satisfies? cc/Complex x)
       (satisfies? p/Polar x)
       (satisfies? Complex x)))
