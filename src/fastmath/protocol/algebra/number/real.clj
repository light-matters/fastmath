(ns fastmath.protocol.algebra.number.real
  (:require
   [fastmath.core :as m]
   [fastmath.protocol.algebra.field :as field]))

(extend-type Double
  field/Field
  (negate [x] (m/- x))
  (inverse [x] (m// x))
  (add [x y] (m/+ x y))
  (multiply [x y] (m/* x y))
  (one [x] 1.0)
  (zero [x] 0.0))

(defn ? [x]
  (satisfies? Double x))

(field/zero 35.0)

