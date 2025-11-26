(ns fastmath.protocol.algebra.coordinate.polar)

(defprotocol PolarCoordinate
  (angle [z])
  (magnitude [z])
  (polar-values [z]))

(defn ? [x]
  (satisfies? PolarCoordinate x))
