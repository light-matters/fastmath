(ns fastmath.protocol.algebra.structure.coordinate.polar)

(defprotocol PolarCoordinate
  (angle [z])
  (magnitude [z])
  (polar-values [z]))

(defn ? [x]
  (satisfies? PolarCoordinate x))
