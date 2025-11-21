(ns fastmath.protocol.algebra.coordinate.polar)

(defprotocol Polar
  (angle [z])
  (magnitude [z]))
