(ns fastmath.protocol.algebra.coordinate.complex)

(defprotocol ComplexCoordinate
  (re [z])
  (im [z])
  (conjugate [z]))

(defn ? [x]
  (satisfies? ComplexCoordinate x))
