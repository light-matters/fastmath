(ns fastmath.protocol.algebra.coordinate.complex)

(defprotocol Complex
  (re [z])
  (im [z])
  (conjugate [z]))
