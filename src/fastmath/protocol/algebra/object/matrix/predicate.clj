(ns fastmath.protocol.algebra.object.matrix.predicate)

(defprotocol  MatrixPredicate
  (square? [m])
  (normal? [m])
  (symmetric? [m])
  (singular? [m])
  (unitary? [m]))
