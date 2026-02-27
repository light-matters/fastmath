(ns fastmath.protocol.algebra.object.matrix.solve
  (:require [fastmath.protocol.algebra.structure.ring :as ring]
            [fastmath.protocol.representation.d2 :as d2]
            [fastmath.protocol.algebra.object.matrix.square :as sq]
            [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]))

(defprotocol MatrixSolve
  (cholesky [m])
  (solve [m v]))

(defn ? [x]
  (and
   (sq/? x)))
