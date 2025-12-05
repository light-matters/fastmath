(ns fastmath.protocol.algebra.object.matrix.square
  (:require [fastmath.protocol.algebra.structure.ring :as ring]
            [fastmath.protocol.algebra.object.matrix.rectangular :as rmatrix]))

(defprotocol SquareMatrix
  (determinant [A])
  (trace [A])

  ;; TODO: Maybe these predicates should be implemented separately (function ns)? They tend to be a composite of required features 
  (normal? [A])
  (symmetric? [A])
  (singular? [A])
  (unitary? [A]))

(def multiply ring/multiply)
(def one ring/one)

(defn ? [x]
  (and
   (ring/? x)
   (satisfies? rmatrix/RectangularMatrix x)
   (apply = (rmatrix/shape x))
   (satisfies? rmatrix/SquareMatrix x)))
