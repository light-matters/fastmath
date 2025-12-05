(ns fastmath.protocol.algebra.object.matrix.square
  (:require [fastmath.protocol.algebra.structure.ring :as ring]
            [fastmath.protocol.algebra.object.matrix.rectangular :as rmatrix]))

(defprotocol SquareMatrix)

(def multiply ring/multiply)
(def one ring/one)

(defn ? [x]
  (and
   (satisfies? rmatrix/RectangularMatrix x)
   (satisfies? rmatrix/SquareMatrix x)))
