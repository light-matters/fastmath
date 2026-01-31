(ns fastmath.protocol.algebra.object.matrix.square
  (:require [fastmath.protocol.algebra.structure.ring :as ring]
            [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]))

(defprotocol SquareMatrix
  (determinant [A])
  (trace [A]))

(def add rmat/add)
(def zero rmat/zero)
(def negate rmat/negate)
(def scale rmat/scale)
(def norm rmat/norm)

(def multiply ring/multiply)
(def one ring/one)

(defn ? [x]
  (and
   (apply = (rmat/shape x))
   (ring/? x)
   (satisfies? rmat/RectangularMatrix x)
   (satisfies? rmat/SquareMatrix x)))
