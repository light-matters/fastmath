(ns fastmath.protocol.algebra.object.matrix.rectangular
  (:require
   [fastmath.protocol.algebra.structure.space.vector :as vspace]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]))

(defprotocol RectangularMatrix)

(def add vspace/add)
(def zero vspace/zero)
(def negate vspace/negate)
(def norm nspace/norm)
(def scale vspace/scale)

(defn ? [x]
  (and
   ;; maths
   (vspace/? x)
   (nspace/? x)
   ;;
   (satisfies? RectangularMatrix x)))
