(ns fastmath.protocol.algebra.object.matrix.general
  "Miscellaneous functions that we expect any useful matrix implementation to have, but which are not mathematically necessary or only partally defined, e.g. only work on certain shapes or are only implemented for efficiency reasons.
  "
  (:require
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.structure.space.vector :as vspace]
   [fastmath.protocol.algebra.structure.space.vector :as nspace]))

(defprotocol GeneralMatrix
  (add--s [A s])
  (inner [A B])
  (outer [A B])
  (map--m [A f])
  ;; TODO: Should `map--m` be part of the D2 protocol?
  (multiply [A B])
  (multiply--e [A B])
  (multiply--v [A v])
  (subtract [A B])
  (transpose [A])

  (square? [A])
  ;; TODO: square? should maybe be part of the function nses
  )
(defn ? [x]
  (and
   ;; maths
   (vspace/? x)
   (nspace/? x)
   ;;
   (satisfies? rmat/RectangularMatrix x)))
