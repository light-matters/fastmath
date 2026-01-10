(ns fastmath.protocol.algebra.object.matrix.friendly
  "Miscellaneous functions that we expect any useful matrix implementation to have, but which are not mathematically necessary or only partally defined, e.g. only work on certain shapes or are only implemented for efficiency reasons.
  "
  (:require
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.structure.space.vector :as vspace]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]))

;; TODO: Add docstrings

(defprotocol GeneralMatrix
  (add--s [A s])
  (inner [A B])
  (outer [A B])

  (multiply [A B])
  ;; TODO: Where to put matrix multiplies? Doesn't always work (dimensions).
  ;; TODO: Should these be replaced with mathematical symbols?
  (multiply--e [A B])
  (multiply--v [A v])

  (subtract [A B])
  ;; NOTE: Available for efficiency reasons.
  (transpose [A])

  (square? [A])
  ;; TODO: square? should maybe be part of the function nses
  )
(defn ? [x]
  (and
   ;;
   (satisfies? rmat/RectangularMatrix x)))

(comment (-> [(into [] (range 10)) (into [] (range 10))]
             flatten
             double-array
             count))
