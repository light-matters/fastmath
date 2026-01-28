(ns fastmath.protocol.algebra.object.matrix.extra
  "Miscellaneous functions that we expect any useful matrix implementation to have, but which are not mathematically necessary or only partally defined, e.g. only work on certain shapes or are only implemented for efficiency reasons.
  "
  (:require
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]))

;; TODO: Add docstrings
;; - Change to 'util'
;; - consider splitting into more categories, e.g. contingent, compatible shape, predicates ...

(defprotocol MatrixExtra
  (add--s [A s])
  (inner [A B])
  (kronecker [A B])

  (multiply [A B])
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
   (satisfies? rmat/RectangularMatrix x)
   (satisfies? MatrixExtra x)))

(comment (-> [(into [] (range 10)) (into [] (range 10))]
             flatten
             double-array
             count))
