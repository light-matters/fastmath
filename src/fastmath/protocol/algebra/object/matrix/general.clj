(ns fastmath.protocol.algebra.object.matrix.general
  "Miscellaneous functions that we expect any useful matrix implementation to have, but which are not mathematically necessary or only partally defined, e.g. only work on certain shapes or are only implemented for efficiency reasons.
  ")

(defprotocol GeneralMatrix
  (add--s [A s])
  (inner [A B])
  (outer [A B])
  (map--m [A f])
  (multiply [A B])
  (multiply--e [A B])
  (multiply--v [A v])
  (subtract [A B])

  (square? [A]))

