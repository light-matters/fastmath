(ns fastmath.protocol.algebra.object.vector.general
  "Miscellaneous functions that we expect any useful matrix implementation to have, but which are not mathematically necessary or only partally defined, e.g. only work on certain shapes or are only implemented for efficiency reasons.")

;; TODO: Separate protocols into 'definition' and 'extras'-like thing

(defprotocol  GeneralVector
  (row? [v]))

(defn ? [v])
