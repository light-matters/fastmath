(ns fastmath.protocol.algebra.structure.module)

(defprotocol Module
  (scale
    ;; "Scalar `a` acts on `x`."
    [x a]))

(defn ? [x]
  (satisfies? Module x))
