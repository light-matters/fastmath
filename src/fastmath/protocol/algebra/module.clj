(ns fastmath.protocol.algebra.module)

(defprotocol Module
  (scale "Scalar `a` acts on `x`."
    [x a]))

(defn ? [x]
  (satisfies? Module x))
