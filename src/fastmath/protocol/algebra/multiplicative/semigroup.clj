(ns fastmath.protocol.algebra.multiplicative.semigroup)

(defprotocol Semigroup
  (multiply [x y]))

(defn ? [x]
  (satisfies? Semigroup x))
