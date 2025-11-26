(ns fastmath.protocol.algebra.additive.semigroup)

(defprotocol Semigroup
  (add [x y]))

(defn ? [x]
  (satisfies? Semigroup x))
