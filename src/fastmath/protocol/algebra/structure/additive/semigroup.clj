(ns fastmath.protocol.algebra.structure.additive.semigroup)

(defprotocol Semigroup
  (add [x y]))

(defn ? [x]
  (satisfies? Semigroup x))
