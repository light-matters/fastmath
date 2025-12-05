(ns fastmath.protocol.algebra.structure.additive.semigroup)

(defprotocol AdditiveSemigroup
  (add [x y]))

(defn ? [x]
  (satisfies? AdditiveSemigroup x))
