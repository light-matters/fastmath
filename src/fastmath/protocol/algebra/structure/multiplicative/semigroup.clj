(ns fastmath.protocol.algebra.structure.multiplicative.semigroup)

(defprotocol Semigroup
  (multiply [x y]))

(defn ? [x]
  (satisfies? Semigroup x))
