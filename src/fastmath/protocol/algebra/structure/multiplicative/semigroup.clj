(ns fastmath.protocol.algebra.structure.multiplicative.semigroup)

(defprotocol MultiplicativeSemigroup
  (multiply [x y]))

(defn ? [x]
  (satisfies? MultiplicativeSemigroup x))
