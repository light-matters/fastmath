(ns fastmath.protocol.algebra.structure.normed-space)

(defprotocol NormedSpace
  (norm [x]))

(defn ? [x]
  (satisfies? NormedSpace x))
