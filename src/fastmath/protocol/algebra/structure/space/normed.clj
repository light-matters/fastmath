(ns fastmath.protocol.algebra.structure.space.normed)

(defprotocol NormedSpace
  (norm [x]))

(defn ? [x]
  (satisfies? NormedSpace x))
