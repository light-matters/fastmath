(ns fastmath.protocol.algebra.space.normed)

(defprotocol NormedSpace
  (norm [x]))

(defn ? [x]
  (satisfies? NormedSpace x))
