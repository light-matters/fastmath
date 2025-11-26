(ns fastmath.protocol.algebra.normed-space)

(defprotocol NormedSpace
  (norm [x]))

(defn ? [x]
  (satisfies? NormedSpace x))
