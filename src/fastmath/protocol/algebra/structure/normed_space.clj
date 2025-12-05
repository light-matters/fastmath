<<<<<<<< HEAD:src/fastmath/protocol/algebra/space/normed.clj
(ns fastmath.protocol.algebra.space.normed)
========
(ns fastmath.protocol.algebra.structure.normed-space)
>>>>>>>> 7d8ae7ecc20c693007f6529b34c4f3a90d495bb1:src/fastmath/protocol/algebra/structure/normed_space.clj

(defprotocol NormedSpace
  (norm [x]))

(defn ? [x]
  (satisfies? NormedSpace x))
