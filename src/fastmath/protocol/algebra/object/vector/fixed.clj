(ns fastmath.protocol.algebra.object.vector.fixed
  ;; TODO: Should only be an implementation rather than a protocol? 
  "These static, fixed-size vector types will be implemented for efficiency reasons. Currently limited to four dimensions to cover the majority of physics problems. The row and transpose versions should be implemented for similar reasons.")

(defprotocol  M2x1)
(defprotocol  M3x1)
(defprotocol  M4x1)

(defprotocol  M1x2)
(defprotocol  M1x3)
(defprotocol  M1x4)
