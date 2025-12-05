(ns fastmath.protocol.representation.d2
  "Access and manipulation methods for '2D' things.")

(defprotocol D2
;; Information
  (shape [m])
  (num-rows [m])
  (num-cols [m])

;; Retrieval
  (element [m row col])

  (column [m id])
  (row [m id])
  (columns [m])
  (rows [m])

  (diagonal [m])

  ;; TODO:
  ;; - below would be useful
  ;; (extract [m [i0 i1] [j0 j1]])

;; Transformation
  ;; Assuming that `clojure.lang.Seqable` is implemented
  ;; (seq<- [m])
  (array<- ;; "Returns a 2D-array"
    [m]))
