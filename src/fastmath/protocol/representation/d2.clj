(ns fastmath.protocol.representation.d2
  "Access and manipulation methods for '2D', e.g. table-like, things.")

(defprotocol D2
;; Information
  (shape [d2])
  (nrows [d2])
  (ncols [d2])

;; Retrieval
  (element [d2 row col])
  (column [d2 id])
  (row [d2 id])
  (columns [d2])
  (rows [d2])

  ;; TODO:
  ;; - below would be useful
  ;; (extract [m [i0 i1] [j0 j1]])

;; Transformation
  ;; Assuming that `clojure.lang.Seqable` is implemented
  (->array ; "Returns a 2D-array"
    [d2])
  (fmap [d2 f])
  ;; WARNING:
  ;; - the first argument of`map--m` is not the same as `core`. This is due to a limitation of protocols.
  ;; - `f` should also be a function of two arguments when applied to a complex matrix
  ;; 
  )
