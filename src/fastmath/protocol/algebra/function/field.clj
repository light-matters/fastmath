(ns fastmath.protocol.algebra.function.field
  "Standard 'extra' field functions that we expect real and complex numbers to implement, particulaly for efficiency.")

(defprotocol FieldFunction
  (subtract [z1 z2])
  (divide [z1 z2])
  (square [z])
  (square-root [z]))

