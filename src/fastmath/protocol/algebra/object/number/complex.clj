(ns fastmath.protocol.algebra.object.number.complex
  (:require
   [fastmath.protocol.algebra.structure.coordinate.complex :as complex-coordinate]
   [fastmath.protocol.algebra.structure.coordinate.polar :as polar-coordinate]
   [fastmath.protocol.algebra.structure.field :as field]
   [fastmath.protocol.algebra.structure.space.normed :as normed-space]))

(defprotocol ComplexNumber
  "A placeholder to indicate that a number has implemented all of the methods necessary to function as a `fastmath` complex number. This can be ensured by calling `?`.

   All of the strictly necessary mathematical properties are implemented conventionally with protocols (if you think like a mathematician), i.e. with appeal to the abstract properties. In this way, a complex number can be defined as a field on a normed space that has complex and polar coordinates. 

Other convenience methods are listed below."
  ;; TODO:
  ;; - (temp) put core (but not strictly necessary) operations in another protocol
  ;; - maybe add a way of getting the polar coordinates in one go

  (subtract [z1 z2])
  (divide [z1 z2])
  (square [z])
  (square-root [z])
  ;; (real? [z])
  ;; (imaginary? [z])
  )

(defn ? [x]
  (and
   ;; maths
   (field/? x)
   (normed-space/? x)
   (complex-coordinate/?  x)
   (polar-coordinate/?  x)
   ;; convenience
   (instance? clojure.lang.Seqable x)
   ;; (satisfies? ComplexNumber x)
   ))
