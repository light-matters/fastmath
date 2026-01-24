(ns fastmath.protocol.algebra.object.number.complex
  (:require
   [fastmath.protocol.algebra.structure.coordinate.complex :as cc]
   [fastmath.protocol.algebra.structure.coordinate.polar :as pc]
   [fastmath.protocol.algebra.structure.field :as field]

   [fastmath.protocol.algebra.structure.space.normed :as nspace]))

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

(def add field/add)
(def inverse field/inverse)
(def multiply field/multiply)
(def zero field/zero)
(def negate field/negate)
(def norm nspace/norm)

(def conjugate cc/conjugate)
(def re cc/re)
(def im cc/im)

(def angle pc/angle)
(def magnitude pc/magnitude)
(def polar pc/polar-values)

(defn ? [x]
  (and
   ;; maths
   (field/? x)
   (nspace/? x)
   (cc/?  x)
   (pc/?  x)
   ;; convenience
   (instance? clojure.lang.Seqable x)
   ;; (satisfies? ComplexNumber x)
   ))
