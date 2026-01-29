(ns fastmath.api.v2.algebra.complex
  "A convenient 'front-end' to working with complex entities. Outside of constructors, like `i`, all functions should work on both numbers and matrices."
  (:require [fastmath.protocol.algebra.object.number.complex :as pc]
            [fastmath.algebra.object.number.complex.create :as cc]))
;; TODO: Add tests!
;; - Should all of the individual methods really be re-exported here?

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;             Constructor             ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def i "Creates a complex number aording to supplied arguments: (0,0), (re,0) or (re,im)."
  cc/i)

(def add pc/add)
(def inverse pc/inverse)
(def multiply pc/multiply)
(def zero pc/zero)
(def negate pc/negate)
(def norm pc/norm)

(def conjugate pc/conjugate)
(def re pc/re)
(def im pc/im)

(def angle pc/angle)
(def magnitude pc/magnitude)
(def polar pc/polar)
