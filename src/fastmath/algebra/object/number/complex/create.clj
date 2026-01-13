(ns fastmath.algebra.object.number.complex.create
  "Currently implemented using EJML."
  (:require [fastmath.algebra.object.number.complex.ejml :as cn])
  (:import
   (java.lang Math)
   (org.ejml.data Complex_F64 ComplexPolar_F64)
   (org.ejml.ops ComplexMath_F64)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 Constructor                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; TODO: Decide whether it's worth keeping 'contstructors' in a separate namespace. 

(defn create
  "Creates a complex number represented as an ejml `ComplexNumber`.
  Takes optional real and imaginary parts."

  ([^double x ^double y] (cn/ComplexNumber. (Complex_F64. x y)))
  ([^double x] (cn/ComplexNumber. (Complex_F64. x 0.0)))
  ([] cn/ZERO))

(defn <-real [^double x] (create x))

(def i "Same as `complex` above, but a 'friendlier' syntax."
  create)
