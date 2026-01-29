(ns fastmath.algebra.object.number.complex.create
  "
  Presents an API for complex number creation, for the sake of implementation-independent reference. Specifically, putting it here allows a referent to be added to the complex number protocol.

  This implementation currently uses **EJML**, but can/needs to be updated when the implementation changes.
  "
  (:require [fastmath.algebra.object.number.complex.ejml :as cn])
  (:import
   (org.ejml.data Complex_F64)))

(defn create
  "Creates a complex number, (currently) represented as an ejml `ComplexNumber`.
  Takes optional real and imaginary parts."

  ([^double x ^double y] (cn/->ComplexNumber (Complex_F64. x y)))
  ([^double x] (cn/->ComplexNumber (Complex_F64. x 0.0)))
  ([] cn/ZERO))

(defn <-real
  "A complex number from a real number, for consistency with matrix creation naming. In general, `i` is to be preferred for consistency."
  [^double x] (create x))

(def i "Creates a complex number according to supplied arguments: (0,0), (re,0) or (re,im)."
  create)
