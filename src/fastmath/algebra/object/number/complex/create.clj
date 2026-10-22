(ns fastmath.algebra.object.number.complex.create
  "
  Presents an API for complex number creation, for the sake of implementation-independent reference. Specifically, putting it here allows a referent to be added to the complex number protocol.

  This implementation uses a backend-independent pair of doubles; see [[fastmath.algebra.object.number.complex.primitive]].
  "
  (:require [fastmath.algebra.object.number.complex.primitive :as cn]))

(defn create
  "Creates a complex number. Takes optional real and imaginary parts."

  ([^double x ^double y] (cn/->ComplexNumber x y))
  ([^double x] (cn/->ComplexNumber x 0.0))
  ([] cn/ZERO))

(defn <-real
  "A complex number from a real number, for consistency with matrix creation naming. In general, `i` is to be preferred for consistency."
  [^double x] (create x))

(def i "Creates a complex number according to supplied arguments: (0,0), (re,0) or (re,im)."
  create)
