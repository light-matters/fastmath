(ns fastmath.algebra.function.field
  (:require
   [fastmath.protocol.algebra.additive.group :as ag]
   [fastmath.protocol.algebra.additive.semigroup :as asg]
   [fastmath.protocol.algebra.multiplicative.group :as mg]
   [fastmath.protocol.algebra.multiplicative.semigroup :as msg]))

(defn subtract
  "x - y = x + (-y) or -x."
  ([x] (ag/negate x))
  ([x y]
   (asg/add x (ag/negate y))))

(defn divide
  "x / y = x * y⁻¹
   Assumes x, y lie in a field domain and y != 0."
  [x y]
  (msg/multiply x (mg/inv y)))

(defn square
  "x² = x * x"
  [x]
  (msg/multiply x x))
