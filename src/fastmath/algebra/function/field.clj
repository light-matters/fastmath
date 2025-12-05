(ns fastmath.algebra.function.field
  (:require
   [fastmath.protocol.algebra.additive.group :as ag]
   [fastmath.protocol.algebra.multiplicative.group :as mg]))

(defn subtract
  "x - y = x + (-y) or -x."
  ([x] (ag/negate x))
  ([x y]
   (ag/add x (ag/negate y))))

(defn divide
  "x / y = x * y⁻¹
   Assumes x, y lie in a field domain and y != 0."
  [x y]
  (mg/multiply x (mg/inverse y)))

(defn square
  "x² = x * x"
  [x]
  (mg/multiply x x))
