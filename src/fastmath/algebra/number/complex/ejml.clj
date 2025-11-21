(ns fastmath.algebra.number.complex.ejml
  "An EJML implementation of complex numbers.

   TODO: Finish migrating definitons into formal maths-based protocols.
  "
  (:require
   [fastmath.default :as default]
   [fastmath.protocol.algebra.number.complex :as pC]
   [fastmath.protocol.algebra.additive.semigroup :as as]
   [fastmath.protocol.algebra.additive.monoid :as am]
   [fastmath.protocol.algebra.additive.group :as ag]
   [fastmath.protocol.algebra.multiplicative.semigroup :as ms]
   [fastmath.protocol.algebra.multiplicative.monoid :as mm]
   [fastmath.protocol.algebra.multiplicative.group :as mg]
   [fastmath.protocol.algebra.field :as f]
   [fastmath.protocol.algebra.normed-space :as ns]
   [fastmath.protocol.algebra.coordinate.complex :as cc]
   [fastmath.protocol.algebra.coordinate.polar :as p])
  (:import
   (java.lang Math)
   (org.ejml.data Complex_F64 ComplexPolar_F64)
   (org.ejml.ops ComplexMath_F64)))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                   Utility                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn  format--complex ^String [^double re ^double im ^long p ^double eps]
  (let [rz (if (< (Math/abs re) eps) 0.0 re)
        iz (if (< (Math/abs im) eps) 0.0 im)]
    (str (format (str "%." p "f") rz)
         (format (str "%+" p "f") iz)
         "i")))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype ComplexNumber [^Complex_F64 z]
  pC/Complex

  as/Semigroup
  (add [_ z2]
    (let [out (Complex_F64.)]
      (ComplexMath_F64/plus z (.-z ^ComplexNumber z2) out)
      (ComplexNumber. out)))
  am/Monoid
  ;; TODO: Cache the zero?
  (zero [_] (-> (Complex_F64. 0.0 0.0)
                ComplexNumber.))
  ag/Group
  (negate [_]
    (-> (Complex_F64. (.-real z) (.-imaginary z))
        ComplexNumber.))

  ms/Semigroup
  (multiply [_ z2]
    (let [out (Complex_F64.)]
      (ComplexMath_F64/multiply z ^Complex_F64 (.-z ^ComplexNumber z2) out)
      (ComplexNumber. out)))
  mm/Monoid
  (one [_] (-> (Complex_F64. 1.0 0.0)
               ComplexNumber.))
  mg/Group
  (inverse [_]
    (let [out (Complex_F64.)]
      (ComplexMath_F64/divide 1.0 ^Complex_F64 z out)
      (ComplexNumber. out)))

  f/Field

  ns/NormedSpace
  (norm [_]
    (let [polar (ComplexPolar_F64. 0.0 0.0)]
      (ComplexMath_F64/convert z polar)
      (.-r polar)))

  cc/Complex
  (re [_] (.-real z))
  (im [_] (.-imaginary z))
  (conjugate [_] (ComplexNumber. (Complex_F64. (.-real z) (* -1.0 (.-imaginary z)))))

  p/Polar
  (angle [_] (let [polar (ComplexPolar_F64. 0.0 0.0)]
               (ComplexMath_F64/convert z polar)
               (.-theta polar)))
  (magnitude [_]
    (let [polar (ComplexPolar_F64. 0.0 0.0)]
      (ComplexMath_F64/convert z polar)
      (.-r polar)))

;; (->seq [_] (seq [(.-real z) (.-imaginary z)]))
  ;; (->array [_] (double-array [(.-real z) (.-imaginary z)]))

  ;; (sub [_ z2]
  ;;   (let [out (Complex_F64.)]
  ;;     (ComplexMath_F64/minus z (.-z ^ComplexNumber z2) out)
  ;;     (ComplexNumber. out)))

  ;; (div [_ z2]
  ;;   (let [out (Complex_F64.)]
  ;;     (ComplexMath_F64/divide z ^Complex_F64 (.-z ^ComplexNumber z2) out)
  ;;     (ComplexNumber. out)))

  Object
  (toString [_]
    (format--complex ^double (.-real z) ^double (.-imaginary z) 5 default/tolerance)))
