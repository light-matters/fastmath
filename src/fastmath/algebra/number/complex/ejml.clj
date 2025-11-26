(ns fastmath.algebra.number.complex.ejml
  "An EJML implementation of complex numbers.

   TODO: Finish migrating definitons into formal maths-based protocols.
  "
  (:require
   [clojure.math :as math]
   [fastmath.default :as default]
   [fastmath.protocol.algebra.number.complex :as protoComplex]
   [fastmath.protocol.algebra.additive.semigroup :as as]
   [fastmath.protocol.algebra.additive.monoid :as am]
   [fastmath.protocol.algebra.additive.group :as ag]
   [fastmath.protocol.algebra.multiplicative.semigroup :as ms]
   [fastmath.protocol.algebra.multiplicative.monoid :as mm]
   [fastmath.protocol.algebra.multiplicative.group :as mg]
   [fastmath.protocol.algebra.ring :as r]
   [fastmath.protocol.algebra.field :as f]
   [fastmath.protocol.algebra.normed-space :as ns]
   [fastmath.protocol.algebra.coordinate.complex :as cc]
   [fastmath.protocol.algebra.coordinate.polar :as p])
  (:import
   (java.lang Math)
   (org.ejml.data Complex_F64 ComplexPolar_F64)
   (org.ejml.ops ComplexMath_F64)))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                    Cache                                    ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def ^:private zeroF64
  (Complex_F64. 0.0 0.0))
(def ^:private oneF64 (Complex_F64. 1.0 0.0))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                   Utility                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn  format--complex ^String [^double re ^double im ^long p ^double eps]
  (let [rz (if (< (Math/abs re) eps) 0.0 re)
        iz (if (< (Math/abs im) eps) 0.0 im)]
    (str (format (str "%." p "f") rz)
         (format (str "%+" p "f") iz)
         "i")))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                   Type                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype ComplexNumber [^Complex_F64 z]
  protoComplex/ComplexNumber

  as/Semigroup
  (add [_ z2]
    (let [out (Complex_F64.)]
      (ComplexMath_F64/plus z (.-z ^ComplexNumber z2) out)
      (ComplexNumber. out)))
  am/Monoid
  (zero [_] (ComplexNumber. zeroF64))
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
  (one [_] (ComplexNumber. oneF64))
  mg/Group
  (inverse [_]
    (let [out (Complex_F64.)]
      (ComplexMath_F64/divide 1.0 ^Complex_F64 z out)
      (ComplexNumber. out)))

  r/Ring

  f/Field

  ns/NormedSpace
  (norm ^double [_]
    (let [^ComplexPolar_F64 polar zeroF64]
      (ComplexMath_F64/convert z polar)
      (.-r polar)))

  cc/ComplexCoordinate
  (re ^double  [_] (.-real z))
  (im  ^double [_] (.-imaginary z))
  (conjugate [_]
    (ComplexNumber. (Complex_F64. (.-real z) (* -1.0 (.-imaginary z)))))

  p/PolarCoordinate
  (angle ^double [_]
    (let [^ComplexPolar_F64  polar zeroF64]
      (ComplexMath_F64/convert z polar)
      (.-theta polar)))
  (magnitude ^double [_]
    (let [^ComplexPolar_F64 polar zeroF64]
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

  clojure.lang.Seqable
  (seq [_] (list (.-real z) (.-imaginary z)))

  Object
  (toString [_]
    (format--complex ^double (.-real z) ^double (.-imaginary z) 5 default/tolerance)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                  Constants                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; TODO: Should the functions be in the protocol or API ns?
(def ^{:doc "z=0+i"} I (ComplexNumber. (Complex_F64. 0.0 1.0)))
(def ^{:doc "z=0-i"} I- (ComplexNumber. (Complex_F64. 0.0 -1.0)))
(def ^{:doc "z=0-i"} -I (ComplexNumber. (Complex_F64. 0.0 -1.0)))
(def ^{:doc "z=1+0i"} ONE (ComplexNumber. oneF64))
(def ^{:doc "z=2+0i"} TWO (ComplexNumber. (Complex_F64. 2.0 0.0)))
(def ^{:doc "z=0+0i"} ZERO (ComplexNumber. zeroF64))
(def ^{:doc "z=pi+0i"} PI (ComplexNumber. (Complex_F64. math/PI 0.0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 Constructor                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; TODO: Should constructors be moved to some kind of API ns?

(defn complex
  "Creates a complex number represented as an ejml `ComplexNumber`.
  Takes optional real and imaginary parts."

  ([^double x ^double y] (ComplexNumber. (Complex_F64. x y)))
  ([^double x] (ComplexNumber. (Complex_F64. x 0.0)))
  ([] ZERO))
(defn i
  "Same as `complex` above, but in a 'friendlier' syntax."
  ([^double x ^double y] (ComplexNumber. (Complex_F64. x y)))
  ([^double x] (ComplexNumber. (Complex_F64. x 0.0)))
  ([] ZERO))
