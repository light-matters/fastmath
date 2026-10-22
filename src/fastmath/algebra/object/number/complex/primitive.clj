(ns fastmath.algebra.object.number.complex.primitive
  "A primitive, backend-independent implementation of complex numbers as a pair of doubles (re, im). This is the canonical complex number representation: matrix backends (EJML, ojAlgo, ...) convert to/from this at their own boundary rather than this type depending on any of them.
  "
  (:require
   [fastmath.default :as default]
   [fastmath.vector :as v]
   [fastmath.protocol.algebra.function.field :as ff]
   [fastmath.protocol.algebra.object.number.complex :as C]
   [fastmath.protocol.algebra.structure.additive.semigroup :as as]
   [fastmath.protocol.algebra.structure.additive.monoid :as am]
   [fastmath.protocol.algebra.structure.additive.group :as ag]
   [fastmath.protocol.algebra.structure.multiplicative.semigroup :as ms]
   [fastmath.protocol.algebra.structure.multiplicative.monoid :as mm]
   [fastmath.protocol.algebra.structure.multiplicative.group :as mg]
   [fastmath.protocol.algebra.structure.ring :as r]
   [fastmath.protocol.algebra.structure.field :as f]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]
   [fastmath.protocol.algebra.structure.coordinate.complex :as cc]
   [fastmath.protocol.algebra.structure.coordinate.polar :as polar])
  (:import
   (java.lang Math)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

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

(deftype ComplexNumber [^double re ^double im]

  as/AdditiveSemigroup
  (add [_ z2]
    (ComplexNumber. (+ re (.-re ^ComplexNumber z2)) (+ im (.-im ^ComplexNumber z2))))
  am/AdditiveMonoid
  (zero [_] (ComplexNumber. 0.0 0.0))
  ag/AdditiveGroup
  (negate [_] (ComplexNumber. (- re) (- im)))

  ms/MultiplicativeSemigroup
  (multiply [_ z2]
    (let [re2 (.-re ^ComplexNumber z2)
          im2 (.-im ^ComplexNumber z2)]
      (ComplexNumber. (- (* re re2) (* im im2))
                       (+ (* re im2) (* im re2)))))
  mm/MultiplicativeMonoid
  (one [_] (ComplexNumber. 1.0 0.0))
  mg/MultiplicativeGroup
  (inverse [_]
    (let [denom (+ (* re re) (* im im))]
      (ComplexNumber. (/ re denom) (/ (- im) denom))))

  r/Ring

  f/Field

  nspace/NormedSpace
  (norm ^double [_] (Math/hypot re im))

  cc/ComplexCoordinate
  (re ^double [_] re)
  (im ^double [_] im)
  (conjugate [_] (ComplexNumber. re (- im)))

  polar/PolarCoordinate
  (angle ^double [_] (Math/atan2 im re))
  (magnitude ^double [_] (Math/hypot re im))
  (polar-values [_] [(Math/hypot re im) (Math/atan2 im re)])

  ff/FieldFunction
  (subtract [_ z2]
    (ComplexNumber. (- re (.-re ^ComplexNumber z2)) (- im (.-im ^ComplexNumber z2))))
  (divide [_ z2]
    (let [re2 (.-re ^ComplexNumber z2)
          im2 (.-im ^ComplexNumber z2)
          denom (+ (* re2 re2) (* im2 im2))]
      (ComplexNumber. (/ (+ (* re re2) (* im im2)) denom)
                       (/ (- (* im re2) (* re im2)) denom))))
  (square [_]
    (ComplexNumber. (- (* re re) (* im im)) (* 2.0 re im)))
  (square-root [_]
    ;; principal branch: non-negative real part
    (let [magnitude (Math/hypot re im)
          sign      (if (neg? im) -1.0 1.0)]
      (ComplexNumber. (Math/sqrt (/ (+ magnitude re) 2.0))
                       (* sign (Math/sqrt (/ (- magnitude re) 2.0))))))

;; ==================================================
  clojure.lang.Seqable
;; ==================================================
  (seq [_] (list re im))
;; ==================================================
  clojure.lang.Sequential
;; ==================================================

  Object
  (toString [_]
    (format--complex re im 5 default/tolerance))
  (equals [_ z2]
    (and (instance? ComplexNumber z2)
         (let [^ComplexNumber z2 z2]
           (and (== (.-re z2) re)
                (== (.-im z2) im)))))
  (hashCode [_]
    (mix-collection-hash (unchecked-int (v/dhash-code (v/dhash-code re) im)) 2))

  clojure.lang.IHashEq
  (hasheq [_]
    (mix-collection-hash (unchecked-int (v/dhash-code (v/dhash-code re) im)) 2)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                  Constants                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; TODO: Should the functions be in the protocol or API ns?
(def ^{:doc "z=0+i"} I (ComplexNumber. 0.0 1.0))
(def ^{:doc "z=0-i"} I- (ComplexNumber. 0.0 -1.0))
(def ^{:doc "z=0-i"} -I (ComplexNumber. 0.0 -1.0))
(def ^{:doc "z=1+0i"} ONE (ComplexNumber. 1.0 0.0))
(def ^{:doc "z=2+0i"} TWO (ComplexNumber. 2.0 0.0))
(def ^{:doc "z=0+0i"} ZERO (ComplexNumber. 0.0 0.0))
(def ^{:doc "z=pi+0i"} PI (ComplexNumber. Math/PI 0.0))

(comment

  (println "testing")
  (C/norm I)

  (instance? clojure.lang.Seqable (ComplexNumber. 1.0 1.0))
  (f/inverse (ComplexNumber. 3.0 4.0))
  (f/negate (ComplexNumber. 1.0 -2.0))
  (polar/angle I-)
  (C/square (ComplexNumber. 5.0 1.0)))
