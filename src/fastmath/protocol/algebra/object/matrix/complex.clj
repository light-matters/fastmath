(ns fastmath.protocol.algebra.object.matrix.complex
  "Holds definitions that *only* apply to complex, rectangular matrices."
  (:require [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
            [fastmath.protocol.algebra.structure.coordinate.complex :as cc]))

(defprotocol ComplexMatrix
  (adjoint [A])
  ;; (hermitian? [A])
  ;; - TODO: move hermitian? to function or square ns
  (real?
    ;; "Checks if the imaginary components of all of the matrix elements are close to zero, using the default tolerance.
    ;; N.B. This is a protocol method due to its low-level implementation (for efficiency).
    ;;")
    [A]))

(def add rmat/add)
(def zero rmat/zero)
(def negate rmat/negate)
(def scale rmat/scale)
(def norm rmat/norm)

(defn ? [x]
  (and
   (satisfies? rmat/RectangularMatrix x)
   (cc/? x)
   (satisfies? ComplexMatrix x)))
