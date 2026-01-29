(ns fastmath.algebra.object.type
  (:require
   [fastmath.api.v2.algebra.predicate :as pred]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]))

(derive ::matrix ::type)
(derive ::vector ::type)
(derive ::scalar ::type)

(derive ::matrix--real ::matrix)
(derive ::matrix--real ::real)
(derive ::matrix--complex ::matrix)
(derive ::matrix--complex ::complex)

(derive ::vector ::matrix)
(derive ::vector--real ::vector)
(derive ::vector--real ::real)
(derive ::vector--complex ::vector)
(derive ::vector--complex ::complex)

(derive ::scalar--real ::scalar)
(derive ::scalar--real ::real)
(derive ::scalar--complex ::scalar)
(derive ::scalar--complex ::complex)

(defn ?
  "Classify an argument so the arithmetic multimethods can dispatch on it."
  [x]
  (if (rmat/? x)
    (if (pred/linear-shape? x)
      (if (cmat/? x)
        ::vector--complex
        ::vector--real)
      (if (cmat/? x)
        ::matrix--complex
        ::matrix--real))
    (cond
      (number? x) ::scalar--real
      (pred/complex-number? x) ::scalar--complex)))
