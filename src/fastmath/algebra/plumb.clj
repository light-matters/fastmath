(ns fastmath.algebra.plumb
  "The necessary plumbing: it needs to exist, but we don't want to see it..."
  (:refer-clojure :exclude [type])
  (:require
   [fastmath.algebra.object.type :as tp]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.algebra.object.matrix.create :as mat]))

(defn rank--domain [x]
  (cond
    (isa? x ::tp/real) 0
    (isa? x ::tp/complex) 1
    :else (throw (ex-info "No recognised domain for " {:value x}))))

(defn rank--type [x]
  (cond
    (isa? x ::tp/scalar) 0
    (isa? x ::tp/vector) 1
    (isa? x ::tp/matrix) 2
    :else (throw (ex-info "No recognised algebraic type for " {:value x}))))

(defn promote-domain
  "Increases the set type of the element (if possible), e.g. a real number becomes a complex number or a real matrix becomes a complex matrix."
  [x]
  (case (tp/? x)
    ::tp/scalar--real (C/<-real (double x))
    ::tp/scalar--complex x
    ::tp/matrix--real (mat/<-real x)
    ::tp/matrix--complex x

    (throw (ex-info "No promotion rule for " {:value x}))))

(defn ensure-domain-match
  "Takes a pair of arguments and promotes arguments where necessary to ensure compatible domains.

   e.g. when adding a real number to a complex number, the `double` becomes a complex number with 0 for the 'imaginary' part."
  [a b]
  (let [ta (tp/? a)
        tb (tp/? b)
        ra (rank--domain ta)
        rb (rank--domain tb)]
    (cond
      (= ta tb) [a b]
      (= ra rb) [a b]
      (< ra rb) [(promote-domain a) b]
      :else [a (promote-domain b)])))

