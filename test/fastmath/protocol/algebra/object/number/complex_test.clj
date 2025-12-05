(ns fastmath.protocol.algebra.object.number.complex-test
  ;; TODO:
  ;; - Check that all protocol methods are implemented.
  (:require
   [clojure.test :refer [is testing]]
   [fastmath.core :as m]
   [fastmath.protocol.algebra.structure.coordinate.complex :as coor]
   [fastmath.protocol.algebra.structure.coordinate.polar :as polar]
   [fastmath.protocol.algebra.structure.field :as field]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]
   [fastmath.protocol.algebra.object.number.complex :as Z]))

(defn protocol-tests
  "Takes a complex number constructor, `i`, and uses it to check conformity with the protocol."
  [i]
  (let [;; TODO: Should check that the following constants are implemented
        ONE (i 1.0 0.0)
        I (i 0.0 1.0)
        I- (i 0.0 -1.0)]

    (testing "Conveniences"
      (let [z (i 1.5 3.2)]
        (is (= (list 1.5 3.2) (seq z)))
        (is (= (i 1.5 3.2) z))
        (is (= (i 1.5) (i 1.5 0.0)))
        (is (= (i) (i 0.0 0.0)))

        (is (= 3.2 (second z)))
        (is (= [2.5 4.2] (mapv inc z)))))

    (testing "Predicates"
      (let [z I]
        (is (=  (Z/? z) true))))

    (testing "Essential operations"
      (is (= (i 1.0 1.0) (field/add I ONE)))
      (is (= (i 1.0 1.0) (Z/subtract ONE I-)))
      (is (= -90.0 (m/degrees (polar/angle I-))))
      (is (= (i 0.0 -1.0) (coor/conjugate I)))
      (is (= (i 0.44 0.08) (Z/divide (i 1 2) (i 3 4))))
      (is (= (i 0.12 -0.16) (field/inverse (i 3 4))))
      (is (= (i -5.0 10.0) (field/multiply (i 1 2) (i 3 4))))
      (is (= (i -1.0 2.0) (field/negate (i 1.0 -2.0))))
      (is (= (i -3.0 4.0) (Z/square (i 1.0 2.0))))
      (is (= (i 25.0 0.0) (Z/square (i 5.0 0.0))))
      (is (= (m/sqrt 2.0) (nspace/norm (i 1.0 1.0)))))))

(comment
  (require '[fastmath.algebra.object.number.complex.ejml :as sut])
  (require '[fastmath.protocol.algebra.structure.ring :as ring])
  (require '[fastmath.protocol.algebra.structure.additive.semigroup :as asg])
  (require '[fastmath.protocol.algebra.structure.additive.monoid :as am])
  (require '[fastmath.protocol.algebra.structure.additive.group :as ag])
  (require '[fastmath.protocol.algebra.structure.multiplicative.semigroup :as msg])
  (require '[fastmath.protocol.algebra.structure.multiplicative.monoid :as mmonoid])
  (require '[fastmath.protocol.algebra.structure.multiplicative.group :as mgroup])
  (def I (sut/i 0.0 1.0))
  (satisfies? msg/MultiplicativeSemigroup I))
