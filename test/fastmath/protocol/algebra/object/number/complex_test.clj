(ns fastmath.protocol.algebra.object.number.complex-test
  ;; TODO:
  ;; - Check that all protocol methods are implemented.
  (:require
   [clojure.test :refer [is testing]]
   [fastmath.core :as m]
<<<<<<< HEAD:test/fastmath/protocol/algebra/number/complex_test.clj
   [fastmath.protocol.algebra.coordinate.complex :as coor]
   [fastmath.protocol.algebra.coordinate.polar :as polar]
   [fastmath.protocol.algebra.field :as field]
   [fastmath.protocol.algebra.space.normed.clj :as nspace]
   [fastmath.protocol.algebra.number.complex :as Z]))
=======
   [fastmath.protocol.algebra.structure.coordinate.complex :as coor]
   [fastmath.protocol.algebra.structure.coordinate.polar :as polar]
   [fastmath.protocol.algebra.structure.field :as field]
   [fastmath.protocol.algebra.structure.normed-space :as nspace]
   [fastmath.protocol.algebra.object.number.complex :as Z]))
>>>>>>> 7d8ae7ecc20c693007f6529b34c4f3a90d495bb1:test/fastmath/protocol/algebra/object/number/complex_test.clj

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



