(ns fastmath.protocol.algebra.number.complex
  (:require [fastmath.protocol.algebra.number.complex :as sut]
            [fastmath.protocol.algebra.field :as field]
            [fastmath.protocol.algebra.normed-space :as nspace]
            [fastmath.protocol.algebra.coordinate.complex :as ccoord]
            [clojure.test :refer [deftest is are testing]]))

(defn complex-number-tests [ns]
  (let [complex (ns-resolve ns 'complex)
        i (ns-resolve ns 'i)
        ZERO (ns-resolve ns 'ZERO)
        I (ns-resolve ns 'I)
        I- (ns-resolve ns 'I-)
        ONE (ns-resolve ns 'ONE)]

    (testing "Conveniences"
      (let [z (complex 1.5 3.2)]
        (is (= (list 1.5 3.2) z))
        (is (= (i 1.5 3.2) z))
        (is (= (i 1.5) (complex 1.5 0.0)))
        (is (= (i) (complex 0.0 0.0)))
        (is (= ZERO (i 0.0 0.0)))

        (is (= 3.2 (second z)))
        (is (= [2.5 4.2] (mapv inc z)))))

    (testing "Essential operations"
      (is (= (i 1.0 1.0) (field/add I ONE)))
      (is (= (i 1.0 1.0) (c/sub c/ONE I-)))
      (is (= -90.0 (m/degrees (ccoord/angle I-))))
      (is (= (i 0.0 -1.0) (c/conjugate c/I)))
      (is (= (i 0.44 0.08) (c/div (i 1 2) (i 3 4))))
      (is (= (i 0.12 -0.16) (c/reciprocal (i 3 4))))
      (is (= (i -5.0 10.0) (c/mult (i 1 2) (i 3 4))))
      (is (= (i -1.0 2.0) (c/neg (i 1.0 -2.0))))
      (is (= (i -3.0 4.0) (c/sq (i 1.0 2.0))))
      (is (= (i 25.0 0.0) (c/sq (i 5.0 0.0))))
      (is (= (m/sqrt 2.0) (nspace/magnitude (i 1.0 1.0)))))

    (testing "Predicates"
      (let [z (sut/ZERO)]
        (is (=  (sut/? z) true))))))




