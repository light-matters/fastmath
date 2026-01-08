(ns fastmath.protocol.algebra.object.matrix.complex-test
  (:require [fastmath.protocol.algebra.object.matrix.complex :as sut]
            [fastmath.algebra.object.matrix.create :as mcreate]
            [clojure.test :as t]))

(t/deftest basic
  (let [mat2x2 (mcreate/<-coll 2 2 [1.0 2.0 3.0 -4.0 5.0 6.0 7.0 0.0])
        mat2x2--2 (mcreate/<-coll 2 2 [45 21 34 -456 -45 67 56 0.0])
        mat2x2--real (mcreate/<-coll 2 2 [45 0.0 34 0.0 -45 0.0 56 0.0])]

    (t/testing "Fixed-size, square matrices."
      (t/are [q a] (= q a)
        (sut/add mat2x2 mat2x2--2)
        (mcreate/<-coll 2 2 [46.0 23.0 37.0 -460.0 -40.0 73.0 63.0 0.0])

        (sut/zero mat2x2)
        (mcreate/<-coll 2 2 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]))

      (sut/negate mat2x2)
      (mcreate/<-coll 2 2 [-1.0 -2.0 -3.0 +4.0 -5.0 -6.0 -7.0 -0.0])

      (sut/scale mat2x2 5.0)
      (mcreate/<-coll 2 2 [5.0 10.0 15.0 -20.0 25.0 30.0 35.0 0.0])

      (sut/norm mat2x2)
      11.8321595662)

    (sut/adjoint mat2x2)
    (mcreate/<-coll 2 2 [1.0 -2.0 5.0 -6.0
                         3.0 +4.0 7.0 -0.0])

    (t/testing "Predicates"
      (t/are [q a] (= q a)
        (sut/real? mat2x2) false
        (sut/real? mat2x2--real) true
        (sut/real? (mcreate/<-coll 2 2 [45 34  -45  56])) true
        (sut/? (mcreate/<-coll 2 2 [45 34  -45  56])) false
        (sut/? mat2x2) true))))


