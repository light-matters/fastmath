(ns fastmath.protocol.algebra.object.matrix.extra-test
  (:require [fastmath.protocol.algebra.object.matrix.extra :as sut]
            [clojure.test :as t]
            [fastmath.algebra.object.matrix.create :as mcreate]))

(def mat2x2 (mcreate/<-coll 2 2
                            [1.0 2.0 3.0 -4.0
                             5.0 6.0 7.0 0.0]))
(def mat2x2--m3 (mcreate/<-coll 2 2
                                [-2.0 -1.0 0.0 -7.0
                                 2.0 3.0 4.0 -3.0]))
(def mat2x2--real (mcreate/<-coll 2 2
                                  [45 0.0 34 0.0
                                   -45 0.0 56 0.0]))
(def mat2x2--real--m3 (mcreate/<-coll 2 2
                                      [42 -3.0 31 -3.0
                                       -48 -3.0 53 -3.0]))

(def mat2x3 (mcreate/<-coll 2 3
                            [1.0 2.0 3.0 -4.0 5.0 6.0 7.0 0.0 -1 2 3 4]))

(t/deftest basic
  (t/testing "Fixed-size, square matrices."
    (t/are [q a] (= q a)
      (sut/add--s mat2x2 -3) mat2x2--m3
      (sut/add--s mat2x2--real -3) mat2x2--real--m3

      ;; (sut/inner mat2x2 mat2x2)
      ;; (sut/inner mat2x2--real mat2x2--real)

      (sut/outer mat2x2 mat2x2)
      (mcreate/<-coll 4 4
                      [-3 4 , 11 2
                       -7 , 16  , 7 , 14  , 11 , 2  , -7 , 24
                       39 , 2  , 21 , 28
                       -7 , 16  , 39 , 2
                       -11 , 60  , 35 , 42  , 7 , 14  , 21 , 28
                       35 , 42  , 49 , 0])
      (sut/outer mat2x2--real mat2x2--real)
      (mcreate/<-rows
       [[2025 1530 1530 1156]
        [-2025 2520 -1530 1904]
        [-2025 -1530 2520 1904]
        [2025 -2520 -2520 3136]])))

  (t/testing "Predicates"
    (t/are [q a] (= q a)
      (sut/? mat2x2) true
      (sut/? mat2x2--real) true)
    (sut/square? mat2x2) true
    (sut/square? mat2x2--real) true))
