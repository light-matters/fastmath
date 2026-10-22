(ns fastmath.protocol.algebra.object.matrix.extra-test
  (:require [fastmath.protocol.algebra.object.matrix.extra :as sut]
            [fastmath.protocol.representation.d2 :as d2]
            [clojure.test :as t]
            [fastmath.algebra.object.number.complex.create :as C]
            [fastmath.algebra.object.matrix.create :as mat]))

(def mat2x2 (mat/<-coll 2 2
                        [1.0 2.0 3.0 -4.0
                         5.0 6.0 7.0 0.0]))

(def mat2x2--m3 (mat/<-coll 2 2
                            [-2.0 2.0 0.0 -4.0
                             2.0 6.0 4.0 0.0]))
(def mat2x2--real (mat/<-coll 2 2
                              [45 0.0 34 0.0
                               -45 0.0 56 0.0]))

(def mat2x2--real--m3 (mat/<-coll 2 2
                                  [42 0.0 31 0.0
                                   -48 0.0 53 0.0]))

(def mat2x3 (mat/<-coll 2 3
                        [1.0 2.0 3.0 -4.0 5.0 6.0 7.0 0.0 -1 2 3 4]))

(t/deftest basic
  (t/testing "Fixed-size, square matrices."
    (t/are [q a] (= q a)

      (sut/diagonal mat2x2) (mat/<-coll 2 1 [1.0 2.0 7 0])

      (sut/add--s mat2x2 (C/i -3)) mat2x2--m3
      (sut/add--s mat2x2--real (C/i -3)) mat2x2--real--m3

      (sut/inner mat2x2 mat2x2) (C/i 140)
      (sut/inner mat2x2--real mat2x2--real) (C/i 8342)

      (sut/kronecker mat2x2 mat2x2)
      (mat/<-rows
       [[[-3.0 4.0] [11.0 2.0] [11.0 2.0] [-7.0 -24.0]]
        [[-7.0 16.0] [7.0 14.0] [39.0 -2.0] [21.0 -28.0]]
        [[-7.0 16.0] [39.0 -2.0] [7.0 14.0] [21.0 -28.0]]
        [[-11.0 60.0] [35.0 42.0] [35.0 42.0] [49.0 0.0]]])

      (sut/kronecker mat2x2--real mat2x2--real)
      (mat/<-real (mat/<-rows
                   [[2025 1530 1530 1156]
                    [-2025 2520 -1530 1904]
                    [-2025 -1530 2520 1904]
                    [2025 -2520 -2520 3136]])))))

;; -------- Real (genuine RealDense) fixtures --------
(def real2x2 (mat/<-rows [[1.0 2.0] [3.0 4.0]]))
(def real2x2-b (mat/<-rows [[5.0 6.0] [7.0 8.0]]))
(def real2x1 (mat/<-rows [[1.0] [2.0]]))

(t/deftest real-matrix-extra
  (t/testing "diagonal, add--s, inner, kronecker"
    (t/are [e a] (= e a)
      (sut/diagonal real2x2) (mat/<-rows [[1.0] [4.0]])
      (sut/add--s real2x2 3.0) (mat/<-rows [[4.0 5.0] [6.0 7.0]])
      (sut/inner real2x2 real2x2) 30.0
      (sut/kronecker real2x2 real2x2)
      (mat/<-rows [[1.0 2.0 2.0 4.0]
                   [3.0 4.0 6.0 8.0]
                   [3.0 6.0 4.0 8.0]
                   [9.0 12.0 12.0 16.0]])))
  (t/testing "multiply, multiply--s, multiply--e, multiply--v, subtract, transpose"
    (t/are [e a] (= e a)
      (sut/multiply real2x2 real2x2-b) (mat/<-rows [[19.0 22.0] [43.0 50.0]])
      (sut/multiply--s real2x2 2.0) (mat/<-rows [[2.0 4.0] [6.0 8.0]])
      (sut/multiply--e real2x2 real2x2-b) (mat/<-rows [[5.0 12.0] [21.0 32.0]])
      (sut/multiply--v real2x2 real2x1) (mat/<-rows [[5.0] [11.0]])
      (sut/subtract real2x2-b real2x2) (mat/<-rows [[4.0 4.0] [4.0 4.0]])
      (sut/transpose real2x2) (mat/<-rows [[1.0 3.0] [2.0 4.0]]))))

;; -------- Complex (ComplexDense) fixtures --------
(def c2x2 (mat/<-rows [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]]))
(def c2x2-b (mat/<-rows [[[5.0 0.0] [6.0 0.0]] [[7.0 0.0] [8.0 0.0]]]))
(def c2x1 (mat/<-rows [[[1.0 0.0]] [[2.0 0.0]]]))

(t/deftest complex-matrix-extra
  (t/testing "multiply, multiply--s, multiply--e, multiply--v, subtract, transpose"
    (t/are [e a] (= e a)
      (sut/multiply c2x2 c2x2-b)
      (mat/<-rows [[[19.0 5.0] [22.0 6.0]] [[43.0 -7.0] [50.0 -8.0]]])

      (sut/multiply--s c2x2 [2.0 0.0])
      (mat/<-rows [[[2.0 2.0] [4.0 0.0]] [[6.0 0.0] [8.0 -2.0]]])

      (sut/multiply--e c2x2 c2x2-b)
      (mat/<-rows [[[5.0 5.0] [12.0 0.0]] [[21.0 0.0] [32.0 -8.0]]])

      (sut/multiply--v c2x2 c2x1)
      (mat/<-rows [[[5.0 1.0]] [[11.0 -2.0]]])

      (sut/subtract c2x2-b c2x2)
      (mat/<-rows [[[4.0 -1.0] [4.0 0.0]] [[4.0 0.0] [4.0 1.0]]])

      (sut/transpose c2x2)
      (mat/<-rows [[[1.0 1.0] [3.0 0.0]] [[2.0 0.0] [4.0 -1.0]]]))))

(t/deftest fmap-test
  (t/testing "RealDense fmap - works directly, no complex promotion needed"
    (t/is (= (mat/<-rows [[2.0 3.0] [4.0 5.0]])
             (d2/fmap real2x2 inc)))
    (t/testing "does not mutate the original"
      (d2/fmap real2x2 inc)
      (t/is (= (mat/<-rows [[1.0 2.0] [3.0 4.0]]) real2x2))))
  (t/testing "ComplexDense fmap"
    (t/is (= (mat/<-rows [[[2.0 1.0] [3.0 0.0]] [[4.0 0.0] [5.0 -1.0]]])
             (d2/fmap c2x2 (fn [re im] [(inc re) im]))))
    (t/testing "does not mutate the original"
      (d2/fmap c2x2 (fn [re im] [(inc re) im]))
      (t/is (= (mat/<-rows [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]]) c2x2)))))
