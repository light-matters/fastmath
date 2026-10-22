(ns fastmath.algebra.object.matrix.rectangular.real.ojalgo-test
  "RealDense is not wired into `fastmath.algebra.object.matrix.create`'s
  backend dispatcher (EJML is the only one it builds), so this exercises
  the ojAlgo backend directly rather than through the protocol-level
  tests the EJML backend already has."
  (:require [fastmath.algebra.object.matrix.rectangular.real.ojalgo :as sut]
            [fastmath.protocol.representation.d2 :as d2]
            [fastmath.protocol.algebra.structure.additive.semigroup :as asg]
            [fastmath.protocol.algebra.structure.additive.monoid :as am]
            [fastmath.protocol.algebra.structure.additive.group :as ag]
            [fastmath.protocol.algebra.structure.space.normed :as nspace]
            [fastmath.protocol.algebra.object.matrix.extra :as emat]
            [fastmath.protocol.algebra.object.matrix.square :as sq]
            [fastmath.protocol.algebra.object.matrix.solve :as solve]
            [fastmath.protocol.algebra.object.matrix.predicate :as pred]
            [clojure.test :as t :refer [deftest testing are is]]))

(def A (sut/<-rows [[1.0 2.0] [3.0 4.0]]))
(def B (sut/<-rows [[5.0 6.0] [7.0 8.0]]))
(def v (sut/<-rows [[1.0] [2.0]]))

(deftest construction-and-retrieval
  (are [e a] (= e a)
    [2 2] (d2/shape A)
    2.0 (d2/element A 0 1)
    [[1.0 2.0] [3.0 4.0]] (d2/->vectors A)
    (sut/<-rows [[1.0] [3.0]]) (d2/column A 0)
    (sut/<-rows [[1.0 2.0]]) (d2/row A 0)
    '(1.0 2.0 3.0 4.0) (seq A)))

(deftest additive-and-module
  (are [e a] (= e a)
    (sut/<-rows [[6.0 8.0] [10.0 12.0]]) (asg/add A B)
    (sut/<-rows [[0.0 0.0] [0.0 0.0]]) (am/zero A)
    (sut/<-rows [[-1.0 -2.0] [-3.0 -4.0]]) (ag/negate A)
    (sut/<-rows [[2.0 4.0] [6.0 8.0]]) (fastmath.protocol.algebra.structure.module/scale A 2.0))
  (is (< (Math/abs (- (nspace/norm (sut/<-rows [[1.0 1.0]])) (Math/sqrt 2.0))) 1.0e-9)))

(deftest matrix-extra
  (are [e a] (= e a)
    (sut/<-rows [[1.0] [4.0]]) (emat/diagonal A)
    (sut/<-rows [[4.0 5.0] [6.0 7.0]]) (emat/add--s A 3.0)
    30.0 (emat/inner A A)
    (sut/<-rows [[1.0 2.0 2.0 4.0]
                 [3.0 4.0 6.0 8.0]
                 [3.0 6.0 4.0 8.0]
                 [9.0 12.0 12.0 16.0]])
    (emat/kronecker A A)

    (sut/<-rows [[19.0 22.0] [43.0 50.0]]) (emat/multiply A B)
    (sut/<-rows [[2.0 4.0] [6.0 8.0]]) (emat/multiply--s A 2.0)
    (sut/<-rows [[5.0 12.0] [21.0 32.0]]) (emat/multiply--e A B)
    (sut/<-rows [[5.0] [11.0]]) (emat/multiply--v A v)
    (sut/<-rows [[4.0 4.0] [4.0 4.0]]) (emat/subtract B A)
    (sut/<-rows [[1.0 3.0] [2.0 4.0]]) (emat/transpose A)
    (sut/<-rows [[2.0 3.0] [4.0 5.0]]) (d2/fmap A inc))
  (testing "fmap does not mutate the original"
    (d2/fmap A inc)
    (is (= (sut/<-rows [[1.0 2.0] [3.0 4.0]]) A))))

(deftest square-and-solve
  (are [e a] (= e a)
    -2.0 (sq/determinant A)
    5.0 (sq/trace A))
  (let [{:keys [L lower? spd?]} (solve/cholesky (sut/<-rows [[4.0 2.0] [2.0 3.0]]))]
    (is (true? lower?))
    (is (true? spd?))
    (is (= (sut/<-rows [[4.0 2.0] [2.0 3.0]]) (emat/multiply L (emat/transpose L)))))
  (is (thrown? clojure.lang.ExceptionInfo
               (solve/cholesky (sut/<-rows [[1.0 2.0] [2.0 1.0]]))))
  (is (= (sut/<-rows [[2.0] [3.0]])
         (solve/solve (sut/<-rows [[2.0 0.0] [0.0 3.0]]) (sut/<-rows [[4.0] [9.0]])))))

(deftest predicates
  (are [e a] (= e a)
    true  (pred/square? A)
    false (pred/square? (sut/<-rows [[1.0 2.0 3.0] [4.0 5.0 6.0]]))
    true  (pred/symmetric? (sut/<-rows [[2.0 1.0] [1.0 2.0]]))
    false (pred/symmetric? A)
    true  (pred/normal? (sut/<-rows [[0.0 -1.0] [1.0 0.0]]))
    false (pred/normal? (sut/<-rows [[1.0 1.0] [0.0 1.0]]))
    true  (pred/singular? (sut/<-rows [[1.0 2.0] [2.0 4.0]]))
    false (pred/singular? A)
    true  (pred/unitary? (sut/<-rows [[1.0 0.0] [0.0 1.0]]))
    true  (pred/unitary? (sut/<-rows [[0.0 -1.0] [1.0 0.0]]))
    false (pred/unitary? A)))

(deftest equality
  (is (= A (sut/<-rows [[1.0 2.0] [3.0 4.0]])))
  (is (not= A B)))
