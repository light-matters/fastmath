(ns fastmath.algebra.object.matrix.rectangular.complex.ojalgo-test
  "ComplexDense is not wired into `fastmath.algebra.object.matrix.create`'s
  backend dispatcher (EJML is the only one it builds), so this exercises
  the ojAlgo backend directly rather than through the protocol-level
  tests the EJML backend already has."
  (:require [fastmath.algebra.object.matrix.rectangular.complex.ojalgo :as sut]
            [fastmath.algebra.object.matrix.rectangular.real.ojalgo :as realdense]
            [fastmath.algebra.object.number.complex.create :as C]
            [fastmath.protocol.representation.d2 :as d2]
            [fastmath.protocol.algebra.structure.additive.semigroup :as asg]
            [fastmath.protocol.algebra.structure.additive.monoid :as am]
            [fastmath.protocol.algebra.structure.additive.group :as ag]
            [fastmath.protocol.algebra.structure.module :as module]
            [fastmath.protocol.algebra.structure.coordinate.complex :as cc]
            [fastmath.protocol.algebra.object.number.complex :as pc]
            [fastmath.protocol.algebra.object.matrix.extra :as emat]
            [fastmath.protocol.algebra.object.matrix.square :as sq]
            [fastmath.protocol.algebra.object.matrix.solve :as solve]
            [fastmath.protocol.algebra.object.matrix.complex :as cmat]
            [fastmath.protocol.algebra.object.matrix.complex-predicate :as cpred]
            [fastmath.protocol.algebra.object.matrix.predicate :as pred]
            [clojure.test :as t :refer [deftest testing are is]]))

(def A (sut/<-rows [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]]))
(def B (sut/<-rows [[[5.0 0.0] [6.0 0.0]] [[7.0 0.0] [8.0 0.0]]]))
(def v (sut/<-rows [[[1.0 0.0]] [[2.0 0.0]]]))

(deftest construction-and-retrieval
  (are [e a] (= e a)
    [2 2] (d2/shape A)
    (C/i 2.0 0.0) (d2/element A 0 1)
    [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]] (d2/->vectors A)
    (sut/<-rows [[[1.0 1.0]] [[3.0 0.0]]]) (d2/column A 0)
    (sut/<-rows [[[1.0 1.0] [2.0 0.0]]]) (d2/row A 0)
    (list (C/i 1.0 1.0) (C/i 2.0 0.0) (C/i 3.0 0.0) (C/i 4.0 -1.0)) (seq A)))

(deftest additive-and-module
  (are [e a] (= e a)
    (sut/<-rows [[[6.0 1.0] [8.0 0.0]] [[10.0 0.0] [12.0 -1.0]]]) (asg/add A B)
    (sut/<-rows [[[0.0 0.0] [0.0 0.0]] [[0.0 0.0] [0.0 0.0]]]) (am/zero A)
    (sut/<-rows [[[-1.0 -1.0] [-2.0 0.0]] [[-3.0 0.0] [-4.0 1.0]]]) (ag/negate A)
    (sut/<-rows [[[2.0 2.0] [4.0 0.0]] [[6.0 0.0] [8.0 -2.0]]])
    (module/scale A [2.0 0.0])))

(deftest matrix-extra
  (are [e a] (= e a)
    (sut/<-rows [[[4.0 1.0] [5.0 0.0]] [[6.0 0.0] [7.0 -1.0]]]) (emat/add--s A [3.0 0.0])

    (sut/<-rows [[[19.0 5.0] [22.0 6.0]] [[43.0 -7.0] [50.0 -8.0]]]) (emat/multiply A B)
    (sut/<-rows [[[2.0 2.0] [4.0 0.0]] [[6.0 0.0] [8.0 -2.0]]]) (emat/multiply--s A [2.0 0.0])
    (sut/<-rows [[[5.0 5.0] [12.0 0.0]] [[21.0 0.0] [32.0 -8.0]]]) (emat/multiply--e A B)
    (sut/<-rows [[[5.0 1.0]] [[11.0 -2.0]]]) (emat/multiply--v A v)
    (sut/<-rows [[[4.0 -1.0] [4.0 0.0]] [[4.0 0.0] [4.0 1.0]]]) (emat/subtract B A)
    (sut/<-rows [[[1.0 1.0] [3.0 0.0]] [[2.0 0.0] [4.0 -1.0]]]) (emat/transpose A)
    (C/i 32.0 0.0) (emat/inner A A)
    (sut/<-rows [[[2.0 1.0] [3.0 0.0]] [[4.0 0.0] [5.0 -1.0]]])
    (d2/fmap A (fn [re im] [(inc re) im])))
  (testing "fmap does not mutate the original"
    (d2/fmap A (fn [re im] [(inc re) im]))
    (is (= (sut/<-rows [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]]) A))))

(deftest complex-coordinate-and-matrix
  (are [e a] (= e a)
    (realdense/<-rows [[1.0 2.0] [3.0 4.0]]) (cc/re A)
    (realdense/<-rows [[1.0 0.0] [0.0 -1.0]]) (cc/im A)
    (sut/<-rows [[[1.0 -1.0] [2.0 0.0]] [[3.0 0.0] [4.0 1.0]]]) (cc/conjugate A)
    (sut/<-rows [[[1.0 -1.0] [3.0 0.0]] [[2.0 0.0] [4.0 1.0]]]) (cmat/adjoint A))
  (is (true? (cmat/real? (sut/<-rows [[[1.0 0.0] [2.0 0.0]] [[3.0 0.0] [4.0 0.0]]]))))
  (is (false? (cmat/real? A))))

(deftest square-and-solve
  (let [eps 1.0e-9
        approx= (fn [expected actual]
                  (< (pc/norm
                      (pc/subtract expected actual))
                     eps))]
    (is (approx= (C/i -1.0 3.0) (sq/determinant A)))
    (is (= (C/i 5.0 0.0) (sq/trace A))))
  (let [herm-spd (sut/<-rows [[[4.0 0.0] [2.0 1.0]] [[2.0 -1.0] [3.0 0.0]]])
        {:keys [L lower? spd?]} (solve/cholesky herm-spd)]
    (is (true? lower?))
    (is (true? spd?))
    (is (= herm-spd (emat/multiply L (cmat/adjoint L)))))
  (is (= (sut/<-rows [[[2.0 0.0]] [[3.0 0.0]]])
         (solve/solve (sut/<-rows [[[2.0 0.0] [0.0 0.0]] [[0.0 0.0] [3.0 0.0]]])
                      (sut/<-rows [[[4.0 0.0]] [[9.0 0.0]]])))))

(deftest predicates
  (let [identity-2 (sut/<-rows [[[1.0 0.0] [0.0 0.0]] [[0.0 0.0] [1.0 0.0]]])
        symmetric-not-hermitian (sut/<-rows [[[1.0 2.0] [3.0 4.0]] [[3.0 4.0] [1.0 2.0]]])
        hermitian-not-symmetric (sut/<-rows [[[1.0 0.0] [2.0 3.0]] [[2.0 -3.0] [1.0 0.0]]])
        singular (sut/<-rows [[[1.0 0.0] [2.0 0.0]] [[2.0 0.0] [4.0 0.0]]])
        phase-unitary (sut/<-rows [[[0.0 1.0] [0.0 0.0]] [[0.0 0.0] [0.0 1.0]]])]
    (are [e a] (= e a)
      true  (pred/square? A)
      true  (pred/symmetric? symmetric-not-hermitian)
      false (pred/symmetric? hermitian-not-symmetric)
      true  (cpred/hermitian? hermitian-not-symmetric)
      false (cpred/hermitian? symmetric-not-hermitian)
      true  (pred/normal? phase-unitary)
      true  (pred/singular? singular)
      false (pred/singular? identity-2)
      true  (pred/unitary? identity-2)
      true  (pred/unitary? phase-unitary)
      false (pred/unitary? A))))

(deftest equality
  (is (= A (sut/<-rows [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]])))
  (is (not= A B)))
