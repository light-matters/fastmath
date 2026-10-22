(ns fastmath.protocol.algebra.object.matrix.predicate-test
  (:require [fastmath.protocol.algebra.object.matrix.predicate :as sut]
            [fastmath.protocol.algebra.object.matrix.complex-predicate :as csut]
            [clojure.test :as t :refer [deftest testing are]]
            [fastmath.algebra.object.matrix.create :as mat]))

;; -------- Real (RealDense) fixtures --------
(def r-identity (mat/<-rows [[1.0 0.0] [0.0 1.0]]))
(def r-symmetric (mat/<-rows [[2.0 1.0] [1.0 2.0]]))
(def r-nonsymmetric (mat/<-rows [[1.0 2.0] [3.0 4.0]]))
(def r-singular (mat/<-rows [[1.0 2.0] [2.0 4.0]]))
(def r-rotation (mat/<-rows [[0.0 -1.0] [1.0 0.0]])) ; orthogonal, not symmetric
(def r-nonnormal (mat/<-rows [[1.0 1.0] [0.0 1.0]]))
(def r-nonsquare (mat/<-rows [[1.0 2.0 3.0] [4.0 5.0 6.0]]))

;; -------- Complex (ComplexDense) fixtures --------
(def c-identity (mat/<-rows [[[1.0 0.0] [0.0 0.0]] [[0.0 0.0] [1.0 0.0]]]))
(def c-symmetric-not-hermitian
  (mat/<-rows [[[1.0 2.0] [3.0 4.0]] [[3.0 4.0] [1.0 2.0]]]))
(def c-hermitian-not-symmetric
  (mat/<-rows [[[1.0 0.0] [2.0 3.0]] [[2.0 -3.0] [1.0 0.0]]]))
(def c-singular
  (mat/<-rows [[[1.0 0.0] [2.0 0.0]] [[2.0 0.0] [4.0 0.0]]]))
(def c-phase-unitary ; i * identity - unitary, not real-orthogonal
  (mat/<-rows [[[0.0 1.0] [0.0 0.0]] [[0.0 0.0] [0.0 1.0]]]))
(def c-nonnormal
  (mat/<-rows [[[1.0 0.0] [1.0 0.0]] [[0.0 0.0] [1.0 0.0]]]))
(def c-nonsquare
  (mat/<-rows [[[1.0 0.0] [2.0 0.0] [3.0 0.0]] [[4.0 0.0] [5.0 0.0] [6.0 0.0]]]))

(deftest real-predicates
  (testing "square?"
    (are [e a] (= e a)
      true  (sut/square? r-identity)
      false (sut/square? r-nonsquare)))
  (testing "symmetric?"
    (are [e a] (= e a)
      true  (sut/symmetric? r-identity)
      true  (sut/symmetric? r-symmetric)
      false (sut/symmetric? r-nonsymmetric)
      false (sut/symmetric? r-rotation)))
  (testing "normal?"
    (are [e a] (= e a)
      true  (sut/normal? r-identity)
      true  (sut/normal? r-symmetric)
      true  (sut/normal? r-rotation)
      false (sut/normal? r-nonnormal)))
  (testing "singular?"
    (are [e a] (= e a)
      true  (sut/singular? r-singular)
      false (sut/singular? r-identity)
      false (sut/singular? r-symmetric)))
  (testing "unitary? (orthogonal)"
    (are [e a] (= e a)
      true  (sut/unitary? r-identity)
      true  (sut/unitary? r-rotation)
      false (sut/unitary? r-nonsymmetric))))

(deftest complex-predicates
  (testing "square?"
    (are [e a] (= e a)
      true  (sut/square? c-identity)
      false (sut/square? c-nonsquare)))
  (testing "symmetric? (not the same as hermitian?)"
    (are [e a] (= e a)
      true  (sut/symmetric? c-identity)
      true  (sut/symmetric? c-symmetric-not-hermitian)
      false (sut/symmetric? c-hermitian-not-symmetric)))
  (testing "normal?"
    (are [e a] (= e a)
      true  (sut/normal? c-identity)
      true  (sut/normal? c-phase-unitary)
      false (sut/normal? c-nonnormal)))
  (testing "singular?"
    (are [e a] (= e a)
      true  (sut/singular? c-singular)
      false (sut/singular? c-identity)))
  (testing "unitary?"
    (are [e a] (= e a)
      true  (sut/unitary? c-identity)
      true  (sut/unitary? c-phase-unitary)
      false (sut/unitary? c-nonnormal)))
  (testing "hermitian? (not the same as symmetric?)"
    (are [e a] (= e a)
      true  (csut/hermitian? c-identity)
      true  (csut/hermitian? c-hermitian-not-symmetric)
      false (csut/hermitian? c-symmetric-not-hermitian))))
