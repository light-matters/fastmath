(ns fastmath.protocol.algebra.object.matrix.square-test
  (:require [fastmath.protocol.algebra.object.matrix.square :as sut]
            [fastmath.protocol.algebra.object.number.complex :as pc]
            [clojure.test :as t :refer [deftest testing are is]]
            [fastmath.algebra.object.number.complex.create :as C]
            [fastmath.algebra.object.matrix.create :as mat]))

(def real2x2 (mat/<-rows [[1.0 2.0] [3.0 4.0]]))
(def complex2x2 (mat/<-rows [[[1.0 1.0] [2.0 0.0]] [[3.0 0.0] [4.0 -1.0]]]))

(defn- approx=
  "The EJML determinant is computed via floating-point multiplication, so
  compare within a small tolerance rather than exactly."
  [expected actual]
  (< (pc/norm (pc/subtract expected actual)) 1.0e-9))

(deftest real-determinant-and-trace
  (testing "determinant/trace return plain doubles for RealDense"
    (are [e a] (= e a)
      -2.0 (sut/determinant real2x2)
      5.0  (sut/trace real2x2))))

(deftest complex-determinant-and-trace
  (testing "determinant/trace return the canonical ComplexNumber, not a raw EJML type"
    (is (approx= (C/i -1.0 3.0) (sut/determinant complex2x2)))
    (is (= (C/i 5.0 0.0) (sut/trace complex2x2)))))
