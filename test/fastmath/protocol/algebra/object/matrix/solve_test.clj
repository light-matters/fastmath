(ns fastmath.protocol.algebra.object.matrix.solve-test
  (:require [fastmath.protocol.algebra.object.matrix.solve :as sut]
            [fastmath.protocol.algebra.object.matrix.extra :as emat]
            [fastmath.protocol.algebra.object.matrix.complex :as cmat]
            [clojure.test :as t :refer [deftest testing is are]]
            [fastmath.algebra.object.matrix.create :as mat]))

;; Cholesky factors involve irrational numbers (sqrt), so instead of
;; hardcoding L we check the defining property L*L^T = A (or L*L* = A for
;; complex), which is also a stronger test of correctness than a fixed value.

(deftest real-cholesky
  (let [A    (mat/<-rows [[4.0 2.0] [2.0 3.0]])
        {:keys [L lower? spd?]} (sut/cholesky A)]
    (testing "reports lower-triangular SPD factorization"
      (is (true? lower?))
      (is (true? spd?)))
    (testing "L * L^T reconstructs A"
      (is (= A (emat/multiply L (emat/transpose L)))))
    (testing "does not mutate the original matrix"
      (is (= (mat/<-rows [[4.0 2.0] [2.0 3.0]]) A))))
  (testing "throws on a non-square matrix"
    (is (thrown? clojure.lang.ExceptionInfo
                 (sut/cholesky (mat/<-rows [[1.0 2.0 3.0] [4.0 5.0 6.0]])))))
  (testing "throws on a non-SPD matrix"
    (is (thrown? clojure.lang.ExceptionInfo
                 (sut/cholesky (mat/<-rows [[1.0 2.0] [2.0 1.0]]))))))

(deftest complex-cholesky
  (let [A    (mat/<-rows [[[4.0 0.0] [2.0 1.0]] [[2.0 -1.0] [3.0 0.0]]])
        {:keys [L lower? spd?]} (sut/cholesky A)]
    (testing "reports lower-triangular SPD factorization"
      (is (true? lower?))
      (is (true? spd?)))
    (testing "L * L* (conjugate transpose) reconstructs A"
      (is (= A (emat/multiply L (cmat/adjoint L)))))
    (testing "does not mutate the original matrix"
      (is (= (mat/<-rows [[[4.0 0.0] [2.0 1.0]] [[2.0 -1.0] [3.0 0.0]]]) A)))))

(deftest real-solve
  (testing "A x = b, exact diagonal case"
    (is (= (mat/<-rows [[2.0] [3.0]])
           (sut/solve (mat/<-rows [[2.0 0.0] [0.0 3.0]])
                      (mat/<-rows [[4.0] [9.0]])))))
  (testing "A x = b round-trip for a non-trivial matrix"
    (let [A (mat/<-rows [[4.0 2.0] [2.0 3.0]])
          b (mat/<-rows [[6.0] [5.0]])
          x (sut/solve A b)]
      (is (= b (emat/multiply A x)))))
  (testing "does not mutate the original matrix"
    (let [A (mat/<-rows [[2.0 0.0] [0.0 3.0]])]
      (sut/solve A (mat/<-rows [[4.0] [9.0]]))
      (is (= (mat/<-rows [[2.0 0.0] [0.0 3.0]]) A)))))

(deftest complex-solve
  (testing "A x = b, exact diagonal case"
    (is (= (mat/<-rows [[[2.0 0.0]] [[3.0 0.0]]])
           (sut/solve (mat/<-rows [[[2.0 0.0] [0.0 0.0]] [[0.0 0.0] [3.0 0.0]]])
                      (mat/<-rows [[[4.0 0.0]] [[9.0 0.0]]])))))
  (testing "A x = b round-trip for a non-trivial matrix"
    (let [A (mat/<-rows [[[4.0 0.0] [2.0 1.0]] [[2.0 -1.0] [3.0 0.0]]])
          b (mat/<-rows [[[6.0 0.0]] [[5.0 0.0]]])
          x (sut/solve A b)]
      (is (= b (emat/multiply A x))))))
