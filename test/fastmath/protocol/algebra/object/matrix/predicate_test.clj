(ns fastmath.protocol.algebra.object.matrix.predicate-test
  (:require [fastmath.protocol.algebra.object.matrix.extra :as sut]
            [clojure.test :as t :refer [deftest are]]
            [fastmath.algebra.object.number.complex.create :as C]
            [fastmath.algebra.object.matrix.create :as matrix]))

;; (t/testing "Predicates"
;;   (t/are [q a] (= q a)
;;     (sut/? mat2x2) true
;;     (sut/? mat2x2--real) true)
;;   (sut/square? mat2x2) true
;;   (sut/square? mat2x2--real) true)

;; (deftest square
;;   (are [e a] (= e a)
;;     (sut/square?)))
