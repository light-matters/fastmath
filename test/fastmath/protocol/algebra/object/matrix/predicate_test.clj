(ns fastmath.protocol.algebra.object.matrix.predicate-test
  (:require [fastmath.protocol.algebra.object.matrix.extra :as sut]
            [clojure.test :as t :refer [deftest are]]
            [fastmath.algebra.object.number.complex.create :as C]
            [fastmath.algebra.object.matrix.create :as matrix]))

(deftest square
  (are [e a] (= e a)

    (sut/square?)))
