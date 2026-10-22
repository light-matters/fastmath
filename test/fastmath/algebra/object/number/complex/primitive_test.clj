(ns fastmath.algebra.object.number.complex.primitive-test
  (:require
   [clojure.test :refer :all]
   [fastmath.algebra.object.number.complex.create :as sut]
   [fastmath.protocol.algebra.object.number.complex-test :as ctest]))

(deftest complex-number--primitive
  (ctest/protocol-tests sut/i))

