(ns fastmath.algebra.object.number.complex.ejml-test
  (:require
   [clojure.test :refer :all]
   [fastmath.algebra.object.number.complex.ejml :as sut]
   [fastmath.protocol.algebra.object.number.complex-test :as ctest]))

(deftest complex-number--ejml
  (ctest/protocol-tests sut/i))

