(ns fastmath.algebra.number.complex.ejml-test
  (:require
   [clojure.test :refer :all]
   [fastmath.algebra.number.complex.ejml :as sut]
   [fastmath.protocol.algebra.number.complex-test :as ctest]))

(deftest complex-number--ejml
  (ctest/protocol-tests sut/i))

