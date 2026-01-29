(ns fastmath.algebra.object.type-test
  (:require [fastmath.algebra.object.type :as sut]
            [fastmath.api.v2.algebra.complex :as C]
            [fastmath.algebra.test-object :as t]
            [clojure.test :refer [deftest are]]))

(deftest type-test
  (are [q a] (= q a)
    (sut/? t/m--c) ::sut/matrix--complex
    (sut/? t/m--vc) ::sut/vector--complex
    (sut/? t/m--vrand) ::sut/vector--complex
    (sut/? t/m--r) ::sut/matrix--real
    (sut/? 1.0) ::sut/scalar--real
    (sut/? (C/i 1.0)) ::sut/scalar--complex))
