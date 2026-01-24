(ns fastmath.api.v2.algebra.general-test
  (:require [fastmath.api.v2.algebra.general :as sut]
            [fastmath.algebra.object.matrix.create :as mat]
            [fastmath.algebra.object.number.complex.create :as C]
            [clojure.test :refer [deftest is are]]))

(def m--c
  (mat/<-coll 3 3 (partition 2 (range 18))))
(def m--vc
  (mat/<-coll 1 3 (partition 2 (range 6))))
(def m--vrand
  (let [n (+ 1 (rand-int 10))]
    (mat/<-coll n 1 (range (* 2 n)))))
(def m--r
  (mat/<-coll 2 2 [1 0 0 1]))

(deftest linear-shape?-test
  (are [q a] (= q a)
    (sut/linear-shape? m--c) false
    (sut/linear-shape? m--vc) true
    (sut/linear-shape? m--vrand) true))

(deftest scalar?-test
  (are [q a] (= q a)
    (sut/scalar? 1.0) true
    (sut/scalar? (C/i 1.0 0)) true
    (sut/scalar? (C/i 1.0 -10.5)) true
    (sut/scalar? "1.0") false
    (sut/scalar? m--c) false
    (sut/scalar? m--vc) false
    (sut/scalar? m--vrand) false))

(deftest matrix?-test
  (are [q a] (= q a)
    (sut/matrix? m--c) true
    (sut/matrix? m--vc) true
    (sut/matrix? m--vrand) true
    (sut/matrix? "thng") false
    (sut/matrix? 1.0) false
    (sut/matrix? (C/i 1.0)) false))

(deftest vector?-test
  (are [q a] (= q a)
    (sut/vector? m--c) false
    (sut/vector? m--vc) true
    (sut/vector? m--vrand) true
    (sut/vector? "thng") false
    (sut/vector? 1.0) false
    (sut/vector? (C/i 1.0)) false))

(deftest real?-test
  (are [q a] (= q a)
    (sut/real? m--c) false
    (sut/real? m--vc) false
    (sut/real? m--vrand) false
    (sut/real? m--r) true
    ;; (sut/real? "thng") (is (thrown? java.lang.Exception))
    (sut/real? 1.0) true
    (sut/real? (C/i 1.0)) false))

(deftest type-test
  (are [q a] (= q a)
    (sut/type m--c) ::sut/matrix--complex
    (sut/type m--vc) ::sut/vector--complex
    (sut/type m--vrand) ::sut/vector--complex
    (sut/type m--r) ::sut/matrix--real
    (sut/type 1.0) ::sut/scalar--real
    (sut/type (C/i 1.0)) ::sut/scalar--complex))

(#'sut/promote--domain m--r)

(comment
  (sut/scalar? (C/i 1.0))
  (sut/matrix? m--c))
