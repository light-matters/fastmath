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
  (mat/<-coll 3 3 [1 0 0 0 1 0 0 0 1]))

(def m--r-promoted
  (mat/<-coll 3 3 (interleave [1 0 0 0 1 0 0 0 1]
                              (repeat 9 0))))

(def cn (C/i 1.0 -3.0))

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

(deftest ensure-domain-match-test
  (are [q a] (= q a)
    (#'sut/ensure-domain-match cn 1.0) [(C/i 1 -3) (C/i 1.0)]

    (#'sut/ensure-domain-match m--c m--r)
    [m--c (mat/<-rows
           [[(C/i 1) (C/i 0) (C/i 0)]
            [(C/i 0) (C/i 1) (C/i 0)]
            [(C/i 0) (C/i 0) (C/i 1)]])]

    (#'sut/ensure-domain-match m--c m--c)
    [m--c m--c]

    (#'sut/ensure-domain-match m--c (C/i))
    [m--c (C/i)]

    (#'sut/ensure-domain-match m--r (C/i))
    [m--r-promoted (C/i)]

    (#'sut/ensure-domain-match m--r 1)
    [m--r 1]))

(deftest add*-matrix-test
  (are [q a] (= q a)
    (sut/add* m--c m--r) (mat/<-coll 3 3 [1 1 2 3 4 5 6 7 9 9 10 11 12 13 14 15 17 17])
    (sut/add* m--r m--c) (mat/<-coll 3 3 [1 1 2 3 4 5 6 7 9 9 10 11 12 13 14 15 17 17])))

(deftest add*-m-s-test
  (let [cp5 (mat/<-coll 3 3 [5 1 7 3 9 5 11 7 13 9 15 11 17 13 19 15 21 17])
        rp7 (mat/<-coll 3 3 [8 7 7 7 8 7 7 7 8])
        cp7 (mat/<-real rp7)]
    (are [q a] (= q a)
      (sut/add* m--c 5.0) cp5
      (sut/add*  5.0 m--c) cp5

      (sut/add* m--c (C/i 5.0)) cp5
      (sut/add* (C/i 5.0) m--c) cp5

      (sut/add* m--r 7.0) rp7
      (sut/add*  7.0 m--r) rp7

      (sut/add* m--r (C/i 7.0)) cp7
      (sut/add*  (C/i 7.0) m--r) cp7)))

(deftest add*-s-s
  (are [q a] (= q a)
    (sut/add* 1.0 5) 6.0
    (sut/add* -1.0 57.65432) 56.65432

    (sut/add* 1.0 (C/i 6.0)) (C/i 7 0)
    (sut/add* (C/i 6.0) 1.0) (C/i 7 0)
    (sut/add* (C/i 1.0 5) (C/i 6.0)) (C/i 7 5)))

(deftest +-test
  (are [q a] (= q a)
    (sut/+ m--c m--r m--r m--c)
    (mat/<-coll 3 3 [2 2 4 6 8 10 12 14 18 18 20 22 24 26 28 30 34 34])))

(comment
  (sut/scalar? (C/i 1.0))
  (sut/matrix? m--c))
