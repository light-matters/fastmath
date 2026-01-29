(ns fastmath.algebra.plumb-test
  (:require
   [clojure.test :as t :refer [are deftest]]
   [fastmath.algebra.object.matrix.create :as mat]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.algebra.plumb :as sut]))

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

(deftest ensure-domain-match-test
  (are [q a] (= q a)
    (sut/ensure-domain-match (C/i 1 -3) 1.0)
    [(C/i 1 -3) (C/i 1.0)]

    (sut/ensure-domain-match m--c m--r)
    [m--c (mat/<-rows
           [[(C/i 1) (C/i 0) (C/i 0)]
            [(C/i 0) (C/i 1) (C/i 0)]
            [(C/i 0) (C/i 0) (C/i 1)]])]

    (sut/ensure-domain-match m--c m--c)
    [m--c m--c]

    (sut/ensure-domain-match m--c (C/i))
    [m--c (C/i)]

    (sut/ensure-domain-match m--r (C/i))
    [m--r-promoted (C/i)]

    (sut/ensure-domain-match m--r 1)
    [m--r 1]))
