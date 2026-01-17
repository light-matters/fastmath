(ns fastmath.api.v2.algebra.general-test
  (:require [fastmath.api.v2.algebra.general :as sut]
            [fastmath.algebra.object.matrix.create :as mat]
            [fastmath.algebra.object.number.complex.create :as cn]
            [clojure.test :refer [deftest is are]]))

(def m--c
  (mat/<-coll 3 3 (partition 2 (range 18))))
(def m--vc
  (mat/<-coll 1 3 (partition 2 (range 6))))
(def m--vr
  (let [n (rand-int 10)]
    (mat/<-coll n 1 (range (* 2 n)))))

(deftest linear-shape?-test
  (are [q a] (= q a)
    (sut/linear-shape? m--c) false
    (sut/linear-shape? m--vc) true
    (sut/linear-shape? m--vr) true))





