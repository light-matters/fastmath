(ns fastmath.algebra.test-object
  "For holding useful entities for comparison."
  (:require
   [fastmath.algebra.object.matrix.create :as mat]
   [clojure.test :as t]))

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
