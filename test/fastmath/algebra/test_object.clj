(ns fastmath.algebra.test-object
  "For holding useful entities for comparison."
  (:require
   [fastmath.algebra.object.matrix.create :as mat]
   [clojure.test :as t]))

(def m--c
  (mat/<-coll 3 3 (partition 2 (range 18))))
;; :shape [3 3] :type complex128
;; [[0.000+1.000000i   2.000+3.000000i   4.000+5.000000i  ]
;;  [6.000+7.000000i   8.000+9.000000i   10.000+11.000000i]
;;  [12.000+13.000000i 14.000+15.000000i 16.000+17.000000i]]

(def m--vc
  (mat/<-coll 1 3 (partition 2 (range 6))))
;; :shape [1 3] :type complex128
;; [[0.000+1.000000i 2.000+3.000000i 4.000+5.000000i]]

(def m--vrand
  (let [n (+ 1 (rand-int 10))]
    (mat/<-coll n 1 (range (* 2 n)))))
;; Something like...
;; :shape [5 1] :type complex128
;; [[0.000+1.000000i]
;;  [2.000+3.000000i]
;;  [4.000+5.000000i]
;;  [6.000+7.000000i]
;;  [8.000+9.000000i]]

(def m--r
  (mat/<-coll 3 3 [1 0 0 0 1 0 0 0 1]))
;; :shape [3 3] :type float64
;; [[1.000 0.000 0.000]
;;  [0.000 1.000 0.000]
;;  [0.000 0.000 1.000]]

(def m--r-promoted
  (mat/<-coll 3 3 (interleave [1 0 0 0 1 0 0 0 1]
                              (repeat 9 0))))
;; :shape [3 3] :type float64
;; [[1.000 0.000 0.000]
;;  [0.000 1.000 0.000]
;;  [0.000 0.000 1.000]]

