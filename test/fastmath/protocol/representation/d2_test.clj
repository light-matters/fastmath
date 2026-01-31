(ns fastmath.protocol.representation.d2-test
  (:require [fastmath.protocol.representation.d2 :as sut]
            [fastmath.protocol.algebra.object.matrix.rectangular :as rmat] [fastmath.algebra.object.number.complex.create :as C]
            [fastmath.algebra.test-object :as to :refer [m--c, m--r m--vc m--vrand]]
            [clojure.test :as t :refer [deftest are testing]]
            [fastmath.algebra.object.matrix.create :as mat]))

(deftest information
  (are [e a] (= a e)
    [3 3] (sut/shape m--c)
    [3 3] (sut/shape m--r)
    [1 3] (sut/shape m--vc)
    [(count (sut/rows m--vrand))
     (count (sut/columns m--vrand))]
    (sut/shape m--vrand)))

(deftest retrieval
  (are [a e] (= a e)
    (C/i 0 1) (sut/element m--c 0 0)
    0.0 (sut/element m--r 1 2)
    (C/i 2 3) (sut/element m--vc 0 1)
    ;; column and row
    (mat/<-coll 3 1 [2 3 8 9 14 15]) (sut/column m--c 1)
    (mat/<-coll 1 3 [12 13 14 15 16 17]) (sut/row m--c 2)

    ;; columns and rows
    [(mat/<-coll 3 1 [0 1 6 7 12 13])
     (mat/<-coll 3 1 [2 3 8 9 14 15])
     (mat/<-coll 3 1 [4 5 10 11 16 17])]
    (sut/columns m--c)

    (mapv #(mat/<-coll 1 3 %)
          (partition (* 2 3) (range 18)))
    (sut/rows m--c)))

(deftest transformation
  (are [a e] (= a e)
    ;; ->array
    ;; (sut/->array m--c)

    ;; map

    (mat/<-coll 3 3
                (flatten (map (fn [[x y]] [(inc x) y])
                              (partition 2 (range 18)))))

    (sut/fmap m--c #(apply rmat/add [%1 1]))))

(comment
  m--c
  (sut/fmap m--c #(rmat/add % (C/i 100.0))))
