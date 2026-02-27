(ns fastmath.api.v2.algebra.general-test
  "
  Note: In general, the test expectations were derived from Wolframite. 
  "
  (:require
   [clojure.test :refer [are deftest is]]
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.algebra.object.matrix.create :as mat]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.algebra.test-object :as t]
   [fastmath.api.v2.algebra.general :as sut]
   [fastmath.api.v2.algebra.predicate :as pred]
   [fastmath.protocol.algebra.object.matrix.rectangular :as prot-mat]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;              Constants              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;                Tests                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftest same-shape?-test
  (are [q a] (= q a)
    (pred/same-shape? t/m--r t/m--c) true
    (pred/same-shape? t/m--r t/m--r) true
    (pred/same-shape? t/m--c t/m--c) true
    (pred/same-shape? t/m--c t/m--vc) false
    (pred/same-shape? t/m--c t/m--vrand) false))

(deftest add*-matrix-test
  (are [q a] (= q a)
    (sut/add* t/m--c t/m--r) (mat/<-coll 3 3 [1 1 2 3 4 5 6 7 9 9 10 11 12 13 14 15 17 17])
    (sut/add* t/m--r t/m--c) (mat/<-coll 3 3 [1 1 2 3 4 5 6 7 9 9 10 11 12 13 14 15 17 17]))

  (is (thrown? clojure.lang.ExceptionInfo
               (sut/add* t/m--r t/m--vrand))))

(deftest add*-m-s-test
  (let [cp5 (mat/<-coll 3 3 [5 1 7 3 9 5 11 7 13 9 15 11 17 13 19 15 21 17])
        rp7 (mat/<-coll 3 3 [8 7 7 7 8 7 7 7 8])
        cp7 (mat/<-real rp7)]
    (are [q a] (= q a)
      (sut/add* t/m--c 5.0) cp5
      (sut/add*  5.0 t/m--c) cp5

      (sut/add* t/m--c (C/i 5.0)) cp5
      (sut/add* (C/i 5.0) t/m--c) cp5

      (sut/add* t/m--r 7.0) rp7
      (sut/add*  7.0 t/m--r) rp7

      (sut/add* t/m--r (C/i 7.0)) cp7
      (sut/add*  (C/i 7.0) t/m--r) cp7)))

(deftest add*-s-s
  (are [q a] (= q a)
    (sut/add* 1.0 5) 6.0
    (sut/add* -1.0 57.65432) 56.65432

    (sut/add* 1.0 (C/i 6.0)) (C/i 7 0)
    (sut/add* (C/i 6.0) 1.0) (C/i 7 0)
    (sut/add* (C/i 1.0 5) (C/i 6.0)) (C/i 7 5)))

(deftest +-test
  (are [q a] (= q a)
    (sut/+ t/m--c) t/m--c
    (sut/+ t/m--r) t/m--r
    (sut/+ (C/i)) (C/i)
    (sut/+ 5) 5

    (sut/+ t/m--c t/m--r t/m--r t/m--c)
    (mat/<-coll 3 3 [2 2 4 6 8 10 12 14 18 18 20 22 24 26 28 30 34 34])))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;               Subtract              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftest subtract*-matrix-test
  (is (thrown? clojure.lang.ExceptionInfo
               (sut/subtract* t/m--r t/m--vrand)))
  (are [q a] (= q a)
    (sut/subtract* t/m--c t/m--c) (mat/<-real (mat/zero 3 3))
    (sut/subtract* t/m--r t/m--r) (mat/zero 3 3)

    (sut/subtract* t/m--c t/m--r) (mat/<-coll 3 3 '(-1.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 7.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 15.0 17.0))
    (sut/subtract* t/m--r t/m--c) (prot-mat/negate (mat/<-coll 3 3 '(-1.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 7.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 15.0 17.0)))))

(deftest subtract*-m-s-test
  (let [c5 (mat/<-coll 3 3 [-5.0 1.0 -3.0 3.0 -1.0 5.0 1.0 7.0 3.0 9.0 5.0 11.0 7.0 13.0 9.0 15.0 11.0 17.0])
        r7 (mat/<-coll 3 3 [-6.0 -7.0 -7.0 -7.0 -6.0 -7.0 -7.0 -7.0 -6.0])
        c7 (mat/<-real r7)]
    (are [q a] (= q a)
      (sut/subtract* c5) (prot-mat/negate c5)
      (sut/subtract* t/m--c 5.0) c5
      (sut/subtract*  5.0 t/m--c) (prot-mat/negate c5)

      (sut/subtract* t/m--c (C/i 5.0)) c5
      (sut/subtract* (C/i 5.0) t/m--c) (prot-mat/negate c5)

      (sut/subtract* t/m--r 7.0) r7
      (sut/subtract*  7.0 t/m--r) (prot-mat/negate r7)

      (sut/subtract* t/m--r (C/i 7.0)) c7
      (sut/subtract*  (C/i 7.0) t/m--r) (prot-mat/negate c7))))

(comment (mat/<-real
          (mat/<-coll 3 3 [-6.0 -7.0 -7.0 -7.0 -6.0 -7.0 -7.0 -7.0 -6.0])))

(deftest subtract*-s-s
  (are [q a] (= q a)
    (sut/subtract* 1.0 5) -4.0
    (sut/subtract* -1.0 57.65432) -58.65432

    (sut/subtract* 1.0 (C/i 6.0)) (C/i -5 0)
    (sut/subtract* (C/i 6.0) 1.0) (C/i 5 0)
    (sut/subtract* (C/i 1.0 5) (C/i 6.0)) (C/i -5 5)))

(deftest --test
  (are [q a] (= q a)
    (sut/- t/m--c t/m--c) (prot-mat/zero t/m--c)

    (sut/- t/m--c t/m--r t/m--r t/m--c)
    (mat/<-coll 3 3 [-2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0])

    (sut/- t/m--c t/m--r t/m--r t/m--c -2)
    (mat/<-real (mat/<-coll 3 3 [0 2 2 2 0 2 2 2 0]))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;            Multiplication           ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(deftest multiply*-matrix-test
  ;; (is (thrown? clojure.lang.ExceptionInfo
  ;;              (sut/multiply* t/m--c t/m--vc)))
  (are [q a] (= q a)
    (d2/->vectors (sut/multiply* t/m--c t/m--vc))
    [[[-1.0 0.0] [-3.0 2.0] [-5.0 4.0] [-3.0 2.0] [-5.0 12.0] [-7.0 22.0] [-5.0 4.0] [-7.0 22.0] [-9.0 40.0]] [[-7.0 6.0] [-9.0 32.0] [-11.0 58.0] [-9.0 8.0] [-11.0 42.0] [-13.0 76.0] [-11.0 10.0] [-13.0 52.0] [-15.0 94.0]] [[-13.0 12.0] [-15.0 62.0] [-17.0 112.0] [-15.0 14.0] [-17.0 72.0] [-19.0 130.0] [-17.0 16.0] [-19.0 82.0] [-21.0 148.0]]]

    (sut/multiply* (mat/<-coll 3 1 [4 5 6])
                   (mat/<-coll 1 3 [1 2 3]))
    (mat/<-coll 3 3 [4 8 12 5 10 15 6 12 18]))

  (are [q a] (= q a)
    (sut/multiply* t/m--c t/m--c) (mat/<-coll 3 3 [-27.0 144.0 -33.0 174.0 -39.0 204.0 -45.0 378.0 -51.0 480.0 -57.0 582.0 -63.0 612.0 -69.0 786.0 -75.0 960.0])
    (sut/multiply* t/m--r t/m--r) (mat/identity 3)
    (sut/multiply* t/m--c t/m--r) t/m--c
    (sut/multiply* t/m--r t/m--c) t/m--c))

(deftest multiply*-m-s-test
  (let [mc5 (mat/<-coll 3 3 [0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0 80.0 85.0])
        mr7 (mat/diagonal (repeat 3 7))]
    (are [a q] (= q a)
      (sut/multiply* t/m--c 5.0) mc5
      (sut/multiply*  5.0 t/m--c) mc5

      (sut/multiply* t/m--c (C/i 5.0)) mc5
      (sut/multiply* (C/i 5.0) t/m--c) mc5

      (sut/multiply* t/m--r 7.0) mr7
      (sut/multiply*  7.0 t/m--r) mr7

      (sut/multiply* t/m--r (C/i 7.0)) (mat/<-real mr7)
      (sut/multiply*  (C/i 7.0) t/m--r) (mat/<-real mr7))

    ;; TODO: outer products
    ))

(deftest *-test
  (are [q a] (= q a)
    (sut/* t/m--c t/m--c) (mat/<-coll 3 3 [-27.0 144.0 -33.0 174.0 -39.0 204.0 -45.0 378.0 -51.0 480.0 -57.0 582.0 -63.0 612.0 -69.0 786.0 -75.0 960.0])

    (sut/* t/m--c t/m--r t/m--r t/m--c)
    (mat/<-coll 3 3 [-27.0 144.0 -33.0 174.0 -39.0 204.0 -45.0 378.0 -51.0 480.0 -57.0 582.0 -63.0 612.0 -69.0 786.0 -75.0 960.0])

    (sut/* t/m--c t/m--r t/m--r t/m--c -2)
    (mat/<-coll 3 3 [54.0 -288.0 66.0 -348.0 78.0 -408.0 90.0 -756.0 102.0 -960.0 114.0 -1164.0 126.0 -1224.0 138.0 -1572.0 150.0 -1920.0])))

