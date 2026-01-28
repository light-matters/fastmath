(ns fastmath.api.v2.algebra.general-test
  "
  Note: In general, the test expectations were derived from Wolframite. 
  "
  (:require
   [clojure.test :refer [are deftest is]]
   [fastmath.protocol.algebra.object.matrix.rectangular :as prot-mat]
   [fastmath.algebra.object.matrix.create :as mat]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.api.v2.algebra.general :as sut]))

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

(deftest same-shape?-test
  (are [q a] (= q a)
    (#'sut/same-shape? m--r m--c) true
    (#'sut/same-shape? m--r m--r) true
    (#'sut/same-shape? m--c m--c) true
    (#'sut/same-shape? m--c m--vc) false
    (#'sut/same-shape? m--c m--vrand) false))

(deftest add*-matrix-test
  (are [q a] (= q a)
    (sut/add* m--c m--r) (mat/<-coll 3 3 [1 1 2 3 4 5 6 7 9 9 10 11 12 13 14 15 17 17])
    (sut/add* m--r m--c) (mat/<-coll 3 3 [1 1 2 3 4 5 6 7 9 9 10 11 12 13 14 15 17 17]))

  (is (thrown? clojure.lang.ExceptionInfo
               (sut/add* m--r m--vrand))))

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
    (sut/+ m--c) m--c
    (sut/+ m--r) m--r
    (sut/+ (C/i)) (C/i)
    (sut/+ 5) 5

    (sut/+ m--c m--r m--r m--c)
    (mat/<-coll 3 3 [2 2 4 6 8 10 12 14 18 18 20 22 24 26 28 30 34 34])))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;               Subtract              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftest subtract*-matrix-test
  (is (thrown? clojure.lang.ExceptionInfo
               (sut/subtract* m--r m--vrand)))
  (are [q a] (= q a)
    (sut/subtract* m--c m--c) (mat/<-real (mat/zero 3 3))
    (sut/subtract* m--r m--r) (mat/zero 3 3)

    (sut/subtract* m--c m--r) (mat/<-coll 3 3 '(-1.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 7.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 15.0 17.0))
    (sut/subtract* m--r m--c) (prot-mat/negate (mat/<-coll 3 3 '(-1.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 7.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 15.0 17.0)))))

(deftest subtract*-m-s-test
  (let [c5 (mat/<-coll 3 3 [-5.0 1.0 -3.0 3.0 -1.0 5.0 1.0 7.0 3.0 9.0 5.0 11.0 7.0 13.0 9.0 15.0 11.0 17.0])
        r7 (mat/<-coll 3 3 [-6.0 -7.0 -7.0 -7.0 -6.0 -7.0 -7.0 -7.0 -6.0])
        c7 (mat/<-real r7)]
    (are [q a] (= q a)
      (sut/subtract* c5) (prot-mat/negate c5)
      (sut/subtract* m--c 5.0) c5
      (sut/subtract*  5.0 m--c) (prot-mat/negate c5)

      (sut/subtract* m--c (C/i 5.0)) c5
      (sut/subtract* (C/i 5.0) m--c) (prot-mat/negate c5)

      (sut/subtract* m--r 7.0) r7
      (sut/subtract*  7.0 m--r) (prot-mat/negate r7)

      (sut/subtract* m--r (C/i 7.0)) c7
      (sut/subtract*  (C/i 7.0) m--r) (prot-mat/negate c7))))

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
    (sut/- m--c m--c) (prot-mat/zero m--c)

    (sut/- m--c m--r m--r m--c)
    (mat/<-coll 3 3 [-2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0])

    (sut/- m--c m--r m--r m--c -2)
    (mat/<-real (mat/<-coll 3 3 [0 2 2 2 0 2 2 2 0]))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;            Multiplication           ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(deftest multiply*-matrix-test
  (is (thrown? clojure.lang.ExceptionInfo
               (sut/multiply* m--c m--vc)))
  (are [q a] (= q a)
    (sut/multiply* m--c m--c) (mat/<-coll 3 3 [-27.0 144.0 -33.0 174.0 -39.0 204.0 -45.0 378.0 -51.0 480.0 -57.0 582.0 -63.0 612.0 -69.0 786.0 -75.0 960.0])
    (sut/multiply* m--r m--r) (mat/identity 3)
    (sut/multiply* m--c m--r) m--c
    (sut/multiply* m--r m--c) m--c))

(deftest multiply*-m-s-test
  (let [mc5 (mat/<-coll 3 3 [0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0 80.0 85.0])
        mr7 (mat/diagonal (repeat 3 7))]
    (are [a q] (= q a)
      (sut/multiply* m--c 5.0) mc5
      (sut/multiply*  5.0 m--c) mc5

      (sut/multiply* m--c (C/i 5.0)) mc5
      (sut/multiply* (C/i 5.0) m--c) mc5

      (sut/multiply* m--r 7.0) mr7
      (sut/multiply*  7.0 m--r) mr7

      (sut/multiply* m--r (C/i 7.0)) (mat/<-real mr7)
      (sut/multiply*  (C/i 7.0) m--r) (mat/<-real mr7))))

(deftest *-test
  (are [q a] (= q a)
    (sut/* m--c m--c) (mat/<-coll 3 3 [-27.0 144.0 -33.0 174.0 -39.0 204.0 -45.0 378.0 -51.0 480.0 -57.0 582.0 -63.0 612.0 -69.0 786.0 -75.0 960.0])

    (sut/* m--c m--r m--r m--c)
    (mat/<-coll 3 3 [-27.0 144.0 -33.0 174.0 -39.0 204.0 -45.0 378.0 -51.0 480.0 -57.0 582.0 -63.0 612.0 -69.0 786.0 -75.0 960.0])

    (sut/* m--c m--r m--r m--c -2)
    (mat/<-coll 3 3 [54.0 -288.0 66.0 -348.0 78.0 -408.0 90.0 -756.0 102.0 -960.0 114.0 -1164.0 126.0 -1224.0 138.0 -1572.0 150.0 -1920.0])))
