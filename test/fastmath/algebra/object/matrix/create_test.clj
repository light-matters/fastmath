(ns fastmath.algebra.object.matrix.create-test
  (:require
   [clojure.test :refer [deftest is are]]
   [fastmath.algebra.object.number.complex.ejml :as C :refer [i]]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as rejml]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as cejml]

   [fastmath.algebra.object.matrix.create :as sut])
  (:import
   (org.ejml.data DMatrixRMaj ZMatrixRMaj Complex_F64)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;              Utilities              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn- rdense<-vs
  "A 1D Java array of complex numbers (real and imaginary separated) from a 3D Clojure vector system of double-valued pairs, e.g.

[[[1 0] [2 1]]
 [[3 -1] [4 0]]].
  "
  ;; TODO: Add something like this to the official create method
  ([rows]
   (rdense<-vs rows true))
  ([rows rows?]
   (let [sh [(count rows) (count (first rows))]
         [nrows ncols] (if rows? sh (reverse sh))
         data (->> rows
       ;; (into [] (comp cat))
                   flatten
                   double-array)]
     (rejml/->RealDense (DMatrixRMaj. nrows ncols rows? data)))))

(defn- cdense<-vs
  "A 1D Java array of complex numbers (real and imaginary separated) from a 3D Clojure vector system of double-valued pairs, e.g.

[[[1 0] [2 1]]
 [[3 -1] [4 0]]].
  "
  ;; TODO: Add something like this to the official create method
  ([rows]
   (cdense<-vs rows true))
  ([rows rows?]
   (let [sh [(count rows) (count (first rows))]
         [nrows ncols] (if rows? sh (reverse sh))
         data (->> rows
       ;; (into [] (comp cat))
                   flatten
                   double-array)]
     (cejml/->ComplexDense (ZMatrixRMaj. nrows ncols rows? data)))))

(def RI2
  (rdense<-vs [[1.0 0.0]
               [0.0 1.0]]))
(def CI2
  (cdense<-vs [[[1.0 0.0] [0.0 0.0]]
               [[0.0 0.0] [1.0 0.0]]]))
(def real2x3
  (rdense<-vs [[10.0 1.0  -5.0]
               [0.0 4.0  2.789]]))
(def real3x2
  (rdense<-vs [[10.0 1.0  -5.0]
               [0.0 4.0  2.789]]
              false))

(def real5x5
  (rdense<-vs [[-1.1235 0.0  0.0 0 0.0]
               [0.0 4.0  0.0 0.0 0.0]
               [0 0 8 0 0]
               [0  0 0 23 0]
               [0 0 0 0 5e-3]]))

(def real5x3
  (rdense<-vs [[10.0 1.0  -5.0]
               [0.0 4.0  2.789]
               [-5 31 8]
               [5e10 23 -18e11]
               [2 13 56]]))

(def real3x5
  (rdense<-vs [[10.0 1.0  -5.0]
               [0.0 4.0  2.789]
               [-5 31 8]
               [5e10 23 -18e11]
               [2 13 56]] false))

(def complex2x3
  (cdense<-vs [[(i 78 0.0) (i 0.0 77) (i 16.13456 56)]
               [(i 9134 -341) (i 24 2341) (i 10 -56)]]))
(def complex3x2
  (cdense<-vs [[(i 78 0.0) (i 0.0 77) (i 16.13456 56)]
               [(i 9134 -341) (i 24 2341) (i 10 -56)]]
              false))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;                                   Tests                                                                       ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(deftest identity-test
  (are [q a] (= q a)
    (sut/identity 2 {:domain :complex})
    (cdense<-vs [[[1.0 0.0] [0.0 0.0]]
                 [[0.0 0.0] [1.0 0.0]]])

    (sut/identity 3 {:domain :real})
    (rdense<-vs [[1.0 0.0 0.0]
                 [0.0 1.0 0.0]
                 [0.0 0.0 1.0]])))

(deftest zero-test
  (are [q a] (= q a)
    (sut/zero 2)
    (rdense<-vs [[0.0  0.0]
                 [0.0  0.0]])
    (sut/zero 2 3)
    (rdense<-vs [[0.0 0.0  0.0]
                 [0.0 0.0  0.0]])

    (sut/zero 2 2 {:domain :complex})
    (cdense<-vs [[[0.0 0.0]  [0.0 0.0]]
                 [[0.0 0.0]  [0.0 0.0]]])))

(deftest <-coll-test
  (are [q a] (= q a)
    (sut/<-coll 2 3 [10.0 1.0 -5.0 0.0 4.0 2.789])
    real2x3

    (sut/<-coll 2 2 [(i 78 0.0) (i 0.0 77) (i 9134 -341) (i 24 2341)])
    (cdense<-vs [[[78 0.0] [0.0 77]]
                 [[9134 -341] [24 2341]]])))

(deftest <-rows-test
  (are [q a] (= q a)
    (sut/<-rows [[10.0 1.0 -5.0] [0.0 4.0 2.789]])
    real2x3

    (sut/<-rows [[(i 78 0.0) (i 0.0 77) (i 16.13456 56)]
                 [(i 9134 -341) (i 24 2341) (i 10 -56)]])
    complex2x3

    (sut/<-rows [[10.0 1.0  -5.0]
                 [0.0 4.0  2.789]
                 [-5 31 8]
                 [5e10 23 -18e11]
                 [2 13 56]])
    real5x3))

(deftest <-cols-test
  (are [q a] (= q a)
    (sut/<-cols [[10.0 1.0 -5.0] [0.0 4.0 2.789]])
    real3x2

    (sut/<-cols [[(i 78 0.0) (i 0.0 77) (i 16.13456 56)]
                 [(i 9134 -341) (i 24 2341) (i 10 -56)]])
    complex3x2

    (sut/<-cols [[10.0 1.0  -5.0]
                 [0.0 4.0  2.789]
                 [-5 31 8]
                 [5e10 23 -18e11]
                 [2 13 56]])
    real3x5))

(deftest <-diagonal-test
  (are [q a] (= q a)
    (sut/diagonal [1.0 1.0])
    RI2)

  (sut/diagonal [(i 1.0 0.0) (i 1.0 0.0)])
  CI2

  (sut/diagonal [-1.1235 4.0 8 23 5e-3])
  real5x5)

(deftest <-real-test
  (are [q a] (= q a)
    (sut/<-real RI2)
    CI2

    (sut/<-real real3x2)
    (cdense<-vs [[[10.0 0.0] [1.0 0.0]  [-5.0 0.0]]
                 [[0.0 0.0] [4.0 0.0]  [2.789 0.0]]]
                false)))

(println CI2)
(println complex3x2)

(comment
  (-> (sut/zero 2 3)
      type)

  (flatten [(i 78 0.0) (i 0.0 77) (i 9134 -341) (i 24 2341)])

  (sut/<-coll 2 2 [(i 78 0.0) (i 0.0 77) (i 9134 -341) (i 24 2341)])

  (println (-> (ZMatrixRMaj. 2 2 true)
               cejml/->ComplexDense))

  (-> (double-array [10.0 1.0 -5.0 0.0 4.0 2.789])
      println))

