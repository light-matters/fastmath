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
  [rows]
  (let [[nrows ncols] [(count rows) (count (first rows))]
        data (->> rows
       ;; (into [] (comp cat))
                  flatten
                  double-array)]
    (rejml/->RealDense (DMatrixRMaj. nrows ncols true data))))
(defn- cdense<-vs
  "A 1D Java array of complex numbers (real and imaginary separated) from a 3D Clojure vector system of double-valued pairs, e.g.

[[[1 0] [2 1]]
 [[3 -1] [4 0]]].
  "
  ;; TODO: Add something like this to the official create method
  [rows]
  (let [[nrows ncols] [(count rows) (count (first rows))]
        data (->> rows
       ;; (into [] (comp cat))
                  flatten
                  double-array)]
    (cejml/->ComplexDense (ZMatrixRMaj. nrows ncols true data))))

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
    (rdense<-vs [[10.0 1.0  -5.0]
                 [0.0 4.0  2.789]])

    (sut/<-coll 2 2 [(i 78 0.0) (i 0.0 77) (i 9134 -341) (i 24 2341)])
    (cdense<-vs [[[78 0.0] [0.0 77]]
                 [[9134 -341] [24 2341]]])))

(comment
  (-> (sut/zero 2 3)
      type)

  (flatten [(i 78 0.0) (i 0.0 77) (i 9134 -341) (i 24 2341)])

  (sut/<-coll 2 2 [(i 78 0.0) (i 0.0 77) (i 9134 -341) (i 24 2341)])

  (println (-> (ZMatrixRMaj. 2 2 true)
               cejml/->ComplexDense))

  (-> (double-array [10.0 1.0 -5.0 0.0 4.0 2.789])
      println))

