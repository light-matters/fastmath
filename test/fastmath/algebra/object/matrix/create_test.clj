(ns fastmath.algebra.object.matrix.create-test
  (:require
   [clojure.test :refer [deftest is]]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as cejml]

   [fastmath.algebra.object.matrix.create :as sut])
  (:import
   (org.ejml.data ZMatrixRMaj Complex_F64)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;              Utilities              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
    (println [nrows ncols])
    (cejml/->ComplexDense (ZMatrixRMaj. nrows ncols true data))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;                                   Tests                                                                       ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftest identity-test
  (is (= (sut/identity 2) (cdense<-vs [[[1.0 0.0] [0.0 0.0]]
                                       [[0.0 0.0] [1.0 0.0]]]))))

(comment
  (println (cejml/->ComplexDense (ZMatrixRMaj. 2 2 false (double-array [0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0])))))
