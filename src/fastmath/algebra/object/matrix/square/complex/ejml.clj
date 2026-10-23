(ns fastmath.algebra.object.matrix.square.complex.ejml
  (:refer-clojure :exclude [identity])
  (:require
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as complexdense])
  (:import
   (java.lang Math)
   (org.ejml.data ZMatrixRMaj)
   (org.ejml.dense.row CommonOps_ZDRM)))
;; TODO: WIP placeholder for the future squarematrix type. 

(defn complexdense
  ;; zero
  ([^long n ^long o]
   (complexdense/complexdense n o))
  ;; diagonal, interleaved
  ([data]
   (complexdense/complexdense data))
  ;; elements, row-major interleaved
  ([^long n ^long o data]
   (complexdense/complexdense n o data)))

(defn identity [^long n]
  (complexdense/->ComplexDense (CommonOps_ZDRM/identity n)))

(defn <-real  [A]
  (let [M (.-M A)
        Z (ZMatrixRMaj. ^long (.numRows M) ^long (.numCols M))]
    (CommonOps_ZDRM/convert M Z)
    (complexdense/->ComplexDense Z)))

(defn <-rows  [rows]
  (complexdense/<-rows rows))

(defn <-cols  [cols]
  (let [^ZMatrixRMaj A (.M (<-rows cols))
        ^ZMatrixRMaj out (ZMatrixRMaj. (.numCols A) (.numRows A))]
    (CommonOps_ZDRM/transpose A out)

    (complexdense/->ComplexDense out)))
