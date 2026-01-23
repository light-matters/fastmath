(ns fastmath.algebra.object.matrix.square.complex.ejml
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
   (complexdense/->ComplexDense (ZMatrixRMaj. n o)))

  ;; diagonal
  ([^doubles data]
   (-> (CommonOps_ZDRM/diag (double-array data))
       complexdense/->ComplexDense))

;; elements
  ([^long n ^long o ^doubles data]
   (complexdense/->ComplexDense (ZMatrixRMaj. n o true (double-array (flatten data))))))

(defn identity [^long n]
  (complexdense/->ComplexDense (CommonOps_ZDRM/identity n)))

(defn <-real  [^fastmath.algebra.object.matrix.rectangular.real.ejml.RealDense A]
  (let [M (.-M A)
        Z (ZMatrixRMaj. ^long (.numRows M) ^long (.numCols M))]
    (CommonOps_ZDRM/convert M Z)
    (complexdense/->ComplexDense Z)))

(defn <-rows  [rows]
  (let [nrows (count rows)
        ncols (count (first rows))
        data  (double-array (flatten rows))]
    (println "data: ")
    (println data)
    (complexdense/->ComplexDense (ZMatrixRMaj. nrows ncols true data))))

(defn <-cols  [cols]
  (let [^ZMatrixRMaj A (.M (<-rows cols))
        ^ZMatrixRMaj out (ZMatrixRMaj. (.numCols A) (.numRows A))]
    (CommonOps_ZDRM/transpose A out)

    (complexdense/->ComplexDense out)))
