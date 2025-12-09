(ns fastmath.algebra.object.matrix.square.real.ejml
  (:require
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as realdense])
  (:import
   (org.ejml.data DMatrixRMaj)
   (org.ejml.dense.row CommonOps_DDRM)))
;; TODO: WIP placeholder for the future squarematrix type. 

(defn realdense

  ;; zero
  ([^long n ^long o]
   (realdense/->RealDense (DMatrixRMaj. n o)))

  ;; diagonal
  ([^doubles data]
   (-> (CommonOps_DDRM/diag (double-array data))
       realdense/->RealDense))

;; elements
  ([^long n ^long o ^doubles data]
   (realdense/->RealDense (DMatrixRMaj. n o false data))))

(defn identity [^long n]
  (realdense/->RealDense (CommonOps_DDRM/identity n)))

(defn <-rows  [rows]
  (let [nrows (count rows)
        ncols (count (first rows))
        data  (double-array (apply concat rows))]
    (realdense/->RealDense (DMatrixRMaj. nrows ncols true data))))

(defn <-cols  [cols]
  (let [^DMatrixRMaj A (.M (<-rows cols))
        ^DMatrixRMaj out (DMatrixRMaj. (.numCols A) (.numRows A))]
    (CommonOps_DDRM/transpose A out)

    (realdense/->RealDense out)))
