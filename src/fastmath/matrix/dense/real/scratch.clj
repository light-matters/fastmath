(ns fastmath.matrix.dense.real.scratch
  (:import
   (org.ejml.dense.row
    CommonOps_DDRM)
   (org.ejml.data DMatrixRMaj)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def mat (DMatrixRMaj. 4 4 true (double-array (flatten [[1 2 3 10] [4 5 6 11] [7 8 9 12] (range 13 17)]))))

(defn ^DMatrixRMaj extract
  "Extracts a submatrix m[i0:i1, j0:j1) into a new DMatrixRMaj.
   Indices are zero-based and end-exclusive."
  [^DMatrixRMaj m
   ^longs [i0 i1]
   ^longs [j0 j1]]
  (let [^DMatrixRMaj out (DMatrixRMaj. (- i1 i0) (- j1 j0))]
    (CommonOps_DDRM/extract m i0 i1 j0 j1 out 0 0)
    out))

(println mat)
(-> mat
    (extract [0 4] [0 1])
    println)
