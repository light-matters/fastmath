(ns fastmath.matrix.dense.complex.ojalg
  "Implementing complex matrices using ojAlgo as a backend [https://github.com/optimatika/ojAlgo/tree/bd85bf0fe485016817fad945a5975eba8164b369?tab=readme-ov-file].

[https://javadoc.io/doc/org.ojalgo/ojalgo/latest/ojalgo/module-summary.html]
   TODO:
   - Understand the API!

   Legend:
   m - matrix type
   i - row index
   j - column index
  "

  (:import
   (java.lang.reflect Array)
   (org.ojalgo.matrix.store GenericStore PhysicalStore)
   (org.ojalgo.scalar ComplexNumber)
   (org.ojalgo.structure Access1D)
   (org.ojalgo.array Array1D)
   (org.ojalgo.matrix
    MatrixR064
    MatrixC128)
   (org.ejml.dense.row
    CommonOps_DDRM)
   (org.ejml.data DMatrixRMaj)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; Type hint the return value
(defn complex ^ComplexNumber
  ([real] (ComplexNumber/of (double real) 0.0))
  ([real imag] (ComplexNumber/of (double real) (double imag))))

;; Create a mutable complex matrix
(defn create-complex-matrix ^PhysicalStore [^long rows ^long cols]
  (.make (GenericStore/C128) rows cols))

;; Example usage
(defn demo []
  (let [^PhysicalStore A (create-complex-matrix 2 2)]
    (.set A (long 0) (long 0) ^ComplexNumber (complex 1.0 0.0))
    (.set A (long 0) (long 1) ^ComplexNumber (complex 0.0 1.0))
    (.set A (long 1) (long 0) ^ComplexNumber (complex 2.0 3.0))
    (.set A (long 1) (long 1) ^ComplexNumber (complex -1.0 -2.0))
    A))

;; Print function - use double cast instead of .doubleValue
(defn print-complex-matrix [^PhysicalStore matrix]
  (dotimes [i (.countRows matrix)]
    (dotimes [j (.countColumns matrix)]
      (let [^ComplexNumber c (.get matrix (long i) (long j))
            real (double (.getReal c))
            imag (double (.getImaginary c))]
        (print (format "%7.2f%+.2fi  " real imag))))
    (println)))

;; Matrix operations
(defn matrix-multiply ^PhysicalStore [^PhysicalStore m1 ^PhysicalStore m2]
  (.multiply m1 m2))

(defn matrix-add ^PhysicalStore [^PhysicalStore m1 ^PhysicalStore m2]
  (.add m1 m2))

(defn matrix-conjugate ^PhysicalStore [^PhysicalStore m]
  (.conjugate m))

(defn matrix-transpose ^PhysicalStore [^PhysicalStore m]
  (.transpose m))

;; Get element
(defn get-element ^ComplexNumber [^PhysicalStore matrix ^long row ^long col]
  (.get matrix row col))

;; Create from nested vectors
(defn matrix-from-vectors ^PhysicalStore [vectors]
  (let [rows (long (count vectors))
        cols (long (count (first vectors)))
        ^PhysicalStore matrix (create-complex-matrix rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (let [[real imag] (get-in vectors [i j])]
          (.set matrix (long i) (long j) ^ComplexNumber (complex real imag)))))
    matrix))

;; Create identity matrix
(defn identity-matrix ^PhysicalStore [^long n]
  (.makeIdentity (GenericStore/C128) n))

;; Convert matrix to Clojure data
(defn matrix->vectors [^PhysicalStore matrix]
  (let [rows (.countRows matrix)
        cols (.countColumns matrix)]
    (vec (for [i (range rows)]
           (vec (for [j (range cols)]
                  (let [^ComplexNumber c (.get matrix (long i) (long j))]
                    [(double (.getReal c))
                     (double (.getImaginary c))])))))))

;; Helper to extract real and imaginary parts
(defn complex->pair [^ComplexNumber c]
  [(double (.getReal c)) (double (.getImaginary c))])

;; Examples
(comment
  ;; Create and print a matrix
  (def A (demo))
  (print-complex-matrix A)

  ;; Create from vectors
  (def B (matrix-from-vectors
          [[[1 0] [0 -1]]
           [[0 1] [1 0]]]))
  (print-complex-matrix B)

  ;; Matrix multiplication
  (def C (matrix-multiply A B))
  (print-complex-matrix C)

  ;; Identity matrix
  (def I (identity-matrix 3))
  (print-complex-matrix I)

  ;; Conjugate
  (def A-conj (matrix-conjugate A))
  (print-complex-matrix A-conj)

  ;; Get individual element
  (let [elem (get-element A 0 1)
        [r i] (complex->pair elem)]
    (println "Element at (0,1):" r "+" i "i"))

  ;; Convert to Clojure data
  (matrix->vectors A))
