(ns fastmath.matrix.dense.scratch
  "A playground for working through ideas.

  Currently playing with the complex implementation of ojalgo matrices.

https://gist.github.com/apete/b3278dc2f8c2db6a00369c211ba321db
  Reading through this ^^^ (CommonMistake.java) helped to get an idea of the API design. It's stil pretty implicit though and I'm sure I'm not using it right.
  "
  (:import
   (org.ojalgo.matrix BasicMatrix MatrixC128 MatrixC128$DenseReceiver)
   (org.ojalgo.matrix.store GenericStore)
   (org.ojalgo.scalar ComplexNumber)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(println "===start===")

;; ============================================================================
;; Complex number helpers
;; ============================================================================

(defn ^ComplexNumber complex
  "Create a complex number."
  ([real] (ComplexNumber/of (double real) 0.0))
  ([real imag] (ComplexNumber/of (double real) (double imag))))

(defn complex->pair
  "Extract [real imag] from a ComplexNumber."
  [^ComplexNumber c]
  [(double (.getReal c)) (double (.getImaginary c))])

;; ============================================================================
;; Efficient matrix creation from arrays
;; ============================================================================

(defn ^GenericStore create-store-from-array
  "Create a GenericStore from a flat array of ComplexNumbers.
  Data is in row-major order.
  
  Examples:
  (create-store-from-array 2 2 [(complex 1 0) (complex 0 1)
                                 (complex 2 3) (complex -1 -2)])
  Creates:
  [1+0i   0+1i  ]
  [2+3i  -1-2i  ]"
  [^long rows ^long cols data]
  (let [arr (into-array ComplexNumber data)
        factory (GenericStore/C128)]
    (.rows factory (into-array (Class/forName "[Lorg.ojalgo.scalar.ComplexNumber;")
                               [(into-array ComplexNumber (take cols arr))
                                (into-array ComplexNumber (drop cols arr))]))))

;; Actually, let's use a simpler approach - just create and populate
(defn ^GenericStore matrix-from-vectors-fast
  "Fast creation of GenericStore from nested vectors of [real imag] pairs."
  [vectors]
  (let [rows (long (count vectors))
        cols (long (count (first vectors)))
        store (create-store rows cols)]
    ;; Populate using array copying would be faster, but this is still reasonably fast
    (dotimes [i rows]
      (dotimes [j cols]
        (let [[real imag] (get-in vectors [i j])
              c (complex real imag)]
          (.set ^GenericStore store (long i) (long j) ^ComplexNumber c))))
    store))

(defn ^MatrixC128 matrix-from-vectors
  "Create a MatrixC128 from nested vectors of [real imag] pairs."
  [vectors]
  (store->matrix (matrix-from-vectors-fast vectors)))

;; ============================================================================
;; Alternative: Use factory's rows/columns methods
;; ============================================================================

(defn ^MatrixC128 matrix-from-rows
  "Create a MatrixC128 from rows of ComplexNumbers.
  Each row is a sequence of ComplexNumber objects."
  [rows-data]
  (let [factory MatrixC128/FACTORY
        rows-arrays (into-array (Class/forName "[Lorg.ojalgo.scalar.ComplexNumber;")
                                (map #(into-array ComplexNumber %) rows-data))]
    (.rows factory rows-arrays)))

(defn ^MatrixC128 matrix-from-vectors-via-rows
  "Create a MatrixC128 from nested vectors using the rows factory method."
  [vectors]
  (let [rows-data (for [row vectors]
                    (for [[real imag] row]
                      (complex real imag)))]
    (matrix-from-rows rows-data)))

;; ============================================================================
;; Store operations
;; ============================================================================

(defn ^GenericStore create-store
  "Create an empty mutable GenericStore for complex matrices."
  [^long rows ^long cols]
  (.make (GenericStore/C128) rows cols))

(defn ^MatrixC128 store->matrix
  "Wrap a GenericStore in an immutable MatrixC128."
  [^GenericStore store]
  (.makeWrapper MatrixC128/FACTORY store))

(defn ^GenericStore store-from-vectors
  "Create a GenericStore from nested vectors (element-by-element)."
  [vectors]
  (let [rows (long (count vectors))
        cols (long (count (first vectors)))
        store (create-store rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (let [[real imag] (get-in vectors [i j])
              c (complex real imag)]
          (.set ^GenericStore store (long i) (long j) ^ComplexNumber c))))
    store))

(defn ^MatrixC128 identity-matrix
  "Create an identity matrix of size n x n."
  [^long n]
  (.makeIdentity MatrixC128/FACTORY (int n)))

;; ============================================================================
;; Common operations
;; ============================================================================

(defn ^Comparable get-element
  "Get an element from a matrix."
  [matrix ^long row ^long col]
  (.get ^BasicMatrix matrix row col))

(defn ^BasicMatrix matrix-multiply
  "Multiply two matrices."
  [^BasicMatrix m1 ^BasicMatrix m2]
  (.multiply m1 m2))

(defn ^BasicMatrix matrix-add
  "Add two matrices."
  [^BasicMatrix m1 ^BasicMatrix m2]
  (.add m1 m2))

(defn ^BasicMatrix matrix-conjugate
  "Conjugate of a matrix."
  [^BasicMatrix m]
  (.conjugate m))

(defn ^BasicMatrix matrix-transpose
  "Transpose of a matrix."
  [^BasicMatrix m]
  (.transpose m))

(defn ^BasicMatrix matrix-invert
  "Inverse of a matrix."
  [^BasicMatrix m]
  (.invert m))

(defn ^BasicMatrix matrix-solve
  "Solve Ax = B for x."
  [^BasicMatrix A ^BasicMatrix B]
  (.solve A B))

(defn ^BasicMatrix matrix-scale
  "Multiply matrix by a scalar."
  [^BasicMatrix m scalar]
  (.multiply m ^Comparable scalar))

(defn ^BasicMatrix matrix-scale-real
  "Multiply matrix by a real scalar."
  [^BasicMatrix m ^double scalar]
  (.multiply m scalar))

(defn matrix-shape
  "Get [rows cols] dimensions of a matrix."
  [^BasicMatrix matrix]
  [(.countRows matrix) (.countColumns matrix)])

;; ============================================================================
;; Conversion and display
;; ============================================================================

(defn matrix->vectors
  "Convert a matrix to nested vectors of [real imag] pairs."
  [^BasicMatrix matrix]
  (let [rows (.countRows matrix)
        cols (.countColumns matrix)]
    (vec (for [i (range rows)]
           (vec (for [j (range cols)]
                  (let [^ComplexNumber c (.get matrix (long i) (long j))]
                    [(double (.getReal c))
                     (double (.getImaginary c))])))))))

(defn print-matrix
  "Print a matrix in a readable format."
  [^BasicMatrix matrix]
  (dotimes [i (.countRows matrix)]
    (dotimes [j (.countColumns matrix)]
      (let [^ComplexNumber c (.get matrix (long i) (long j))
            real (double (.getReal c))
            imag (double (.getImaginary c))]
        (print (format "%7.2f%+.2fi  " real imag))))
    (println)))

;; ============================================================================
;; Examples
;; ============================================================================

(comment
  ;; METHOD 1: Using matrix-from-vectors (via GenericStore)
  (def A (matrix-from-vectors
          [[[1 0] [0 1]]
           [[2 3] [-1 -2]]]))
  (print-matrix A)

  (type A)

  ;; METHOD 2: Using matrix-from-vectors-via-rows (via MatrixC128 factory)
  (def B (matrix-from-vectors-via-rows
          [[[1 0] [0 -1]]
           [[0 1] [1 0]]]))
  (print-matrix B)

  ;; Matrix operations

  (def A-inv (matrix-invert A))
  (print-matrix A-inv)

  (def C (matrix-multiply A A-inv))
  (print-matrix C)

  (def I (identity-matrix 3))
  (print-matrix I))
