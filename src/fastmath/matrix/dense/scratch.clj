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

(ns fastmath.matrix.complex
  (:import [org.ojalgo.matrix MatrixC128 MatrixC128$Factory MatrixC128$DenseReceiver BasicMatrix]
           [org.ojalgo.scalar ComplexNumber]
           [org.ojalgo.matrix.store GenericStore]))

(set! *warn-on-reflection* true)

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
;; High-level API using MatrixC128 (recommended for most use cases)
;; ============================================================================

(defn ^MatrixC128$DenseReceiver create-matrix-builder
  "Create a mutable builder for constructing a MatrixC128.
  Use set! to populate it, then call build to get the immutable matrix."
  [^long rows ^long cols]
  (let [size (* rows cols)]
    (.makeDense MatrixC128/FACTORY (int size))))

(defn ^MatrixC128 build-matrix
  "Build an immutable MatrixC128 from a builder."
  [^MatrixC128$DenseReceiver builder]
  (.build builder))

(defn set-element!
  "Set an element in a matrix builder. Returns the builder for chaining."
  [^MatrixC128$DenseReceiver builder ^long row ^long col value]
  ;; Use the specific overload: set(long, long, Comparable)
  (.set builder row col ^Comparable value)
  builder)

(defn ^MatrixC128 matrix-from-vectors
  "Create a MatrixC128 from nested vectors of [real imag] pairs."
  [vectors]
  (let [rows (long (count vectors))
        cols (long (count (first vectors)))
        ^MatrixC128$DenseReceiver builder (create-matrix-builder rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (let [[real imag] (get-in vectors [i j])
              ^Comparable c (complex real imag)]
          (.set builder (long i) (long j) c))))
    (build-matrix builder)))

(defn ^MatrixC128 identity-matrix
  "Create an identity matrix of size n x n."
  [^long n]
  (.makeIdentity MatrixC128/FACTORY (int n)))

;; ============================================================================
;; Low-level API using GenericStore (for maximum performance and mutability)
;; ============================================================================

(defn ^GenericStore create-store
  "Create a mutable GenericStore for complex matrices.
  This is the low-level, high-performance option."
  [^long rows ^long cols]
  (.make (GenericStore/C128) rows cols))

(defn ^MatrixC128 store->matrix
  "Wrap a GenericStore in an immutable MatrixC128."
  [^GenericStore store]
  (.makeWrapper MatrixC128/FACTORY store))

(defn ^GenericStore store-from-vectors
  "Create a GenericStore from nested vectors of [real imag] pairs."
  [vectors]
  (let [rows (long (count vectors))
        cols (long (count (first vectors)))
        ^GenericStore store (create-store rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (let [[real imag] (get-in vectors [i j])
              ^ComplexNumber c (complex real imag)]
          (.set store (long i) (long j) c))))
    store))

;; ============================================================================
;; Common operations (work with both MatrixC128 and GenericStore)
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
  ;; HIGH-LEVEL APPROACH (Recommended - immutable, functional)
  (def A (matrix-from-vectors
          [[[1 0] [0 1]]
           [[2 3] [-1 -2]]]))

  (print-matrix A)

  (def B (matrix-from-vectors
          [[[1 0] [0 -1]]
           [[0 1] [1 0]]]))

  ;; Matrix operations
  (def C (matrix-multiply A B))
  (print-matrix C)

  (def A-inv (matrix-invert A))
  (print-matrix A-inv)

  (def I (identity-matrix 3))
  (print-matrix I)

  ;; LOW-LEVEL APPROACH (For performance)
  (def store (create-store 2 2))
  (.set store (long 0) (long 0) (complex 1.0 0.0))
  (.set store (long 0) (long 1) (complex 0.0 1.0))
  (print-matrix (store->matrix store)))

(comment
  (defn ->complex-array
    [rows3d]
    (into-array
     (mapv (fn [plane]
             (into-array
              (mapv #(apply complex %) plane)))
           rows3d)))

  (def carray (->complex-array [[[1 2]  [1 0]]
                                [[0 2]  [1 0]]]))

  (.rows MatrixC128$Factory carray))




