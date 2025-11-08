(ns fastmath.matrix.dense.scratch2
  "A playground for working through ideas.

  Currently playing with the complex implementation of ojalgo matrices.

https://gist.github.com/apete/b3278dc2f8c2db6a00369c211ba321db
  Reading through this ^^^ (CommonMistake.java) helped to get an idea of the API design. It's stil pretty implicit though and I'm sure I'm not using it right.
  "
  (:require [criterium.core :as crit])
  (:import
   (org.ojalgo.matrix BasicMatrix MatrixC128 MatrixC128$DenseReceiver)
   (org.ojalgo.matrix.store GenericStore GenericStore$Factory)
   (org.ojalgo.scalar ComplexNumber)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(println "===start===")

(defn cn ^ComplexNumber [r i]
  (ComplexNumber/of (double r) (double i)))

(defn ->store
  [rows--complex]
  (let [nrows (count rows--complex)
        data (->carray rows--complex)]
    (print data)
    (GenericStore/wrap (.array GenericStore$Factory) data nrows)))

(def store
  (GenericStore/wrap GenericStore/C128 (->carray vcs) 2))
(->store vcs)

(println store)

(def vcs [[[1 0] [2 1]]
          [[3 -1] [4 0]]])

(defn ->carray
  "TODO: keep this code, regardless!
  I learned about transducers again..."
  [rows--complex]
  (->> rows--complex
       (into [] (comp cat))
       (mapv #(apply cn %))
       (into-array ComplexNumber)))

(defn ^GenericStore wrap-array-reflection
  "Use reflection to call GenericStore.wrap, bypassing access restrictions."
  [rows cols data]
  (let [arr (into-array ComplexNumber data)
        factory (GenericStore/C128)
        ;; Find the wrap method with 3 parameters
        wrap-method (->> (.getMethods GenericStore)
                         (filter #(= "wrap" (.getName ^java.lang.reflect.Method %)))
                         (filter #(= 3 (alength (.getParameterTypes ^java.lang.reflect.Method %))))
                         first)]
    (when wrap-method
      (.setAccessible ^java.lang.reflect.Method wrap-method true)
      ;; invoke takes: (method, object-instance, arg-array)
      ;; For static methods, object-instance is nil
      (.invoke ^java.lang.reflect.Method wrap-method nil (into-array Object [factory arr (int cols)])))))

;; Or the practical solution - populate manually
(defn ^GenericStore create-from-nested-vectors
  "Efficiently create GenericStore from nested vectors."
  [vectors]
  (let [rows (long (count vectors))
        cols (long (count (first vectors)))
        arr (->carray vectors)
        store (.make (GenericStore/C128) rows cols)]
    ;; Populate from the flat array
    (dotimes [i rows]
      (dotimes [j cols]
        (let [idx (+ (* i cols) j)]
          (.set store (long i) (long j) ^ComplexNumber (aget arr idx)))))
    store))

;; Print function (was missing)
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

;; Test
(def vcs [[[1 0] [2 1]]
          [[3 -1] [4 0]]])

(def result (wrap-array-reflection 2 2 (->carray vcs)))

;; Quick bench
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                    bench                                    ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Test data
(def test-vcs-small [[[1 0] [2 1]]
                     [[3 -1] [4 0]]])

(def test-vcs-medium
  (vec (for [i (range 10)]
         (vec (for [j (range 10)]
                [(double i) (double j)])))))

(def test-vcs-large
  (vec (for [i (range 100)]
         (vec (for [j (range 100)]
                [(double i) (double j)])))))

;; Method 1: Reflection-based wrap
(defn bench-wrap-reflection [vcs]
  (let [rows (count vcs)
        cols (count (first vcs))]
    (wrap-array-reflection rows cols (->carray vcs))))

;; Method 2: Manual population
(defn bench-manual-populate [vcs]
  (create-from-nested-vectors vcs))

;; Method 3: Element-by-element (baseline)
(defn bench-element-by-element [vcs]
  (let [rows (long (count vcs))
        cols (long (count (first vcs)))
        store (.make (GenericStore/C128) rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (let [[real imag] (get-in vcs [i j])]
          (.set store (long i) (long j) ^ComplexNumber (cn real imag)))))
    store))

;; Run benchmarks
(println "=== 2x2 Matrix ===")
(println "Reflection wrap:")
(crit/quick-bench (bench-wrap-reflection test-vcs-small))

(println "\nManual populate:")
(crit/quick-bench (bench-manual-populate test-vcs-small))

(println "\nElement-by-element:")
(crit/quick-bench (bench-element-by-element test-vcs-small))

(println "\n=== 10x10 Matrix ===")
(println "Reflection wrap:")
(crit/quick-bench (bench-wrap-reflection test-vcs-medium))

(println "\nManual populate:")
(crit/quick-bench (bench-manual-populate test-vcs-medium))

(println "\nElement-by-element:")
(crit/quick-bench (bench-element-by-element test-vcs-medium))

(println "\n=== 100x100 Matrix ===")
(println "Reflection wrap:")
(crit/quick-bench (bench-wrap-reflection test-vcs-large))

(println "\nManual populate:")
(crit/quick-bench (bench-manual-populate test-vcs-large))

(println "\nElement-by-element:")
(crit/quick-bench (bench-element-by-element test-vcs-large))
