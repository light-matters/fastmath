
(ns benchmark.apache-ejml--minimum
  (:import
   [java.util Random]
   [org.ejml.data DMatrixRMaj]
   [org.ejml.dense.row CommonOps_DDRM]
   [org.apache.commons.math3.linear Array2DRowRealMatrix RealMatrix
    LUDecomposition DecompositionSolver]))

(set! *warn-on-reflection* true)

;; -------- stable RNG (repeatable runs) --------
(def ^Random rng (Random. 42))
(defn rnd ^double [] (- (.nextDouble rng) 0.5))

;; -------- tiny helpers --------
(defn gflops [flops seconds] (/ flops (* seconds 1.0e9)))
(defn gemm-flops [m n k] (* 2.0 m n k))
(defn lu-flops   [n]     (* (/ 2.0 3.0) n n n))

(defn min-ms-for
  "Heuristic: longer min-measure times for tiny sizes to reduce noise."
  [m]
  (cond
    (<= m 8)   1000
    (<= m 16)   500
    (<= m 64)   300
    :else       200))

(defn measure
  "Median per-call seconds for thunk f. Each trial runs ≥ :min-ms."
  [f {:keys [min-ms trials warmups] :or {min-ms 200 trials 7 warmups 3}}]
  ;; warm-up first
  (dotimes [_ warmups] (f))
  (let [arr    (double-array trials)
        min-ns (long (* 1.0e6 (double min-ms)))]
    (dotimes [t trials]
      (let [start    (System/nanoTime)
            deadline (+ start min-ns)]
        (loop [it 0]
          (if (< (System/nanoTime) deadline)
            (do (f) (recur (unchecked-inc it)))
            (let [elapsed   (- (System/nanoTime) start)
                  seconds   (/ (double elapsed) 1.0e9)
                  per-op    (/ seconds (double (max 1 it)))]
              (aset-double arr t per-op))))))
    ;; median
    (aget (double-array (sort (seq arr))) (quot trials 2))))

;; -------- EJML (real, dense) --------
(defn ejml-matrix ^DMatrixRMaj [rows cols]
  ;; Avoid varargs constructor; fill via setters.
  (let [r (int rows) c (int cols) M (DMatrixRMaj. r c)]
    (dotimes [i rows]
      (dotimes [j cols]
        (.set M (int i) (int j) (rnd))))
    M))

(defn ejml-mmul ^DMatrixRMaj [^DMatrixRMaj A ^DMatrixRMaj B]
  (let [C (DMatrixRMaj. (.numRows A) (.numCols B))]
    (CommonOps_DDRM/mult A B C)
    C))

(defn ejml-solve ^DMatrixRMaj [^DMatrixRMaj A]
  ;; Solve A x = b with random b; no varargs.
  (let [n (.numRows A)
        b (DMatrixRMaj. n 1)
        x (DMatrixRMaj. n 1)]
    (dotimes [i n] (.set b (int i) 0 (rnd)))
    (CommonOps_DDRM/solve A b x)
    x))

;; -------- Apache Commons Math (real, dense) --------
(defn commons-matrix ^RealMatrix [rows cols]
  ;; Avoid double[][] ctor; use (rows, cols) then setEntry.
  (let [r (int rows) c (int cols) M (Array2DRowRealMatrix. r c)]
    (dotimes [i rows]
      (dotimes [j cols]
        (.setEntry M (int i) (int j) (rnd))))
    M))

(defn commons-mmul ^RealMatrix [^RealMatrix A ^RealMatrix B]
  (.multiply A B))

(defn commons-solve ^RealMatrix [^RealMatrix A]
  (let [n (.getRowDimension A)
        b (double-array n)]
    (dotimes [i n] (aset-double b (int i) (rnd)))
    (.solve (.getSolver (LUDecomposition. A)) b)))

;; -------- sizes & printing --------
(def sizes-gemm
  [[4 4 4] [8 8 8] [16 16 16] [32 32 32]
   [64 64 64] [128 128 128] [256 256 256] [512 512 512]
   [1024 1024 1024] [2048 64 2048] [64 2048 64]])

(def sizes-solve [4 8 16 32 64 128 256 512 1024])

(def header
  (str (format "%8s  %-6s  %5s  %5s  %5s   %17s   %8s"
               "LIB" "OP" "M" "K" "N" "MED(s) (µs)" "GF/s")
       "\n" (apply str (repeat 78 "-"))))

(defn fmt-times [t]
  (format "%.9f s (%.2f µs)" t (* 1e6 t)))

(defn print-row [{:keys [lib op shape time-s gflops]}]
  (let [[m k n] (concat shape (repeat nil))]
    (println (format "%8s  %-6s  %5s  %5s  %5s   %17s   %8.2f"
                     (name lib) (name op)
                     (or m "-") (or k "-") (or n "-")
                     (fmt-times time-s) (double (or gflops 0.0))))))

;; -------- benches --------
(defn bench-gemm [lib m k n]
  (case lib
    :ejml
    (let [A (ejml-matrix m k)
          B (ejml-matrix k n)
          t (measure #(ejml-mmul A B) {:min-ms (min-ms-for m)})]
      {:lib lib :op :gemm :shape [m k n]
       :time-s t :gflops (gflops (gemm-flops m n k) t)})

    :commons
    (let [A (commons-matrix m k)
          B (commons-matrix k n)
          t (measure #(commons-mmul A B) {:min-ms (min-ms-for m)})]
      {:lib lib :op :gemm :shape [m k n]
       :time-s t :gflops (gflops (gemm-flops m n k) t)})))

(defn bench-solve [lib n]
  (case lib
    :ejml
    (let [A (ejml-matrix n n)
          t (measure #(ejml-solve A) {:min-ms (max 300 (min-ms-for n))})]
      {:lib lib :op :solve :shape [n n 1]
       :time-s t :gflops (gflops (lu-flops n) t)})

    :commons
    (let [A (commons-matrix n n)
          t (measure #(commons-solve A) {:min-ms (max 300 (min-ms-for n))})]
      {:lib lib :op :solve :shape [n n 1]
       :time-s t :gflops (gflops (lu-flops n) t)})))

;; -------- entrypoint --------
(defn -main [& _]
  (println header)
  ;; GEMM
  (doseq [[m k n] sizes-gemm, lib [:ejml :commons]]
    (try
      (print-row (bench-gemm lib m k n))
      (catch Throwable t
        (println (format "%8s  %-6s  %5d  %5d  %5d   ERROR: %s"
                         (name lib) "gemm" m k n (.getMessage t))))))
  ;; SOLVE
  (doseq [n sizes-solve, lib [:ejml :commons]]
    (try
      (print-row (bench-solve lib n))
      (catch Throwable t
        (println (format "%8s  %-6s  %5d  %5s  %5s   ERROR: %s"
                         (name lib) "solve" n "-" "-" (.getMessage t)))))))
