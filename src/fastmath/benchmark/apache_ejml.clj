(ns fastmath.benchmark.apache-ejml
  (:import
   [org.ejml.data DMatrixRMaj]
   [org.ejml.dense.row CommonOps_DDRM]
   [org.apache.commons.math3.linear Array2DRowRealMatrix RealMatrix
    LUDecomposition DecompositionSolver]))

;; ---------- helpers ----------
(defn rnd ^double [] (unchecked-double (- (rand) 0.5)))

(defn gflops [^double flops ^double sec]
  (/ flops (* sec 1.0e9)))

(defn gemm-flops [m n k] (* 2.0 m n k))
(defn lu-flops   [n]     (* (/ 2.0 3.0) n n n))

(defn measure
  "Median per-call seconds for thunk f. Each trial runs ≥ :min-ms."
  [f {:keys [min-ms trials warmups] :or {min-ms 200 trials 7 warmups 3}}]
  ;; warm up
  (dotimes [_ warmups] (f))
  (let [arr (double-array trials)
        min-ns (long (* 1e6 min-ms))]
    (dotimes [t trials]
      (let [deadline (+ (System/nanoTime) min-ns)
            start    (System/nanoTime)]
        (loop [it 0]
          (if (< (System/nanoTime) deadline)
            (do (f) (recur (inc it)))
            (let [elapsed (- (System/nanoTime) start)]
              (aset-double arr t (/ (double elapsed) 1.0e9 (max 1 it))))))))
    (->> arr seq sort (nth (quot trials 2)))))

;; ---------- EJML (real, dense) ----------
(defn ejml-matrix ^DMatrixRMaj [rows cols]
  ;; Avoid varargs ctor; fill via setters.
  (let [M (DMatrixRMaj. rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (.set M i j (rnd))))
    M))

(defn ejml-gemm [^DMatrixRMaj A ^DMatrixRMaj B]
  (let [C (DMatrixRMaj. (.numRows A) (.numCols B))]
    (CommonOps_DDRM/mult A B C)
    C))

(defn ejml-solve [^DMatrixRMaj A]
  ;; Solve A x = b with random b; avoid any array-based ctor.
  (let [n (.numRows A)
        b  (DMatrixRMaj. n 1)
        x  (DMatrixRMaj. n 1)]
    (dotimes [i n] (.set b i 0 (rnd)))
    (CommonOps_DDRM/solve A b x)
    x))

;; ---------- Apache Commons Math (real, dense) ----------
(defn commons-matrix ^RealMatrix [rows cols]
  ;; Avoid (double[][] ...) ctors; allocate and fill via setEntry.
  (let [M (Array2DRowRealMatrix. rows cols)]
    (dotimes [i rows]
      (dotimes [j cols]
        (.setEntry M i j (rnd))))
    M))

(defn commons-gemm [^RealMatrix A ^RealMatrix B]
  (.multiply A B))

(defn commons-solve [^RealMatrix A]
  (let [n (.getRowDimension A)
        b (double-array n)]
    (dotimes [i n] (aset-double b i (rnd)))
    (let [solver ^DecompositionSolver (.getSolver (LUDecomposition. A))]
      (.solve solver b))))

;; ---------- harness & I/O ----------
(def sizes-gemm
  [[4 4 4] [8 8 8] [16 16 16] [32 32 32]
   [64 64 64] [128 128 128] [256 256 256] [512 512 512]
   [1024 1024 1024] [2048 64 2048] [64 2048 64]])

(def sizes-solve [4 8 16 32 64 128 256 512 1024])

(def header
  (str (format "%8s  %-6s  %5s  %5s  %5s   %8s   %6s"
               "LIB" "OP" "M" "K" "N" "MED(s)" "GF/s")
       "\n" (apply str (repeat 66 "-"))))

(defn fmt2 [x] (format "%.2f" (double x)))
(defn fmt6 [x] (format "%.6f" (double x)))

(defn row [{:keys [lib op shape time-s gflops]}]
  (let [[m k n] (concat shape (repeat nil))]
    (println (format "%8s  %-6s  %5s  %5s  %5s   %8s   %6s"
                     (name lib) (name op)
                     (or m "-") (or k "-") (or n "-")
                     (fmt6 time-s) (fmt2 (or gflops 0.0))))))

(defn bench-gemm [lib m k n]
  (case lib
    :ejml
    (let [A (ejml-matrix m k)
          B (ejml-matrix k n)
          t (measure #(ejml-gemm A B) {})]
      {:lib lib :op :gemm :shape [m k n]
       :time-s t :gflops (gflops (gemm-flops m n k) t)})
    :commons
    (let [A (commons-matrix m k)
          B (commons-matrix k n)
          t (measure #(commons-gemm A B) {})]
      {:lib lib :op :gemm :shape [m k n]
       :time-s t :gflops (gflops (gemm-flops m n k) t)})))

(defn bench-solve [lib n]
  (case lib
    :ejml
    (let [A (ejml-matrix n n)
          t (measure #(ejml-solve A) {:min-ms 300})]
      {:lib lib :op :solve :shape [n n 1]
       :time-s t :gflops (gflops (lu-flops n) t)})
    :commons
    (let [A (commons-matrix n n)
          t (measure #(commons-solve A) {:min-ms 300})]
      {:lib lib :op :solve :shape [n n 1]
       :time-s t :gflops (gflops (lu-flops n) t)})))

(defn run! []
  (println header)
  (doseq [[m k n] sizes-gemm, lib [:ejml :commons]]
    (try
      (row (bench-gemm lib m k n))
      (catch Throwable t
        (println (format "%8s  %-6s  %5d  %5d  %5d   ERROR: %s"
                         (name lib) "gemm" m k n (.getMessage t))))))
  (doseq [n sizes-solve, lib [:ejml :commons]]
    (try
      (row (bench-solve lib n))
      (catch Throwable t
        (println (format "%8s  %-6s  %5d  %5s  %5s   ERROR: %s"
                         (name lib) "solve" n "-" "-" (.getMessage t)))))))

(defn -main [& _] (run!))
