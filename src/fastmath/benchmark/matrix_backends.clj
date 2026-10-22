(ns fastmath.benchmark.matrix-backends
  "Compares the EJML and ojAlgo matrix backends against each other, across
  the real and complex domains and a range of sizes.

  Measures through fastmath's own protocol-level API (`rmat/add`,
  `emat/multiply`, `sq/determinant`, `solve/cholesky`, `solve/solve`, ...)
  rather than calling either library directly - this is what the choice
  of backend actually costs/saves a fastmath user, including protocol
  dispatch and the wrapper deftype, not just the raw library's own speed.

  `determinant`/`cholesky`/`solve` each come in two variants: the plain
  one builds a fresh matrix per trial (the cost of computing on data the
  backend hasn't seen before), and the `-repeated` one queries the SAME
  matrix repeatedly on purpose, to surface ojAlgo's immutable-matrix
  caching (confirmed empirically: 10-400x faster from the second call on
  an unchanged matrix) as a real, reportable advantage rather than
  something the benchmark papers over.

  Run with `(run!)` for the default sweep, or `(-main)` from the CLI.
  "
  (:require
   [fastmath.benchmark.apache-ejml :refer [measure gflops gemm-flops lu-flops]]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as real-ejml]
   [fastmath.algebra.object.matrix.rectangular.real.ojalgo :as real-ojalgo]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as complex-ejml]
   [fastmath.algebra.object.matrix.rectangular.complex.ojalgo :as complex-ojalgo]
   [fastmath.protocol.algebra.object.matrix.rectangular :as rmat]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.square :as sq]
   [fastmath.protocol.algebra.object.matrix.solve :as solve]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.structure.space.normed :as nspace]
   [fastmath.protocol.representation.d2 :as d2])
  (:import
   (java.util Random)))

(set! *warn-on-reflection* true)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                               Timing harness                                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; `measure`/`gflops`/`gemm-flops`/`lu-flops` are reused from
;; fastmath.benchmark.apache-ejml (median-of-trials, each trial running
;; for at least :min-ms to average out timer noise).

(defn- min-ms-for
  "Longer measurement windows for tiny sizes, where a single call is too
  fast for System/nanoTime noise to average out otherwise."
  ^long [^long n]
  (cond (<= n 8) 600 (<= n 32) 300 :else 150))

(defn measure-fresh
  "Like `measure`, but for operations where calling `f` repeatedly on the
  SAME input would be unfairly cheap: confirmed empirically that ojAlgo's
  immutable BasicMatrix caches determinant/cholesky/solve results on the
  matrix instance itself, so repeated calls after the first are ~10-400x
  faster regardless of size - not representative of computing it once per
  (fresh) matrix, which is what `measure`'s inner repeat-many-times-fast
  loop would otherwise end up mostly measuring for those three ops.
  Builds a fresh input via `(make)` for every trial and times exactly one
  call of `(f input)` on it - this does fold construction cost in with the
  operation's own cost, a reasonable trade for correctness here since O(n^3)
  decomposition cost dominates construction at any size worth reporting."
  [make f {:keys [trials] :or {trials 7}}]
  (let [times (repeatedly trials
                          (fn []
                            (let [input (make)
                                  start (System/nanoTime)]
                              (f input)
                              (/ (double (- (System/nanoTime) start)) 1.0e9))))]
    (nth (sort times) (quot trials 2))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                         Per-domain/backend data generation                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A fixed-seed RNG so repeated runs generate the same matrices and are
;; comparable to each other.

(def ^:private ^Random rng (Random. 42))
(defn- rnd ^double [] (- (.nextDouble rng) 0.5))

(defn- real-rows [^long n ^long m]
  (vec (repeatedly n (fn [] (vec (repeatedly m rnd))))))

(defn- complex-rows [^long n ^long m]
  (vec (repeatedly n (fn [] (vec (repeatedly m (fn [] [(rnd) (rnd)])))))))

(def ^:private backend-info
  {[:real :ejml]      {:<-rows real-ejml/<-rows      :rows real-rows    :adjoint emat/transpose}
   [:real :ojalgo]    {:<-rows real-ojalgo/<-rows    :rows real-rows    :adjoint emat/transpose}
   [:complex :ejml]   {:<-rows complex-ejml/<-rows   :rows complex-rows :adjoint cmat/adjoint}
   [:complex :ojalgo] {:<-rows complex-ojalgo/<-rows :rows complex-rows :adjoint cmat/adjoint}})

(defn- matrix [domain backend ^long n ^long m]
  (let [{:keys [<-rows rows]} (backend-info [domain backend])]
    (<-rows (rows n m))))

(defn- spd
  "A random symmetric/Hermitian positive (semi-)definite matrix: R* R for
  a random square R - needed so cholesky has something valid to chew on."
  [domain backend ^long n]
  (let [{:keys [adjoint]} (backend-info [domain backend])
        R (matrix domain backend n n)]
    (emat/multiply (adjoint R) R)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                  Operations                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Each returns {:time-s ... :gflops ...} (:gflops nil when not meaningful).
;; `n` is always the square size being benchmarked at.

(defn- bench-add [domain backend ^long n]
  (let [A (matrix domain backend n n)
        B (matrix domain backend n n)]
    {:time-s (measure #(rmat/add A B) {:min-ms (min-ms-for n)}) :gflops nil}))

(defn- bench-transpose-view [domain backend ^long n]
  ;; The bare call. ojAlgo's `.transpose()` returns a lazy view over the
  ;; original store rather than eagerly copying, so this measures "build
  ;; a view" (near zero, regardless of n) for ojAlgo but real, eager work
  ;; for EJML - the asymmetry itself is the point of having both this and
  ;; `bench-transpose` below.
  (let [A (matrix domain backend n n)]
    {:time-s (measure #(emat/transpose A) {:min-ms (min-ms-for n)}) :gflops nil}))

(defn- bench-transpose [domain backend ^long n]
  ;; Forces every element of the transposed result to be read, so a lazy
  ;; view and an eager copy end up paying comparable total cost - this is
  ;; what the transpose actually costs once something downstream uses it.
  (let [A (matrix domain backend n n)]
    {:time-s (measure #(d2/->arrays (emat/transpose A)) {:min-ms (min-ms-for n)}) :gflops nil}))

(defn- bench-multiply [domain backend ^long n]
  (let [A (matrix domain backend n n)
        B (matrix domain backend n n)
        t (measure #(emat/multiply A B) {:min-ms (min-ms-for n)})]
    {:time-s t :gflops (gflops (gemm-flops n n n) t)}))

(defn- bench-determinant [domain backend ^long n]
  ;; A fresh matrix per trial, so this measures the cost of computing a
  ;; determinant on data the backend hasn't seen before - see
  ;; `bench-determinant-repeated` below for the complementary number.
  (let [t (measure-fresh #(matrix domain backend n n) sq/determinant {})]
    {:time-s t :gflops (gflops (lu-flops n) t)}))

(defn- bench-determinant-repeated [domain backend ^long n]
  ;; The SAME matrix queried repeatedly, deliberately. ojAlgo's BasicMatrix
  ;; is immutable and caches determinant/cholesky/solve on the instance -
  ;; confirmed empirically (10-400x faster from the second call on). That's
  ;; a genuine advantage of the immutable design for workloads that query
  ;; the same matrix's properties more than once (which is a realistic
  ;; pattern - e.g. checking singular?/rank before solving), so it belongs
  ;; in the comparison rather than being designed away. EJML's API is
  ;; mutable and doesn't (can't safely) do this, so it should look roughly
  ;; the same here as in the fresh-matrix version above.
  (let [A (matrix domain backend n n)
        t (measure #(sq/determinant A) {:min-ms (min-ms-for n)})]
    {:time-s t :gflops (gflops (lu-flops n) t)}))

(defn- bench-cholesky [domain backend ^long n]
  (let [t (measure-fresh #(spd domain backend n) solve/cholesky {})]
    {:time-s t :gflops (gflops (/ (lu-flops n) 2.0) t)}))

(defn- bench-cholesky-repeated [domain backend ^long n]
  ;; See `bench-determinant-repeated` - same idea, for cholesky.
  (let [A (spd domain backend n)
        t (measure #(solve/cholesky A) {:min-ms (min-ms-for n)})]
    {:time-s t :gflops (gflops (/ (lu-flops n) 2.0) t)}))

(defn- bench-solve [domain backend ^long n]
  (let [{:keys [<-rows rows]} (backend-info [domain backend])
        t (measure-fresh #(vector (spd domain backend n) (<-rows (rows n 1)))
                         (fn [[A b]] (solve/solve A b))
                         {})]
    {:time-s t :gflops (gflops (lu-flops n) t)}))

(defn- bench-solve-repeated [domain backend ^long n]
  ;; See `bench-determinant-repeated` - same idea, for solve (same A, b
  ;; reused across trials; only A is what ojAlgo would cache against).
  (let [{:keys [<-rows rows]} (backend-info [domain backend])
        A (spd domain backend n)
        b (<-rows (rows n 1))
        t (measure #(solve/solve A b) {:min-ms (min-ms-for n)})]
    {:time-s t :gflops (gflops (lu-flops n) t)}))

(defn- bench-norm [domain backend ^long n]
  (let [A (matrix domain backend n n)]
    {:time-s (measure #(nspace/norm A) {:min-ms (min-ms-for n)}) :gflops nil}))

(def ops
  {:add                 {:f #'bench-add                 :sizes [8 16 32 64 128 256]}
   :transpose           {:f #'bench-transpose           :sizes [8 16 32 64 128 256]}
   :transpose-view      {:f #'bench-transpose-view      :sizes [8 16 32 64 128 256]}
   :norm                {:f #'bench-norm                :sizes [8 16 32 64 128 256]}
   :multiply            {:f #'bench-multiply            :sizes [8 16 32 64 128]}
   :determinant         {:f #'bench-determinant         :sizes [8 16 32 64 128]}
   :determinant-repeated {:f #'bench-determinant-repeated :sizes [8 16 32 64 128]}
   :cholesky            {:f #'bench-cholesky            :sizes [8 16 32 64 128]}
   :cholesky-repeated   {:f #'bench-cholesky-repeated   :sizes [8 16 32 64 128]}
   :solve               {:f #'bench-solve               :sizes [8 16 32 64 128]}
   :solve-repeated      {:f #'bench-solve-repeated      :sizes [8 16 32 64 128]}})

(def domains [:real :complex])
(def backends [:ejml :ojalgo])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 Reporting                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def header
  (str (format "%8s  %-10s  %-7s  %5s   %10s   %8s" "DOMAIN" "OP" "BACKEND" "N" "MED(s)" "GF/s")
       "\n" (apply str (repeat 56 "-"))))

(defn- fmt6 [x] (format "%.6f" (double x)))
(defn- fmt2 [x] (format "%.2f" (double x)))

(defn- row [{:keys [domain op backend n time-s gflops]}]
  (println (format "%8s  %-10s  %-7s  %5d   %10s   %8s"
                   (name domain) (name op) (name backend) (long n)
                   (fmt6 time-s) (if gflops (fmt2 gflops) "-"))))

(defn- error-row [domain op backend n ^Throwable t]
  (println (format "%8s  %-10s  %-7s  %5d   ERROR: %s"
                   (name domain) (name op) (name backend) (long n) (.getMessage t))))

(defn run!
  "Runs the full sweep and prints a table. Pass `:ops-filter`/`:domains-filter`/
  `:backends-filter` to narrow it, e.g.
  (run! {:ops-filter [:multiply] :domains-filter [:complex]})."
  ([] (run! {}))
  ([{:keys [ops-filter domains-filter backends-filter]
     :or   {ops-filter (keys ops) domains-filter domains backends-filter backends}}]
   (println header)
   (doseq [op ops-filter
           :let [{:keys [f sizes]} (ops op)]
           n sizes
           domain domains-filter
           backend backends-filter]
     (try
       (row (merge {:domain domain :op op :backend backend :n n} (f domain backend n)))
       (catch Throwable t
         (error-row domain op backend n t))))))

(defn -main [& _] (run!))
