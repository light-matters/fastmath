(ns fastmath.api.v2.algebra.fast-test
  (:require
   [clojure.test :refer [deftest is testing]]
   [fastmath.algebra.function.matrix :as gmat]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as cd]
   [fastmath.algebra.object.number.complex.create :as C]
   [fastmath.api.v2.algebra.fast :as sut]))

(defn- rnd-rows [^java.util.Random rng r c]
  (vec (for [_ (range r)] (vec (for [_ (range c)] [(- (.nextDouble rng) 0.5) (- (.nextDouble rng) 0.5)])))))

(defn- close? [a b tol]
  (every? true? (map (fn [[x y] [u v]] (and (<= (Math/abs (- x u)) tol) (<= (Math/abs (- y v)) tol)))
                     (apply concat (sut/->vectors a)) (apply concat (sut/->vectors b)))))

(defn- cmul [[a b] [c d]] [(- (* a c) (* b d)) (+ (* a d) (* b c))])

(deftest constructors
  (let [rows [[[1.0 2.0] 3.0] [(C/i 4.0 -1.0) [0.0 5.0]]]
        A (sut/<-rows rows)]
    (is (= [[[1.0 2.0] [3.0 0.0]] [[4.0 -1.0] [0.0 5.0]]] (sut/->vectors A)))
    (is (= [2 2] (sut/shape A))))
  (testing "complexdense from a double[] (previously failed: flatten drops arrays)"
    (is (= [[[1.0 0.0] [0.0 0.0]] [[0.0 0.0] [1.0 0.0]]]
           (sut/->vectors (cd/complexdense 2 2 (double-array [1 0 0 0 0 0 1 0]))))))
  (testing "interleaved round trip copies, never aliases"
    (let [d (double-array [1 2 3 4 5 6])
          A (sut/<-interleaved 1 3 d)]
      (aset d 0 99.0)
      (is (= [[[1.0 2.0] [3.0 4.0] [5.0 6.0]]] (sut/->vectors A)))
      (let [out (sut/->interleaved A)]
        (aset out 1 -7.0)
        (is (= 2.0 (second (ffirst (sut/->vectors A))))))))
  (is (= [[[1.0 0.0] [0.0 0.0]] [[0.0 0.0] [1.0 0.0]]] (sut/->vectors (sut/identity 2))))
  (is (= [[[2.0 1.0] [0.0 0.0]] [[0.0 0.0] [3.0 0.0]]] (sut/->vectors (sut/diagonal [[2.0 1.0] 3.0])))))

(deftest kronecker-matches-generic-definition
  (let [rng (java.util.Random. 1)]
    (doseq [[r1 c1 r2 c2] [[2 2 2 2] [4 4 4 4] [2 3 3 1] [1 4 3 2]]]
      (let [A (sut/<-rows (rnd-rows rng r1 c1))
            B (sut/<-rows (rnd-rows rng r2 c2))]
        (is (close? (sut/kronecker A B) (cd/<-rows (gmat/kronecker-rows A B)) 1e-15) [r1 c1 r2 c2])
        (is (= [(* r1 r2) (* c1 c2)] (sut/shape (sut/kronecker A B))))))))

(deftest linear-combination-matches-scale-and-add
  (let [rng (java.util.Random. 2)
        ms (repeatedly 4 #(sut/<-rows (rnd-rows rng 5 3)))
        cs [1.0 [0.3 -0.2] [0.0 1.0] -2.5]]
    (is (close? (sut/linear-combination cs ms)
                (reduce sut/add (map sut/scale ms cs)) 1e-14))
    (is (thrown? clojure.lang.ExceptionInfo
                 (sut/linear-combination [1.0 1.0] [(sut/zeros 2 2) (sut/zeros 2 3)])))))

(deftest solve-square-and-least-squares
  (let [rng (java.util.Random. 3)]
    (doseq [n [1 2 4 16 33]]
      (let [A (sut/add (sut/<-rows (rnd-rows rng n n)) (sut/scale (sut/identity n) 4.0))
            B (sut/<-rows (rnd-rows rng n 2))
            X (sut/solve A B)]
        (is (close? (sut/multiply A X) B 1e-12) n)
        (testing "inputs untouched (the solver is reused)"
          (is (close? (sut/solve A B) X 0.0)))))
    (is (thrown? clojure.lang.ExceptionInfo (sut/solve (sut/zeros 3 3) (sut/<-rows (rnd-rows rng 3 1)))))
    (testing "overdetermined: least squares (normal equations hold)"
      (let [A (sut/<-rows (rnd-rows rng 6 3)) B (sut/<-rows (rnd-rows rng 6 1))
            X (sut/solve A B)
            r (sut/subtract (sut/multiply A X) B)]
        (is (close? (sut/multiply (sut/adjoint A) r) (sut/zeros 3 1) 1e-12))))))

(deftest solve-is-thread-safe
  (let [rng (java.util.Random. 4)
        problems (vec (for [n (take 60 (cycle [4 16 7]))]
                        [(sut/add (sut/<-rows (rnd-rows rng n n)) (sut/scale (sut/identity n) 3.0))
                         (sut/<-rows (rnd-rows rng n 1))]))
        seqr (mapv (fn [[A B]] (sut/->vectors (sut/solve A B))) problems)
        parr (vec (pmap (fn [[A B]] (sut/->vectors (sut/solve A B))) (concat problems problems problems)))]
    (is (= (concat seqr seqr seqr) parr))))

(deftest element-and-vectors-agree
  (let [A (sut/<-rows (rnd-rows (java.util.Random. 5) 3 4))]
    (doseq [i (range 3) j (range 4)
            :let [e (sut/element A i j) [re im] (get-in (sut/->vectors A) [i j])]]
      (is (= [re im] [(fastmath.protocol.algebra.structure.coordinate.complex/re e)
                      (fastmath.protocol.algebra.structure.coordinate.complex/im e)])))))

(deftest add-diagonal-matches-linear-combination
  (let [rng (java.util.Random. 6)
        A (sut/<-rows (rnd-rows rng 5 5))
        v (sut/<-rows (mapv vector (first (rnd-rows rng 1 5))))
        c [0.7 -1.3]]
    (is (close? (sut/add-diagonal A c v)
                (sut/linear-combination [1.0 c] [A (sut/diagonal (map first (sut/->vectors v)))]) 1e-15))
    (is (close? (sut/add-diagonal A 2.0 (sut/transpose v))
                (sut/linear-combination [1.0 2.0] [A (sut/diagonal (map first (sut/->vectors v)))]) 1e-15))
    (testing "A unchanged"
      (is (close? A (sut/<-rows (sut/->vectors A)) 0.0)))
    (is (thrown? clojure.lang.ExceptionInfo (sut/add-diagonal A 1.0 (sut/zeros 4 1))))))

(deftest sparse-terms-in-linear-combination
  (let [rng (java.util.Random. 8)
        A (sut/<-rows (rnd-rows rng 6 6))
        S (sut/<-rows (vec (for [i (range 6)] (vec (for [j (range 6)] (if (= (mod (+ i (* 2 j)) 5) 0) [(- i j) 0.5] 0.0))))))
        cs [1.0 [0.2 -0.9] -3.0]]
    (is (close? (sut/linear-combination cs [A (sut/sparse S) (sut/sparse A)])
                (sut/linear-combination cs [A S A]) 1e-14))
    (is (close? (sut/linear-combination [[0.0 1.0]] [(sut/sparse S)]) (sut/scale S [0.0 1.0]) 0.0))
    (is (thrown? clojure.lang.ExceptionInfo (sut/linear-combination [1.0 1.0] [A (sut/sparse (sut/zeros 5 6))])))))

(deftest solve-shifted-matches-solve-of-add-diagonal
  (let [rng (java.util.Random. 9)
        A (sut/add (sut/<-rows (rnd-rows rng 7 7)) (sut/scale (sut/identity 7) 3.0))
        v (sut/<-rows (mapv vector (first (rnd-rows rng 1 7))))
        B (sut/<-rows (rnd-rows rng 7 2))
        A-before (sut/->vectors A)]
    (doseq [c [0.0 1.5 [0.3 -2.0] -7.0]]
      (is (close? (sut/solve-shifted A c v B) (sut/solve (sut/add-diagonal A c v) B) 1e-13) (str c)))
    (is (= A-before (sut/->vectors A)) "A unchanged")
    (testing "thread-safe"
      (let [cs (vec (for [k (range 40)] [(* 0.1 k) (- 1.0 (* 0.05 k))]))
            seqr (mapv #(sut/->vectors (sut/solve-shifted A % v B)) cs)]
        (is (= seqr (vec (pmap #(sut/->vectors (sut/solve-shifted A % v B)) cs))))))
    (is (thrown? clojure.lang.ExceptionInfo
                 (sut/solve-shifted (sut/zeros 3 3) 0.0 (sut/zeros 3 1) (sut/zeros 3 1))))))
