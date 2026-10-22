(ns fastmath.algebra.object.matrix.rectangular.complex.fourier.ejml-test
  "Requiring fastmath.algebra.object.matrix.rectangular.complex.fourier.ejml
  (not just the protocol) is what actually installs the extend-type - a
  bare require of the protocol namespace alone wouldn't dispatch to
  anything."
  (:require [fastmath.algebra.object.matrix.rectangular.complex.fourier.ejml]
            [fastmath.algebra.object.matrix.rectangular.complex.ejml :as sut]
            [fastmath.protocol.analysis.transform.fourier :as dft]
            [fastmath.protocol.representation.d2 :as d2]
            [clojure.test :as t :refer [deftest testing are is]]))

(defn- col [pairs] (sut/<-rows (mapv vector pairs)))

(deftest power-of-two-identities
  (testing "impulse -> constant"
    (is (= (col [[1.0 0.0] [1.0 0.0] [1.0 0.0] [1.0 0.0]])
           (dft/fft (col [[1.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]])))))
  (testing "DC -> single spike"
    (is (= (col [[4.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]])
           (dft/fft (col [[1.0 0.0] [1.0 0.0] [1.0 0.0] [1.0 0.0]])))))
  (testing "genuinely complex input is handled correctly - a purely
            imaginary impulse's DFT is [i, 1, -i, -1], not all zeros
            (this is the exact case that exposed a real bug in ojAlgo's
            own DiscreteFourierTransform, which silently drops the
            imaginary part of complex input - see complex-branch-handover.md)"
    (let [result (dft/fft (col [[0.0 0.0] [0.0 1.0] [0.0 0.0] [0.0 0.0]]))
          [[r0 i0] [r1 i1] [r2 i2] [r3 i3]] (mapv first (d2/->vectors result))]
      (is (< (Math/abs (- r0 0.0)) 1.0e-9)) (is (< (Math/abs (- i0 1.0)) 1.0e-9))
      (is (< (Math/abs (- r1 1.0)) 1.0e-9)) (is (< (Math/abs i1) 1.0e-9))
      (is (< (Math/abs r2) 1.0e-9))         (is (< (Math/abs (- i2 -1.0)) 1.0e-9))
      (is (< (Math/abs (- r3 -1.0)) 1.0e-9)) (is (< (Math/abs i3) 1.0e-9)))))

(deftest non-power-of-two-identities
  (testing "n=5 impulse -> constant, via the O(n^2) DFT fallback"
    (is (= (col (repeat 5 [1.0 0.0]))
           (dft/fft (col [[1.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]])))))
  (testing "n=3 DC -> single spike"
    (is (= (col [[3.0 0.0] [0.0 0.0] [0.0 0.0]])
           (dft/fft (col [[1.0 0.0] [1.0 0.0] [1.0 0.0]]))))))

(deftest round-trips
  (testing "power of 2 (n=8)"
    (let [v (col (mapv (fn [x] [(double x) (double (- x))]) (range 8)))]
      (is (= v (dft/ifft (dft/fft v))))))
  (testing "not a power of 2 (n=5)"
    (let [v (col (mapv (fn [x] [(double x) (double (* 2 x))]) (range 5)))]
      (is (= v (dft/ifft (dft/fft v))))))
  (testing "n=1 (edge case)"
    (let [v (col [[3.0 4.0]])]
      (is (= v (dft/fft v)))
      (is (= v (dft/ifft v))))))

(deftest row-vector-too
  (is (= (sut/<-rows [[[1.0 0.0] [1.0 0.0] [1.0 0.0] [1.0 0.0]]])
         (dft/fft (sut/<-rows [[[1.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]]])))))

(deftest rejects-non-vectors
  (is (thrown? clojure.lang.ExceptionInfo
               (dft/fft (sut/<-rows [[[1.0 0.0] [2.0 0.0]] [[3.0 0.0] [4.0 0.0]]])))))

(deftest fft2d-identities
  (testing "2D impulse -> all-ones (both dimensions power of 2)"
    (is (= (sut/<-rows (vec (repeat 4 (vec (repeat 4 [1.0 0.0])))))
           (dft/fft2d (sut/<-rows [[[1.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]]
                                    [[0.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]]
                                    [[0.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]]
                                    [[0.0 0.0] [0.0 0.0] [0.0 0.0] [0.0 0.0]]])))))
  (testing "round trip, non-square, mixed complex values (power of 2)"
    (let [M (sut/<-rows [[[1.0 2.0] [3.0 -1.0] [0.0 0.0]]
                          [[2.0 0.0] [1.0 1.0] [-1.0 4.0]]
                          [[0.0 1.0] [2.0 2.0] [3.0 0.0]]
                          [[4.0 -2.0] [0.0 0.0] [1.0 1.0]]])]
      (is (= M (dft/ifft2d (dft/fft2d M))))))
  (testing "round trip, square, neither dimension a power of 2 (O(n^2) fallback in both passes)"
    (let [M (sut/<-rows [[[1.0 0.0] [2.0 1.0] [0.0 -1.0]]
                          [[3.0 2.0] [1.0 0.0] [0.0 1.0]]
                          [[0.0 0.0] [1.0 1.0] [2.0 -2.0]]])]
      (is (= M (dft/ifft2d (dft/fft2d M)))))))
