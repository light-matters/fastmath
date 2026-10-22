(ns fastmath.analysis.transform.fourier
  "A backend-independent discrete Fourier transform kernel: an iterative
  radix-2 Cooley-Tukey FFT (O(n log n)) for power-of-2 sizes, and a
  direct O(n^2) summation for any other size - the same dual-path design
  ojAlgo's own `DiscreteFourierTransform` uses, since there's no single
  O(n log n) algorithm that works for arbitrary n without extra machinery
  (Bluestein's algorithm, mixed-radix factoring, ...) that isn't
  warranted here.

  Operates on parallel `re`/`im` double arrays rather than fastmath's own
  `ComplexNumber`, to avoid a complex-object allocation per element in
  what's meant to be a fast, performance-sensitive transform - any
  backend can convert to/from this at its own element-access boundary,
  the same boundary-conversion principle the matrix backends already
  follow for their own element type.
  "
  (:import (java.lang Math)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn- power-of-two? [^long n]
  (and (pos? n) (zero? (bit-and n (dec n)))))

(defn- bit-reverse ^long [^long x ^long bits]
  (loop [x x i 0 r 0]
    (if (>= i bits)
      r
      (recur (unsigned-bit-shift-right x 1) (inc i) (bit-or (bit-shift-left r 1) (bit-and x 1))))))

(defn- fft-radix2!
  "In-place iterative radix-2 Cooley-Tukey FFT on `re`/`im` (length a power
  of 2). `sign` is -1.0 for the forward transform, +1.0 for the inverse
  (the sign of the exponent in exp(sign * i * theta)) - not normalized."
  [^doubles re ^doubles im ^long n ^double sign]
  (let [bits (Long/numberOfTrailingZeros n)]
    ;; bit-reversal permutation
    (dotimes [i n]
      (let [j (bit-reverse i bits)]
        (when (< i j)
          (let [tr (aget re i) ti (aget im i)]
            (aset-double re i (aget re j)) (aset-double im i (aget im j))
            (aset-double re j tr) (aset-double im j ti)))))
    ;; butterfly stages
    (loop [size 2]
      (when (<= size n)
        (let [half (quot size 2)
              theta (/ (* sign 2.0 Math/PI) size)
              wr-step (Math/cos theta)
              wi-step (Math/sin theta)]
          (loop [block 0]
            (when (< block n)
              (loop [j 0 wr 1.0 wi 0.0]
                (when (< j half)
                  (let [i1 (+ block j) i2 (+ i1 half)
                        xr (aget re i2) xi (aget im i2)
                        tr (- (* wr xr) (* wi xi))
                        ti (+ (* wr xi) (* wi xr))
                        ur (aget re i1) ui (aget im i1)]
                    (aset-double re i1 (+ ur tr)) (aset-double im i1 (+ ui ti))
                    (aset-double re i2 (- ur tr)) (aset-double im i2 (- ui ti))
                    (recur (inc j)
                           (- (* wr wr-step) (* wi wi-step))
                           (+ (* wr wi-step) (* wi wr-step))))))
              (recur (+ block size))))
          (recur (* size 2)))))))

(defn- dft
  "Direct O(n^2) DFT/IDFT for arbitrary n - not in place, returns a fresh
  [re im] pair, inputs untouched."
  [^doubles re ^doubles im ^long n ^double sign]
  (let [out-re (double-array n)
        out-im (double-array n)]
    (dotimes [k n]
      (let [base (/ (* sign 2.0 Math/PI (double k)) n)]
        (loop [j 0 sr 0.0 si 0.0]
          (if (< j n)
            (let [theta (* base j)
                  c (Math/cos theta) s (Math/sin theta)
                  xr (aget re j) xi (aget im j)]
              (recur (inc j)
                     (+ sr (- (* xr c) (* xi s)))
                     (+ si (+ (* xr s) (* xi c)))))
            (do (aset-double out-re k sr) (aset-double out-im k si))))))
    [out-re out-im]))

(defn transform
  "Forward (`sign` -1.0) or un-normalized inverse (`sign` +1.0) discrete
  Fourier transform of the complex vector given by parallel `re`/`im`
  arrays. Returns a fresh [re im] pair; never mutates the inputs."
  [^doubles re ^doubles im ^double sign]
  (let [n (alength re)]
    (if (power-of-two? n)
      (let [re2 (aclone re) im2 (aclone im)]
        (fft-radix2! re2 im2 n sign)
        [re2 im2])
      (dft re im n sign))))

(defn fft
  "Forward discrete Fourier transform. Returns a fresh [re im] pair."
  [re im]
  (transform re im -1.0))

(defn ifft
  "Inverse discrete Fourier transform, normalized by 1/n. Returns a
  fresh [re im] pair."
  [re im]
  (let [n (alength ^doubles re)
        [^doubles out-re ^doubles out-im] (transform re im 1.0)]
    (dotimes [i n]
      (aset-double out-re i (/ (aget out-re i) n))
      (aset-double out-im i (/ (aget out-im i) n)))
    [out-re out-im]))
