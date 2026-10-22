(ns fastmath.algebra.object.matrix.rectangular.complex.fourier.ejml
  "Wires the discrete Fourier transform kernel
  (`fastmath.analysis.transform.fourier`) onto EJML's `ComplexDense`.

  Operates directly on `ComplexDense`'s raw interleaved re/im double
  array rather than going through `d2/element` per cell - the same
  choice `multiply--e`/`fmap` already make in the EJML backend itself,
  and doubly so here since an FFT's whole point is O(n log n), which a
  per-element protocol-dispatch-plus-ComplexNumber-allocation inner loop
  would undermine.

  A separate file via `extend-type` rather than more methods added to
  the already-large `complex/ejml.clj` deftype: this is a numerical
  transform, not a matrix-algebra operation (see
  `fastmath.protocol.analysis.transform.fourier`'s docstring), and
  neither `complex/ejml.clj` nor the `FourierTransform` protocol need to
  know about each other for this to work - `ComplexDense`'s backing
  field is public (the Clojure deftype default), so extending a new
  protocol onto it from here needs no changes to the type's own
  namespace at all.
  "
  (:require
   [fastmath.analysis.transform.fourier :as dft]
   [fastmath.protocol.analysis.transform.fourier :as proto]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as ejml])
  (:import
   (fastmath.algebra.object.matrix.rectangular.complex.ejml ComplexDense)
   (org.ejml.data ZMatrixRMaj)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; ==================================================
;; Utilities
;; ==================================================

(defn- vector-length
  "Validates M is a row or column vector and returns its length."
  ^long [^ZMatrixRMaj M]
  (let [rows (.numRows M) cols (.numCols M)]
    (when-not (or (== rows 1) (== cols 1))
      (throw (ex-info "Fourier transform requires a row or column vector"
                      {:shape [rows cols]})))
    (max rows cols)))

(defn- deinterleave
  "M's backing array is already contiguous re/im pairs for a row or
  column vector (only one row or one column, so no stride to account
  for) - split it into parallel re/im double arrays."
  [^ZMatrixRMaj M ^long n]
  (let [^doubles data (.getData M)
        re (double-array n)
        im (double-array n)]
    (dotimes [i n]
      (aset-double re i (aget data (* 2 i)))
      (aset-double im i (aget data (inc (* 2 i)))))
    [re im]))

(defn- ->matrix
  "Builds a new ZMatrixRMaj with the same shape as M, filled with the
  given re/im pair, interleaved back the way M itself stores them."
  ^ZMatrixRMaj [^ZMatrixRMaj M ^doubles re ^doubles im ^long n]
  (let [out (ZMatrixRMaj. (.numRows M) (.numCols M))
        ^doubles data (.getData out)]
    (dotimes [i n]
      (aset-double data (* 2 i) (aget re i))
      (aset-double data (inc (* 2 i)) (aget im i)))
    out))

(defn- transform-rows!
  "Applies `pass` (`dft/fft` or `dft/ifft`) along every row of `data` in
  place. Each row is already contiguous in row-major interleaved storage,
  so no striding is needed here."
  [^doubles data ^long rows ^long cols pass]
  (let [re (double-array cols) im (double-array cols)]
    (dotimes [i rows]
      (let [row-base (* i cols)]
        (dotimes [j cols]
          (let [idx (* 2 (+ row-base j))]
            (aset-double re j (aget data idx))
            (aset-double im j (aget data (inc idx)))))
        (let [[^doubles out-re ^doubles out-im] (pass re im)]
          (dotimes [j cols]
            (let [idx (* 2 (+ row-base j))]
              (aset-double data idx (aget out-re j))
              (aset-double data (inc idx) (aget out-im j)))))))))

(defn- transform-cols!
  "Applies `pass` along every column of `data` in place. Columns are
  strided (stride = cols) in row-major interleaved storage."
  [^doubles data ^long rows ^long cols pass]
  (let [re (double-array rows) im (double-array rows)]
    (dotimes [j cols]
      (dotimes [i rows]
        (let [idx (* 2 (+ (* i cols) j))]
          (aset-double re i (aget data idx))
          (aset-double im i (aget data (inc idx)))))
      (let [[^doubles out-re ^doubles out-im] (pass re im)]
        (dotimes [i rows]
          (let [idx (* 2 (+ (* i cols) j))]
            (aset-double data idx (aget out-re i))
            (aset-double data (inc idx) (aget out-im i))))))))

(defn- transform2d
  "Separable 2D transform: `pass` applied along every row, then along
  every column. For `dft/ifft` this is exactly the correct 1/(rows*cols)
  2D normalization, since each pass independently normalizes by 1/length
  and 1/cols * 1/rows = 1/(rows*cols)."
  ^ZMatrixRMaj [^ZMatrixRMaj M pass]
  (let [rows (.numRows M) cols (.numCols M)
        data (aclone ^doubles (.getData M))]
    (transform-rows! data rows cols pass)
    (transform-cols! data rows cols pass)
    (let [out (ZMatrixRMaj. rows cols)]
      (System/arraycopy data 0 (.getData out) 0 (alength data))
      out)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                  Wiring                                     ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(extend-type ComplexDense
  proto/FourierTransform
  (fft [this]
    (let [^ZMatrixRMaj M (.-M this)
          n (vector-length M)
          [re im] (deinterleave M n)
          [out-re out-im] (dft/fft re im)]
      (ejml/->ComplexDense (->matrix M out-re out-im n))))
  (ifft [this]
    (let [^ZMatrixRMaj M (.-M this)
          n (vector-length M)
          [re im] (deinterleave M n)
          [out-re out-im] (dft/ifft re im)]
      (ejml/->ComplexDense (->matrix M out-re out-im n))))
  (fft2d [this]
    (let [^ZMatrixRMaj M (.-M this)]
      (ejml/->ComplexDense (transform2d M dft/fft))))
  (ifft2d [this]
    (let [^ZMatrixRMaj M (.-M this)]
      (ejml/->ComplexDense (transform2d M dft/ifft)))))
