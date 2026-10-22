(ns fastmath.protocol.analysis.transform.fourier
  "Discrete Fourier transform on a complex vector (a row or column dense
  matrix, per this API's convention - there's no separate vector type).

  Lives under `analysis` rather than `algebra`, alongside the (currently
  empty, but already anticipated) `fastmath.protocol.analysis.decomposition.*`
  tree: a Fourier transform isn't a ring/field/group operation the way
  add/multiply/transpose are - it's a numerical-analysis transform that
  happens to consume and produce matrix-shaped data, same category as a
  matrix decomposition.")

(defprotocol FourierTransform
  (fft [v] "Forward discrete Fourier transform.")
  (ifft [v] "Inverse discrete Fourier transform, normalized by 1/n.")
  (fft2d [M] "Forward 2D discrete Fourier transform of a general matrix (row pass then column pass).")
  (ifft2d [M] "Inverse 2D discrete Fourier transform, normalized by 1/(rows*cols)."))

(defn ? [x]
  (satisfies? FourierTransform x))
