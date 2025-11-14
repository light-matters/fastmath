(ns fastmath.protocols.linear-algebra.complex.matrix)

(defprotocol MatrixComplex
  (adjoint [m])
  (conjugate [m])
  (im
    ;; TODO: Should this return a RealDense or a ComplexDense?
    [m])
  (re
    ;; TODO: Should this return a RealDense or a ComplexDense?
    [m])
  (hermitian? [m])
  (real?
    ;; "Whether or not the matrix only has Real elements, subject to tolerance."
    [m]))
