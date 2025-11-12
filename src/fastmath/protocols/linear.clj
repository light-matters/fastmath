(ns fastmath.protocols.linear)

;; TODO:
;; - Should the protocols be 'Real' and 'Complex' or 'Matrix'Real etc?
;; ===================================================
;; Notes on notation 
;; ===================================================
;; A matrix of doubles is considered the default. Other types are specified according to the legend below. 
;; 
;; --e element
;; --s scalar
;; --m matrix
;; --v vector
;; --vt vector-transpose
;; ===================================================

(defprotocol MatrixReal
  ;; Transformations
  (->seq [m])

  (->array [m])
  (->array--2d [m])
  (->array--float [m])
  (->array--2dfloat [m])

;; Retrieval
  (columns [m])
  (rows [m])

  (diagonal [m])
  (element [m row col])

  (column [m id])
  (row [m id])

  (num-rows [m])
  (num-cols [m])

  (shape [m])

  ;; TODO:
  ;; Subset or slicing??

  ;; Operations
  (add [m1 m2])
  (add--s [m s])

  (sub [m] [m1 m2])

  (map--m [f m]) ;; (Keep predictable argument order)

  (mul
    [m1 m2]
    [m1 t1? m2 t2?])
  (mul--e [m1 m2])
  (mul--v [m v])
  (mul--vt [m v])
  (mul--s [m s])

  (cholesky [m])
  (determinant [m])
  (inverse [m])
  (norm [m t])
  (solve [m v])
  (trace [m])
  (transpose [m])

  (normal? [m])
  (singular? [m])
  (square? [m])
  (symmetric? [m] [m tol])
  (unitary? [m]))

(defprotocol MatrixComplex
  (adjoint [m])
  (conjugate [m])
  (imag
    ;; TODO: Should this return a RealDense or a ComplexDense?
    [m])
  (real
    ;; TODO: Should this return a RealDense or a ComplexDense?
    [m])
  (hermitian? [m])
  (real?
    ;; "Whether or not the matrix only has Real elements, subject to tolerance."
    [m]))

(defprotocol MatrixDecomposition
  (component [m c]))

