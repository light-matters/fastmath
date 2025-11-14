(ns fastmath.protocols.linear-algebra.complex.number)

;; ===================================================
;; Notes on notation 
;; ===================================================
;; Doubles are considered the default. Other types are specified according to the legend below. 
;; 
;; --e element
;; --s scalar
;; --z complex scalar
;; --m matrix
;; --v vector
;; --vt vector-transpose
;; ===================================================
(defprotocol NumberComplex
  (->seq [z])
  (->array [z])

  (conjugate [z])
  (re [z])
  (im [z])
  (phase [z])
  (mag [z])

  (add [z1 z2])
  (sub [z1 z2])
  (mul [z1 z2])
  (div [z1 z2])
  (sqrt [z1])

  (real? [z]))
