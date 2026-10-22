(ns fastmath.algebra.function.matrix
  "Matrix operations derived purely from other protocol methods (shape,
  element access, transpose, multiply, determinant, adjoint) rather than
  from any backend's internals.

  Because these only call other protocol methods, the same implementation
  works unmodified for every matrix backend (EJML, ojAlgo, ...) and for
  both the real and complex element domains - there is exactly one place
  to get them right, instead of one per backend per domain.

  These are NOT a replacement for operations whose whole cost is a single
  optimized backend call (add, multiply, transpose, the decompositions,
  unitary?/singular? where a backend offers a dedicated, numerically
  robust primitive) - those stay backend-native for speed and are left
  alone here."
  (:require
   [fastmath.default :as default]
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.protocol.algebra.object.matrix.extra :as emat]
   [fastmath.protocol.algebra.object.matrix.square :as sq]
   [fastmath.protocol.algebra.object.matrix.complex :as cmat]
   [fastmath.protocol.algebra.object.number.complex :as pc]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                      Element-level arithmetic, either domain                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Real matrix elements are plain numbers, which don't satisfy the `Field`
;; protocol (only `ComplexNumber` does) - so these can't just call
;; `field/add`/`field/multiply` directly and need to branch once here,
;; the same way `fastmath.api.v2.algebra.predicate` already does.

(defn- g-add [a b]
  (if (number? a) (+ (double a) (double b)) (pc/add a b)))

(defn- g-multiply [a b]
  (if (number? a) (* (double a) (double b)) (pc/multiply a b)))

(defn- g-conjugate [a]
  (if (pc/? a) (pc/conjugate a) a))

(defn- scalar-magnitude [a]
  (if (pc/? a) (pc/norm a) (Math/abs (double a))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                  Predicates                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn square? [A]
  (apply = (d2/shape A)))

(defn- adjoint-or-transpose
  "Conjugate transpose for a complex matrix, plain transpose for a real one -
  so `normal?` below works unmodified for either domain."
  [A]
  (if (cmat/? A) (cmat/adjoint A) (emat/transpose A)))

(defn symmetric? [A]
  (and (square? A)
       (= A (emat/transpose A))))

(defn hermitian?
  "Only meaningful for complex matrices; a backend's own native check may
  exist, but verify against this before trusting it blindly - ojAlgo
  56.1.0's `BasicMatrix.isHermitian()` was confirmed (empirically, not
  from its docs) to not check that diagonal entries are real, so e.g. a
  diagonal matrix with a non-zero imaginary diagonal passes it although
  it plainly isn't Hermitian."
  [A]
  (and (square? A)
       (= A (cmat/adjoint A))))

(defn normal? [A]
  (and (square? A)
       (let [At (adjoint-or-transpose A)]
         (= (emat/multiply A At) (emat/multiply At A)))))

(defn singular?
  "Whether `determinant` is (approximately) zero. As naive as the
  backend-native versions this replaces - not a substitute for a proper
  rank/condition-number check on ill-conditioned matrices."
  [A]
  (and (square? A)
       (< (scalar-magnitude (sq/determinant A)) ^double default/tolerance)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                           Bilinear / tensor products                        ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn inner
  "Sum over i,j of conj(A[i,j]) * B[i,j]. Returns a plain double for real
  matrices, a ComplexNumber for complex ones - whatever `d2/element`
  itself returns."
  [A B]
  (let [shape-A (d2/shape A)
        shape-B (d2/shape B)]
    (when-not (= shape-A shape-B)
      (throw (ex-info "Shape mismatch for inner product"
                      {:shape-A shape-A :shape-B shape-B})))
    (let [[rows cols] shape-A]
      (reduce g-add
              (for [i (range rows) j (range cols)]
                (g-multiply (g-conjugate (d2/element A i j)) (d2/element B i j)))))))

(defn kronecker-rows
  "A⊗B's entries as nested rows: (A⊗B)[i*pB+p, j*qB+q] = A[i,j]*B[p,q].
  Returns raw nested rows, not a wrapped matrix - the caller already has
  its own zero-dependency-risk `<-rows` constructor to finish the job,
  and building the result here would mean depending on the matrix
  constructor namespace, which depends on every backend in turn."
  [A B]
  (let [[a-rows a-cols] (d2/shape A)
        [b-rows b-cols] (d2/shape B)]
    (for [i (range a-rows) p (range b-rows)]
      (for [j (range a-cols) q (range b-cols)]
        (g-multiply (d2/element A i j) (d2/element B p q))))))
