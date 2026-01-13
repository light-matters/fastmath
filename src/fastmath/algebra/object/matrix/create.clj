(ns fastmath.algebra.object.matrix.create
  "For making matrices, by weird and wonderful ways.

  You can change the implementation choice in the require statement and, assuming you followed the conventions, everything else should just work.
  "
  ;; TODO: Make SquareMatrices where appropriate (currently just placeholders)
  ;; TODO: Make fixed sizes matrices where appropriate
  ;; TODO: Make implementation choice optional?
  (:require
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.protocol.algebra.object.number.complex :as C]
   [fastmath.algebra.object.matrix.rectangular.real.ejml :as rmat]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as cmat]
   [fastmath.algebra.object.matrix.square.real.ejml :as sqrmat]
   [fastmath.algebra.object.matrix.square.complex.ejml :as sqcmat]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;              Identities             ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def domains (set #{:real :complex}))
(def domain--default :real)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;              Constructors             ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn identity
  "A square matrix, with unity on the leading diagonal and all other elements zero. "
  ([nrows] (identity nrows {:domain domain--default}))
  ([nrows {:keys [domain]
           :or {domain domain--default}}]
   (case domain
     :real    (sqrmat/realdense (take nrows (repeat 1)))
     :complex (sqcmat/complexdense  (interleave (take nrows (repeat 1)) (take nrows (repeat 0))))
     (throw (ex-info "Unknown domain" {:domain domain})))))

(defn zero
  "A matrix with all elements zeroed, as specified by the number of rows, rows and columns or rows, columns and optional info. - currently domain."
  ([nrows]
   (^:long zero nrows nrows {:domain domain--default}))
  ([^:long nrows ^:long ncols]
   (zero nrows ncols {:domain domain--default}))
  ([^:long nrows ^:long ncols {:keys [domain] :or {domain domain--default}}]
   (if (= nrows ncols)
     ;; TODO: Implement proper square and fixed-size matrices
     (case domain
       :real    (sqrmat/realdense    nrows ncols)
       :complex (sqcmat/complexdense nrows ncols)
       (throw (ex-info "Unknown domain" {:domain domain})))
     (case domain
       :real    (rmat/realdense    nrows ncols)
       :complex (cmat/complexdense nrows ncols)
       (throw (ex-info "Unknown domain" {:domain domain}))))))

(defn- complex-elements?
  "Check for the different ways that complex elements could be represented."
  [[nrows ncols] coll]
  (let [coll0 (first coll)]
    (or (C/? coll0)
        (and (coll? coll0) (= (count coll0) 2))
        (= (count coll) (* 2 nrows ncols)))))

(defn <-coll
  "From 1-D collection. Chooses domain based on the first element. Don't mix real and complex numbers!"
  [nrows ncols coll]
  (let [domain (if-not (complex-elements? [nrows ncols] coll)
                 :real :complex)]
    (if (= nrows ncols)
      (case domain
        :real    (sqrmat/realdense    nrows ncols coll)
        :complex (sqcmat/complexdense nrows ncols coll)
        (throw (ex-info "Unknown domain" {:domain domain})))
      (case domain
        :real    (rmat/realdense    nrows ncols coll)
        :complex (cmat/complexdense nrows ncols coll)
        (throw (ex-info "Unknown domain" {:domain domain}))))))

(defn <-rows [rows]
  (let [domain (if-not (C/? (first (first rows))) :real :complex)]
    (if (= (count rows) (count (first rows)))
      (case domain
        :real    (sqrmat/<-rows rows)
        :complex (sqcmat/<-rows rows)
        (throw (ex-info "Unknown domain" {:domain domain})))
      (case domain
        :real    (rmat/<-rows rows)
        :complex (cmat/<-rows rows)
        (throw (ex-info "Unknown domain" {:domain domain}))))))

(defn <-cols [cols]
  (let [domain (if-not (C/? (first (first cols))) :real :complex)]
    (if (= (count cols) (count (first cols)))
      (case domain
        :real    (sqrmat/<-cols cols)
        :complex (sqcmat/<-cols cols)
        (throw (ex-info "Unknown domain" {:domain domain})))
      (case domain
        :real    (rmat/<-cols cols)
        :complex (cmat/<-cols cols)
        (throw (ex-info "Unknown domain" {:domain domain}))))))

(defn diagonal
  "From 1-D collection, creates a matrix with zero-values apart from on the leading diagonal. Chooses domain based on the first element. Don't mix real and complex numbers!"
  [coll]
  (let [domain (if-not (C/? (first coll)) :real :complex)]
    ;; TODO: Use square matrix types
    (case domain
      :real    (rmat/realdense     coll)
      :complex (cmat/complexdense  coll)
      (throw (ex-info "Unknown domain" {:domain domain})))))

;; ==================================================
;; COMPLEX-only
;; ==================================================
(defn <-real
  "A complex matrix from an existing, real matrix."
  [M]
  (if-not (instance? fastmath.algebra.object.matrix.rectangular.complex.ejml.ComplexDense M)
    (let [[nrows ncols] (d2/shape M)]
      (if (= nrows ncols)
        (sqcmat/<-real M)
        (cmat/<-real M)))
    M))

