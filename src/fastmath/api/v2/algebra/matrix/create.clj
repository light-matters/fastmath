(ns fastmath.api.v2.algebra.matrix.create
  "For making matrices, by weird and wonderful ways."
  ;; TODO: Make some real matrices!
  ;; TODO: Make SquareMatrices where appropriate
  ;; TODO: Make implementation choice optional?
  (:require
   [fastmath.protocol.representation.D2 :as d2]
   [fastmath.algebra.object.matrix.rectangular.complex.ejml :as cmat]))

(def domains (set #{:real :complex}))
(def domain--default :complex)

(defn identity [num-rows {:keys [domain]
                          :or {domain domain--default}}])

(defn zero
  ([nrows]
   (zero [nrows nrows] {:domain domain--default}))
  ([[nrows ncols]]
   (zero [nrows ncols] {:domain domain--default}))
  ([[nrows ncols] {:keys [domain] :or {domain domain--default}}]
   (if (not= nrows ncols)
     ;; TODO: implement square options
     (case domain
       :complex (cmat/complexdense n m)
       :real    (rmat/realdense    n m)
       (throw (ex-info "Unknown domain" {:domain domain})))
     (case domain
       :complex (cmat/complexdense n m)
       :real    (rmat/realdense    n m)
       (throw (ex-info "Unknown domain" {:domain domain}))))))

(defn <-coll [[num-rows num-cols] coll {:keys [domain]
                                        :or {domain domain--default}}])

(defn <-rows [rows {:keys [domain]
                    :or {domain domain--default}}])
(defn <-cols [cols {:keys [domain]
                    :or {domain domain--default}}])
(defn <-diagonal [coll {:keys [domain]
                        :or {domain domain--default}}])

;; ==================================================
;; COMPLEX-only
;; ==================================================
(defn <-real
  "A complex matrix from an existing, real matrix."
  [M])
