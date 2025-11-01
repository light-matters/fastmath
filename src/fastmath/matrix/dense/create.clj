(ns fastmath.matrix.dense.create
  "All of your favourite ways to make dense matrices."
  (:require
   [fastmath.matrix.dense.complex.ejml :as mc]
   [fastmath.matrix.dense.real.ejml :as mr]
   [fastmath.protocols.matrix2 :as mat]))

;; ===================================================
(def ^:private dense-matrix '{:return-type #{;; mc/ComplexDense
                                             ;; mr/RealDense
                                             }})

(defmacro defn-dense
  "Type-hints the `defn` with `dense-matrix`."
  [name & declaration]
  (let [m (merge (meta name)
                 dense-matrix)]
    `(defn ~(with-meta name m) ~@declaration)))

;; ===================================================

(comment "WIP!")

(defn zero
  [num-rows  num-cols &
   {:keys [complex?]
    :or {complex? false}}]
  (if-not complex?
    (mr/realdense num-rows num-cols)
    (println "not implemented yet!")))

(defn <-elements
  ;; TODO: Make this work with vectors
  [^long num-rows ^long num-cols ^doubles elements]
  (mr/realdense num-rows num-cols elements))

(defn <-rows [rows]
  (mr/realdense<-rows rows))
(defn <-cols [cols]
  (mr/realdense<-cols cols))

(defn diagonal
  [^doubles ss]
  (mr/realdense ss))

(defn I
  [^long num-diag]
  (diagonal (take num-diag (repeat 1.0))))

(defn rotation--2d [rad]
  (println "NIY!"))

(comment (zero 3 3)
         (<-elements 3 3 (double-array [0 0 0 0 0 0 0 0 0.]))
         (<-rows [[1 0] [0 1]])
         (<-cols [[1 0] [0 1]])
         (diagonal [1 2 3])
         (I 2)

         (= (I 3)
            (mat/mul (I 3) (I 3))))

