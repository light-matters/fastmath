(ns fastmath.matrix.dense.create
  "All of your favourite ways to make dense matrices."
  (:require
   [fastmath.matrix.dense.complex.ejml :as mc] [fastmath.matrix.dense.real.ejml :as mr]))

;; ===================================================
(def ^:private dense-matrix '{:return-type #{mc/ComplexDense mr/RealDense}})

(defmacro defn-dense
  "Type-hints the `defn` with `dense-matrix`."
  [name & declaration]
  (let [m (merge (meta name)
                 dense-matrix)]
    (println (meta name))
    (println m)
    `(defn ~(with-meta name m) ~@declaration)))

;; ===================================================

(comment "WIP!")

(defn-dense zero
  [num-rows num-cols])

(defn-dense diagonal
  ([^floats ss])
  ([^long num-diagonal ^float s]))

(defn-dense I
  [^long num-diagonal]
  (diagonal num-diagonal 1.0))

(defn-dense rotation--2d [])
(defn-dense rotation--3d [])

(defn-dense <-rows [])
(defn-dense <-cols [])
