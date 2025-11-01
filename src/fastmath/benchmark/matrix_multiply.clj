(ns fastmath.benchmark.matrix-multiply
  (:require [fastmath.matrix :as m]
            [fastmath.vector :as v]))

(def thing
  (m/rows->mat [1 0 0 0]
               [0 1 0 0]
               [0 0 1 0]
               [0 0 0 1]))

(m/mulv (m/zero 4 3 true) (v/vec->RealVector [0 0 1]))

(comment
  (println thing))
