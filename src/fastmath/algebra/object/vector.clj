(ns fastmath.algebra.object.vector
  (:require
   [fastmath.protocol.algebra.object.matrix.extra :as extra]
   [fastmath.protocol.representation.d2 :as d2]
   [fastmath.algebra.object.matrix.create :as M]))

(def orientations #{:row :column})

(defn create
  "Creates a zero vector of the given length."
  ([coll]
   (create coll :column))
  ([coll type]
   (case type
     :column (M/<-coll (count coll) 1  coll)
     :row (M/<-coll  1 (count coll) coll))))

;; (defn- index-of-first-one [xs]
;;   (loop [xs xs, i 0]
;;     (cond
;;       (empty? xs) nil
;;       (= 1 (first xs)) i
;;       :else (recur (next xs) (unchecked-inc i)))))

(defn orientation [v]
  (let [s (d2/shape v)
        a (s 0)
        b (s 1)]
    (cond
      (= b 1) :column
      (= a 1) :row
      :else   nil)))

