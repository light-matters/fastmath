(ns fastmath.algebra.function.complex
  (:require [fastmath.default :as default]
            [fastmath.protocol.algebra.coordinate.complex :as cc])
  (:import [java.lang Math]))

(defn real? [z]
  (< (Math/abs ^double (cc/im z)) ^double default/tolerance))


