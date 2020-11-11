(ns quadrature.trapezoid
  "Trapezoid method."
  (:require [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.riemann :as qr]
            [quadrature.interpolate.richardson :as ir]
            [sicmutils.function :as f]
            [sicmutils.generic :as g]
            [quadrature.util :as u]
            [quadrature.util.aggregate :as ua]
            [quadrature.util.stream :as us]))
