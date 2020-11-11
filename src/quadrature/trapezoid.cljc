(ns quadrature.trapezoid
  "Trapezoid method."
  (:require [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.riemann :as qr]
            [quadrature.interpolate.richardson :as ir]
            [sicmutils.function :as f]
            [sicmutils.generic :as g]
            [sicmutils.util :as u]
            [quadrature.util.aggregate :as ua]
            [quadrature.util.stream :as us]))

(defn single-trapezoid [f xl xr]
  (g// (g/* (g/- xr xl)
            (g/+ (f xl) (f xr)))
       2))

(defn- trapezoid-sum* [f a b]
  (qr/windowed-sum (partial single-trapezoid f)
                   a b))

(def ^:private pi-estimator*
  (let [f (fn [x] (/ 4 (+ 1 (* x x))))]
    (trapezoid-sum* f 0.0 1.0)))

(defn trapezoid-sum
  "Returns a function of `n`, some number of slices of the total integration
  range, that returns an estimate for the definite integral of $f$ over the
  range $(a, b)$ using the trapezoid method."
  [f a b]
  (let [width (- b a)]
    (fn [n]
      (let [h  (/ width n)
            fx (fn [i] (f (+ a (* i h))))]
        (* h (+ (/ (+ (f a) (f b)) 2)
                (ua/sum fx 1 n)))))))

(def ^:private pi-estimator
  (let [f (fn [x] (/ 4 (+ 1 (* x x))))]
    (trapezoid-sum* f 0.0 1.0)))

(defn trapezoid-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the open interval $(a, b)$ using the Trapezoid method.

  ## Optional arguments:

  `:n`: If `:n` is a number, returns estimates with $n, 2n, 4n, ...$ slices,
  geometrically increasing by a factor of 2 with each estimate.

  If `:n` is a sequence, the resulting sequence will hold an estimate for each
  integer number of slices in that sequence.

  `:accelerate?`: if supplied (and `n` is a number), attempts to accelerate
  convergence using Richardson extrapolation. If `n` is a sequence this option
  is ignored."
  ([f a b] (trapezoid-sequence f a b {:n 1}))
  ([f a b {:keys [n accelerate?] :or {n 1}}]
   (let [S      (trapezoid-sum f a b)
         next-S (qr/Sn->S2n f a b)
         xs     (qr/incrementalize S next-S 2 n)]
     (if (and accelerate? (number? n))
       (ir/richardson-sequence xs 2 2 2)
       xs))))

(qc/defintegrator integral
  "Returns an estimate of the integral of `f` over the closed interval $[a, b]$
  using the Trapezoid method with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `trapezoid-sequence` for information on the optional args in `opts` that
  customize this function's behavior."
  :area-fn single-trapezoid
  :seq-fn trapezoid-sequence)
