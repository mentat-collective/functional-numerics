(ns quadrature.midpoint
  (:require [quadrature.interpolate.richardson :as ir]
            [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.riemann :as qr]
            [sicmutils.generic :as g]
            [quadrature.util :as u]
            [quadrature.util.aggregate :as ua]
            [quadrature.util.stream :as us]))

(defn single-midpoint [f a b]
  (let [width      (g/- b a)
        half-width (g// width 2)
        midpoint   (g/+ a half-width)]
    (g/* width (f midpoint))))

(defn- midpoint-sum* [f a b]
  (let [area-fn (partial single-midpoint f)]
    (qr/windowed-sum area-fn a b)))

(defn- Sn->S3n [f a b]
  (let [width (- b a)]
    (fn [Sn n]
      (let [h        (/ width n)
            delta    (/ h 6)
            l-offset (+ a delta)
            r-offset (+ a (* 5 delta))
            fx (fn [i]
                 (let [ih (* i h)]
                   (+ (f (+ l-offset ih))
                      (f (+ r-offset ih)))))]
        (-> (+ Sn (* h (ua/sum fx 0 n)))
            (/ 3.0))))))

(defn midpoint-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the open interval $(a, b)$ using the Midpoint method.

  ## Optional arguments:

  `:n`: If `:n` is a number, returns estimates with $n, 3n, 9n, ...$ slices,
  geometrically increasing by a factor of 3 with each estimate.

  If `:n` is a sequence, the resulting sequence will hold an estimate for each
  integer number of slices in that sequence.

  `:accelerate?`: if supplied (and `n` is a number), attempts to accelerate
  convergence using Richardson extrapolation. If `n` is a sequence this option
  is ignored."
  ([f a b] (midpoint-sequence f a b {:n 1}))
  ([f a b {:keys [n accelerate?] :or {n 1}}]
   (let [S      (qr/midpoint-sum f a b)
         next-S (Sn->S3n f a b)
         xs     (qr/incrementalize S next-S 3 n)]
     (if (and accelerate? (number? n))
       (ir/richardson-sequence xs 3 2 2)
       xs))))

(qc/defintegrator integral
  "Returns an estimate of the integral of `f` over the open interval $(a, b)$
  using the Midpoint method with $1, 3, 9 ... 3^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `midpoint-sequence` for information on the optional args in `opts` that
  customize this function's behavior."
  :area-fn single-midpoint
  :seq-fn midpoint-sequence)
