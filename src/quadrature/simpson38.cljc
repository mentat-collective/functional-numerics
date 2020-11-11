(ns quadrature.simpson38
  (:require [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.trapezoid :as qt]
            [quadrature.interpolate.richardson :as ir]
            [quadrature.util.stream :as us]))

(defn simpson38-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed interval $[a, b]$ using Simpson's 3/8 rule.

  Simpson's 3/8 rule is equivalent to the trapezoid method subject to:

  - one refinement of Richardson extrapolation, and

  - a geometric increase of integration slices by a factor of 3 for each
  sequence element. (the Trapezoid method increases by a factor of 2 by
  default.)

  The trapezoid method fits a line to each integration slice. Simpson's 3/8 rule
  fits a cubic to each slice.

  Returns estimates with $n, 3n, 9n, ...n3^i$ slices, geometrically increasing by a
  factor of 3 with each estimate.

    ## Optional arguments:

  If supplied, `:n` (default 1) specifies the initial number of slices to use.

  NOTE: the Trapezoid method is able to reuse function evaluations as its
  windows narrow /only/ when increasing the number of integration slices by 2.
  Simpson's 3/8 rule increases the number of slices geometrically by a factor of
  2 each time, so it will never hit the incremental path. You may want to
  memoize your function before calling `simpson38-sequence`."
  ([f a b] (simpson38-sequence f a b {:n 1}))
  ([f a b {:keys [n] :or {n 1}}]
   {:pre [(number? n)]}
   (-> (qt/trapezoid-sequence f a b (us/powers 3 n))
       (ir/richardson-column 1 3 2 2))))

(qc/defintegrator integral
  "Returns an estimate of the integral of `f` over the closed interval $[a, b]$
  using Simpson's 3/8 rule with $1, 3, 9 ... 3^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `simpson38-sequence` for more information about Simpson's 3/8 rule, caveats
  that might apply when using this integration method and information on the
  optional args in `opts` that customize this function's behavior."
  :area-fn (comp first simpson38-sequence)
  :seq-fn simpson38-sequence)
