(ns quadrature.simpson
  (:require [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.trapezoid :as qt]
            [quadrature.interpolate.richardson :as ir]))

(defn simpson-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed interval $[a, b]$ using Simpson's rule.

  Simpson's rule is equivalent to the trapezoid method subject to one refinement
  of Richardson extrapolation. The trapezoid method fits a line to each
  integration slice. Simpson's rule fits a quadratic to each slice.

  Returns estimates with $n, 2n, 4n, ...$ slices, geometrically increasing by a
  factor of 2 with each estimate.

  ## Optional arguments:

  If supplied, `:n` (default 1) specifies the initial number of slices to use."
  ([f a b] (simpson-sequence f a b {:n 1}))
  ([f a b {:keys [n] :or {n 1}}]
   {:pre [(number? n)]}
   (-> (qt/trapezoid-sequence f a b n)
       (ir/richardson-column 1 2 2 2))))

(qc/defintegrator integral
  "Returns an estimate of the integral of `f` over the closed interval $[a, b]$
  using Simpson's rule with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `simpson-sequence` for more information about Simpson's rule, caveats that
  might apply when using this integration method and information on the optional
  args in `opts` that customize this function's behavior."
  :area-fn (comp first simpson-sequence)
  :seq-fn simpson-sequence)
