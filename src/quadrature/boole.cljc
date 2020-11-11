(ns quadrature.boole
  (:require [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.trapezoid :as qt]
            [quadrature.interpolate.richardson :as ir]))

(defn boole-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed interval $[a, b]$ using Boole's rule.

  Boole's rule is equivalent to the trapezoid method subject to two refinements
  of Richardson extrapolation. The trapezoid method fits a line to each
  integration slice. Boole's rule fits a quartic to each slice.

  Returns estimates with $n, 2n, 4n, ...$ slices, geometrically increasing by a
  factor of 2 with each estimate.

  ## Optional arguments:

  If supplied, `:n` (default 1) specifies the initial number of slices to use."
  ([f a b] (boole-sequence f a b {:n 1}))
  ([f a b {:keys [n] :or {n 1}}]
   {:pre [(number? n)]}
   (-> (qt/trapezoid-sequence f a b n)
       (ir/richardson-column 2 2 2 2))))

(qc/defintegrator integral
  "Returns an estimate of the integral of `f` over the closed interval $[a, b]$
  using Boole's rule with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `boole-sequence` for more information about Boole's rule, caveats that
  might apply when using this integration method and information on the optional
  args in `opts` that customize this function's behavior."
  :area-fn (comp first boole-sequence)
  :seq-fn boole-sequence)
