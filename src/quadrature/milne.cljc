(ns quadrature.milne
  (:require [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [quadrature.midpoint :as qm]
            [quadrature.interpolate.richardson :as ir]
            [quadrature.util.stream :as us]))

(defn milne-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the open interval $(a, b)$ using Milne's rule.

  Milne's rule is equivalent to the midpoint method subject to one refinement of
  Richardson extrapolation.

  Returns estimates with $n, 2n, 4n, ...$ slices, geometrically increasing by a
  factor of 2 with each estimate.

  ## Optional arguments:

  If supplied, `:n` (default 1) specifies the initial number of slices to use.

  NOTE: the Midpoint method is able to reuse function evaluations as its windows
  narrow /only/ when increasing the number of integration slices by 3. Milne's
  method increases the number of slices geometrically by a factor of 2 each
  time, so it will never hit the incremental path. You may want to memoize your
  function before calling `milne-sequence`."
  ([f a b] (milne-sequence f a b {:n 1}))
  ([f a b {:keys [n] :or {n 1} :as opts}]
   {:pre [(number? n)]}
   (-> (qm/midpoint-sequence f a b (assoc opts :n (us/powers 2 n)))
       (ir/richardson-column 1 2 2 2))))

(qc/defintegrator integral
  "Returns an estimate of the integral of `f` over the open interval $(a, b)$
  using Milne's rule with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `milne-sequence` for more information about Milne's rule, caveats that
  might apply when using this integration method and information on the optional
  args in `opts` that customize this function's behavior."
  :area-fn (comp first milne-sequence)
  :seq-fn milne-sequence)
