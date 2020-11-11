(ns quadrature.infinite
  (:require [clojure.core.match :refer [match]]
            [quadrature.common :as qc]
            [quadrature.substitute :as qs]))

(defn improper
  "Accepts:

  - An `integrator` (function of `f`, `a`, `b` and `opts`)
  - `a` and `b`, the endpoints of an integration interval, and
  - (optionally) `opts`, a dict of integrator-configuring options

  And returns a new integrator that's able to handle infinite endpoints. (If you
  don't specify `##-Inf` or `##Inf`, the returned integrator will fall through
  to the original `integrator` implementation.)

  All `opts` will be passed through to the supplied `integrator`.

  ## Optional arguments relevant to `improper`:

  `:infinite-breakpoint`: If either `a` or `b` is equal to `##Inf` or `##-Inf`,
  this function will internally perform a change of variables on the regions
  from:

  `(:infinite-breakpoint opts) => ##Inf`

  or

  `##-Inf => (- (:infinite-breakpoint opts))`

  using $u(t) = {1 \\over t}$, as described in the `infinitize` method of
  `substitute.cljc`. This has the effect of mapping the infinite endpoint to an
  open interval endpoint of 0.

  Where should you choose the breakpoint? According to Press in Numerical
  Recipes, section 4.4: \"At a sufficiently large positive value so that the
  function funk is at least beginning to approach its asymptotic decrease to
  zero value at infinity.\"

  References:

  - Press, Numerical Recipes (p138), Section 4.4: http://phys.uri.edu/nigh/NumRec/bookfpdf/f4-4.pdf"
  [integrator]
  (fn rec
    ([f a b] (rec f a b {}))
    ([f a b opts]
     (let [{:keys [infinite-breakpoint] :as opts} (-> {:infinite-breakpoint 1}
                                                      (merge opts))
           call (fn [integrate l r interval]
                  (let [m (qc/with-interval opts interval)]
                    (let [result (integrate f l r m)]
                      (:result result))))
           ab-interval   (qc/interval opts)
           integrate     (partial call integrator)
           inf-integrate (partial call (qs/infinitize integrator))
           r-break       (Math/abs infinite-breakpoint)
           l-break       (- r-break)]
       (match [[a b]]
              [(:or [##-Inf ##-Inf] [##Inf ##Inf])]
              {:converged? true
               :terms-checked 0
               :result 0.0}

              [(:or [_ ##-Inf] [##Inf _])]
              (-> (rec f b a opts)
                  (update-in [:result] -))





              [[##-Inf ##Inf]]
              (let [-inf->l (inf-integrate a l-break qc/open-closed)
                    l->r    (integrate     l-break r-break qc/closed)
                    r->+inf (inf-integrate r-break b qc/closed-open)]
                {:converged? true
                 :result (+ -inf->l l->r r->+inf)})




              [[##-Inf _]]
              (if (<= b l-break)
                (inf-integrate a b ab-interval)
                (let [-inf->l (inf-integrate a l-break qc/open-closed)
                      l->b    (integrate     l-break b (qc/close-l ab-interval))]
                  {:converged? true
                   :result (+ -inf->l l->b)}))




              [[_ ##Inf]]
              (if (>= a r-break)
                (inf-integrate a b ab-interval)
                (let [a->r    (integrate     a r-break (qc/close-r ab-interval))
                      r->+inf (inf-integrate r-break b qc/closed-open)]
                  {:converged? true
                   :result (+ a->r r->+inf)}))





              :else (integrator f a b opts))))))
