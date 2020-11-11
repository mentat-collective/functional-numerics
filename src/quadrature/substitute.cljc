(ns quadrature.substitute
  "## U Substitution and Variable Changes

  This namespace provides implementations of functions that accept an
  `integrator` and perform a variable change to address some singularity, like
  an infinite endpoint, in the definite integral.

  The strategies currently implemented were each described by Press, et al. in
  section 4.4 of ['Numerical
  Recipes'](http://phys.uri.edu/nigh/NumRec/bookfpdf/f4-4.pdf)."
  (:require [clojure.core.match :refer [match]]
            [quadrature.common :as qc]))

(defn infinitize
  "Performs a variable substitution targeted at turning a single infinite endpoint
  of an improper integral evaluation an (open) endpoint at 0 by applying the
  following substitution:

  $$u(t) = {1 \\over t}$$ $$du = {-1 \\over t^2}$$

  This works when the integrand `f` falls off at least as fast as $1 \\over t^2$
  as it approaches the infinite limit.

  The returned function requires that `a` and `b` have the same sign, ie:

  $$ab > 0$$

  Transform the bounds with $u(t)$, and cancel the negative sign by changing
  their order:

  $$\\int_{a}^{b} f(x) d x=\\int_{1 / b}^{1 / a} \\frac{1}{t^{2}} f\\left(\\frac{1}{t}\\right) dt$$

  References:

  - Mathworld, \"Improper Integral\": https://mathworld.wolfram.com/ImproperIntegral.html
  - Press, Numerical Recipes, Section 4.4: http://phys.uri.edu/nigh/NumRec/bookfpdf/f4-4.pdf"
  [integrate]
  (fn call
    ([f a b] (call f a b {}))
    ([f a b opts]
     {:pre [(not
             (and (qc/infinite? a)
                  (qc/infinite? b)))]}
     (let [f' (fn [t]
                (/ (f (/ 1.0 t))
                   (* t t)))
           a' (if (qc/infinite? b) 0.0 (/ 1.0 b))
           b' (if (qc/infinite? a) 0.0 (/ 1.0 a))
           opts (qc/update-interval opts qc/flip)]
       (integrate f' a' b' opts)))))

(defn- inverse-power-law
  "Implements a change of variables to address a power law singularity at the
  lower or upper integration endpoint.

  An \"inverse power law singularity\" means that the integrand diverges as

  $$(x - a)^{-\\gamma}$$

  near $x=a$. Passing true for `lower?` to specify a singularity at the lower
  endpoint, false to signal an upper-endpoint singularity.

  References:

  - Mathworld, \"Improper Integral\": https://mathworld.wolfram.com/ImproperIntegral.html
  - Press, Numerical Recipes, Section 4.4: http://phys.uri.edu/nigh/NumRec/bookfpdf/f4-4.pdf
  - Wikipedia, \"Finite-time Singularity\": https://en.wikipedia.org/wiki/Singularity_(mathematics)#Finite-time_singularity
  "
  [integrate gamma lower?]
  {:pre [(<= 0 gamma 1)]}
  (fn call
    ([f a b] (call f a b {}))
    ([f a b opts]
     (let [inner-pow (/ 1 (- 1 gamma))
           gamma-pow (* gamma inner-pow)
           a' 0
           b' (Math/pow (- b a) (- 1 gamma))
           t->t' (if lower?
                   (fn [t] (+ a (Math/pow t inner-pow)))
                   (fn [t] (- b (Math/pow t inner-pow))))
           f' (fn [t] (* (Math/pow t gamma-pow)
                        (f (t->t' t))))]
       (-> (integrate f' a' b' opts)
           (update-in [:result] (partial * inner-pow)))))))

(defn inverse-power-law-lower [integrate gamma]
  (inverse-power-law integrate gamma true))

(defn inverse-power-law-upper [integrate gamma]
  (inverse-power-law integrate gamma false))

(defn inverse-sqrt-lower
  "Implements a change of variables to address an inverse square root singularity
  at the lower integration endpoint. Use this when the integrand diverges as

  $$1 \\over {\\sqrt{x - a}}$$

  near the lower endpoint $a$."
  [integrate]
  (fn call
    ([f a b] (call f a b {}))
    ([f a b opts]
     (let [f' (fn [t] (* t (f (+ a (* t t)))))]
       (-> (integrate f' 0 (Math/sqrt (- b a)) opts)
           (update-in [:result] (partial * 2)))))))

(defn inverse-sqrt-upper
  "Implements a change of variables to address an inverse square root singularity
  at the upper integration endpoint. Use this when the integrand diverges as

  $$1 \\over {\\sqrt{x - b}}$$

  near the upper endpoint $b$."
  [integrate]
  (fn call
    ([f a b] (call f a b {}))
    ([f a b opts]
     (let [f' (fn [t] (* t (f (- b (* t t)))))]
       (-> (integrate f' 0 (Math/sqrt (- b a)) opts)
           (update-in [:result] (partial * 2)))))))

(defn exponential-upper
  "Implements a change of variables to address an exponentially diverging upper
  integration endpoint. Use this when the integrand diverges as $\\exp{x}$ near
  the upper endpoint $b$."
  [integrate]
  (fn call
    ([f a b] (call f a b {}))
    ([f a b opts]
     {:pre [(qc/infinite? b)]}
     (let [f' (fn [t] (* (- (Math/log t))
                        (/ 1 t)))
           opts (qc/update-interval opts qc/flip)]
       (integrate f 0 (Math/exp (- a)) opts)))))
