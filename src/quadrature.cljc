(ns quadrature
  (:require [sicmutils.numerical.compile :as c]
            [quadrature.adaptive :as qa]
            [quadrature.boole :as boole]
            [quadrature.common :as qc]
            [quadrature.bulirsch-stoer :as bs]
            [quadrature.infinite :as qi]
            [quadrature.midpoint :as mid]
            [quadrature.milne :as milne]
            [quadrature.riemann :as riemann]
            [quadrature.romberg :as romberg]
            [quadrature.simpson :as simp]
            [quadrature.simpson38 :as simp38]
            [quadrature.trapezoid :as trap]
            [quadrature.util :as u]))

(def ^:private quadrature-methods
  {:open                    {:method :adaptive-bulirsch-stoer
                             :interval qc/open}
   :closed                  {:method :adaptive-bulirsch-stoer
                             :interval qc/closed}
   :closed-open             {:method :adaptive-bulirsch-stoer
                             :interval qc/closed-open}
   :open-closed             {:method :adaptive-bulirsch-stoer
                             :interval qc/open-closed}
   :bulirsch-stoer-open     bs/open-integral
   :bulirsch-stoer-closed   bs/closed-integral
   :adaptive-bulirsch-stoer (qa/adaptive bs/open-integral bs/closed-integral)
   :left-riemann            riemann/left-integral
   :right-riemann           riemann/right-integral
   :lower-riemann           riemann/lower-integral
   :upper-riemann           riemann/upper-integral
   :midpoint                mid/integral
   :trapezoid               trap/integral
   :boole                   boole/integral
   :milne                   milne/integral
   :simpson                 simp/integral
   :simpson38               simp38/integral
   :romberg                 romberg/closed-integral
   :romberg-open            romberg/open-integral})

(def available-methods
  (into #{} (keys quadrature-methods)))

(defn- extract-method
  "Attempts to turn the supplied argument into an integration method; returns nil
  if method doesn't exist."
  [method]
  (cond (fn? method)
        [method {}]

        (keyword? method)
        (extract-method
         (quadrature-methods method))

        (map? method)
        (let [[f m] (extract-method
                     (:method method))]
          [f (merge (dissoc method :method) m)])))

(defn get-integrator
  "Takes:

  - An integration method, specified as either:
    - a keyword naming one of the available methods in `available-methods`
    - a function with the proper integrator signature
    - a dictionary of integrator options with a `:method` key

  - `a` and `b` integration endpoints
  - an optional dictionary of options `m`

  And returns a pair of an integrator function and a possibly-enhanced options
  dictionary.

  (Some integration functions require extra options, so the returned dictionary
  may have more entries than the `m` you pass in.)

  If either endpoint is infinite, the returned integrator is wrapped in
  `qi/improper` and able to handle infinite endpoints (as well as non-infinite
  endpoints by passing through directly to the underlying integrator)."
  ([method a b] (get-integrator method a b {}))
  ([method a b m]
   (when-let [[integrate opts] (extract-method method)]
     (let [integrate (if (or (qc/infinite? a)
                             (qc/infinite? b))
                       (qi/improper integrate)
                       integrate)]
       [integrate (dissoc (merge opts m) :method)]))))

(defn definite-integral
  "Evaluates the definite integral of integrand `f` across the interval $a, b$.
  Optionally accepts a dictionary `opts` of customizing options; All `opts` will
  be passed through to the supplied `integrate` functions.

  If you'd like more control, or to retrieve the integration function directly
  without looking it up via `:method` each time, see `get-integrator`.

  All supplied options are passed through to the underlying integrator; see the
  specific integrator for information on what options are available.

  ## Keyword arguments:

  `:method`: Specifies the integration method used. Must be

  - a keyword naming one of the available methods in `available-methods`
  - a function with the proper integrator signature
  - a dictionary of integrator options with a `:method` key

  Defaults to `:open`, which specifies an adaptive bulirsch-stoer quadrature method.

  `:compile?` If true, the generic function will be simplified and compiled
  before execution. (Clojure only for now.) Defaults to false.

  `:info?` If true, `definite-integral` will return a map of integration
  information returned by the underlying integrator. Else, returns an estimate
  of the definite integral."
  ([f a b] (definite-integral f a b {}))
  ([f a b {:keys [method compile? info?]
           :or {method :open
                compile? false
                info? false}
           :as opts}]
   (if-let [[integrate m] (get-integrator method a b opts)]
     (let [f      #?(:clj (if compile? (c/compile-univariate-function f) f)
                     :cljs f)
           result (integrate f a b m)]
       (if info? result (:result result)))
     (u/illegal (str "Unknown method: " method
                     ". Try one of: "
                     available-methods)))))
