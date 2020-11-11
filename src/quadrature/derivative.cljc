(ns quadrature.derivative
  "Different numerical derivative implementations."
  (:require [sicmutils.calculus.derivative :as d]
            [quadrature.interpolate.richardson :as r]
            [sicmutils.function :as f]
            [sicmutils.generic :as g]
            [sicmutils.infix :as if]
            [quadrature.util :as u]
            [quadrature.util.stream :as us]
            [sicmutils.value :as v]))

(defn forward-difference
  "Returns a single-variable function of a step size `h` that calculates the
  forward-difference estimate of the the first derivative of `f` at point `x`:

  f'(x) = [f(x + h) - f(x)] / h

  Optionally accepts a third argument `fx == (f x)`, in case you've already
  calculated it elsewhere and would like to save a function evaluation."
  ([f x] (forward-difference f x (f x)))
  ([f x fx]
   (fn [h]
     (/ (- (f (+ x h)) fx) h))))

(def ^:private fx-h :tangle no
  (->> (d/taylor-series-terms func 'x (g/negate 'h))
       (take 5)
       (reduce g/+)))

(defn backward-difference
  "Returns a single-variable function of a step size `h` that calculates the
  backward-difference estimate of the first derivative of `f` at point `x`:

  f'(x) = [f(x) - f(x - h)] / h

  Optionally accepts a third argument `fx == (f x)`, in case you've already
  calculated it elsewhere and would like to save a function evaluation."
  ([f x] (backward-difference f x (f x)))
  ([f x fx]
   (fn [h]
     (/ (- fx (f (- x h))) h))))

(defn central-difference
  "Returns a single-variable function of a step size `h` that calculates the
  central-difference estimate of the first derivative of `f` at point `x`:

  f'(x) = [f(x + h) - f(x - h)] / 2h"
  [f x]
  (fn [h]
    (/ (- (f (+ x h)) (f (- x h)))
       (* 2 h))))

(defn central-difference-d2
  "Returns a single-variable function of a step size `h` that calculates the
  central-difference estimate of the second derivative of `f` at point `x`:

  f''(x) = [f(x + h) - 2f(x) + f(x - h)] / h^2

  Optionally accepts a third argument `fx == (f x)`, in case you've already
  calculated it elsewhere and would like to save a function evaluation."
  ([f x] (central-difference-d2 f x (f x)))
  ([f x fx]
   (let [fx*2 (* 2 fx)]
     (fn [h]
       (/ (- (+ (f (+ x h))
                (f (- x h)))
             fx*2)
          (* h h))))))

(defn- make-derivative-fn
  [f]
  (fn [x]
    (let [h 1e-5]
      ((central-difference f x) h))))

(defn- central-diff-stream [f x h]
  (map (central-difference f x)
       (us/zeno 2 h)))

(defn- roundoff-units
  "Returns the number of 'roundoff units', ie, multiples of the machine epsilon,
  that roundoff error contributes to the total relative error, given a relative
  error percentage estimated for some initial step size $h$."
  [rel-error-ratio]
  (inc
   (Math/floor
    (Math/abs
     (double rel-error-ratio)))))

(defn- max-iterations
  "Solution for `n`, in:

  `initial-error` * 2^n <= `tolerance`"
  [units tolerance]
  (let [initial-error (* v/machine-epsilon units)]
    (Math/floor
     (/ (Math/log (/ tolerance initial-error))
        (Math/log 2)))))

(defn- terms-before-roundoff
  "Generates a default max number of terms, based on roundoff error estimates."
  [ratio tolerance]
  (inc
   (max-iterations (roundoff-units ratio)
                   tolerance)))

(def valid-methods
  #{:central :central-d2 :forward :backward})

(defn- configs [method f x fx]
  (case method
    :forward
    {:p 1
     :q 1
     :function (forward-difference f x fx)
     :ratio-fn (fn [h] (/ fx (- (f (+ x h)) fx)))}

    :central
    {:p 2
     :q 2
     :function (central-difference f x)
     :ratio-fn (fn [h]
                 (/ fx (- (f (+ x h))
                          (f (- x h)))))}

    :backward
    {:p 1
     :q 1
     :function (backward-difference f x fx)
     :ratio-fn (fn [h]
                 (/ fx (- fx (f (- x h)))))}

    :central-d2
    {:p 2
     :q 2
     :function (central-difference-d2 f x fx)
     :ratio-fn (fn [h]
                 (- (+ (f (+ x h))
                       (f (- x h)))
                    (* 2 fx)))}

    (u/illegal
     (str "Invalid method: " method ". Please try one of " valid-methods))))

(defn- fill-defaults
  "Fills in default values required by `D-numeric`. Any option not used by
  `D-numeric` gets passed on to `us/seq-limit`."
  [m]
  (let [defaults {:tolerance (Math/sqrt v/machine-epsilon)
                  :method    :central}
        {:keys [method] :as opts} (merge defaults m)]
    (assert (contains? valid-methods method)
            (str method " is not a valid method. Please try one of: " valid-methods))
    opts))

(defn D-numeric
  "Takes a function `f: R => R` (function of a single real variable), and returns
  a new function of `x` that approximates the derivative $Df(x)$ (or $D^2f(x)$
  if you pass `:method :central-d2`).

  Returns the estimated value of the derivative at `x`. If you pass `:info?
  true`, the fn returns a dictionary of the results of `us/seq-limit`:

  {:converged? <boolean>
   :terms-checked <int>
   :result <derivative estimate>}

  Make sure to visit `sicmutils.calculus.derivative/D` if you want symbolic or
  automatic differentiation.

  ## Roundoff Estimate

  The returned function will attempt to estimate how many times it can halve the
  step size used to estimate the derivative before roundoff error swamps the
  calculation, and force the function to return (with `:converged? false`, if
  you pass `:info?`)

  ## Optional Arguments

  `D-numeric` takes optional args as its second param. Any of these can be
  overridden by passing a second argument to the function returned by
  `D-numeric`; helpful for setting defaults and then overriding them later.

  The returned function passes through these and any other options to
  `us/seq-limit`, where they control the sequence of richardson
  extrapolation-accelerated estimates.

  Options:

  - `:method`: one of `:central`, `:central-d2`, `:forward` or `:backward`.
  `:central-d2` forces a second derivative estimate; the other methods configure
  a first derivative estimator.

  - `:info?` if false (default), returns the estimated value of `x`. If true,
  returns a dictionary with more information (see `D-numeric`'s docstring for
  more info.)

  - `:initial-h`: the initial `h` to use for derivative estimates before $h \to
  0$. Defaults to 0.1 * abs(x).

  - `:tolerance`: see `us/stream-limit` for a discussion of how this value
  handles relative vs absolute tolerance. $\\sqrt(\\epsilon)$ by default, where
  $\\epsilon$ = machine tolerance.

  - `:maxterms`: the maximum number of terms to consider when hunting for a
  derivative estimate. This defaults to an estimate generated internally,
  designed to prevent roundoff error from swamping the result. If you want to
  disable this feature, set `:maxterms` to something moderately large, like
  `:maxterms 100`. But do so carefully! See the surrounding namespace for a
  larger discussion."
  ([f] (D-numeric f {}))
  ([f opts]
   (let [opts (fill-defaults opts)]
     (fn df
       ([x] (df x {}))
       ([x overrides]
        (let [{:keys [maxterms tolerance initial-h method info?] :as opts} (merge opts overrides)
              {:keys [ratio-fn function p q]} (configs method f x (f x))
              h (or initial-h (* 0.1 (g/abs x)))
              n (or maxterms (terms-before-roundoff
                              (ratio-fn h)
                              tolerance))
              estimates (map function (us/zeno 2 h))
              result    (-> (r/richardson-sequence estimates 2 p q)
                            (us/seq-limit (assoc opts :maxterms n)))]
          (if info? result (:result result))))))))
