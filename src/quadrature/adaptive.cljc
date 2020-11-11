(ns quadrature.adaptive
  (:require [quadrature.common :as qc]
            [quadrature.util.aggregate :as ua]))

(def ^:dynamic *adaptive-maxterms* 10)

(def ^:dynamic *neighborhood-width* 0.05)

(defn- split-point
  "Returns a point within`fuzz-factor` of the midpoint of the interval $[a, b]$.
  `fuzz-factor` defaults to 0 (ie, `split-point` returns the midpoint)."
  ([a b] (split-point a b 0))
  ([a b fuzz-factor]
   {:pre [(>= fuzz-factor 0)
          (< fuzz-factor 1)]}
   (let [width  (- b a)
         offset (if (zero? fuzz-factor)
                  0.5
                  (+ 0.5 (* fuzz-factor (dec (rand 2.0)))))]
     (+ a (* offset width)))))

(defn- fill-defaults
  "Populates the supplied `opts` dictionary with defaults required by `adaptive`.
  Two of these have values controlled by dynamic variables in `adaptive.cljc`."
  [opts]
  (merge {:maxterms *adaptive-maxterms*
          :adaptive-neighborhood-width *neighborhood-width*
          :interval qc/open}
         opts))

(defn adaptive
  "Accepts one or two 'integrators', ie, functions of:

  - `f`: some integrand
  - `a` and `b`: the lower and upper endpoints of integration
  - `opts`, a dictionary of configuration options

  And returns a new integrator that adaptively subdivides the region $a, b$ into
  intervals if integration fails to converge. If two integrators are supplied,
  the first is applied to any interval that's not explicitly closed on both
  ends, and the second integrator is applied to explicitly closed intervals.
  This behavior depends on the interval set in the supplied `opts` dictionary.

  All `opts` will be passed through to the supplied `integrate` functions.

  ## Optional arguments relevant to `adaptive`:

  `:maxterms`: defaults to `*adaptive-maxterms*`. This is passed to the
  underlying integrators, and determines how long each interval attempts to
  converge before a subdivision.

  `:adaptive-neighborhood-width`: When dividing an interval, the split point
  will be within this factor of the midpoint. Set `:adaptive-neighborhood-width`
  to 0 for deterministic splitting."
  ([integrator] (adaptive integrator integrator))
  ([open-integrator closed-integrator]
   (fn rec
     ([f a b] (rec f a b {}))
     ([f a b opts]
      (let [opts      (fill-defaults opts)
            integrate (fn [l r interval]
                        (if (qc/closed? interval)
                          (closed-integrator f l r opts)
                          (open-integrator f l r opts)))]
        (loop [stack [[a b (:interval opts)]]
               sum   (ua/kahan-sum)
               iteration 0]
          (if (empty? stack)
            {:converged? true
             :iterations iteration
             :result (first sum)}
            (let [[l r interval] (peek stack)
                  remaining      (pop stack)
                  {:keys [converged? result]} (integrate l r interval)]
              (if converged?
                (recur remaining
                       (ua/kahan-sum sum result)
                       (inc iteration))
                (let [midpoint (split-point l r (:adaptive-neighborhood-width opts))]
                  (recur (conj remaining
                               [midpoint r (qc/close-l interval)]
                               [l midpoint (qc/close-r interval)])
                         sum
                         (inc iteration))))))))))))
