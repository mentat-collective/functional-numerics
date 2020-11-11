(ns quadrature.common
  "Implements utilities shared by all integrators, for example:

  - code to wrap a sequence of progressively better estimates in a common `integrator` interface
  - data structures implementing various integration intervals."
  #?(:cljs (:refer-clojure :exclude [infinite?]
                           :rename {infinite? core-infinite?}))
  (:require [quadrature.util.stream :as us]
            [taoensso.timbre :as log]))

(def ^:dynamic *roundoff-cutoff* 1e-14)

(def open        [::open ::open])
(def closed      [::closed ::closed])
(def open-closed [::open ::closed])
(def closed-open [::closed ::open])
(def infinities #{##Inf ##-Inf})
(def infinite?
  #?(:cljs core-infinite?
     :clj (comp boolean infinities)))

(defn closed?
  "Returns true if the argument represents an explicit `closed` interval, false
  otherwise."
  [x] (= x closed))
(def open? (complement closed?))

(defn close-l [[_ r]] [::closed r])
(defn close-r [[l _]] [l ::closed])
(defn open-l [[_ r]] [::open r])
(defn open-r [[l _]] [l ::open])
(defn flip [[l r]] [r l])

(defn interval
  "Extracts the interval (or `open` as a default) from the supplied integration
  options dict."
  [opts]
  (get opts :interval open))

(defn with-interval
  "Sets the specified interval to a key inside the suppled `opts` map of arbitrary
  integration options."
  [opts interval]
  (assoc opts :interval interval))

(defn update-interval
  "Accepts:

  - a dictionary of arbitrary options
  - one of the 4 interval modification functions

  and returns a dict of options with `f` applied to the contained interval (or
  `open` if no interval is set).
  "
  [opts f]
  (let [k :interval]
    (assoc opts k (f (interval opts)))))

(defn- narrow-slice?
  "Returns true if the range $[a, b]$ is strip narrow enough to pass the following
  test:

  |b - a| / |a| + |b| <= `cutoff`

  False otherwise. This inequality measures how close the two floating point
  values are, scaled by the sum of their magnitudes."
  [a b cutoff]
  (let [sum (+ (Math/abs a)
               (Math/abs b))]
    (or (<= sum cutoff)
        (<= (Math/abs (- b a))
            (* cutoff sum)))))

(defn make-integrator-fn
  "Generates an `integrator` function from two functions with the following
  signatures and descriptions:

  - `(area-fn f a b)` estimates the integral of `f` over the interval `(a, b)`
  with no subdivision, nothing clever at all.

  - `(seq-fn f a b opts)` returns a sequence of successively refined estimates
  of the integral of `f` over `(a, b)`. `opts` can contain kv pairs that
  configure the behavior of the sequence function (a sequence of the number of
  integration slices to use, for example.)

  The returned function has the signature:

  `(f a b opts)`

  All `opts` are passed on to `seq-fn`, /and/ to `us/seq-limit` internally,
  where the options configure the checks on sequence convergence."
  [area-fn seq-fn]
  (fn call
    ([f a b] (call f a b {}))
    ([f a b {:keys [roundoff-cutoff]
             :or {roundoff-cutoff *roundoff-cutoff*}
             :as opts}]
     (if (narrow-slice? a b roundoff-cutoff)
       (do (log/info "Integrating narrow slice: " a b)
           {:converged? true
            :terms-checked 1
            :result (area-fn f a b)})
       (-> (seq-fn f a b opts)
           (us/seq-limit opts))))))

(defn- name-with-attributes
  "Taken from `clojure.tools.macro/name-with-attributes`.

  Handles optional docstrings and attribute maps for a name to be defined in a
  list of macro arguments. If the first macro argument is a string, it is added
  as a docstring to name and removed from the macro argument list. If afterwards
  the first macro argument is a map, its entries are added to the name's
  metadata map and the map is removed from the macro argument list. The return
  value is a vector containing the name with its extended metadata map and the
  list of unprocessed macro arguments."
  ([name body] (name-with-attributes name body {}))
  ([name body meta]
   (let [[docstring body] (if (string? (first body))
                            [(first body) (next body)]
                            [nil body])
         [attr body]      (if (map? (first body))
                            [(first body) (next body)]
                            [{} body])
         attr             (merge meta attr)
         attr             (if docstring
                            (assoc attr :doc docstring)
                            attr)
         attr             (if (meta name)
                            (conj (meta name) attr)
                            attr)]
     [(with-meta name attr) body])))

(defmacro defintegrator
  "Helper macro for defining integrators."
  [sym & body]
  (let [meta       {:arglists (list 'quote '([f a b] [f a b opts]))}
        [sym body] (name-with-attributes sym body meta)
        {:keys [area-fn seq-fn]} (apply hash-map body)]
    (assert seq-fn (str "defintegrator " sym ": seq-fn cannot be nil"))
    (assert area-fn (str "defintegrator " sym ": area-fn cannot be nil"))
    `(def ~sym
       (make-integrator-fn ~area-fn ~seq-fn))))
