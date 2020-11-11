(ns quadrature.util.stream
  "This namespace contains various standard sequences, as well as utilities for
  working with strict and lazy sequences."
  (:require [clojure.pprint :as pp]
            [sicmutils.generic :as g]
            [sicmutils.value :as v]))

(defn pprint
  "Realizes and pretty-prints `n` elements from the supplied sequence `xs`."
  [n xs]
  (doseq [x (take n xs)]
    (pp/pprint x)))

(defn powers
  "Returns an infinite sequence of `x * n^i`, starting with i == 0. `x` defaults
  to 1."
  ([n] (powers n 1))
  ([n x] (iterate #(* n %) x)))

(defn zeno
  "Returns an infinite sequence of x / n^i, starting with i == 0. `x` defaults to
  1."
  ([n] (zeno n 1))
  ([n x] (iterate #(/ % n) x)))

(defn scan
  "Returns a function that accepts a sequence `xs`, and performs a scan by:

  - Aggregating each element of `xs` into
  - the initial value `init`
  - transforming each intermediate result by `present` (defaults to `identity`).

  The returned sequence contains every intermediate result, after passing it
  through `present` (not `init`, though).

  This is what distinguishes a scan from a 'fold'; a fold would return the final
  result. A fold isn't appropriate for aggregating infinite sequences, while a
  scan works great.

  Arities:

  - the three-arity version takes `init` value, `f` fold function and `present`.
  - the two arity version drops `init`, and instead calls `f` with no arguments.
  The return value is the initial aggregator.
  - The 1-arity version only takes the `merge` fn and defaults `present` to
  `identity`.
    "
  [f & {:keys [present init]
        :or {present identity}}]
  (let [init (or init (f))]
    (fn [xs]
      (->> (reductions f init xs)
           (map present)
           (rest)))))

(defn- close-enuf?
  "relative closeness, transitioning to absolute closeness when we get
  significantly smaller than 1."
  [tolerance]
  (fn [h1 h2]
    (<= (g/abs (- h1 h2))
        (* 0.5 tolerance (+ 2 (g/abs h1) (g/abs h2))))))

(defn seq-limit
  "Accepts a sequence, iterates through it and returns a dictionary of this form:

  {:converged? <boolean>
   :terms-checked <int>
   :result <sequence element>}

  `:converged?` is true if the sequence reached convergence by passing the tests
  described below, false otherwise.

  `:terms-checked` will be equal to the number of items examined in the
  sequence.

  `:result` holds the final item examined in the sequence.

  ## Optional keyword args:

  `:convergence-fn` user-supplied function of two successive elements in `xs`
  that stops iteration and signals convergence if it returns true.

  `:fail-fn` user-supplied function of two successive elements in `xs` that
  stops iteration and signals NO convergence (failure!) if it returns true.

  `:minterms` `seq-limit` won't return until at least this many terms from the
  sequence have been processed.

  `:maxterms` `seq-limit` will return (with `:converged? false`) after
  processing this many elements without passing any other checks.

  `:tolerance` A combination of relative and absolute tolerance. defaults to
  `sqrt(machine epsilon)`."
  ([xs] (seq-limit xs {}))
  ([xs {:keys [minterms
               maxterms
               tolerance
               convergence-fn
               fail-fn]
        :or {minterms       2
             tolerance      (Math/sqrt v/machine-epsilon)
             convergence-fn (close-enuf? tolerance)}}]
   (if (empty? xs)
     {:converged? false
      :terms-checked 0
      :result        nil}
     (let [stop? (if maxterms
                   (fn [i] (>= i maxterms))
                   (constantly false))]
       (loop [[x1 & [x2 :as more]] xs
              terms-checked 1]
         (if (empty? more)
           {:converged?    false
            :terms-checked terms-checked
            :result        x1}
           (let [terms-checked (inc terms-checked)
                 converged?    (convergence-fn x1 x2)]
             (if (and (>= terms-checked minterms)
                      (or converged?
                          (stop? terms-checked)))
               {:converged?    converged?
                :terms-checked terms-checked
                :result        x2}
               (recur more terms-checked)))))))))
