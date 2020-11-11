(ns quadrature.util.aggregate
  "Utilities for aggregating sequences.")

(defn kahan-sum
  "Implements a fold that tracks the summation of a sequence of floating point
  numbers, using Kahan's trick for maintaining stability in the face of
  accumulating floating point errors."
  ([] [0.0 0.0])
  ([[sum c] x]
   (let [y (- x c)
         t (+ sum y)]
     [t (- (- t sum) y)])))

(defn sum
  "Sums either:

  - a series `xs` of numbers, or
  - the result of mapping function `f` to `(range low high)`

  Using Kahan's summation trick behind the scenes to keep floating point errors
  under control."
  ([xs]
   (first
    (reduce kahan-sum [0.0 0.0] xs)))
  ([f low high]
   (sum (map f (range low high)))))

(defn scanning-sum
  "Returns every intermediate summation from summing either:

  - a series `xs` of numbers, or
  - the result of mapping function `f` to `(range low high)`

  Using Kahan's summation trick behind the scenes to keep floating point errors
  under control."
  ([xs]
   (->> (reductions kahan-sum [0.0 0.0] xs)
        (map first)
        (rest)))
  ([f low high]
   (scanning-sum
    (map f (range low high)))))
