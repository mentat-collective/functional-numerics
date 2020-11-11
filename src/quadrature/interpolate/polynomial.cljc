(ns quadrature.interpolate.polynomial
  "This namespace contains a discussion of polynomial interpolation, and different
  methods for fitting a polynomial of degree N-1 to N points and evaluating that
  polynomial at some different `x`."
  (:require [sicmutils.generic :as g]
            [quadrature.util.aggregate :as ua]
            [quadrature.util.stream :as us]))

(defn lagrange
  "Generates a lagrange interpolating polynomial that fits every point in the
  supplied sequence `points` (of form `[x (f x)]`) and returns the value of the
  polynomial evaluated at `x`.

  The Lagrange polynomial has this form:

  g(x) =  (f(a) * [(x-b)(x-c)...] / [(a-b)(a-c)...])
        + (f(b) * [(x-a)(x-c)...] / [(b-a)(b-c)...])
        + ...

  for points `[a f(a)], [b f(b)], [c f(c)]` etc.

  This particular method of interpolating `x` into the polynomial is
  inefficient; any new calculation requires fully recomputing. Takes O(n^2)
  operations in the number of points.
  "
  [points x]
  (let [points     (vec points)
        n          (count points)
        build-term (fn [i [a fa]]
                     (let [others (for [j (range n) :when (not= i j)]
                                    (get-in points [j 0]))
                           p (reduce g/* (map #(g/- x %) others))
                           q (reduce g/* (map #(g/- a %) others))]
                       (g// (g/* fa p) q)))]
    (transduce (map-indexed build-term)
               g/+
               points)))

(defn neville-recursive
  "Top-down implementation of Neville's algorithm.

  Returns the value of `P(x)`, where `P` is a polynomial fit (using Neville's
  algorithm) to every point in the supplied sequence `points` (of form `[x (f
  x)]`)

  The efficiency and results should be identical to
  `quadrature.interpolate/lagrange`. This function represents a step on
  the journey toward more incremental methods of polynomial interpolation.

  References:

  - Press's Numerical Recipes (p103), chapter 3: http://phys.uri.edu/nigh/NumRec/bookfpdf/f3-1.pdf
  - Wikipedia: https://en.wikipedia.org/wiki/Neville%27s_algorithm"
  [points x]
  (letfn [(evaluate [points]
            (if (= 1 (count points))
              (let [[[_ y]] points]
                y)
              (let [l-branch (pop points)
                    r-branch (subvec points 1)
                    [xl]     (first points)
                    [xr]     (peek points)]
                (g// (g/+ (g/* (g/- x xr) (evaluate l-branch))
                          (g/* (g/- xl x) (evaluate r-branch)))
                     (g/- xl xr)))))]
    (evaluate (vec points))))

(defn- neville-prepare
  "Processes each point of the form [x, (f x)] into:

  $$[x_l, x_r, p]$$

  where $p$ is the polynomial that spans all points from $l$ to $r$. The
  recursion starts with $p = f(x)$.
  "
  [[x fx]]
  [x x fx])

(defn- neville-combine-fn
  "Given some value $x$, returns a function that combines $l$ and $r$ entries in
  the tableau, arranged like this:

  l -- return
     /
    /
   /
  r

  generates the `return` entry of the form

  $$[x_l, x_r, p]$$."
  [x]
  (fn [[xl _ pl] [_ xr pr]]
    (let [plr (g// (g/+ (g/* (g/- x xr) pl)
                        (g/* (g/- xl x) pr))
                   (g/- xl xr))]
      [xl xr plr])))

(defn- neville-next-column
  "This function takes some point $x$, and returns a new function that takes some
  column in the tableau and generates the next column."
  [x]
  (fn [prev-column]
    (map (neville-combine-fn x)
         prev-column
         (rest prev-column))))

(defn- neville-tableau [points x]
  (->> (map neville-prepare points)
       (iterate (neville-next-column x))
       (take-while seq)))

(defn first-terms [tableau]
  (map first tableau))

(defn- neville-present [row]
  (map (fn [[_ _ p]] p) row))

(defn neville-incremental*
  "Takes a potentially lazy sequence of `points` and a point `x` and generates a
  lazy sequence of approximations of P(x).

  entry N in the returned sequence is the estimate using a polynomial generated
  from the first N points of the input sequence."
  [points x]
  (neville-present
   (first-terms
    (neville-tableau points x))))

(defn tableau-fn
  "Returns a Newton-style approximation tableau, given:

  - `prepare`: a fn that processes each element of the supplied `points` into
  the state necessary to calculate future tableau entries.

  - `merge`: a fn of `l`and `r` the tableau entries:

  l -- return
     /
    /
   /
  r

  the inputs are of the same form returned by `prepare`. `merge` should return a
  new structure of the same form.

  - `points`: the (potentially lazy) sequence of points used to generate the
  first column of the tableau.
  "
  [prepare merge points]
  (let [next-col (fn [previous-col]
                   (map merge
                        previous-col
                        (rest previous-col)))]
    (->> (map prepare points)
         (iterate next-col)
         (take-while seq))))

(defn- neville-merge
  "Returns a tableau merge function. Identical to `neville-combine-fn` but uses
  native operations instead of generic operations."
  [x]
  (fn [[xl _ pl] [_ xr pr]]
    (let [p (/ (+ (* (- x xr) pl)
                  (* (- xl x) pr))
               (- xl xr))]
      [xl xr p])))

(defn neville
  "Takes:

  - a (potentially lazy) sequence of `points` of the form `[x (f x)]` and
  - a point `x` to interpolate

  and generates a lazy sequence of approximations of P(x). Each entry in the
  return sequence incorporates one more point from `points` into the P(x)
  estimate.

  Said another way: the Nth in the returned sequence is the estimate using a
  polynomial generated from the first N points of the input sequence:

  p0 p01 p012 p0123 p01234

  This function generates each estimate using Neville's algorithm:

  $$P(x) = [(x - x_r) P_l(x) - (x - x_l) P_r(x)] / [x_l - x_r]$$

  ## Column

  If you supply an integer for the third `column` argument, `neville` will
  return that /column/ of the interpolation tableau instead of the first row.
  This will give you a sequence of nth-order polynomial approximations taken
  between point `i` and the next `n` points.

  As a reminder, this is the shape of the tableau:

   p0 p01 p012 p0123 p01234
   p1 p12 p123 p1234 .
   p2 p23 p234 .     .
   p3 p34 .    .     .
   p4 .   .    .     .

  So supplying a `column` of `1` gives a sequence of linear approximations
  between pairs of points; `2` gives quadratic approximations between successive
  triplets, etc.

  References:

  - Press's Numerical Recipes (p103), chapter 3: http://phys.uri.edu/nigh/NumRec/bookfpdf/f3-1.pdf
  - Wikipedia: https://en.wikipedia.org/wiki/Neville%27s_algorithm
  "
  ([points x]
   (neville-present
    (first-terms
     (tableau-fn neville-prepare
                 (neville-merge x)
                 points))))
  ([points x column]
   (-> (tableau-fn neville-prepare
                   (neville-merge x)
                   points)
       (nth column)
       (neville-present))))

(defn- mn-prepare
  "Processes an initial point [x (f x)] into the required state:

  [x_l, x_r, C, D]

  The recursion starts with $C = D = f(x)$."
  [[x fx]]
  [x x fx fx])

(defn- mn-merge
  "Implements the recursion rules described above to generate x_l, x_r, C and D
  for a tableau node, given the usual left and left-up tableau entries."
  [x]
  (fn [[xl _ _ dl] [_ xr cr _]]
    (let [diff   (- cr dl)
          den    (- xl xr)
          factor (/ diff den)
          c      (* factor (- xl x))
          d      (* factor (- xr x))]
      [xl xr c d])))

(defn mn-present
  "Returns a (lazy) sequence of estimates by successively adding C values from the
  first entry of each tableau column. Each C value is the delta from the
  previous estimate."
  [row]
  (ua/scanning-sum
   (map (fn [[_ _ c _]] c) row)))

(defn modified-neville
  "Similar to `neville` (the interface is identical) but slightly more efficient.
  Internally this builds up its estimates by tracking the delta from the
  previous estimate.

  This non-obvious change lets us swap an addition in for a multiplication,
  making the algorithm slightly more efficient.

  See the `neville` docstring for usage information, and info about the required
  structure of the arguments.

  The structure of the `modified-neville` algorithm makes it difficult to select
  a particular column. See `neville` if you'd like to generate polynomial
  approximations between successive sequences of points.

  References:

  - \"A comparison of algorithms for polynomial interpolation\", A. Macleod,
    https://www.sciencedirect.com/science/article/pii/0771050X82900511
  - Press's Numerical Recipes (p103), chapter 3: http://phys.uri.edu/nigh/NumRec/bookfpdf/f3-1.pdf
  "
  [points x]
  (mn-present
   (first-terms
    (tableau-fn mn-prepare
                (mn-merge x)
                points))))

(defn- generate-new-row* [prepare merge]
  (fn [prev-row point]


    (reduce merge (prepare point) prev-row)))

(defn- generate-new-row [prepare merge]
  (fn [prev-row point]
    (reductions merge (prepare point) prev-row)))

(defn tableau-fold-fn
  "Transforms the supplied `prepare` and `merge` functions into a new function
  that can merge a new point into a tableau row (generating the next tableau
  row).

  More detail on the arguments:

  - `prepare`: a fn that processes each element of the supplied `points` into
  the state necessary to calculate future tableau entries.

  - `merge`: a fn of `l`and `r` the tableau entries:

  l -- return
     /
    /
   /
  r

  the inputs are of the same form returned by `prepare`. `merge` should return a
  new structure of the same form."
  [prepare merge]
  (fn [prev-row point]
    (reductions merge (prepare point) prev-row)))

(defn- neville-fold-fn
  "Returns a function that accepts:

  - `previous-row`: previous row of an interpolation tableau
  - a new point of the form `[x (f x)]`

  and returns the next row of the tableau using the algorithm described in
  `neville`."
  [x]
  (tableau-fold-fn neville-prepare
                   (neville-merge x)))

(defn- modified-neville-fold-fn
  "Returns a function that accepts:

  - `previous-row`: previous row of an interpolation tableau
  - a new point of the form `[x (f x)]`

  and returns the next row of the tableau using the algorithm described in
  `modified-neville`."
  [x]
  (tableau-fold-fn mn-prepare
                   (mn-merge x)))

(defn tableau-fold
  "Returns a function that accepts a sequence of points and processes them into a
  tableau by generating successive rows, one at a time.

  The final row is passed into `present-fn`, which generates the final return
  value.

  This is NOT appropriate for lazy sequences! Fully consumes the input."
  [fold-fn present-fn]
  (fn [points]
    (present-fn
     (reduce fold-fn [] points))))

(defn tableau-scan
  "Takes a folding function and a final presentation function (of accumulator type
  => return value) and returns a NEW function that:

  - accepts a sequence of incoming points
  - returns the result of calling `present` on each successive row."
  [fold-fn present-fn]
  (fn [xs]
    (->> (reductions fold-fn [] xs)
         (map present-fn)
         (rest))))

(defn neville-fold
  "Returns a function that consumes an entire sequence `xs` of points, and returns
  a sequence of successive approximations of `x` using polynomials fitted to the
  points in reverse order.

  This function uses the `neville` algorithm internally."
  [x]
  (tableau-fold (neville-fold-fn x)
                neville-present))

(defn neville-scan
  "Returns a function that consumes an entire sequence `xs` of points, and returns
  a sequence of SEQUENCES of successive polynomial approximations of `x`; one
  for each of the supplied points.

  For a sequence a, b, c... you'll see:

  [(neville [a] x)
   (neville [b a] x)
   (neville [c b a] x)
   ...]"
  [x]
  (tableau-scan (neville-fold-fn x)
                neville-present))

(defn modified-neville-fold
  "Returns a function that consumes an entire sequence `xs` of points, and returns
  a sequence of successive approximations of `x` using polynomials fitted to the
  points in reverse order.

  This function uses the `modified-neville` algorithm internally."
  [x]
  (tableau-fold (modified-neville-fold-fn x)
                mn-present))

(defn modified-neville-scan
  "Returns a function that consumes an entire sequence `xs` of points, and returns
  a sequence of SEQUENCES of successive polynomial approximations of `x`; one
  for each of the supplied points.

  For a sequence a, b, c... you'll see:

  [(modified-neville [a] x)
   (modified-neville [b a] x)
   (modified-neville [c b a] x)
   ...]"
  [x]
  (tableau-scan (modified-neville-fold-fn x)
                mn-present))
