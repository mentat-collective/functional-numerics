(ns quadrature.riemann
  (:require [quadrature.interpolate.richardson :as ir]
            [quadrature.common :as qc
             #?@(:cljs [:include-macros true])]
            [sicmutils.generic :as g]
            [sicmutils.util :as u]
            [quadrature.util.aggregate :as ua]
            [quadrature.util.stream :as us]
            [sicmutils.numsymb]))

(defn windowed-sum
  "Takes:

  - `area-fn`, a function of the left and right endpoints of some integration
  slice
  - definite integration bounds `a` and `b`

  and returns a function of `n`, the number of slices to use for an integration
  estimate.

  `area-fn` should return an estimate of the area under some curve between the
  `l` and `r` bounds it receives."
  [area-fn a b]
  (fn [n]
    (let [width       (/ (- b a) n)
          grid-points (concat (range a b width) [b])]
      (ua/sum
       (map area-fn grid-points (rest grid-points))))))

(defn- left-sum* [f a b]
  (-> (fn [l r] (* (f l) (- r l)))
      (windowed-sum a b)))

(defn- left-sum
  "Returns a function of `n`, some number of slices of the total integration
  range, that returns an estimate for the definite integral of $f$ over the
  range $[a, b)$ using a left Riemann sum."
  [f a b]
  (let [width (- b a)]
    (fn [n]
      (let [h  (/ width n)
            fx (fn [i] (f (+ a (* i h))))]
        (* h (ua/sum fx 0 n))))))

(defn- right-sum* [f a b]
  (-> (fn [l r] (* (f r) (- r l)))
      (windowed-sum a b)))

(defn- right-sum
  "Returns a function of `n`, some number of slices of the total integration
  range, that returns an estimate for the definite integral of $f$ over the
  range $(a, b]$ using a right Riemann sum."
  [f a b]
  (let [width (- b a)]
    (fn [n]
      (let [h     (/ width n)
            start (+ a h)
            fx    (fn [i] (f (+ start (* i h))))]
        (* h (ua/sum fx 0 n))))))

(defn- upper-sum
  "Returns an estimate for the definite integral of $f$ over the range $[a, b]$
  using an upper Riemann sum.

  This function may or may not make an evaluation at the endpoints $a$ or $b$,
  depending on whether or not the function is increasing or decreasing at the
  endpoints."
  [f a b]
  (-> (fn [l r] (* (- r l)
                  (max (f l) (f r))))
      (windowed-sum a b)))

(defn- lower-sum
  "Returns an estimate for the definite integral of $f$ over the range $[a, b]$
  using a lower Riemann sum.

  This function may or may not make an evaluation at the endpoints $a$ or $b$,
  depending on whether or not the function is increasing or decreasing at the
  endpoints."
  [f a b]
  (-> (fn [l r] (* (- r l)
                  (min (f l) (f r))))
      (windowed-sum a b)))

(defn- accelerate
  "NOTE - this is only appropriate for Richardson-accelerating sequences with t=2,
  p=q=1.

  This only applies to the Riemann sequences in this namespace!"
  [estimate-seq {:keys [n accelerate?] :or {n 1}}]
  (if (and accelerate? (number? n))
    (ir/richardson-sequence estimate-seq 2 1 1)
    estimate-seq))

(defn midpoint-sum
  "Returns a function of `n`, some number of slices of the total integration
  range, that returns an estimate for the definite integral of $f$ over the
  range $(a, b)$ using midpoint estimates."
  [f a b]
  (let [width (- b a)]
    (fn [n]
      (let [h      (/ width n)
            offset (+ a (/ h 2.0))
            fx     (fn [i] (f (+ offset (* i h))))]
        (* h (ua/sum fx 0 n))))))

(defn Sn->S2n
  "Returns a function of:

  - `Sn`: a sum estimate for `n` partitions, and
  - `n`: the number of partitions

  And returns a new estimate for $S_{2n}$ by sampling the midpoints of each
  slice. This incremental update rule is valid for left and right Riemann sums,
  as well as the midpoint method."
  [f a b]
  (let [midpoints (midpoint-sum f a b)]
    (fn [Sn n]
      (-> (+ Sn (midpoints n))
          (/ 2.0)))))

(defn- left-sequence* [f a b n0]
  (let [first-S ((left-sum f a b) n0)
        steps   (us/powers 2 n0)]
    (reductions (Sn->S2n f a b) first-S steps)))

(defn geometric-estimate-seq
  "Accepts:

  - `S-fn`: a function of `n` that generates a numerical integral estimate from
  `n` slices of some region, and
  - `next-S-fn`: a function of (previous estimate, previous `n`) => new estimate
  - `factor`: the factor by which `n` increases for successive estimates
  - `n0`: the initial `n` to pass to `S-fn`

  The new estimate returned b `next-S-fn` should be of `factor * n` slices."
  [S-fn next-S-fn factor n0]
  (let [first-S (S-fn n0)
        steps   (us/powers factor n0)]
    (reductions next-S-fn first-S steps)))

(defn left-sequence**
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed-open interval $a, b$ by taking left-Riemann sums with

  n0, 2n0, 4n0, ...

  slices."
  ([f a b] (left-sequence** f a b 1))
  ([f a b n0]
   (geometric-estimate-seq (left-sum f a b)
                           (Sn->S2n f a b)
                           2
                           n0)))

(defn- general-estimate-seq
  "Accepts:

  - `S-fn`: a function of `n` that generates a numerical integral estimate from
  `n` slices of some region, and
  - `next-S-fn`: a function of (previous estimate, previous `n`) => new estimate
  - `factor`: the factor by which `next-S-fn` increases `n` in its returned estimate
  - `n-seq`: a monotonically increasing sequence of `n` slices to use.

  Returns a sequence of estimates of returned by either function for each `n` in
  `n-seq`. Internally decides whether or not to use `S-fn` or `next-S-fn` to
  generate successive estimates."
  [S-fn next-S-fn factor n-seq]
  (let [f (fn [[cache _] n]
            (let [Sn (if (zero? (rem n factor))
                       (let [prev (quot n factor)]
                         (if-let [S-prev (get cache prev)]
                           (next-S-fn S-prev prev)
                           (S-fn n)))
                       (S-fn n))]
              [(assoc cache n Sn) Sn]))]
    (->> (reductions f [{} nil] n-seq)
         (map second)
         (rest))))

(defn incrementalize
  "Function that generalizes the ability to create successively-refined estimates
  of an integral, given:

  - `S-fn`: a function of `n` that generates a numerical integral estimate from
  `n` slices of some region, and
  - `next-S-fn`: a function of (previous estimate, previous `n`) => new estimate
  - `factor`: the factor by which `next-S-fn` increases `n` in its returned estimate
  - `n`: EITHER a number, or a monotonically increasing sequence of `n` slices to use.

  If `n` is a sequence, returns a (lazy) sequence of estimates generated for
  each entry in `n`.

  If `n` is a number, returns a lazy sequence of estimates generated for each
  entry in a geometrically increasing series of inputs $n, n(factor),
  n(factor^2), ....$

  Internally decides whether or not to use `S-fn` or `next-S-fn` to generate
  successive estimates."
  [S-fn next-S-fn factor n]
  (let [f (if (number? n)
            geometric-estimate-seq
            general-estimate-seq)]
    (f S-fn next-S-fn factor n)))

(defn left-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed-open interval $a, b$ by taking left-Riemann sums.

  ## Optional Arguments

  `:n`: If `n` is a number, returns estimates with $n, 2n, 4n, ...$ slices,
  geometrically increasing by a factor of 2 with each estimate.

  If `n` is a sequence, the resulting sequence will hold an estimate for each
  integer number of slices in that sequence.

  `:accelerate?`: if supplied (and `n` is a number), attempts to accelerate
  convergence using Richardson extrapolation. If `n` is a sequence this option
  is ignored."
  ([f a b] (left-sequence f a b {}))
  ([f a b opts]
   (let [S      (left-sum f a b)
         next-S (Sn->S2n f a b)]
     (-> (incrementalize S next-S 2 (:n opts 1))
         (accelerate opts)))))

(defn right-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed-open interval $a, b$ by taking right-Riemann sums.

  ## Optional Arguments

  `:n`: If `n` is a number, returns estimates with $n, 2n, 4n, ...$ slices,
  geometrically increasing by a factor of 2 with each estimate.

  If `n` is a sequence, the resulting sequence will hold an estimate for each
  integer number of slices in that sequence.

  `:accelerate?`: if supplied (and `n` is a number), attempts to accelerate
  convergence using Richardson extrapolation. If `n` is a sequence this option
  is ignored."
  ([f a b] (right-sequence f a b {}))
  ([f a b opts]
   (let [S      (right-sum f a b)
         next-S (Sn->S2n f a b)]
     (-> (incrementalize S next-S 2 (:n opts 1))
         (accelerate opts)))))

(defn lower-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed interval $(a, b)$ by taking lower-Riemann sums.

  ## Optional Arguments

  `:n`: If `n` is a number, returns estimates with $n, 2n, 4n, ...$ slices,
  geometrically increasing by a factor of 2 with each estimate.

  If `n` is a sequence, the resulting sequence will hold an estimate for each
  integer number of slices in that sequence.

  `:accelerate?`: if supplied (and `n` is a number), attempts to accelerate
  convergence using Richardson extrapolation. If `n` is a sequence this option
  is ignored."
  ([f a b] (lower-sequence f a b {}))
  ([f a b {:keys [n] :or {n 1} :as opts}]
   (let [n-seq (if (number? n)
                 (us/powers 2 n)
                 n)]
     (-> (map (lower-sum f a b) n-seq)
         (accelerate opts)))))

(defn upper-sequence
  "Returns a (lazy) sequence of successively refined estimates of the integral of
  `f` over the closed interval $(a, b)$ by taking upper-Riemann sums.

  ## Optional Arguments

  `:n`: If `n` is a number, returns estimates with $n, 2n, 4n, ...$ slices,
  geometrically increasing by a factor of 2 with each estimate.

  If `n` is a sequence, the resulting sequence will hold an estimate for each
  integer number of slices in that sequence.

  `:accelerate?`: if supplied (and `n` is a number), attempts to accelerate
  convergence using Richardson extrapolation. If `n` is a sequence this option
  is ignored."
  ([f a b] (upper-sequence f a b {}))
  ([f a b {:keys [n] :or {n 1} :as opts}]
   (let [n-seq (if (number? n)
                 (us/powers 2 n)
                 n)]
     (-> (map (upper-sum f a b) n-seq)
         (accelerate opts)))))

(qc/defintegrator left-integral
  "Returns an estimate of the integral of `f` across the closed-open interval $a,
  b$ using a left-Riemann sum with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `left-sequence` for information on the optional args in `opts` that
  customize this function's behavior."
  :area-fn (fn [f a b] (* (f a) (- b a)))
  :seq-fn left-sequence)

(qc/defintegrator right-integral
  "Returns an estimate of the integral of `f` across the closed-open interval $a,
  b$ using a right-Riemann sum with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `right-sequence` for information on the optional args in `opts` that
  customize this function's behavior."
  :area-fn (fn [f a b] (* (f b) (- b a)))
  :seq-fn right-sequence)

(qc/defintegrator lower-integral
  "Returns an estimate of the integral of `f` across the closed-open interval $a,
  b$ using a lower-Riemann sum with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `lower-sequence` for information on the optional args in `opts` that
  customize this function's behavior."
  :area-fn (fn [f a b] (* (min (f a) (f b)) (- b a)))
  :seq-fn lower-sequence)

(qc/defintegrator upper-integral
  "Returns an estimate of the integral of `f` across the closed-open interval $a,
  b$ using an upper-Riemann sum with $1, 2, 4 ... 2^n$ windows for each estimate.

  Optionally accepts `opts`, a dict of optional arguments. All of these get
  passed on to `us/seq-limit` to configure convergence checking.

  See `upper-sequence` for information on the optional args in `opts` that
  customize this function's behavior."
  :area-fn (fn [f a b] (* (max (f a) (f b)) (- b a)))
  :seq-fn upper-sequence)
