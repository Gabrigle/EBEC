# EBEC 0.1.1

* `ebec()` is now the tail-loss test: total loss, worst-case sums at several tail sizes and false-positive components, with a Monte Carlo p-value, ratio observed/expected, percentage excess with bootstrap interval, and a built-in guide.
* The previous share-of-loss index is kept as `ebec_share()`, `ebec_share_fn()` and `ebec_share_fp()` (descriptive only).
