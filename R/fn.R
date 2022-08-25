fn_table <- function(fn) {
    x <- seq(from = -5, to = 5, length.out = 10000)
    cbind(x, fn(x), (fn(x) - fn(x - 0.001)) * 1000)
}