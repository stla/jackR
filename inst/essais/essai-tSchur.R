R <- function(a, i, j) {
  aa <- a
  aa[i] <- a[i] + 1L
  aa[j] <- a[j] - 1L
  aa
}
