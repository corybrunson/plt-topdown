# multiply by an indicator function
pd <- as_persistence(data.frame(dim = 1, birth = c(0, 1), death = c(2, 4)))
f <- list(c(.5, 3), c(1.5, 2), c(1.75, 2))
# exact landscape
pl <- landscape(pd, degree = 1, exact = TRUE)
par(mfcol = c(2L, 2L), mar = c(0, 2, 0, 2))
plot(pl, xlim = c(0, 5))
# plf <- pl$indicator(f, 0)
plf <- pl_indicator(pl, f)
plot(plf, xlim = c(0, 5))
# discrete landscape
pl <- landscape(pd, degree = 1, min_x = 0, max_x = 5, by = .2)
plot(pl)
# plf <- pl$indicator(f, 0)
plf <- pl_indicator(pl, f)
plot(plf, xlim = c(0, 5))
par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))



\dontrun{
# benchmark one- and two-step computation of indicator-based linear form
bench::mark(
  pl$indicator_form(f, 0, p = 1),
  pl$indicator(f, 0)$integral(p = 1)
)
bench::mark(
  pl$indicator_form(f, r = 1, p = 1),
  pl$indicator(f, r = 1)$integral(p = 1)
)
}
