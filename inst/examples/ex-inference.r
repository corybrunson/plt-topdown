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
pl <- landscape(pd, degree = 1, xmin = 0, xmax = 5, by = .2)
plot(pl)
# plf <- pl$indicator(f, 0)
plf <- pl_indicator(pl, f)
plot(plf, xlim = c(0, 5))
par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))

# apply z-test to two sets of landscapes
set.seed(711018L)
circlescapes <- replicate(
  6,
  tdaunif::sample_circle(n = rpois(n = 1, lambda = 24)) |>
    ripserr::vietoris_rips(max_dim = 2L, threshold = 2) |>
    as_persistence() |>
    landscape(degree = 1, exact = TRUE)
)
toruscapes <- replicate(
  6,
  tdaunif::sample_torus_tube(n = rpois(n = 1, lambda = 24)) |>
    ripserr::vietoris_rips(max_dim = 2L, threshold = 2) |>
    as_persistence() |>
    landscape(degree = 1, exact = TRUE)
)
pl_z_test(circlescapes, toruscapes)
pl_perm_test(circlescapes, toruscapes)

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
