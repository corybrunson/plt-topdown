# scale a landscape
x <- tdaunif::sample_torus_tube(4, 2.5)
pd <- as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2))
pl <- landscape(pd$pairs[[1]], exact = FALSE, max_x=2.5, by=0.1)
print(pl$getInternal()[2, , ])
print(pl$scale(0.5)$getInternal()[2, , ])

# create two landscapes from the same sampling distribution
x <- tdaunif::sample_torus_tube(4, 2.5)
y <- tdaunif::sample_torus_tube(4, 2.5)
x_pd <- as_persistence(ripserr::vietoris_rips(x, dim = 1, threshold = 2))
y_pd <- as_persistence(ripserr::vietoris_rips(y, dim = 1, threshold = 2))
x_pl <- landscape(x_pd$pairs[[1]], exact = FALSE, max_x = 2.5, by = 0.1)
y_pl <- landscape(y_pd$pairs[[1]], exact = FALSE, max_x = 2.5, by = 0.1)

# compare landscapes and calculate mean landscape
print(x_pl$getInternal())
print(y_pl$getInternal())
print(pl_mean(list(x_pl, y_pl))$getInternal())

\dontrun{
set.seed(492869L)

# compute landscape for a large sample
pt <- tdaunif::sample_torus_tube(1000, 5)
pd <- as_persistence(ripserr::vietoris_rips(pt, dim = 2, threshold = 2))
pl <- landscape(pd, degree = 1, exact = FALSE, by = 0.1, min_x = 0, max_x = 2)

# compute landscapes for a large sample of small samples
pl_list <- c()
for (i in seq(100)) {
  pti <- tdaunif::sample_torus_tube(100, 5)
  pdi <- as_persistence(ripserr::vietoris_rips(pti, dim = 2, threshold = 2))
  pli <- landscape(pdi, degree = 1, exact = FALSE, by = 0.1, min_x = 0, max_x = 2)
  pl_list <- c(pl_list, pli)
}

# compute the mean landscape
pl_avg <- pl_mean(pl_list)

# compute the distance between the landscapes
pl_diff <- pl$add(pl_avg$scale(-1))
print(pl_inner(pl_diff, pl_diff))
}

# multiply by an indicator function
pd <- as_persistence(data.frame(dim = 1, birth = c(0, 1), death = c(2, 4)))
f <- list(c(.5, 3), c(1.5, 2), c(1.75, 2))
# exact landscape
pl <- landscape(pd, degree = 1, exact = TRUE)
par(mfcol = c(2L, 2L), mar = c(0, 2, 0, 2))
plot(pl)
plf <- pl$indicator_form(f)
plot(plf, xlim = c(0, 4))
# discrete landscape
pl <- landscape(pd, degree = 1, min_x = 0, max_x = 5, by = .2)
plot(pl)
plf <- pl$indicator_form(f)
plot(plf, xlim = c(0, 5))
par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))
