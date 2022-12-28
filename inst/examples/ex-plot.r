# sample points
set.seed(650637L)
rpp <- tdaunif::sample_projective_plane(48L)

# compute persistence data, retaining parameters
pd <- as_persistence(ripserr::vietoris_rips(rpp, max_dim = 2L, threshold = 4))

# plot landscapes
par(mfrow = c(3L, 1L), mar = c(0, 2, 0, 2))
# palette name
plot(landscape(pd, degree = 0L, min_x = 0, max_x = 1.5, by = 0.01),
     palette = "Cork", lwd = 1)
# palette function name
plot(landscape(pd, degree = 1L, min_x = 0, max_x = 1.5, by = 0.01),
     palette = "topo.colors", lwd = 1, rev = TRUE)
# custom color ramp
plot(landscape(pd, degree = 2L, min_x = 0, max_x = 1.5, by = 0.01),
     palette = c("red", "green", "blue"), lwd = 1)
par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))
