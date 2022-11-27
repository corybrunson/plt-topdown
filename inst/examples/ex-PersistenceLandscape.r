# constructor
pl1 <- new(PersistenceLandscape,
           diagram = matrix(c(0,1, 0,2), nrow = 2L, ncol = 2L, byrow = TRUE),
           exact = TRUE,
           min_b = 0, max_b = 10,
           by = 0.1, max_d = 100)
pl2 <- new(PersistenceLandscape,
           diagram = matrix(c(0,1, 0,2), nrow = 2L, ncol = 2L, byrow = TRUE),
           exact = FALSE,
           min_b = 0, max_b = 10,
           by = 0.1, max_d = 100)

# getters
pl1$getInternal()
pl1$getExact()
pl1$getDiscrete()
pl2$getInternal()
try(pl2$getExact())
pl2$getDiscrete()
