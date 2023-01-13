# constructor
pl1 <- new(PersistenceLandscape,
           diagram = matrix(c(0,1, 0,2), nrow = 2L, ncol = 2L, byrow = TRUE),
           exact = TRUE,
           xmin = 0, xmax = 3,
           by = 0.1, ymax = 100)
summary(pl1)
pl2 <- new(PersistenceLandscape,
           diagram = matrix(c(0,1, 0,2), nrow = 2L, ncol = 2L, byrow = TRUE),
           exact = FALSE,
           xmin = -.25, xmax = 5,
           by = 0.1, ymax = 100)
summary(pl2)

# getters
pl1$getInternal()
pl1$toExact()
pl1$toDiscrete()
pl2$getInternal()
try(pl2$toExact())
pl2$toDiscrete()
