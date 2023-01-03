# constructor
pl1 <- new(PersistenceLandscape,
           diagram = matrix(c(0,1, 0,2), nrow = 2L, ncol = 2L, byrow = TRUE),
           exact = TRUE,
           min_x = 0, max_x = 3,
           by = 0.1, max_y = 100)
summary(pl1)
pl2 <- new(PersistenceLandscape,
           diagram = matrix(c(0,1, 0,2), nrow = 2L, ncol = 2L, byrow = TRUE),
           exact = FALSE,
           min_x = -.25, max_x = 5,
           by = 0.1, max_y = 100)
summary(pl2)

# getters
pl1$getInternal()
pl1$getExact()
pl1$getDiscrete()
pl2$getInternal()
try(pl2$getExact())
pl2$getDiscrete()
