test_that("PL is correct for one persistence pair.", {
  p1 <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pd <- p1
  pl <- landscape(pd, exact = TRUE)
  expected <- list(matrix(c(-Inf, 0, 1, 2, Inf, 0, 0, 1, 0, 0),
                          nrow = 5, ncol = 2))
  
  expect_equal(pl$getInternal(), expected, check.attributes=FALSE)
})

scale <- function(h, pl_exact){
  out <- list()
  
  for(i in 1:length(pl_exact)){
    level <- pl_exact[[i]]
    x <- level[,1]
    y <- h*level[,2]
    out <- list(out,c(x,y))
  }
  
  return(out)
}

test_that("PL sum is correct for simple case.", {
  pd1 <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pd2 <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  pl1 <- landscape(pd1, exact = TRUE)
  pl2 <- landscape(pd2, exact = TRUE)
  
  pl <- pl1$add(pl2)
  expected <- list(matrix(c(-Inf, 0, 1, 1.5, 2, Inf, 0, 0, 1, 1, 0, 0),
                          nrow = 6, ncol = 2))
  
  expect_equal(pl$getInternal(), expected, check.attributes=FALSE)
})

test_that("add PL is correct for simple case.", {
  pd1 <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pd2 <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  pl1 <- landscape(pd1, exact = TRUE)
  pl2 <- landscape(pd2, exact = TRUE)
  
  pl <- pl1$add(pl2)
  expected <- list(matrix(c(-Inf, 0, 1, 1.5, 2, Inf, 0, 0, 1, 1, 0, 0),
                          nrow = 6, ncol = 2))
  
  expect_equal(pl$getInternal(), expected, check.attributes=FALSE)
})

x <- tdaunif::sample_circle(100)
pd <- as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2))

test_that("toExact from exact is correct", {
  pl <- landscape(pd$pairs[[1]], exact=TRUE)
  
  expect_error(pl$toExact(), NA)
})

test_that("toDiscrete from exact is correct", {
  pl <- landscape(pd$pairs[[1]], exact=TRUE)
  
  expect_error(pl$toDiscrete(), NA)
})

test_that("getInternal from discrete is correct", {
  pl <- landscape(pd$pairs[[1]], exact=TRUE)
  
  expect_equal(pl$getInternal(), pl$toExact())
})

test_that("getInternal from discrete is correct from diagram", {
  pd <- as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2))
  pl <- landscape(pd, degree = 1L, exact = TRUE,
                  xmax = 2.5, by = 0.1)
  
  pdref <- as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2))
  plref <- landscape(pdref, degree = 1L, exact = TRUE,
                     xmax = 2.5, by = 0.1, ymax = 2)
  
  expect_equal(pl$getInternal(), plref$getInternal())
})
