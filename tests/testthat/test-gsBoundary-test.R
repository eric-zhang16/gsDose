# Standard inputs used across tests
make_gs_inputs <- function() {
  list(
    w = matrix(
      c(
        1,         0,         0,         0,
        sqrt(2/3), sqrt(1/3), 0,         0,
        sqrt(2/4), sqrt(1/4), sqrt(1/4), 0,
        sqrt(2/5), sqrt(1/5), sqrt(1/5), sqrt(1/5)
      ),
      nrow = 4,
      ncol = 4,
      byrow = TRUE
    ),

    h = matrix(
      rep(
        c(sqrt(60 / 700), sqrt(1 - 60 / 700)),
        4
      ),
      nrow = 4,
      ncol = 2,
      byrow = TRUE
    ),

    d = c(250, 350, 459, 530),
    d2 = c(230, 325, 430, 500),
    s = 4,
    planD = 520,
    alpha = 0.025
  )
}


test_that("gsBoundary returns the expected OF upper boundaries", {
  x <- make_gs_inputs()

  result <- gsBoundary(
    w = x$w,
    h = x$h,
    d = x$d,
    d2 = x$d2,
    s = x$s,
    planD = x$planD,
    alpha = x$alpha,
    sf = "OF",
    bt = "upper"
  )

  expected_boundaries <- c(
    3.029029,
    2.518903,
    2.170391,
    2.081603
  )

  expect_equal(result$Z$stage, 1:4)

  expect_equal(
    result$Z$Zbound,
    expected_boundaries,
    tolerance = 1e-5
  )
})
