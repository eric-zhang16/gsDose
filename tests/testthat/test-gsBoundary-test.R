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
    2.518902,
    2.170359,
    2.081565
  )

  expect_equal(result$Z$stage, 1:4)

  expect_equal(
    result$Z$Zbound,
    expected_boundaries,
    tolerance = 1e-5
  )
})


test_that("gsBoundary returns the expected HSD upper boundaries", {
  x <- make_gs_inputs()

  result <- gsBoundary(
    w = x$w, h = x$h, d = x$d, d2 = x$d2, s = x$s,
    planD = x$planD, alpha = x$alpha,
    sf = "HSD", sfpar = -4, bt = "upper"
  )

  expect_equal(
    result$Z$Zbound,
    c(2.779170, 2.571472, 2.230244, 2.060614),
    tolerance = 1e-5
  )
})


test_that("gsBoundary returns the expected lower (futility) boundaries", {
  x <- make_gs_inputs()

  result <- gsBoundary(
    w = x$w, h = x$h, d = x$d, d2 = x$d2, s = x$s,
    planD = x$planD, alpha = 1 - x$alpha,
    sf = "HSD", sfpar = -6, bt = "lower"
  )

  expect_equal(
    result$Z$Zbound,
    c(-1.739914, -1.129203, -0.051731, 1.959961),
    tolerance = 1e-4
  )
})


test_that("gsBoundary reproduces the illustrated example of Zhang and Jin (2025), Section 3", {
  # Stage 2 (IA) at 388 total events, stage 3 (FA) at 473; cohort-2 events 298 and 371.
  w <- matrix(c(1, 0,
                sqrt(95/107), sqrt(12/107)), nrow = 2, byrow = TRUE)
  h <- matrix(c(sqrt(95/388),  sqrt(293/388),
                sqrt(107/473), sqrt(366/473)), nrow = 2, byrow = TRUE)

  result <- gsBoundary(w, h, d = c(388, 473), d2 = c(298, 371), s = 2,
                       planD = 473, alpha = 0.025, sf = "OF", bt = "upper")

  p_bounds <- pnorm(result$Z$Zbound, lower.tail = FALSE)

  # Paper reports alpha2 = 0.0133 spent at the IA and a p-value boundary of 0.021 at the FA
  expect_equal(p_bounds, c(0.013332, 0.021192), tolerance = 1e-4)
})


test_that("gsBoundary solves boundaries beyond +/-5 at early information fractions", {
  # Regression test: uniroot previously failed when the boundary exceeded 5
  w <- matrix(c(1, 0, sqrt(1/2), sqrt(1/2)), nrow = 2, byrow = TRUE)
  h <- matrix(rep(c(sqrt(60/700), sqrt(640/700)), 2), nrow = 2, byrow = TRUE)

  result <- gsBoundary(w, h, d = c(50, 80), d2 = c(45, 70), s = 4,
                       planD = 520, alpha = 0.025, sf = "OF", bt = "upper")

  expect_true(all(is.finite(result$Z$Zbound)))
  expect_equal(result$Z$Zbound, c(7.133522, 5.595428), tolerance = 1e-5)
})


test_that("gsBoundary warns when interim events reach the planned total", {
  x <- make_gs_inputs()

  expect_warning(
    gsBoundary(
      w = x$w[1:3, 1:3], h = x$h[1:3, ], d = c(250, 350, 530), d2 = c(230, 325, 500),
      s = 4, planD = x$planD, alpha = x$alpha, sf = "OF", bt = "upper"
    ),
    regexp = "reach or exceed 'planD'"
  )
})


test_that("gsBoundary validates its inputs", {
  x <- make_gs_inputs()

  # weights whose squares do not sum to 1
  bad_w <- x$w
  bad_w[2, ] <- c(0.5, 0.5, 0, 0)
  expect_error(
    gsBoundary(bad_w, x$h, x$d, x$d2, x$s, x$planD, x$alpha, sf = "OF", bt = "upper"),
    regexp = "'w' must sum to 1"
  )

  bad_h <- x$h
  bad_h[1, ] <- c(0.9, 0.9)
  expect_error(
    gsBoundary(x$w, bad_h, x$d, x$d2, x$s, x$planD, x$alpha, sf = "OF", bt = "upper"),
    regexp = "'h' must sum to 1"
  )

  # non-increasing total events
  expect_error(
    gsBoundary(x$w, x$h, c(250, 250, 459, 530), x$d2, x$s, x$planD, x$alpha,
               sf = "OF", bt = "upper"),
    regexp = "strictly increasing"
  )

  # cohort-2 events exceeding total events
  expect_error(
    gsBoundary(x$w, x$h, x$d, c(260, 325, 430, 500), x$s, x$planD, x$alpha,
               sf = "OF", bt = "upper"),
    regexp = "cannot exceed 'd'"
  )

  # more observed stages than planned
  expect_error(
    gsBoundary(x$w, x$h, x$d, x$d2, s = 3, x$planD, x$alpha, sf = "OF", bt = "upper"),
    regexp = "cannot exceed the planned number of stages"
  )

  # missing HSD parameter
  expect_error(
    gsBoundary(x$w, x$h, x$d, x$d2, x$s, x$planD, x$alpha, sf = "HSD", bt = "upper"),
    regexp = "'sfpar' must be provided"
  )

  # invalid spending function / bound type
  expect_error(
    gsBoundary(x$w, x$h, x$d, x$d2, x$s, x$planD, x$alpha, sf = "POCOCK", bt = "upper"),
    regexp = "'sf' must be"
  )
  expect_error(
    gsBoundary(x$w, x$h, x$d, x$d2, x$s, x$planD, x$alpha, sf = "OF", bt = "two-sided"),
    regexp = "'bt' must be"
  )

  # invalid alpha
  expect_error(
    gsBoundary(x$w, x$h, x$d, x$d2, x$s, x$planD, alpha = 1.2, sf = "OF", bt = "upper"),
    regexp = "'alpha' must be"
  )
})
