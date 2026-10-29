test_that("combDunnett reproduces the documented example under equal allocation", {
  w <- c(sqrt(2/3), sqrt(1/6), sqrt(1/6))
  z <- matrix(c(1.1, 0.8, 1.2,
                1.5, 1.3, 1.3,
                1.8, 1.6, 1.6), nrow = 3, byrow = TRUE)
  d <- matrix(c(20, 18, 19,
                22, 21, 23,
                25, 23, 26), nrow = 3, byrow = TRUE)

  expect_equal(combDunnett(w, z, d), 1.262167, tolerance = 1e-5)
})


test_that("combDunnett uses the Dunnett correlation implied by the randomization ratio", {
  w <- c(sqrt(2/3), sqrt(1/6), sqrt(1/6))
  z <- matrix(c(1.1, 0.8, 1.2,
                1.5, 1.3, 1.3,
                1.8, 1.6, 1.6), nrow = 3, byrow = TRUE)
  d <- matrix(c(20, 18, 19,
                22, 21, 23,
                25, 23, 26), nrow = 3, byrow = TRUE)

  # ratio = 1 must agree with the historical hard-coded 0.5 correlation
  expect_equal(combDunnett(w, z, d, ratio = 1), combDunnett(w, z, d), tolerance = 1e-5)

  # 2:1 dose-to-control allocation increases the correlation (lambda = 2/3)
  expect_equal(combDunnett(w, z, d, ratio = 2), 1.396802, tolerance = 1e-5)

  # higher correlation implies a less severe multiplicity penalty
  expect_gt(combDunnett(w, z, d, ratio = 2), combDunnett(w, z, d, ratio = 1))
})


test_that("combDunnett drops zero-increment arms from the stage Dunnett test with a warning", {
  w <- c(sqrt(1/2), sqrt(1/2))
  z <- matrix(c(1.1, 0.8, 1.2,
                1.5, 0.8, 1.3), nrow = 2, byrow = TRUE)
  d <- matrix(c(20, 21, 19,
                30, 21, 28), nrow = 2, byrow = TRUE)  # arm 2 accrues no events

  expect_warning(
    res <- combDunnett(w, z, d),
    regexp = "accrued no incremental events"
  )

  # manual recomputation: stage 1 Dunnett over all 3 arms,
  # stage 2 Dunnett over arms 1 and 3 only
  corr2 <- matrix(0.5, 2, 2); diag(corr2) <- 1
  corr3 <- matrix(0.5, 3, 3); diag(corr3) <- 1
  alg <- mvtnorm::GenzBretz(abseps = 1e-7, maxpts = 1e6)

  p1max <- 1 - mvtnorm::pmvnorm(mean = rep(0, 3), corr = corr3,
                                lower = rep(-Inf, 3), upper = rep(max(z[1, ]), 3),
                                algorithm = alg)[1]

  z_incr <- (sqrt(c(30, 28)) * c(1.5, 1.3) - sqrt(c(20, 19)) * c(1.1, 1.2)) /
    sqrt(c(30, 28) - c(20, 19))
  p2max <- 1 - mvtnorm::pmvnorm(mean = rep(0, 2), corr = corr2,
                                lower = rep(-Inf, 2), upper = rep(max(z_incr), 2),
                                algorithm = alg)[1]

  manual <- sum(w * qnorm(c(p1max, p2max), lower.tail = FALSE))
  expect_equal(res, manual, tolerance = 1e-5)
})


test_that("combDunnett validates its inputs", {
  w <- c(sqrt(1/2), sqrt(1/2))
  z <- matrix(c(1.1, 0.8, 1.5, 1.3), nrow = 2, byrow = TRUE)
  d <- matrix(c(20, 18, 30, 25), nrow = 2, byrow = TRUE)

  # weights whose squares do not sum to 1
  expect_error(combDunnett(c(0.5, 0.5), z, d), regexp = "'w' must sum to 1")

  # decreasing cumulative events are a data error
  d_bad <- matrix(c(20, 18, 15, 25), nrow = 2, byrow = TRUE)
  expect_error(combDunnett(w, z, d_bad), regexp = "non-decreasing")

  # dimension mismatches
  expect_error(combDunnett(w, z, d[, 1:1]), regexp = "same dimensions")
  expect_error(combDunnett(c(1), z, d), regexp = "one row per stage")

  # invalid ratio
  expect_error(combDunnett(w, z, d, ratio = -1), regexp = "positive randomization ratios")
  expect_error(combDunnett(w, z, d, ratio = c(1, 1, 1)), regexp = "length 1 or one entry per dose arm")

  # no arm accrues incremental events
  d_flat <- matrix(c(20, 18, 20, 18), nrow = 2, byrow = TRUE)
  expect_error(combDunnett(w, z, d_flat), regexp = "no dose arm accrued")
})
