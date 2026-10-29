make_ct_inputs <- function() {
  list(
    w = c(sqrt(0.85), sqrt(0.15)),
    h = c(sqrt(60/700), sqrt(1 - 60/700)),
    z1 = matrix(c(0.8, 1,   0.7,
                  0.9, 1.2, 0.81), nrow = 2, byrow = TRUE),
    d1 = matrix(c(40, 36, 38,
                  51, 42, 49), nrow = 2, byrow = TRUE),
    z2 = 2.35,
    sid = 2,
    zbound = 2.22
  )
}


test_that("combTest evaluates all intersection hypotheses containing the selected dose", {
  x <- make_ct_inputs()

  res <- combTest(x$w, x$h, x$z1, x$d1, x$z2, x$sid, x$zbound)

  # 2^(k-1) intersection hypotheses contain the selected dose
  expect_equal(nrow(res), 4)
  expect_true(all(res[, "dose 2"] == 1))

  expect_equal(res$ComboZ,
               c(2.388778, 2.466033, 2.466033, 2.599204),
               tolerance = 1e-5)

  # all hypotheses rejected at zbound = 2.22, so H02 is rejected overall
  expect_true(all(res$Reject == "Yes"))
})


test_that("combTest passes the randomization ratio through to combDunnett", {
  x <- make_ct_inputs()

  res <- combTest(x$w, x$h, x$z1, x$d1, x$z2, x$sid, x$zbound, ratio = 2)

  expect_equal(res$ComboZ,
               c(2.426710, 2.488215, 2.488215, 2.599204),
               tolerance = 1e-5)
})


test_that("combTest works at the first testing stage with scalar weight", {
  x <- make_ct_inputs()

  res <- combTest(w = 1, h = x$h, z1 = x$z1[1, , drop = FALSE],
                  d1 = x$d1[1, , drop = FALSE], z2 = 2.01, sid = 2, zbound = 2.9)

  expect_equal(nrow(res), 4)
  expect_true(all(res$Reject == "No"))
})


test_that("combTest validates its inputs", {
  x <- make_ct_inputs()

  expect_error(combTest(c(0.5, 0.5), x$h, x$z1, x$d1, x$z2, x$sid, x$zbound),
               regexp = "'w' must sum to 1")
  expect_error(combTest(x$w, c(0.9, 0.9), x$z1, x$d1, x$z2, x$sid, x$zbound),
               regexp = "'h' must sum to 1")
  expect_error(combTest(x$w, x$h, x$z1, x$d1[, 1:2], x$z2, x$sid, x$zbound),
               regexp = "same dimensions")
  expect_error(combTest(x$w, x$h, x$z1[1, , drop = FALSE], x$d1[1, , drop = FALSE],
                        x$z2, x$sid, x$zbound),
               regexp = "one row per stage")
  expect_error(combTest(x$w, x$h, x$z1, x$d1, x$z2, sid = 5, x$zbound),
               regexp = "'sid' must be")
  expect_error(combTest(x$w, x$h, x$z1, x$d1, x$z2, x$sid, zbound = Inf),
               regexp = "'zbound' must be")
  expect_error(combTest(x$w, x$h, x$z1, x$d1, x$z2, x$sid, x$zbound, ratio = c(1, 2)),
               regexp = "length 1 or one entry per dose arm")
})
