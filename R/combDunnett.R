#' @title  Calculate the inverse normal combination z statistic based on Dunnett's test
#' @description For an intersection hypothesis, calculate the inverse normal combination test based on Dunnett's test up to current stage for cohort 1.
#' @param w A vector of pre-specified weights for incremental test statistics up to the current stage. The squared weights must sum to 1.
#' @param z A matrix of cumulative log-rank test statistics across stages (up to current stage) for individual hypotheses that are in the intersection hypothesis. Columns represent dose levels, rows represent stages.
#' @param d A matrix of cumulative number of events across stages (up to current stage) for dose levels that are included in the intersection hypothesis. Columns represent dose levels, rows represent stages. Event counts must be non-decreasing across stages. An arm with no incremental events at a stage is excluded from the Dunnett test at that stage, with a warning.
#' @param ratio A vector of randomization ratios of each dose arm relative to the control arm (n_dose/n_control). It determines the correlation of the Dunnett test statistics: corr(i,j) = sqrt(lambda_i*lambda_j) with lambda = ratio/(1+ratio). A scalar is recycled across arms. The default ratio = 1 (equal allocation) gives the familiar correlation of 0.5.
#'
#' @return The inverse normal Z statistic for the Dunnett's test
#' @export
#' @import mvtnorm
#' @importFrom stats qnorm
#' @examples
#' # Global intersection H with 3 dose arms. Perform testing at stage 3 post dose selection
#' w <- c(sqrt(2/3),sqrt(1/6),sqrt(1/6))
#' z <- matrix(c(1.1, 0.8,1.2,
#'              1.5, 1.3,1.3,
#'              1.8, 1.6,1.6),nrow = 3, ncol = 3, byrow = TRUE)
#' d <- matrix(c(20, 18,19,
#'              22, 21,23,
#'              25, 23,26),nrow = 3, ncol = 3, byrow = TRUE)
#' combDunnett(w,z,d)
#'
#' # Unequal randomization: dose arms randomized 2:1 against control
#' combDunnett(w,z,d,ratio=2)
#'
combDunnett<- function(w,z,d,ratio=1){

  # ---- input validation ----
  if (!is.numeric(w) || length(w) < 1)
    stop("'w' must be a numeric vector of weights")
  if (abs(sum(w^2) - 1) > 1e-6)
    stop("the squared weights in 'w' must sum to 1")

  s <- length(w) # number of stages
  if (!is.matrix(z)) z <- matrix(z, nrow = 1)
  if (!is.matrix(d)) d <- matrix(d, nrow = 1)
  if (!all(dim(z) == dim(d)))
    stop("'z' and 'd' must have the same dimensions")
  if (nrow(z) != s)
    stop("'z' and 'd' must have one row per stage, matching the length of 'w'")

  k <- ncol(z) # number of dose levels
  if (!is.numeric(d) || any(d <= 0))
    stop("'d' must contain positive event counts")
  if (s > 1 && any(d[-1, , drop = FALSE] < d[-s, , drop = FALSE]))
    stop("cumulative event counts in 'd' must be non-decreasing across stages")
  if (!is.numeric(ratio) || any(ratio <= 0))
    stop("'ratio' must contain positive randomization ratios")
  if (length(ratio) == 1) ratio <- rep(ratio, k)
  if (length(ratio) != k)
    stop("'ratio' must have length 1 or one entry per dose arm")

  # Dunnett correlation: corr(i,j) = sqrt(lambda_i*lambda_j), lambda = ratio/(1+ratio).
  # Under equal allocation (ratio = 1) this reduces to 0.5.
  lambda <- ratio/(1+ratio)

  zdunn <- rep(NA_real_, s)
  for(i in 1:s){
    zi <- z[i,]
    if(i==1){
      incl <- 1:k
      z.incr <- zi
    } else {
      zj <- z[(i-1), ] # i=j+1
      dj <- d[(i-1), ]
      di <- d[i, ]

      incl <- which(di>dj)
      if (length(incl) == 0)
        stop("no dose arm accrued incremental events between stage ", i-1, " and stage ", i)
      if (length(incl) < k)
        warning("dose arm(s) ", paste(setdiff(1:k, incl), collapse = ", "),
                " accrued no incremental events between stage ", i-1, " and stage ", i,
                " and are excluded from the Dunnett test at that stage")

      z.incr <- (sqrt(di[incl])*zi[incl] - sqrt(dj[incl])*zj[incl]) / sqrt(di[incl]-dj[incl])
    }

    ki <- length(z.incr)
    if(ki==1){
      zdunn[i] <- z.incr
    } else {
      zmax <- max(z.incr)
      # define covariance matrix for Dunnett's test
      lam <- lambda[incl]
      cov.z <- sqrt(outer(lam, lam))
      diag(cov.z) <- 1
      # tightened quasi-Monte-Carlo settings so results are reproducible across runs
      alg <- mvtnorm::GenzBretz(abseps = 1e-7, maxpts = 1e6)
      pmax <-   1-mvtnorm::pmvnorm(mean=rep(0,ki ), corr=cov.z, lower=rep(-Inf,ki ), upper=rep(zmax,ki ), algorithm=alg)[1]
      zdunn[i] <- qnorm(pmax, lower.tail = FALSE)
    }


  }

  zDunnett <- sum(w*zdunn)
  return(zDunnett)
}
