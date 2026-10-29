#' @title Calculate efficacy boundary for inverse normal combination test
#' @description For an intersection hypothesis, calculate the efficacy boundary at current stage for the inverse normal combination test. Currently, the function use an alpha spending function to allocate alpha across stages
#' @param w A matrix of pre-specified weights for the incremental test statistics up to current stage based on cohort 1. Each row contains the weights used at the corresponding stage, and the squared weights in each row must sum to 1. A vector is treated as a single-row matrix.
#' @param h A matrix of pre-specified weights for cohort 1 and cohort 2 up to current stage. Each row contains the two weights used at the corresponding stage, and the squared weights in each row must sum to 1. A vector of length 2 is treated as a single-row matrix.
#' @param d A vector of observed number of events from both cohorts across stages (up to current stage). Must be strictly increasing.
#' @param d2 A vector of observed number of events from cohort 2 across stages (up to current stage). Must be non-decreasing and cannot exceed d.
#' @param s The planned total number of stages for efficacy testing (excluding dose selection stage)
#' @param planD The planned total number of events from both cohorts at FA.
#' @param alpha The assigned alpha level
#' @param sf A spending function. 'OF' for O'Brien-Fleming spending function. 'HSD' for Hwang-Shih-DeCani spending function
#' @param sfpar The parameter for the Hwang-Shih-DeCani spending function. It is required when sf='HSD' and ignored when sf='OF'.
#' @param bt Bound type. bt='upper' for efficacy upper bound. bt='lower' for efficacy lower bound.
#' @import mvtnorm
#' @importFrom gsDesign sfLDOF
#' @return A list with efficacy boundaries up to current stage, and the covariance matrix
#' @export
#' @examples
#' # 3 doses at dose selection followed with 4 stages for the selected dose
#' # calculate boundary at 1st stage
#'
#' w <- 1
#' h <- c(sqrt(60/700),sqrt(1-60/700))
#' d <- c(250)
#' d2 <- c(230)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d,d2,s,planD,alpha,sf='OF',bt='upper')
#'
#' # calculate boundary at 2nd stage
#' w <- matrix(c(1,        0,
#'              sqrt(2/3), sqrt(1/3)),nrow = 2, ncol = 2, byrow = TRUE)
#' h <- matrix(c(sqrt(60/700), sqrt(1-60/700),
#'              sqrt(60/700),  sqrt(1-60/700)),nrow = 2, ncol = 2, byrow = TRUE)
#'
#' d <- c(250,350)
#' d2 <- c(230,325)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d,d2,s,planD,alpha,sf='OF',bt='upper')
#'
#' # calculate boundary at 3rd stage
#' w <- matrix(c(1,        0,         0,
#'               sqrt(2/3),sqrt(1/3), 0,
#'               sqrt(2/4),sqrt(1/4), sqrt(1/4)),nrow = 3, ncol = 3, byrow = TRUE)
#' h <- matrix(c(sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700)),nrow = 3, ncol = 2, byrow = TRUE)
#'
#' d <- c(250,350,459)
#' d2 <- c(230,325,430)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d,d2,s,planD,alpha,sf='OF',bt='upper')
#'
#' # calculate boundary at 4th stage
#' w <- matrix(c(1,        0,         0,         0,
#'               sqrt(2/3),sqrt(1/3), 0,         0,
#'               sqrt(2/4),sqrt(1/4), sqrt(1/4), 0,
#'               sqrt(2/5),sqrt(1/5), sqrt(1/5), sqrt(1/5)), nrow = 4, ncol = 4, byrow = TRUE)
#' h <- matrix(c(sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700)),nrow = 4, ncol = 2, byrow = TRUE)
#'
#' d <- c(250,350,459,530)
#' d2 <- c(230,325,430,500)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d,d2,s,planD,alpha,sf='OF',bt='upper')
#'
#' # Add futility (lower) bounds. By convention, the futility boundary is derived
#' # by spending the complement of alpha (i.e., 1 - alpha) across stages with a
#' # conservative Hwang-Shih-DeCani spending function.
#' gsBoundary(w,h,d,d2,s,planD,alpha=1-alpha,sf='HSD',sfpar=-6,bt='lower')
#'
gsBoundary <- function(w,h,d,d2,s,planD,alpha,sf,sfpar=NULL,bt){

  # ---- input validation ----
  if (!is.numeric(d) || length(d) < 1 || any(d <= 0))
    stop("'d' must be a vector of positive event counts")
  if (!is.numeric(d2) || length(d2) != length(d) || any(d2 <= 0))
    stop("'d2' must be a vector of positive event counts with the same length as 'd'")
  sc <- length(d)
  if (sc > 1 && any(diff(d) <= 0))
    stop("'d' must be strictly increasing across stages")
  if (sc > 1 && any(diff(d2) < 0))
    stop("'d2' must be non-decreasing across stages")
  if (any(d2 > d))
    stop("'d2' cannot exceed 'd' at any stage")
  if (!is.numeric(s) || length(s) != 1 || s < 1 || s != round(s))
    stop("'s' must be a positive integer")
  if (sc > s)
    stop("length of 'd' (", sc, ") cannot exceed the planned number of stages 's' (", s, ")")
  if (!is.numeric(planD) || length(planD) != 1 || planD <= 0)
    stop("'planD' must be a positive scalar")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a scalar in (0, 1)")
  if (length(sf) != 1 || !sf %in% c("OF", "HSD"))
    stop("'sf' must be 'OF' or 'HSD'")
  if (sf == "HSD" && is.null(sfpar))
    stop("'sfpar' must be provided when sf='HSD'")
  if (length(bt) != 1 || !bt %in% c("upper", "lower"))
    stop("'bt' must be 'upper' or 'lower'")

  if (!is.matrix(w)) w <- matrix(w, nrow = 1)
  if (!is.matrix(h)) h <- matrix(h, nrow = 1)
  if (nrow(w) < sc)
    stop("'w' must have one row of weights per observed stage")
  if (nrow(h) < sc || ncol(h) != 2)
    stop("'h' must have one row per observed stage and two columns")
  if (any(abs(rowSums(w[1:sc, , drop = FALSE]^2) - 1) > 1e-6))
    stop("the squared weights in each row of 'w' must sum to 1")
  if (any(abs(rowSums(h[1:sc, , drop = FALSE]^2) - 1) > 1e-6))
    stop("the squared weights in each row of 'h' must sum to 1")

  # warn when events overrun the planned total at an interim stage
  n.int <- min(sc, s - 1)
  if (n.int >= 1 && any(d[1:n.int] >= planD))
    warning("observed events reach or exceed 'planD' at an interim stage: the information fraction is capped at 1, all remaining alpha is spent at that analysis, and boundaries for later stages will be infinite")

  b.lst <- rep(NA,sc)

  if(s==1){
    if(bt=='upper'){
      b.lst <- qnorm(alpha, lower.tail = FALSE)
    } else if(bt=='lower'){
      b.lst <- qnorm(alpha, lower.tail = TRUE)
    }

    cov.m <- NULL

  } else {
    for(ss in 1:sc){

      if(ss==1){
        if(sf=='OF'){
          alpha.spent <- sfLDOF(alpha, t=c( d[1])/planD )$spend
        } else if(sf=='HSD'){
          alpha.spent <- sfHSD(alpha, t=c( d[1])/planD, sfpar)$spend
        }

        if(bt=='upper'){
          b.lst[ss] <- qnorm(alpha.spent,lower.tail = FALSE)
        } else if(bt=='lower'){
          b.lst[ss] <- qnorm(alpha.spent,lower.tail = TRUE)
        }

        cov.m <- NULL

      } else {
        if(ss==s){
          if(sf=='OF'){
            alpha.cum <- sfLDOF(alpha, t=c( d[1:(ss-1)],planD)/planD )$spend
          } else if(sf=='HSD'){
            alpha.cum <- sfHSD(alpha, t=c( d[1:(ss-1)],planD)/planD, sfpar)$spend
          }

          alpha.spent <- alpha - alpha.cum[(ss-1)]

        } else {
          if(sf=='OF'){
            alpha.cum <- sfLDOF(alpha, t=c( d[1:ss] )/planD )$spend
          } else if(sf=='HSD'){
            alpha.cum <- sfHSD(alpha, t=c( d[1:ss] )/planD, sfpar )$spend
          }
          alpha.spent <- alpha.cum[ss]-alpha.cum[(ss-1)]
        }

        cov.m <- matrix(0, nrow = ss, ncol = ss)

        for (i in 1:ss) { # loop over rows
          for (j in i:ss) { # loop over columns
            if (i == j) {
              cov.m[i, j] <- 1  # Set diagonal elements to 1
            } else {
              par1.wsum<- 0
              for(v in 1:i){
                par1.wsum <- par1.wsum+w[i,v]*w[j,v]
              }
              par1.hprod <- h[i,1]*h[j,1]
              par2.hprod <- h[i,2]*h[j,2]
              cov.m[i, j] <- par1.hprod*par1.wsum + par2.hprod*sqrt(d2[i]/d2[j])

              cov.m[j, i] <- cov.m[i, j]  # Maintain symmetry
            }
          }
        }

        # tightened quasi-Monte-Carlo settings so boundaries are reproducible
        # across runs (the mvtnorm default abseps of 1e-3 lets results wobble)
        alg <- mvtnorm::GenzBretz(abseps = 1e-7, maxpts = 1e6)

        if (alpha.spent < 1e-12) {
          warning("no alpha remains to be spent at stage ", ss,
                  "; boundary set to ", if (bt == "upper") "Inf" else "-Inf")
          b.lst[ss] <- if (bt == "upper") Inf else -Inf
        } else {
          if(bt=='upper'){
            # function to calculate Pr(Zi>ci, Z1<c1,Z2<c2...Zi-1<ci-1)
            inc.prob <- function(x){
              alpha.spent-mvtnorm::pmvnorm(mean=rep(0,ss ), corr=cov.m, lower=c(rep(-Inf,(ss-1)),x ), upper=c(b.lst[!is.na(b.lst)],Inf), algorithm=alg   )[1]
            }
          } else if(bt=='lower'){
            # function to calculate Pr(Zi<=ci, Z1>c1,Z2>c2...Zi-1>ci-1)
            inc.prob <- function(x){
              alpha.spent-mvtnorm::pmvnorm(mean=rep(0,ss ), corr=cov.m, lower=c(b.lst[!is.na(b.lst)],-Inf ), upper=c(rep(Inf,(ss-1)),x), algorithm=alg   )[1]
            }
          }

          b.lst[ss] <- uniroot(inc.prob, interval = c(-15, 15), extendInt = "yes", tol = 1e-7)$root
        }
      }

    }
  }


  res <- data.frame(stage=c(1:sc), Zbound=b.lst)

  return(list(Z=res,Covariance=cov.m)  )
}
