#' @title Calculate efficacy boundary for inverse normal combination test
#' @description For an intersection hypothesis, calculate the efficacy boundary at current stage for the inverse normal combination test.
#' @param w A matrix of pre-specified weights for the incremental test statistics up to current stage based on cohort 1. Columns represent dose levels and rows represent stages
#' @param h A matrix of pre-specified weights for cohort 1 and cohort 2 up to current stage.
#' @param d1 A vector of observed number of events from cohort 1 across stages (up to current stage)
#' @param d2 A vector of observed number of events from cohort 2 across stages (up to current stage)
#' @param s The planned total number of stages for efficacy testing (excluding dose selection stage)
#' @param planD The planned total number of events from cohort 2 at FA.
#' @param alpha: The assigned alpha level
#' @import mvtnorm
#' @importFrom gsDesign sfLDOF
#' @return The efficacy boundaries up to current stage
#' @export
#'
#' @examples
#' # 3 doses at dose selection followed with 4 stages for the selected dose
#' # calculate boundary at 1st stage
#'
#' w <- 1
#' h <- c(sqrt(60/700),sqrt(1-60/700))
#' d1 <- c(20)
#' d2 <- c(230)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d1,d2,s,planD,alpha)
#'
#' # calculate boundary at 2nd stage
#' w <- matrix(c(1,        0,
#'              sqrt(2/3), sqrt(1/3)),nrow = 2, ncol = 2, byrow = TRUE)
#' h <- matrix(c(sqrt(60/700), sqrt(1-60/700),
#'              sqrt(60/700),  sqrt(1-60/700)),nrow = 2, ncol = 2, byrow = TRUE)
#'
#' d1 <- c(20,25)
#' d2 <- c(230,325)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d1,d2,s,planD,alpha)
#'
#' # calculate boundary at 3rd stage
#' w <- matrix(c(1,        0,         0,
#'               sqrt(2/3),sqrt(1/3), 0,
#'               sqrt(2/4),sqrt(1/4), sqrt(1/4)),nrow = 3, ncol = 3, byrow = TRUE)
#' h <- matrix(c(sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700),
#'               sqrt(60/700), sqrt(1-60/700)),nrow = 3, ncol = 2, byrow = TRUE)
#'
#' d1 <- c(20,25,29)
#' d2 <- c(230,325,430)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d1,d2,s,planD,alpha)
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
#' d1 <- c(20,25,29,30)
#' d2 <- c(230,325,430,500)
#' s <- 4
#' planD <- 520
#' alpha <- 0.025
#' gsBoundary(w,h,d1,d2,s,planD,alpha)


gsBoundary <- function(w,h,d1,d2,s,planD,alpha){
  sc <- length(d2)
  d <- d1+d2
  b.lst <- rep(NA,sc)
  for(ss in 1:sc){

    if(ss==1){
      alpha.spent <- sfLDOF(alpha, t=c( d,planD)/planD )$spend[1]
      b.lst[ss] <- qnorm(alpha.spent,lower.tail = F)
    } else {
      if(ss==s){
        alpha.cum <- sfLDOF(alpha, t=c( d[1:(ss-1)],planD)/planD )$spend
        alpha.spent <- alpha.cum[ss]-alpha.cum[(ss-1)]
      } else {
        alpha.cum <- sfLDOF(alpha, t=c( d,planD)/planD )$spend
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

      # function to calculate Pr(Zi>ci, Z1<c1,Z2<c2...Zi-1<ci-1)
      inc.prob <- function(x){
        alpha.spent-mvtnorm::pmvnorm(mean=rep(0,ss ), corr=cov.m, lower=c(rep(-Inf,(ss-1)),x ), upper=c(b.lst[!is.na(b.lst)],Inf)   )[1]
      }

      b.lst[ss] <- uniroot(inc.prob, interval = c(0, 5))$root
    }

  }

  res <- data.frame(stage=c(1:sc), Zbound=b.lst)
  return(res)
}
