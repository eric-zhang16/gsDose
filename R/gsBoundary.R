#' @title Calculate efficacy boundary for inverse normal combination test
#' @description For an intersection hypothesis, calculate the efficacy boundary at current stage for the inverse normal combination test. Currently, the function use an alpha spending function to allocate alpha across stages
#' @param w A matrix of pre-specified weights for the incremental test statistics up to current stage based on cohort 1. Each row contains the weights used at the corresponding stage
#' @param h A matrix of pre-specified weights for cohort 1 and cohort 2 up to current stage.
#' @param d A vector of observed number of events from both cohorts across stages (up to current stage)
#' @param d2 A vector of observed number of events from cohort 2 across stages (up to current stage)
#' @param s The planned total number of stages for efficacy testing (excluding dose selection stage)
#' @param planD The planned total number of events from cohort 2 at FA.
#' @param alpha The assigned alpha level
#' @param sf A spending function. 'OF' for O'Brien-Fleming spending function. 'HSD' for Hwang-Shih-DeCani spending function
#' @param sfpar It specifies the parameter for Hwang-Shih-DeCani spending function . It will be ignored if sfu='OBF'
#' @param bt Bound type. bt='upper' for efficacy upper bound. bt='lower' for efficacy lower bound.
#' @import mvtnorm
#' @importFrom gsDesign sfLDOF
#' @return The efficacy boundaries up to current stage
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
#' # Add futility bounds
#' gsBoundary(w,h,d,d2,s,planD,alpha=1-alpha,sf='HSD',sfpar=-6,bt='lower')
#'
gsBoundary <- function(w,h,d,d2,s,planD,alpha,sf,sfpar=NULL,bt){
  sc <- length(d)
  b.lst <- rep(NA,sc)

  if(s==1){
    if(bt=='upper'){
      b.lst <- qnorm(alpha, lower.tail = FALSE)
    } else if(bt=='lower'){
      b.lst <- qnorm(alpha, lower.tail = TRUE)
    }

  } else {
    for(ss in 1:sc){

      if(ss==1){
        if(sf=='OF'){
          alpha.spent <- sfLDOF(alpha, t=c( d,planD)/planD )$spend[1]
        } else if(sf=='HSD'){
          alpha.spent <- sfHSD(alpha, t=c( d,planD)/planD, sfpar)$spend[1]
        }

        if(bt=='upper'){
          b.lst[ss] <- qnorm(alpha.spent,lower.tail = FALSE)
        } else if(bt=='lower'){
          b.lst[ss] <- qnorm(alpha.spent,lower.tail = TRUE)
        }

      } else {
        if(ss==s){
          if(sf=='OF'){
            alpha.cum <- sfLDOF(alpha, t=c( d[1:(ss-1)],planD)/planD )$spend
          } else if(sf=='HSD'){
            alpha.cum <- sfHSD(alpha, t=c( d[1:(ss-1)],planD)/planD, sfpar)$spend
          }

          alpha.spent <- alpha.cum[ss]-alpha.cum[(ss-1)]

        } else {
          if(sf=='OF'){
            alpha.cum <- sfLDOF(alpha, t=c( d,planD)/planD )$spend
          } else if(sf=='HSD'){
            alpha.cum <- sfHSD(alpha, t=c( d,planD)/planD, sfpar )$spend
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

        if(bt=='upper'){
          # function to calculate Pr(Zi>ci, Z1<c1,Z2<c2...Zi-1<ci-1)
          inc.prob <- function(x){
            alpha.spent-mvtnorm::pmvnorm(mean=rep(0,ss ), corr=cov.m, lower=c(rep(-Inf,(ss-1)),x ), upper=c(b.lst[!is.na(b.lst)],Inf)   )[1]
          }
        } else if(bt=='lower'){
          # function to calculate Pr(Zi<=ci, Z1>c1,Z2>c2...Zi-1>ci-1)
          inc.prob <- function(x){
            alpha.spent-mvtnorm::pmvnorm(mean=rep(0,ss ), corr=cov.m, lower=c(b.lst[!is.na(b.lst)],-Inf ), upper=c(rep(Inf,(ss-1)),x)   )[1]
          }
        }

        b.lst[ss] <- uniroot(inc.prob, interval = c(-5, 5))$root
      }

    }
  }


  res <- data.frame(stage=c(1:sc), Zbound=b.lst)
  return(res)
}
