#' @title  Calculate the inverse normal combination z statistic based on Dunnett's test
#' @description For an intersection hypothesis, calculate the inverse normal combination test based on Dunnett's test up to current stage for cohort 1.
#' @param w A vector of pre-specified weights for incremental test statistics up to the current stage.
#' @param z A matrix of cumulative log-rank test statistics across stages (up to current stage) for individual hypotheses that are in the intersection hypothesis. Columns represent dose levels, rows represent stages.
#' @param d A matrix of cumulative number of events across stages (up to current stage) for dose levels that are included in the intersection hypothesis. Columns represent dose levels, rows represent stages.
#'
#' @return The inverse normal Z statistic for the Dunnett's test
#' @export
#' @import mvtnorm
#' @examples
#' # Global intersection H with 3 dose arms. Perform testing at stage 3 post dose selection
#' w <- c(sqrt(2/3),sqrt(1/6),sqrt(1/6))
#' z <- matrix(c(1.1, 0.8,1.2,
#'              1.5, 1.3,1.3,
#'              1.8, 1.6,1.6),nrow = 3, ncol = 3, byrow = TRUE)
#' d <- matrix(c(20, 18,19,
#'              22, 21,23,
#'              25, 23,23),nrow = 3, ncol = 3, byrow = TRUE)
#'combDunnett(w,z,d)
#'
combDunnett<- function(w,z,d){

  s <- length(w) # number of stages
  k <- ncol(z) # number of dose levels

  zdunn <- rep(NA, s)
  for(i in 1:s){
    zi <- z[i,]
    if(i==1){
      z.incr <- zi
    } else {
      zj <- z[(i-1), ] # i=j+1
      dj <- d[(i-1), ]
      di <- d[i, ]

      z.incr <- rep(0,k)
      popid <- which(di>dj)

      z.incr[popid] <- (sqrt(di[popid])*zi[popid] - sqrt(dj[popid])*zj[popid]) / sqrt(di[popid]-dj[popid])
    }


    if(k==1){
      zdunn[i] <- z.incr
    } else {
      zmax <- max(z.incr)
      # define covariance matrix for Dunnett's test
      cov.z <- matrix(0.5, k, k)
      diag(cov.z) <- 1
      pmax <-   1-mvtnorm::pmvnorm(mean=rep(0,k ), corr=cov.z, lower=rep(-Inf,k ), upper=rep(zmax,k ))[1]
      zdunn[i] <- qnorm(pmax, lower.tail = FALSE)
    }


  }

  zDunnett <- sum(w*zdunn)
  return(zDunnett)
}
