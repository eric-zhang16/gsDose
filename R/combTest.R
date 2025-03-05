#' @title Calculate inverse normal combination z statistics combining cohort 1 and 2
#' @description For the selected dose, calculate inverse normal combination z statistics combining cohort 1 and 2 at current stage for all possible intersection hypotheses that include the selected dose
#' @param w A vector of pre-specified weights for incremental test statistics across stages based on cohort 1. Same for all intersection hypotheses.
#' @param h A vector of pre-specified weights for cohort 1 and cohort 2, respectively, at current stage. Same for all intersection hypotheses.
#' @param z1 A matrix of cumulative log-rank test statistics across stages (up to current stage) for individual hypotheses that are in the intersection hypothesis. Columns represent dose levels, rows represent stages.
#' @param d1 A matrix of cumulative number of events across stages (up to current stage) for dose levels that are included in the intersection hypothesis. Columns represent dose levels, rows represent stages.
#' @param z2 The log-rank test statistic at current stage for the selected dose arm based on cohort 2.
#' @param sid ID of the selected dose
#' @param zbound Standardized efficacy boundary
#'
#' @return The inverse normal combination test Z statistics for all possible intersection hypotheses.
#' @export
#' @examples
#' # One dose is selected from 3 dose arms, and performing testing at stage 3 post dose selection
#' # calculate boundary at 1st stage
#' w <- c(sqrt(2/3),sqrt(1/6),sqrt(1/6))
#' h <- c(sqrt(60/700),sqrt(1-60/700))
#' z1 <- matrix(c(0.7,  1,   0.7,
#'               0.8,  1.1, 0.9,
#'               0.85, 1.2, 1),nrow = 3, ncol = 3, byrow = TRUE)
#' d1 <- matrix(c(20, 18,19,
#'               22, 21,23,
#'               25, 23,23),nrow = 3, ncol = 3, byrow = TRUE)
#' z2 <- 2.1
#' sid <- 2 # dose 2 is selected
#' zbound <- 2.01

#' combTest(w,h,z1,d1, z2, sid,zbound)
#'
combTest<- function(w,h,z1,d1, z2, sid,zbound){

  s <- length(w) # number of stages
  k <- ncol(z1) # number of dose levels
  hset <- rev(expand.grid(rep(list(1:0), k))) # all possible intersection hypotheses
  hset <- hset[hset[,sid]==1,] # all intersection hypotheses including selected dose

  # Perform Dunnett test for each intersection H that includes selected dose
  l <- nrow(hset)
  dose.id <- c(1:k)
  zDunnett <- rep(NA,l)
  zCombtest <- rep(NA,l)

  for(j in 1:l){
    hj <- hset[j,]
    dose.idj <- dose.id[hj==1 ]
    zj <- as.matrix(z1[,dose.idj])
    dj <- as.matrix(d1[,dose.idj])
    z1.dunn  <- combDunnett(w,zj,dj)
    zDunnett[j] <- z1.dunn
    zCombtest[j] <- sum(h*c(z1.dunn,z2))
  }

  generate_dose_labels <- function(k) {
    labels <- paste0("dose ", 1:k)
    return(labels)
  }

  generate_H_labels <- function(l) {
    labels <- paste0("H", 1:l)
    return(labels)
  }

  colnames(hset) <- generate_dose_labels(k)
  z1Dunnett <- data.frame(DunnettZ1=zDunnett)
  z2logrank <- data.frame(LogrankZ2=z2)
  zCombtest <- data.frame(ComboZ=zCombtest)

  zm <- cbind(hset,data.frame(stage=s),z1Dunnett,z2logrank,zCombtest)
  row.names(zm) <- generate_H_labels(l)
  zm$Reject <- 'No'
  zm$Reject[zm$ComboZ>zbound] <- 'Yes'

  return( zm )
}
