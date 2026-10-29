#' @title Calculate inverse normal combination z statistics combining cohort 1 and 2
#' @description For the selected dose, calculate inverse normal combination z statistics combining cohort 1 and 2 at current stage for all possible intersection hypotheses that include the selected dose
#' @param w A vector of pre-specified weights for incremental test statistics across stages based on cohort 1. Same for all intersection hypotheses. The squared weights must sum to 1.
#' @param h A vector of pre-specified weights for cohort 1 and cohort 2, respectively, at current stage. Same for all intersection hypotheses. The squared weights must sum to 1.
#' @param z1 A matrix of cumulative log-rank test statistics across stages (up to current stage) for individual hypotheses that are in the intersection hypothesis. Columns represent dose levels, rows represent stages.
#' @param d1 A matrix of cumulative number of events across stages (up to current stage) for dose levels that are included in the intersection hypothesis. Columns represent dose levels, rows represent stages. Event counts must be non-decreasing across stages.
#' @param z2 The log-rank test statistic at current stage for the selected dose arm based on cohort 2.
#' @param sid ID of the selected dose
#' @param zbound Standardized efficacy boundary
#' @param ratio A vector of randomization ratios of each dose arm relative to the control arm (n_dose/n_control), passed to \code{\link{combDunnett}} to determine the correlation of the Dunnett test statistics. A scalar is recycled across arms. The default ratio = 1 (equal allocation) gives the familiar correlation of 0.5.
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
#'               25, 23,26),nrow = 3, ncol = 3, byrow = TRUE)
#' z2 <- 2.1
#' sid <- 2 # dose 2 is selected
#' zbound <- 2.01

#' combTest(w,h,z1,d1, z2, sid,zbound)
#'
combTest<- function(w,h,z1,d1, z2, sid,zbound,ratio=1){

  # ---- input validation ----
  if (!is.numeric(w) || length(w) < 1)
    stop("'w' must be a numeric vector of weights")
  if (abs(sum(w^2) - 1) > 1e-6)
    stop("the squared weights in 'w' must sum to 1")
  if (!is.numeric(h) || length(h) != 2)
    stop("'h' must contain two weights (cohort 1 and cohort 2)")
  if (abs(sum(h^2) - 1) > 1e-6)
    stop("the squared weights in 'h' must sum to 1")

  s <- length(w) # number of stages
  if (!is.matrix(z1)) z1 <- matrix(z1, nrow = 1)
  if (!is.matrix(d1)) d1 <- matrix(d1, nrow = 1)
  if (!all(dim(z1) == dim(d1)))
    stop("'z1' and 'd1' must have the same dimensions")
  if (nrow(z1) != s)
    stop("'z1' and 'd1' must have one row per stage, matching the length of 'w'")

  k <- ncol(z1) # number of dose levels
  if (!is.numeric(sid) || length(sid) != 1 || !sid %in% seq_len(k))
    stop("'sid' must be a single dose index between 1 and ", k)
  if (!is.numeric(z2) || length(z2) != 1 || !is.finite(z2))
    stop("'z2' must be a finite numeric scalar")
  if (!is.numeric(zbound) || length(zbound) != 1 || !is.finite(zbound))
    stop("'zbound' must be a finite numeric scalar")
  if (!is.numeric(ratio) || any(ratio <= 0))
    stop("'ratio' must contain positive randomization ratios")
  if (length(ratio) == 1) ratio <- rep(ratio, k)
  if (length(ratio) != k)
    stop("'ratio' must have length 1 or one entry per dose arm")

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
    zj <- matrix(NA, nrow=s, ncol=length(dose.idj))
    zj[,1:length(dose.idj)] <- z1[,dose.idj]
    dj <- matrix(NA, nrow=s, ncol=length(dose.idj))
    dj[,1:length(dose.idj)] <- d1[,dose.idj]

    z1.dunn  <- combDunnett(w,zj,dj,ratio=ratio[dose.idj])
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
