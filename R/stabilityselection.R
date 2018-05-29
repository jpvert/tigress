#' Score variable importance with stability selection
#'
#' Score the importance of covariates to predict an output by combining LARS
#' sparse regression with stability selection. Given a matrix of covariates and
#' a vector or outputs to predict, stability selection works by solving
#' repeatedly a sparse regression problem (here we use the LARS method) on a
#' randomly modified covariate matrix (here we subsample the rows and randomly
#' reweight the columns). The stability selection (SS) score then evaluates the
#' importance of each covariates based on how often it is selected. We implement
#' two SS scores, the original one of Meinshausen and Buhlmann which measure of
#' frequency of selection among the top \code{L}L covariates, and the area score
#' of Haury and Vert which combines the frequency of selection among the top
#' \code{L} covariates for different values of \code{L}.
#'
#' @param x The input matrix, each row is a sample, each column a feature.
#' @param y A vector of response variable.
#' @param nsplit The number of splits of the samples into two subsamples
#'   (default \code{100})
#' @param nstepsLARS The maximum number of LARS steps performed at each
#'   iteration (default \code{20})
#' @param alpha The random multiplicative weights of each column are uniformly
#'   sampled in the interval [\code{alpha},1] (default \code{0.2})
#' @param scoring How to score a feature. If \code{"area"} we compute the area
#'   under the stability curve, as proposed by Haury et al. If \code{"max"} we
#'   just compute the stability curve, as propose by Meinshausen and Buhlmann
#'   (default \code{"area"})
#'
#' @return A matrix of SS scores. Each column corresponds to a covariate. Each
#'   row corresponds to a number of LARS steps.
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 40
#' x <- matrix(rnorm(n*p),n,p)
#' beta <- c(rnorm(5),numeric(p-5))
#' y <- x%*%beta + rnorm(n)
#' s <- stabilityselection(x, y, nsplit=500, nstepsLARS=5)
#' matplot(s, type='b', lwd=2, ylab="SS area score", xlab="LARS steps")
stabilityselection <-
  function(x,y,nsplit=100,nstepsLARS=20,alpha=0.2,scoring="area")
  {
    if (!is.numeric(y) || sd(y)==0) stop("y should be a vector of scalars not constant.")
    n <- nrow(x)
    p <- ncol(x)
    halfsize <- as.integer(n/2)
    freq <- matrix(0,nstepsLARS,p)

    i <- 0
    while (i < 2*nsplit) {
      # Randomly reweight each variable
      xs <- t(t(x)*runif(p,alpha,1))

      # Ramdomly split the sample in two sets
      badsplit <- TRUE
      while (badsplit) {
        perm <- sample(n)
        i1 <- perm[1:halfsize]
        i2 <- perm[(halfsize+1):n]
        if (max(sd(y[i1]),sd(y[i2]))>0) {badsplit=FALSE}
      }

      # run LARS on each randomized, sample and check which variables are selected
      if (sd(y[i1]>0)) {
        r <- lars(xs[i1,],y[i1],max.steps=nstepsLARS,normalize=FALSE,trace=FALSE)
        freq<-freq + abs(sign(r$beta[2:(nstepsLARS+1),]))
        i <- i+1
      }
      if (sd(y[i2]>0)) {
        r <- lars(xs[i2,],y[i2],max.steps=nstepsLARS,normalize=FALSE,trace=FALSE)
        freq<-freq + abs(sign(r$beta[2:(nstepsLARS+1),]))
        i <- i+1
      }
    }

    # normalize frequence in [0,1] to get the stability curves
    freq <- freq/i

    # Compute normalized area under the stability curve
    if (scoring=="area")
      score <- apply(freq,2,cumsum)/seq(nstepsLARS)
    else
      score <- apply(freq, 2, cummax)

    invisible(score)
  }
