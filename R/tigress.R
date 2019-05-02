#' Trustful Inference of Gene REgulation using Stability Selection
#'
#' This function implements the TIGRESS model of Haury et al. for gene
#' regulatory network (GRN) inference from gene expression data. Given a matrix
#' of expression data, TIGRESS infers regulation from transcription factor (TFs)
#' to any target gene by performing a regression from the TF expression data to
#' the target gene expression, and scoring the candidate TFs by a stability
#' selection (SS) score.
#'
#' @param expdata A matrix of expression. Each row is an experiment, each column
#'   a gene. The gene names are the column names.
#' @param tflist The list of TFs names. The TF names should match the names in
#'   the expression matrix (default: all genes are considered TFs).
#' @param targetlist The list of targets' names. The targets' names should match
#'   the names in the expression matrix (default: all genes are considered
#'   targets).
#' @param alpha The \code{alpha} parameter for randomization in stability
#'   selection. Should be between 0 and 1 (default \code{0.2}).
#' @param nstepsLARS Number of LARS steps to perform in stability selection
#'   (default \code{5}).
#' @param nsplit Number of sample splits to perform in stability selection
#'   (default \code{100})
#' @param normalizeexp A boolean indicating whether we should mean center and
#'   scale to unit variance the expression data for each gene (default
#'   \code{TRUE}).
#' @param scoring The method for scoring a feature in stability selection. If
#'   \code{"area"}, the score is the area under the stability curve up to
#'   nstepsLARS steps, as proposed by Haury et al. If \code{"max"}, the score is
#'   the maximum value of the stability curve, as proposed by Meinshausen and
#'   Bühlmann in the original paper (default \code{"area"}).
#' @param allsteps A boolean indicating whether we should output the solutions
#'   for all values of LARS steps up to nstepsLARS, or only for nstepsLARS. It
#'   does not cost more computation to compute all solutions (default
#'   \code{TRUE}).
#' @param verb A boolean indicating verbose mode. If \code{TRUE}, print messages
#'   about what we are doing, otherwise remain silent (default \code{FALSE}).
#' @param usemulticore A boolean indicating whether multicore parallel computing
#'   is used. This requires the package \code{parallel} (default \code{FALSE}).
#'
#' @return A list of matrices (or a single matrix if \code{allsteps=FALSE}) with
#'   the scores of each TF x target candidate interaction. Each row corresponds
#'   to a TF, each column to a target.
#'
#' @export
tigress <-
  function(expdata, tflist=colnames(expdata), targetlist=colnames(expdata), alpha=0.2, nstepsLARS=5, nsplit=100, normalizeexp=TRUE, scoring="area", allsteps=TRUE, verb=FALSE, usemulticore=FALSE)
  {
    # Check if we can run multicore
    if (usemulticore) {
      require(parallel)
    }

    # Gene names
    genenames <- colnames(expdata)
    ngenes <- length(genenames)

    # Normalize expression data for each gene
    if (normalizeexp)
      expdata <- scale(expdata)

    # Make sure there are no more steps than variables
    if (nstepsLARS>length(tflist)-1){
      nstepsLARS<-length(tflist)-1
      if (nstepsLARS==0){cat('Too few transcription factors! \n',stderr())}
      if (verb){
        cat(paste('Variable nstepsLARS was changed to: ',nstepsLARS,'\n')) }}

    # Locate TF in gene list by matching their names
    ntf <- length(tflist)
    tfindices <- match(tflist,genenames)
    if (max(is.na(tfindices))) {
      stop('Error: could not find all TF in the gene list!')
    }

    # Locate targets in gene list by matching their names
    ntargets <- length(targetlist)
    targetindices <- match(targetlist,genenames)
    if (max(is.na(targetindices))) {
      stop('Error: could not find all targets in the gene list!')
    }

    # Prepare scoring matrix
    if (allsteps) {
      scorestokeep <- nstepsLARS
    } else {
      scorestokeep <- 1
    }
    score <- list()

    # A small function to score the regulators of a single gene
    stabselonegene <- function(itarget) {
      if (verb) {
        cat('.')
      }

      # Name of the target gene
      targetname <- targetlist[itarget]
      # Find the TF to be used for prediction (all TF except the target if the target is itself a TF)
      predTF <- !match(tflist,targetname,nomatch=0)
      r <- stabilityselection(as.matrix(expdata[,tfindices[predTF]]), as.matrix(expdata[,targetname]), nsplit=nsplit, nsteps=nstepsLARS, alpha=alpha)
      sc <- array(0,dim=c(ntf,scorestokeep),dimnames = list(tflist,seq(scorestokeep)))
      if (allsteps) {
        sc[predTF,] <- t(r)
      } else {
        sc[predTF,] <- t(r[nstepsLARS,])
      }
      invisible(sc)
    }

    # Treat target genes one by one
    if (usemulticore) {
      if (requireNamespace("foreach") && foreach::getDoParRegistered()) {
          `%dopar%` = foreach::`%dopar%`
          score <- foreach::foreach(i=seq(ntargets), .packages=c("tigress", "lars"),
                    .export=c("targetlist", "tflist", "tfindices", "expdata",
                              "verb", "nsplit", "nstepsLARS", "alpha", "ntf",
                              "scorestokeep", "allsteps")) %dopar% stabselonegene(i)
      } else {
        score <- mclapply(seq(ntargets),stabselonegene,mc.cores=detectCores()-1)
      }
    } else {
      score <- lapply(seq(ntargets),stabselonegene)
    }
    # Rank scores
    edgepred <- list()
    for (i in seq(scorestokeep)) {
      # Combine all scores in a single vectors
      edgepred[[i]] <- matrix(unlist(lapply(score,function(x) x[,i,drop=FALSE])), nrow=ntf)
      rownames(edgepred[[i]]) <- tflist
      colnames(edgepred[[i]]) <- targetlist
    }

    # Return the result
    if (allsteps) {
      return(edgepred)
    } else {
      return(edgepred[[1]])
    }
  }

