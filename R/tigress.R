#' Trustful Inference of Gene REgulation using Stability Selection
#'
#' This function implements the TIGRESS model of Haury et al. for gene
#' regulatory network (GRN) inference from gene expression data. Given a matrix
#' of expression data, TIGRESS infers regulation from transcription factor (TFs)
#' to any target gene by performing a regression from the TF expression data to
#' the target gene expression, and scoring the candidate TFs by a stability
#' selection (SS) score. The TF-target pairs are then sorted by decreasing SS
#' score.
#'
#' @param expdata Either a matrix of expression, or the name of a file
#'   containing it. Each row is an experiment, each column a gene. The gene
#'   names are the column names (or are in the first row of the file).
#' @param tflist The list of TFs, or the name of a file containing them. The TF
#'   name should match the names in the gene list of the expression data file.
#'   If NULL, then all genes are considered TF (default \code{NULL}).
#' @param outfile. A file name where we TF-target predictions are written. If
#'   empty, do not write anything (default \code{''}).
#' @param K Total number of TF-target predicted edges to return, by decreasing
#'   SS score. K=0 means that all edges are returned (default \code{0}).
#' @param alpha The \code{alpha) parameter for randomization in stability
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
#'   \code{FALSE}).
#' @param verb A boolean indicating verbose mode. If \code{TRUE}, print messages
#'   about what we are doing, otherwise remain silent (default \code{FALSE}).
#' @param usemulticore A boolean indicating whether multicore parallel computing
#'   is used. This requires the package \code{parallel} (default \code{FALSE}).
#'
#' @return A dataframe (or list of dataframes if \code{allsteps=TRUE}) with the
#'   top \code{K} predicted edges. First column is the TF, second column the
#'   target gene, third column the score. If outfile is provided, the result is
#'   also written to the file OUTFILE (if allsteps=FALSE) or to several files
#'   OUTFILE1, OUTFILE1, ... (if allsteps=TRUE)
#'
#' @export
tigress <-
  function(expdata , tflist=NULL , outfile="prediction.txt" , K=-1 , alpha=0.2 , nstepsLARS=5 , split=100 , normalizeexp=TRUE , scoring="area" , allsteps=FALSE , verb=FALSE , usemulticore=TRUE)
  {
    #
    # OUTPUT
    # A dataframe (or list of dataframes if allsteps=TRUE) with the top K predicted edges. First column is the TF, second column the target gene, third column the score. If outfile is provided, the result is also written to the file OUTFILE (if allsteps=FALSE) or to several files OUTFILE1, OUTFILE1, ... (if allsteps=TRUE)

    # Check if we can run multicore
    if (usemulticore) {
      require(parallel)
    }

    # If needed, load expression data
    if (is.character(expdata))
      expdata <- read.table(expdata, header=1)

    # Gene names
    genenames <- colnames(expdata)
    ngenes <- length(genenames)

    # Normalize expression data for each gene
    if (normalizeexp)
      expdata <- scale(expdata)

    # If needed, load TF list
    if (is.null(tflist)) {
      # No TF list or file provided, we take all genes as TF
      tflist <- genenames
    } else if (length(tflist)==1 && is.na(match(tflist,genenames))) {
      # If this is a single string which is not a gene name, then it should be a file name
      tflist <- read.table(tflist , header=0)
      tflist <- as.matrix(tflist)[,1]
    }

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

    # Number of predictions to return
    Kmax <- ntf*ngenes
    if (K==-1) K<-Kmax
    K <- min(K,Kmax)

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
      targetname <- genenames[itarget]
      # Find the TF to be used for prediction (all TF except the target if the target is itself a TF)
      predTF <- tfindices[!match(tflist,targetname,nomatch=0)]
      r <- stabilityselection(as.matrix(expdata[,predTF]), as.matrix(expdata[,itarget]), nsplit=nsplit, nsteps=nstepsLARS, alpha=alpha)
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
      score <- mclapply(seq(ngenes),stabselonegene,mc.cores=detectCores()-1)
    } else {
      score <- lapply(seq(ngenes),stabselonegene)
    }
    # Rank scores
    edgepred <- list()
    for (i in seq(scorestokeep)) {
      # Combine all scores in a single vectors
      myscore <- unlist(lapply(score,function(x) x[,1,drop=FALSE]))
      ranki <- order(myscore,decreasing=TRUE)[1:K]
      edgepred[[i]] <- data.frame(list(tf=tflist[(ranki-1)%%ntf+1] , target=genenames[(ranki-1)%/%ntf+1] , score=myscore[ranki]))
    }

    # Print and return the result
    if (allsteps) {
      if (nchar(outfile)>0) {
        for (i in seq(length(edgepred))) {
          write.table( edgepred[[i]], file=paste(outfile,i,sep=''), quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
        }
      }
      return(edgepred)
    } else {
      if (nchar(outfile)>0) {
        write.table( edgepred[[1]], file=outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
      }
      return(edgepred[[1]])
    }
  }

