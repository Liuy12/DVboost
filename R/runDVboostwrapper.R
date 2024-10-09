#' Core function for DVboost model training
#'
#' Building training models using subset of total SVs (Only 'DEL' type)
#' @param var.atr.mtx a data.frame containing subset of total SVs (Only 'DEL' type)
#' @param var.ID.vec a vector of characters indicating IDs for each SV in \code{var.atr.mtx}
#' @param is.known.var.vec a vector contains 0 and 1s in the same order as \code{var.atr.mtx} indicating whether a SV is known(1) or novel(0)
#' @param output.DIR.name path to output folder where all figure files and data files will be saved
#' @param input.sample.ID character value used as prefix for output files
#' @param fitting.verbose logical value, if TRUE, it will print out progress and performance indicators. Default to FALSE
#' @param min.N.known.var minimum number of known SVs needed to build training model. Default to 50
#' @param bySVlength logical value, if true, performance will also be evaluated seperately for different categories of SV length. See details section for more
#' @details When bySVlength is set to true, training performance will be evaluated seperately for 5 different categories of SV length. To be specific,
#' the 5 SV length categories are 50-100bp, 100-500bp,500-1k, 1k-10k, >10k. For each SV length category, training performance will be evaluated and exported as
#' figures and text files.
#' @return a list contains DVboost.res which is essentially a \code{\link{gbm.object}} with several additional fields:
#' \itemize{
#' \item fitted.values: converted to probability based on \code{fit} field
#' \item ID: supplied IDs for SVs
#' \item is.known.variant: 0/1 indicating whether the SV is known(1) or novel(0)
#' \item DVboost.Q.score: Q scores for SVs
#' }
#' Besides, this will also save a set of figures and text files to \code{output.DIR.name}. Please see \code{\link{metricView}} for details
#'
#' @examples
#' data(ExampleData, package='DVboost')
#' sample <- 'NA12878'
#' outdir <- getwd()
#' tmp.mtx.DEL <- ExampleData[ExampleData$SVType == 'DEL',]
#' truth.vec <- tmp.mtx.DEL$CNVMAP == 1 | tmp.mtx.DEL$CNVR ==1
#' is.semi.truth.vec <- as.numeric(truth.vec)
#' DVb.res <- runDVboostwrapper( var.atr.mtx = tmp.mtx.DEL, var.ID.vec = rownames(tmp.mtx.DEL),
#'                              is.known.var.vec = is.semi.truth.vec,
#'                              output.DIR.name = outdir,input.sample.ID=sample, bySVlength=FALSE)
#'
#' @seealso
#' \code{\link{metricView}}, \code{\link{fitDVboostmodel}}, \code{\link{DVboostQscore}}
#'
#' @export

runDVboostwrapper <- function(var.atr.mtx, var.ID.vec, is.known.var.vec,
                                output.DIR.name,input.sample.ID,
                                fitting.verbose = FALSE, min.N.known.var = 50, bySVlength=FALSE,caller)
{
  if (!file.exists(output.DIR.name) ){
    cat('The specified output directory \"',output.DIR.name,'\" does not exist, creating directory...\n',sep='')
    dir.create(file.path(output.DIR.name),showWarnings=FALSE)
  }

  DVboost.score.vec <- mat.or.vec(length(var.ID.vec),1) + NA
  names(DVboost.score.vec) <- var.ID.vec

  sel.var.idx <- 1:nrow(var.atr.mtx)
  sel.var.atr.mtx <- var.atr.mtx[sel.var.idx,]
  rownames(sel.var.atr.mtx) <- var.ID.vec[sel.var.idx]
  train.label.vec <- is.known.var.vec[sel.var.idx]

  cat("\n First 10 rows of input variant by attribute matrix... \n")
  print(head(sel.var.atr.mtx,10))

  cat("\n Last 10 rows of input variant by attribute matrix... \n")
  print(tail(sel.var.atr.mtx,10))
  DV.fit.res1 <- fitDVboostmodel(input.mtx=sel.var.atr.mtx, is.known.variant = train.label.vec,
                                   fitting.verbose = fitting.verbose,
                                   min.N.known.var = min.N.known.var,caller=caller)

  DVboost.score.vec[DV.fit.res1$ID] <- DV.fit.res1$DVboost.Q.score
  if(bySVlength){
    lenClass <- c('50_100', '100_500', '500_1000','1000_10000','10000')
    for(i in 1:5){
      plot.filename <- paste(output.DIR.name,"/",input.sample.ID,"_DVboost_summaryPlot_",lenClass[i], '.png',sep = "")
      DV.sampleQC.filename <- paste(output.DIR.name,"/",input.sample.ID,"_DVboost_summaryQC_", lenClass[i], '.txt',sep = "")
      tmp <- strsplit(lenClass[i], '_')
      tmp <- as.numeric(tmp[[1]])
      if(length(tmp)==1) idx <- which(sel.var.atr.mtx$LENGTH >= tmp[1]) else idx <-  which(sel.var.atr.mtx$LENGTH >= tmp[1] & sel.var.atr.mtx$LENGTH < tmp[2])
      if(length(unique(DV.fit.res1$is.known.variant[idx])) == 2)
        metricView(DV.fit.res = DV.fit.res1, input.sample.ID=input.sample.ID, plot.filename = plot.filename, DV.sampleQC.filename = DV.sampleQC.filename, subset = idx)
    }
  }
  ## generate figures for all variants
  plot.filename <- paste(output.DIR.name,"/",input.sample.ID,"_DVboost_summaryPlot.png",sep = "")
  DV.sampleQC.filename <- paste(output.DIR.name,"/",input.sample.ID,"_DVboost_summaryQC.txt",sep = "")
  metricView(DV.fit.res = DV.fit.res1, input.sample.ID=input.sample.ID, plot.filename = plot.filename, DV.sampleQC.filename = DV.sampleQC.filename, subset = 1:length(var.ID.vec))

  return(list(DVboost.res = DV.fit.res1))
}
