#' Fitting a DVboost model based on annotated SVs
#'
#' Internal function called by \code{\link{runDVboostwrapper}}
#'
#' @param input.mtx annotated SVs loaded and formatted via \code{\link{loadVariants}}
#' @param is.known.variant numeric vector of 0/1 in the same order as \code{input.mtx} indicating whether a SV is known(1) or novel(0)
#' @param fitting.verbose logical value, if TRUE, it will print out progress and performance indicators. Default to FALSE
#' @param min.N.known.var minimum number of known SVs needed to build training model. Default to 50
#' @param sampleid same sample name used for annotation procedure
#'
#' @return DV.res which is essentially a \code{\link{gbm.object}} with several additional fields:
#' \itemize{
#' \item fitted.values: converted to probability based on \code{fit} field
#' \item ID: supplied IDs for SVs
#' \item is.known.variant: 0/1 indicating whether the SV is known(1) or novel(0)
#' \item DVboost.Q.score: Q scores for SVs
#' }
#'
#' @examples
#' data(ExampleData, package='DVboost')
#' sample <- gsub('_PE', '', grep('_PE$', colnames(ExampleData), value=TRUE))
#' outdir <- getwd()
#' tmp.mtx.DEL <- ExampleData[ExampleData$SVType == 'DEL',]
#' truth.vec <- tmp.mtx.DEL$CNVMAP == 1 | tmp.mtx.DEL$CNVR ==1
#' is.semi.truth.vec <- as.numeric(truth.vec)
#' DV.fit.res1 <- fitDVboostmodel(input.mtx=tmp.mtx.DEL, is.known.variant = is.semi.truth.vec,
#' sampleid=sample)
#'
#' @seealso
#' \code{\link{runDVboostwrapper}}
#'
#' @export
#' @import gbm

fitDVboostmodel <- function( input.mtx, is.known.variant,
                             fitting.verbose = FALSE, min.N.known.var = 50,sampleid)
{
  cat("\n whether variants are known (1 = known, 0 = others) \n")
  print(table(is.known.variant))
  if( length(which(is.known.variant==1)) < min.N.known.var ){
    stop( paste("######## (DVboost error)\n \t",
                "Too few known variants; less than specified 'min.N.known.var' =", min.N.known.var) )
  }
  cat("\n fitting DVboost model ... \n")
  sel.var.idx <- 1:nrow(input.mtx)
  sel.var.atr.mtx <- input.mtx[sel.var.idx,]
  train_is.known.variant <- is.known.variant[sel.var.idx]
  train.df<-as.data.frame(sel.var.atr.mtx)
  y<-train_is.known.variant

  DV.res=gbm(y~ PE + SR + avgL + avgH + avgZ + avgQ  + Germline,
             distribution="adaboost",
             data=train.df,
             n.trees=2000,cv.folds=5,n.cores=1,
             shrinkage=.01,
             n.minobsinnode=20)

  bestTreeForPrediction= gbm.perf(DV.res, plot.it = F)
  DV.res$fitted.values <-  plogis(2*DV.res$fit)
  # convert it to the 0-1 scale since the adaboost method gives the predictions on logit scale.
  # http://stats.stackexchange.com/questions/37497/how-to-use-r-gbm-with-distribution-adaboost
  DV.res$ID <- rownames(sel.var.atr.mtx)
  DV.res$is.known.variant <- train_is.known.variant
  DVboost.ECDF <- ecdf(DV.res$fitted.values[which(DV.res$is.known.variant==1)])
  DV.res$DVboost.Q.score <- DVboost.ECDF(DV.res$fitted.values)
  return(DV.res)
}
