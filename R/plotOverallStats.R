#' Generate overall performance metrics for multiple samples
#'
#' If you are running DVboost on a set of samples and output to the same directory, the you can use \code{plotOverallStats} to examine the overall
#' performance for all samples
#'
#' @param outputdir full path to output directory as specified in \code{runDVboostwrapper}
#' @param type specify which type of metric is preferred, could be one of 'Fscore', 'AUROC' or 'Both'
#' @param bylen logical value indicating whether to examine the performance by SV length instead. Default to FALSE
#'
#' @return
#' a ggplot object
#'
#' @export
#' @import ggplot2
#'

plotOverallStats <- function(outputdir, type=c('Both', 'Fscore','AUROC'), bylen=FALSE){
  if(!bylen){
    files <- list.files(outputdir, '*_DVboost_summaryQC.txt', full.names = T)
    dataMat <- do.call(rbind, lapply(files, function(i) read.delim(i, header = T)))
    if(!type %in% c('Fscore','AUROC', 'Both')) stop('The specified type is not valid, has to be one of Fscore, AUROC or Both')
    if(type == 'Fsocre')
      df1 <- data.frame(sample=1:nrow(dataMat),
                        values=dataMat[,2],
                        type=rep('Fscore', times=nrow(dataMat)))
    if(type == 'AUROC')
      df1 <- data.frame(sample=1:nrow(dataMat),
                        values=dataMat[,4],
                        type=rep('Fscore', times=nrow(dataMat)))
    if(type == 'Both')
      df1 <- data.frame(sample=rep(1:nrow(dataMat),times=2),
                        values=c(dataMat[,2], dataMat[,4]),
                        type=rep(c('Fscore', 'AUROC'), each=nrow(dataMat)))

    gp1 <- ggplot(df1) + geom_line(aes(x=sample, y=values, color=type)) + theme_classic() + labs(x='Samples', y='Values') + ylim(0,1) + scale_x_continuous(breaks = 1:nrow(dataMat))
    return(gp1)
  }
  else{
    svclass <- c('50_100', '100_500', '500_1000', '1000_10000', '10000')
    filenames <- lapply(svclass, function(i) list.files(outputdir, paste0('*_DVboost_summaryQC_', i, '.txt'), full.names = T))
    nsamples <- max(sapply(filenames, function(i) length(i)))
    svclass <- svclass[sapply(filenames, function(i) length(i)) == nsamples]
    xaxis <- c('50-100', '100-500', '500-1k', '1k-10k', '>10k')[sapply(filenames, function(i) length(i)) == nsamples]
    filenames <- filenames[sapply(filenames, function(i) length(i)) == nsamples]
    dataMat <- lapply(filenames, function(i) do.call(rbind, lapply(i, function(j) read.delim(j, header = T, stringsAsFactors = F))))
    if(type=='Both')
      df1 <- cbind(sapply(dataMat, function(i) i[,2]), sapply(dataMat, function(i) i[,4]))
    tmp <- stack(as.data.frame(df1))
    tmp$idx1 <- rep(c('Fscore', "AUROC"), each = length(svclass)*length(filenames[[2]]))
    xindx <- (length(filenames)+1):(length(filenames)*2)
  }
  if(type=='Fscore'){
    df1 <- sapply(dataMat, function(i) i[,2])
    tmp <- stack(as.data.frame(df1))
    tmp$idx1 <- rep('Fscore', each = length(svclass)*length(filenames[[2]]))
    xindx <- 1:length(filenames)
  }
  if(type == 'AUROC'){
    df1 <- sapply(dataMat, function(i) i[,4])
    tmp <- stack(as.data.frame(df1))
    tmp$idx1 <- rep("AUROC", each = length(svclass)*length(filenames[[2]]))
    xindx <- 1:length(filenames)
  }
  df1 <- do.call(cbind, lapply(dataMat, function(i) i[,5:6]))
  num1 <- apply(df1, 2, median)
  labels <- sapply(1:length(filenames), function(i) paste0(num1[c(2*i-1, 2*i)], collapse = '/'))
  gp1 <- ggplot() + geom_boxplot(aes(x=ind, y = values, color = idx1), data = tmp) +
    theme_classic() + labs(x='SV Length', y='Values') +
    guides(color = guide_legend(title=NULL)) + geom_label(aes(x=xindx, y=0.2, label = labels), size=2) +
    annotate("text", label = '#pos/#neg', x = xindx[length(xindx)], y = 0.15, size = 3) +
    scale_x_discrete(labels=rep(xaxis, times= (xindx[length(xindx)])/length(svclass))) +
    ylim(0,1)
  return(gp1)
}
