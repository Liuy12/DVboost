#' Load annotated SVs file
#'
#' After annotation procedure, SVs are then loaded and formatted properly for model training
#' @aliases ExampleData
#' @param input input file after annotation procedure
#' @return formatted data.frame to be used for model training including following fields:
#' \itemize{\item ChrA: first chromsome number;
#' \item Start: start position of first chromsome;
#' \item ChrB: second chromsome number;
#' \item End: end position of second chromsome;
#' \item SVLen: length of SV;
#' \item Germline: integer indicating germline artifact;
#' \item SVType: character value indicating SV type, could be one of the following 'DEL', 'DUP', 'INV', 'BND';
#' \item PE: number of paired end reads support;
#' \item SR: number of split reads support;
#' \item CNVMAP: 0/1 indicating whether the same SV is a known SV based on polymophic CNV database (https://www.nature.com/articles/nrg3871);
#' \item CNVR: 0/1 indicating whether the same SV is a known SV based on 1k genome database (https://www.nature.com/articles/nature15393);
#' \item avgL: average of strict mask L (depth of coverage is much lower than average). All strict mask averages are calcualted using start (100bp) and end (100bp) regions of the SV
#' \item avgH: average of strict mask H (depth of coverage is much higher than average).
#' \item avgZ: average of strict mask Z (too many reads with zero mapping quality overlap this position).
#' \item: avgQ: average of strict mask Q (average mapping quality is too low).
#' }
#' @examples
#' data(ExampleData, package='DVboost')
#' ### example data loaded and formatted with loadVariants
#' str(ExampleData)
#' @export

loadVariants <- function(input){
  r1 <- read.delim(input, stringsAsFactors = FALSE, header = T, check.names = FALSE)
  ## for delly calls
  if(sum(is.na(r1$SVLen[r1$SVType != 'BND']))/sum(r1$SVType != 'BND') > 0.3)
    r1$SVLen[r1$SVType != 'BND'] <- abs(r1$Start-r1$End)[r1$SVType != 'BND']
  SV.ID.vec <- paste(r1$ChrA,":",r1$Start,"_", r1$End, "@",r1$SVType, sep = "")
  r1$POLYMORPHIC_CNVMAP[is.na(r1$POLYMORPHIC_CNVMAP)] <- 0
  r1$POLYMORPHIC_CNVR[is.na(r1$POLYMORPHIC_CNVR)] <- 0
  avg.L.vec <- (r1$L_START_STRICTMASK + r1$L_END_STRICTMASK)/2
  avg.H.vec <- (r1$H_START_STRICTMASK + r1$H_END_STRICTMASK)/2
  avg.Z.vec <- (r1$Z_START_STRICTMASK + r1$Z_END_STRICTMASK)/2
  avg.Q.vec <- (r1$Q_START_STRICTMASK + r1$Q_END_STRICTMASK)/2


  file_PE = grep('*_PE$', colnames(r1), value = TRUE)
  file_SR = grep('*_SR$', colnames(r1), value = TRUE)

  ## CNVMAP from CNV map publication, CNVR from 1000genomes
  ## see /data2/external_data/Kocher_Jean-Pierre_m026645/s200929.CNV_SV_Strategy_for_WGS/annotation/polymorphic/readme.txt
  tmp.mtx <- data.frame(ChrA=r1$ChrA,
                        Start=r1$Start,
                        ChrB=r1$ChrB,
                        End=r1$End,
                        SVLen=r1$SVLen,
                        Germline=r1$GERMLINE_ARTIFACT,
                        SVType=r1$SVType,
                        PE = r1[,file_PE],
                        SR = r1[,file_SR],
                        CNVMAP=r1$POLYMORPHIC_CNVMAP,
                        CNVR=r1$POLYMORPHIC_CNVR,
                        avgL = avg.L.vec,
                        avgH = avg.H.vec,
                        avgZ = avg.Z.vec,
                        avgQ = avg.Q.vec,stringsAsFactors=F)
  rownames(tmp.mtx) <- SV.ID.vec
  return(tmp.mtx)
}
