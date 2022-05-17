#' Annotate SVs with polymorphic CNVs and genomic features
#'
#' @param input input file after converter scripts for each CNV/SV caller
#' @param genome genome version, e.g. hg19 or hg38
#' @param rocutoff reciprocal overlapping threshold to identify overlaps with polymorphic SV database
#' @param buffersize default buffer size to use if not supplied by caller (e.g. CIPOS is not available)
#' @param seltype select a subset of SVs matching seltype for annotation 
#' @param oneKSV path to .bed file containing polymorphic events from 1000 genome database
#' @param cnvMap path to .bed file containing polymorphic events from CNV map
#' @param dangerTrack path to .bed file containing dangertrack annotations
#' @param strictMask path to .fa file containing strictmask information
#' 
#' @return formatted data.frame to be used for model training including following fields besides the supplied fields:
#' \itemize{
#' \item oneKSV: 0/1 indicating whether the same SV is a known SV based on 1k genome database (https://www.nature.com/articles/nature15393);
#' \item cnvMap: 0/1 indicating whether the same SV is a known SV based on polymophic CNV database (https://www.nature.com/articles/nrg3871);
#' \item dt_score: danger track score from here (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5405793/);
#' \item avgsm_H: average of strict mask H (depth of coverage is much higher than average).
#' \item avgsm_L: average of strict mask L (depth of coverage is much lower than average).
#' \item avgsm_Z: average of strict mask Z (too many reads with zero mapping quality overlap this position).
#' \item avgsm_Q: average of strict mask Q (average mapping quality is too low).
#' }
#' @examples
#' @export

######### with one line per event, CHR1 & CHR2 to incorporate BND or CTX
annotateSV <- function(input,genome='hg38',rocutoff=0.5,buffersize=1000,seltype='DEL',oneKSV,cnvMap,dangerTrack,strictMask){
  t1 <- Sys.time()
  message('loading input data......')
  input <- fread(input,header = T,sep='\t',stringsAsFactors = F,check.names = F)
  setnames(input,c(1:11),c("ID","CHROM1","CHROM2","POS","END","TYPE",'GENOME','LENGTH','CALLER',"INNER_RANGE","OUTER_RANGE"))
  input <- input[TYPE %in% seltype,]
  genome <- SeqinfoForUCSCGenome(genome = genome)
  message('identifying overlaps with 1000 genome and CNVMap polymorphic SVs......')
  idx <- which(input[,apply(.SD,1,function(i) !any(is.na(i))),.SDcols=c("POS","END")])
  if(length(idx) != nrow(input)) stop('NA exists for POS or END columns')
  genomelength <- data.frame(CHROM=names(seqlengths(genome)),size=seqlengths(genome))
  input[,size1:= genomelength$size[match(CHROM1,genomelength$CHROM)]]
  input[,size2:= genomelength$size[match(CHROM2,genomelength$CHROM)]]
  input[INNER_RANGE > size1,INNER_RANGE:= size1]
  input[OUTER_RANGE > size2,OUTER_RANGE:= size2]
  ### only look for overlap for events with breakends on the same chrom
  idx <- which(input$CHROM1 == input$CHROM2)
  gr0 <- GRanges(seqnames = input[idx,CHROM1],ranges = IRanges(start=input[idx,apply(.SD,1,function(i) min(i)),.SDcols=c('POS','END')],end =input[idx,apply(.SD,1,function(i) max(i)),.SDcols=c('POS','END')]),strand = '*')
  gr1 <- import(oneKSV, genome = genome)
  gr2 <- import(cnvMap, genome = genome)
  ol1 <- findOverlaps_ro(gr0,gr1,rocutoff)
  ol2 <- findOverlaps_ro(gr0,gr2,rocutoff)
  input[,oneKSV :=0,]
  input[,cnvMap := 0,]
  input[idx[queryHits(ol1)],oneKSV := 1]
  input[idx[queryHits(ol2)],cnvMap := 1]
  message('calculating danger track score......')
  dangertrack <- fread(dangerTrack,header=F)
  dt <- foreach(i =1:nrow(input),.inorder = TRUE,.combine = c,.packages = 'data.table') %dopar% {
    idx <- dangertrack[V1 == input[i,CHROM1] & V2 <= input[i,POS] & V3 >= input[i,POS],which=TRUE]
    start <- ifelse(length(idx),dangertrack[idx,V5],0)
    idx1 <- dangertrack[V1 == input[i,CHROM2] & V2 <= input[i,END] & V3 >= input[i,END],which=TRUE]
    end <- ifelse(length(idx1),dangertrack[idx1,V5],0)
    mean(start,end)
  }
  message('calculating strick mask score......')
  strictmask <- readBStringSet(strictMask)
  names(strictmask) <- gsub(' \\[.*\\]','',names(strictmask))
  sm <- foreach(i = 1L:nrow(input),.inorder = TRUE,.combine = dplyr::bind_rows,.packages = 'data.table') %dopar% {
    if(! ((input[i,CHROM1] %in% names(strictmask)) & (input[i,CHROM2] %in% names(strictmask))))
      tmp <- c(0,0,0,0) else{
        if(is.na(input[i,INNER_RANGE]) || input[i,INNER_RANGE] >= input[i,POS]) start <- max(0,input[i,POS] - buffersize) else start <- input[i,INNER_RANGE]
        seq1 <- subseq(strictmask[[which(names(strictmask)==input[i,CHROM1])]],start = start,end = input[i,POS])
        if(is.na(input[i,OUTER_RANGE]) || input[i,OUTER_RANGE] <= input[i,END]) end <- min(input[i,END] + buffersize,input[i,size2]) else end <- input[i,OUTER_RANGE]
        seq2 <- subseq(strictmask[[which(names(strictmask)==input[i,CHROM2])]],start = input[i,END],end = end)
        tmp <- as.numeric(strictmaskMetrics(seq1,as.prob=T)+strictmaskMetrics(seq2,as.prob=T)/2)
      }
    names(tmp) <- c('H','L','Z','Q')
    tmp
  }
  input[,dt_score := as.numeric(dt)]
  for(i in 1:4) input[,c('avgsm_H','avgsm_L','avgsm_Z','avgsm_Q')[i] := as.numeric(sm[[i]]),]
  input[,size1:=NULL];input[,size2:=NULL]
  t2 <- Sys.time()
  message(paste0('total time spent: ',as.double(t2-t1),units(t2-t1)))
  return(input)
}

findOverlaps_ro <- function(gr0,gr1,rocutoff){
  ol1 <- findOverlaps(gr0, gr1, ignore.strand = TRUE)
  overlaps <- pintersect(gr0[queryHits(ol1)], gr1[subjectHits(ol1)])
  percentOverlap1 <- width(overlaps) / width(gr0[queryHits(ol1)])
  percentOverlap2 <- width(overlaps) / width(gr1[subjectHits(ol1)])
  ol1 <- ol1[percentOverlap1 > rocutoff & percentOverlap2 > rocutoff]
  ol1
}


strictmaskMetrics <- function(seq1,as.prob=T){
  pH <- letterFrequency(seq1,letters = 'H',as.prob = as.prob)
  pL <- letterFrequency(seq1,letters = 'L',as.prob = as.prob)
  pZ <- letterFrequency(seq1,letters = 'Z',as.prob = as.prob)
  pQ <- letterFrequency(seq1,letters = 'Q',as.prob = as.prob)
  c(pH,pL,pZ,pQ)
}



