df2vcf <- function(df, sampleid, metainfo, output){
  #df <- read.delim(df, header = T, stringsAsFactors = F, check.names = F)
  infofields <- c("CHROM2", "END", "TYPE", "GENOME", "LENGTH", "CALLER", "INNER_RANGE", "OUTER_RANGE", "oneKSV", "cnvMap", "dt_score", "avgsm_H", "avgsm_L", "avgsm_Z", "avgsm_Q", "Qscore")
  formatfields <- c("L2R", "CN", "PE", "SR")
  df[,formatfields][is.na(df[,formatfields])] <- "."
  dfout <- data.frame("#CHROM" = df$CHROM1,
                      POS = df$POS,
                      ID = df$ID, 
                      REF = "N",
                      ALT = paste0("<",df$TYPE, ">"),
                      QUAL = ".",
                      FILTER = ".",
                      INFO = sapply(1:nrow(df), function(i) paste0(infofields, "=",df[i,infofields],collapse = ";")),
                      FORMAT = paste0(formatfields,collapse = ":"),
                      sampleid = sapply(1:nrow(df), function(i) paste0(df[i,formatfields],collapse = ":")),check.names = F)
  colnames(dfout)[10] <- sampleid
  callers <- paste0(unique(unlist(strsplit(df$CALLER, ":"))),collapse = "_")
  headerinfo <- paste0("##fileformat=VCFv4.2\n##fileDate=",gsub("-","", Sys.Date()), "##callers=",callers)
  metainfo <- readLines(metainfo)
  writeLines(c(headerinfo,metainfo), output)
  write.table(dfout, output, append = T, quote = F, sep = "\t", row.names = F, col.names = T)
}


