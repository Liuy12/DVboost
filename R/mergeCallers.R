#' Merging SV events from different callers 
#'
#' @param df a data.frame containing technical annotations and q score from multiple callers (lumpy, manta, e.g.)
#' @param rocutoff threshold for merging events, e.g. events with reciprocal overlap bigger than the set value will be merged in the same cluster
#' @param summethod method to summarize numeric columns for overlapping events, can be either 'mean' or 'median',default to 'median'

#' 
#' @return a data.frame with overlapping events merged, inter-chromosomal events will not be merged
#' @examples
#' @export


mergeCallers <- function(df,rocutoff,summethod='median'){
  df$TYPE[is.na(df$TYPE)] <- 'CNV'
  df <- as.data.frame(df)
  summethod <- get(summethod)
  output <- foreach(i = unique(df$TYPE[which(df$CHROM1 == df$CHROM2)]),.inorder = TRUE,.combine = dplyr::bind_rows) %dopar% {
    message(paste0('merging ',i, ' events'))
    idx1 <- which(df$CHROM1 == df$CHROM2 & df$TYPE == i)
    df1 <- df[idx1,]
    gr0 <- GRanges(seqnames = df1$CHROM1,
                   ranges = IRanges(start = apply(df1[,c("POS","END")],1,min),
                                    end = apply(df1[,c("POS","END")],1,max)),
                   strand = '*',
                   ID = df1$ID)
    hits <- findOverlaps(gr0)
    x <- gr0[queryHits(hits)]
    y <- gr0[subjectHits(hits)]
    relative_overlap <- width(pintersect(x, y)) / pmax(width(x), width(y))
    hits <- hits[relative_overlap >= rocutoff]
    gr1 <- mergeConnectedRanges(gr0, hits)
    output1 <- data.frame(ID=paste0(i,'_clusterevent',1:length(gr1)),
                          CHROM1=seqnames(gr1),CHROM2=seqnames(gr1),POS=start(gr1),END=end(gr1),TYPE=i,GENOME=unique(df1$GENOME),
                          CALLER=sapply(gr1$revmap,function(j) paste0(unique(df1$CALLER[j]),collapse = ':'))
    )
    output2 <- t(sapply(gr1$revmap,function(j) apply(df1[j,c(8,10:ncol(df1))],2,function(k) summethod(k,na.rm = T))))
    output3 <- cbind(output1,output2)
    output3[,match(colnames(df1),colnames(output3))]
  }
  output <- rbind(output,df[df$CHROM1 != df$CHROM2,])
  return(output)
}



mergeConnectedRanges <- function(x, hits)
{
  stopifnot(is(x, "GenomicRanges"))
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  stopifnot(queryLength(hits) == length(x))
  clusters <- extractClustersFromSelfHits(hits)
  ans <- range(extractList(x, clusters))
  if (any(elementNROWS(ans) != 1L))
    stop(wmsg("some connected ranges are not on the same ",
              "chromosome and strand, and thus cannot be ",
              "merged"))
  ans <- unlist(ans)
  mcols(ans)$revmap <- clusters
  ans
}




## Extract clusters from Hits object.
extractClustersFromSelfHits <- function(hits)
{
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  hits <- union(hits, t(hits))
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cid <- seq_len(queryLength(hits))  # cluster ids
  while (TRUE) {
    h <- Hits(qh, cid[sh],
              queryLength(hits), subjectLength(hits))
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
  }
  unname(splitAsList(seq_len(queryLength(hits)), cid))
}

