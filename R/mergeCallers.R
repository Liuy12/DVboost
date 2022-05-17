#' Merging SV events from different callers 
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


mergeCallers <- function(df,rocutoff){
  for(i in unique(df$TYPE[which(df$CHROM1 == df$CHROM2)])){
    idx1 <- which(df$CHROM1 == df$CHROM2 & df$TYPE == i)
    gr0 <- GRanges(seqnames = df$CHROM1,
                   ranges = IRanges(start = apply(df[,c("POS","END")],1,min),
                                    end = apply(df[,c("POS","END")],1,max)),
                   strand = '*',
                   ID = df$ID)
    hits <- findOverlaps(gr0)
    x <- gr0[queryHits(hits)]
    y <- gr0[subjectHits(hits)]
    relative_overlap <- width(pintersect(x, y)) / pmax(width(x), width(y))
    hits <- hits[relative_overlap >= rocutoff]
    gr1 <- mergeConnectedRanges(gr0, hits)
  }
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

