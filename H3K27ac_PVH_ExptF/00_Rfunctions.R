
remove_blacklisted2Granges<-function(peakTable){
  colnames(peakTable)<-c("chrom", "chromStart", "ChromEnd", "name", "score", "strand", "signalValue", "pValue", "qvalue")
  GRangesObject=GenomicRanges::GRanges(seqnames=peakTable[,1],
                                       ranges=IRanges(start=peakTable[,2],
                                                      end=peakTable[,3]),
                                       mcols=peakTable[, 4:4:ncol(peakTable)])
  colnames(GRangesObject@elementMetadata)<-c( "name", "score", "stranded", "signalValue", "pValue", "qvalue")
  
  
  blacklist.df = read_tsv(file="/projects/jph712/CUTnTag/H3K27ac_ARC_scCUTnTAG_B/deepTools/mm10-blacklist.v3.bed" ,col_names = F) 
  blacklist.gr=GenomicRanges::GRanges(seqnames=blacklist.df$X1,
                                      ranges=IRanges(start=as.numeric(blacklist.df$X2),
                                                     end=as.numeric(blacklist.df$X3)))
  simpleGR_diff<-GenomicRanges::setdiff(GRangesObject, blacklist.gr, ignore.strand=TRUE)
  
  binCluster_GR_white<-subsetByOverlaps(GRangesObject,simpleGR_diff, ignore.strand=TRUE)
  return(binCluster_GR_white)
}


bed3_2white_Granges<-function(peakTable){
  colnames(peakTable)<-c("chrom", "chromStart", "ChromEnd")
  GRangesObject=GenomicRanges::GRanges(seqnames=peakTable[,1],
                                       ranges=IRanges(start=peakTable[,2],
                                                      end=peakTable[,3])
                                       )
  blacklist.df = read_tsv(file="/projects/jph712/CUTnTag/H3K27ac_ARC_scCUTnTAG_B/deepTools/mm10-blacklist.v3.bed" ,col_names = F) 
  blacklist.gr=GenomicRanges::GRanges(seqnames=blacklist.df$X1,
                                      ranges=IRanges(start=as.numeric(blacklist.df$X2),
                                                     end=as.numeric(blacklist.df$X3)))
  simpleGR_diff<-GenomicRanges::setdiff(GRangesObject, blacklist.gr, ignore.strand=TRUE)
  
  GR_white<-subsetByOverlaps(GRangesObject,simpleGR_diff, ignore.strand=TRUE)
  return(GR_white)
}

bed6_ToWhite_GRanges<-function(peakTable){
  colnames(peakTable)<-c("chrom", "chromStart", "ChromEnd", "width", "score", "name")
  GRangesObject=GenomicRanges::GRanges(seqnames=peakTable[,1],
                                       ranges=IRanges(start=peakTable[,2],
                                                      end=peakTable[,3]),
                                       mcols=peakTable[, 4:ncol(peakTable)])
  colnames(GRangesObject@elementMetadata)<-c(  "PeakName", "score", "location")
  
  
  blacklist.df = read_tsv(file="/projects/jph712/CUTnTag/H3K27ac_ARC_scCUTnTAG_B/deepTools/mm10-blacklist.v3.bed" ,col_names = F) 
  blacklist.gr=GenomicRanges::GRanges(seqnames=blacklist.df$X1,
                                      ranges=IRanges(start=as.numeric(blacklist.df$X2),
                                                     end=as.numeric(blacklist.df$X3)))
  simpleGR_diff<-GenomicRanges::setdiff(GRangesObject, blacklist.gr, ignore.strand=TRUE)
  
  binCluster_GR_white<-subsetByOverlaps(GRangesObject,simpleGR_diff, ignore.strand=TRUE)
  return(binCluster_GR_white)
}
