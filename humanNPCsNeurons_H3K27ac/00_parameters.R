outsPath="/projects/oliver/SCOaaP/outs/SCOP_2021_0139/scATACseq/03_PipelineOut/210218_kirkeby_cutNtag/cellranger/EndoPool_3/outs"
hashPath="/projects/oliver/SCOaaP/outs/SCOP_2021_0139/scATACseq/03_PipelineOut/210218_kirkeby_cutNtag/kite/HashPool_3/featurecounts"
projPath="/nfsdata/projects/jph712/cutntag_analysis/humanNPCsNeurons_H3K27ac"
ExptName="humanNPCsNeurons_H3K27ac"

#Set thresholds for experiment:
set.seed(1234)

cutoff_reads_min            = 3
cutoff_reads_max            = log10(50000)
cutoff_enhancer_ratio       = 0.35
EnhPromDNAseHyp_ratio       = 0.4
cutoff_fragmentsInRegions = 200

ndim = 30
BinRes=1
PeakRes=1

nproc=30


#Species specific parameters:
##mm10:
# annotations<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/mm10/EnsDb.Mmusculus.v79_getRangesAnnotation.RDS")
# GenomeSeqlengths<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/mm10/EnsDb.Mmusculus.v79_GenomeSeqlengths.RDS")
# blacklist_file="/nfsdata/projects/jph712/cutntag_analysis/mm10/mm10-blacklist.v2.bed"

##hg38:
annotations<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/hg38/EnsDb.Hsapiens.v86_getRangesAnnotation.RDS")
GenomeSeqlengths<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/hg38/EnsDb.Hsapiens.v86_GenomeSeqlengths.RDS")
blacklist_file="/nfsdata/projects/jph712/cutntag_analysis/hg38/hg38-blacklist.v2.bed"

