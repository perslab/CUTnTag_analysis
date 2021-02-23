outsPath="/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/210128_ExptFG_ExptH/cellranger/EndoPool_1/outs"
hashPath="/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/210128_ExptFG_ExptH/kite/HashPool_1_rc/featurecounts"
projPath="/nfsdata/projects/jph712/cutntag_analysis/H3K27ac_PVH_ExptF"
ExptName="H3K27ac_PVH_ExptF"

#Set thresholds for experiment:
set.seed(1234)

cutoff_reads_min            = 3.2
cutoff_reads_max            = log10(50000)
cutoff_enhancer_ratio       = 0.2 

ndim = 20
BinRes=1
PeakRes=1

nproc=30
genome="hg38"


#Species specific parameters:
##mm10:
# annotations<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/mm10/EnsDb.Mmusculus.v79_getRangesAnnotation.RDS")
# GenomeSeqlengths<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/mm10/EnsDb.Mmusculus.v79_GenomeSeqlengths.RDS")
# blacklist_file="/nfsdata/projects/jph712/cutntag_analysis/mm10/mm10-blacklist.v2.bed"

##hg38:
annotations<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/hg38/EnsDb.Hsapiens.v86_getRangesAnnotation.RDS")
GenomeSeqlengths<-readRDS("/nfsdata/projects/jph712/cutntag_analysis/hg38/EnsDb.Hsapiens.v86_GenomeSeqlengths.RDS")
blacklist_file="/nfsdata/projects/jph712/cutntag_analysis/hg38/hg38-blacklist.v2.bed"


