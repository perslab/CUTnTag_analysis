outsPath="/nfsdata/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/201217_VMH-PVH/cellranger/EndoPool_1/outs"
hashPath="/nfsdata/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/201217_VMH-PVH/kite/HashPool_1/featurecounts"
projPath="/nfsdata/projects/jph712/cutntag_analysis/H3K27ac_VMH_ExptD"
ExptName="H3K27ac_VMH_ExptD"

#Set thresholds for experiment:
set.seed(1234)

cutoff_reads_min            = 3.5
cutoff_reads_max            = 5
cutoff_enhancer_ratio       = 0.2  

ndim = 30
BinRes=0.4
PeakRes=0.7
