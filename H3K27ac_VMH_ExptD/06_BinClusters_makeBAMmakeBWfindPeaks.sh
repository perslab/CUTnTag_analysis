#!bin/bash

#06 Bin Clusters: Sinto make bam, deeptools make bigwig, find peaks
outsPath="/nfsdata/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/201217_VMH-PVH/cellranger/EndoPool_1/outs"
hashPath="/nfsdata/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/201217_VMH-PVH/kite/HashPool_1/featurecounts"
projPath="/nfsdata/projects/jph712/cutntag_analysis/H3K27ac_VMH_ExptD"
ExptName="H3K27ac_VMH_ExptD"
nproc=20
GRCm38=2652783500
blacklist_mm10="/projects/jph712/CUTnTag/H3K27ac_ARC_scCUTnTAG_B/deepTools/mm10-blacklist.v2.bed"


#######################
#06a sinto make bam file from filtered cells
#######################
cd $projPath/bam; 
sinto filterbarcodes -b ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam -p $nproc -c $projPath/bam/FilteredCellBarcodes_BinClusters.tsv

#######################
#06b deeptools bigwig
#######################

#Filtered cells:
for i in ${projPath}/bam/BinCluster*.bam; do 
pathName=${i%.bam};
clusterName=${pathName#$projPath/bam/}

samtools index $i
bamCoverage --bam $i --outFileName ${projPath}/bigwig/${clusterName}_nucFree.bigWig --outFileFormat bigwig --numberOfProcessors $nproc   --minFragmentLength 20 --maxFragmentLength 121 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
bamCoverage --bam $i --outFileName ${projPath}/bigwig/${clusterName}_Nucleosomal.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 120 --maxFragmentLength 800 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
bamCoverage --bam $i --outFileName ${projPath}/bigwig/${clusterName}_extendReads.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 20 --maxFragmentLength 800 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
done

#######################
#06c MACS2
#######################

for i in ${projPath}/bam/BinCluster*.bam; do 
pathName=${i%.bam};
clusterName=${pathName#$projPath/bam/}
# narrow
macs2 callpeak -t ${i} -g mm -f BAMPE -n ${ExptName}_$clusterName --outdir ${projPath}/macs2/BinClusters/narrow -q 0.01 -B --SPMR --keep-dup=1 2>&1
# broad
macs2 callpeak -t ${i} -g mm -f BAMPE -n ${ExptName}_$clusterName --outdir ${projPath}/macs2/BinClusters/broad -q 0.01 -B --SPMR --keep-dup=1 --broad-cutoff=0.1 --broad 2>&1
done



#######################
#06d SEACR
#######################
module load samtools
module load bedtools
genome="/home/cbmr/jph712/projects/CUTnTag/H3K27ac_ARC_scCUTnTAG_A/SEACR/tools/mm10GenomeSize.txt"
SEACR="/projects/jph712/CUTnTag/H3K27ac_VMH_bulk/SEACR/SEACR_1.3.sh"
 
mkdir $projPath/seacr -p
mkdir $projPath/seacr/bed -p
mkdir $projPath/seacr/tmpBAM -p

for i in ${projPath}/bam/BinCluster*.bam; do 
pathName=${i%.bam};clusterName=${pathName#$projPath/bam/}
#samtools
samtools view -b -@ $nproc -f 2 -F 1804 $i > ${projPath}/seacr/tmpBAM/${clusterName}_nodup.bam                               
samtools sort ${projPath}/seacr/tmpBAM/${clusterName}_nodup.bam > ${projPath}/seacr/tmpBAM/${clusterName}.sorted.bam
samtools index  ${projPath}/seacr/tmpBAM/${clusterName}.sorted.bam
### better genomecov (from protocols.io):
bedtools genomecov -bg -pc -ibam ${projPath}/seacr/tmpBAM/${clusterName}.sorted.bam > ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph
#SEACR:
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.1 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.1peaks.bed
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.05 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.05peaks.bed
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.01 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.01peaks.bed
done
