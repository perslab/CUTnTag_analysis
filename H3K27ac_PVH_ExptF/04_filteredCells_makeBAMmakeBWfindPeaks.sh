#!/bin/bash

#04 Filtered cells: Sinto make bam, deeptools make bigwig, find peaks
outsPath="/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/210128_ExptFG_ExptH/cellranger/EndoPool_1/outs"
hashPath="/projects/oliver/SCOaaP/outs/SCOP_2020_0113/scATACseq/03_PipelineOut/210128_ExptFG_ExptH/kite/HashPool_1_rc/featurecounts"
projPath="/nfsdata/projects/jph712/cutntag_analysis/H3K27ac_PVH_ExptF"
ExptName="H3K27ac_PVH_ExptF"


nproc=30
GRCm38=2652783500
blacklist_mm10="/projects/jph712/CUTnTag/H3K27ac_ARC_scCUTnTAG_B/deepTools/mm10-blacklist.v2.bed"

#filter, sort and index cellranger possorted.bam file:
mkdir $projPath/bam/tmpUnsorted -p
samtools view -b -f 2 -F 1804 -@ $nproc -o ${projPath}/bam/tmpUnsorted/${ExptName}_NoDupUnsorted.bam $outsPath/possorted_bam.bam
samtools sort ${projPath}/bam/tmpUnsorted/${ExptName}_NoDupUnsorted.bam > ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam
samtools index ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam

#######################
#04a sinto make bam file from filtered cells
#######################
#pip install sinto
mkdir $projPath/bam -p

cd $projPath/bam; sinto filterbarcodes -b ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam -p $nproc -c $projPath/bam/FilteredCellBarcodes.tsv
samtools index $projPath/bam/${ExptName}_FilteredCells.bam

#######################
#04b deeptools bigwig
#######################
mkdir $projPath/bigwig -p

#all cellranger fragments:
bamCoverage --bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam --outFileName ${projPath}/bigwig/${ExptName}_AllCellrangerFragments_nucFree.bigWig --outFileFormat bigwig --numberOfProcessors $nproc   --minFragmentLength 20 --maxFragmentLength 121 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
bamCoverage --bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam --outFileName ${projPath}/bigwig/${ExptName}_AllCellrangerFragments_Nucleosomal.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 120 --maxFragmentLength 800 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
bamCoverage --bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam --outFileName ${projPath}/bigwig/${ExptName}_AllCellrangerFragments_extendReads.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 20 --maxFragmentLength 800 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38

#Filtered cells:
bamCoverage --bam $projPath/bam/${ExptName}_FilteredCells.bam --outFileName ${projPath}/bigwig/${ExptName}_filteredCells_nucFree.bigWig --outFileFormat bigwig --numberOfProcessors $nproc   --minFragmentLength 20 --maxFragmentLength 121 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
bamCoverage --bam $projPath/bam/${ExptName}_FilteredCells.bam --outFileName ${projPath}/bigwig/${ExptName}_filteredCells_Nucleosomal.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 120 --maxFragmentLength 800 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38
bamCoverage --bam $projPath/bam/${ExptName}_FilteredCells.bam --outFileName ${projPath}/bigwig/${ExptName}_filteredCells_extendReads.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 20 --maxFragmentLength 800 --blackListFileName $blacklist_mm10 -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $GRCm38



#######################
#04c MACS2
#######################
mkdir macs2/BinClusters/narrow -p
mkdir macs2/BinClusters/broad -p

####
## For all cellranger fragments:
for i in ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam; do pathName=${i%.bam};
clusterName=${pathName#$projPath/bam/}
# narrow
macs2 callpeak -t ${i} -g mm -f BAMPE -n ${ExptName}_$clusterName --outdir ${projPath}/macs2/BinClusters/narrow -q 0.01 -B --SPMR --keep-dup=1 2>&1
# broad
macs2 callpeak -t ${i} -g mm -f BAMPE -n ${ExptName}_$clusterName --outdir ${projPath}/macs2/BinClusters/broad -q 0.01 -B --SPMR --keep-dup=1 --broad-cutoff=0.1 --broad 2>&1
done

####
##For all filtered cells.

for i in $projPath/bam/${ExptName}_FilteredCells.bam; do 
clusterName=${pathName#$projPath/bam/}
# narrow
macs2 callpeak -t ${i} -g mm -f BAMPE -n ${ExptName}_$clusterName --outdir ${projPath}/macs2/BinClusters/narrow -q 0.01 -B --SPMR --keep-dup=1 2>&1
# broad
macs2 callpeak -t ${i} -g mm -f BAMPE -n ${ExptName}_$clusterName --outdir ${projPath}/macs2/BinClusters/broad -q 0.01 -B --SPMR --keep-dup=1 --broad-cutoff=0.1 --broad 2>&1
done



#######################
#04d SEACR
#######################

module load samtools
module load bedtools
genome="/home/cbmr/jph712/projects/CUTnTag/H3K27ac_ARC_scCUTnTAG_A/SEACR/tools/mm10GenomeSize.txt"
SEACR="/projects/jph712/CUTnTag/H3K27ac_VMH_bulk/SEACR/SEACR_1.3.sh"
 
mkdir $projPath/seacr -p
mkdir $projPath/seacr/bed -p
mkdir $projPath/seacr/tmpBAM -p

#########
## For all cellranger fragments:
##For all filtered cells.
for i in $projPath/bam/${ExptName}_FilteredCells.bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam; do 
pathName=${i%.bam};clusterName=${pathName#$projPath/bam/}
#samtools
samtools view -b -@ $nproc -f 2 -F 1804 $i > ${projPath}/seacr/tmpBAM/${clusterName}_nodup.bam                               
samtools sort ${projPath}/seacr/tmpBAM/${clusterName}_nodup.bam > ${projPath}/seacr/tmpBAM/${clusterName}.sorted.bam
samtools index  ${projPath}/seacr/tmpBAM/${clusterName}.sorted.bam

### better genomecov (from protocols.io):
bedtools genomecov -bg -pc -ibam ${projPath}/seacr/tmpBAM/${clusterName}.sorted.bam > ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph
#SEACR:
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.1 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.1peaks.bed
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.01 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.01peaks.bed
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.001 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.001peaks.bed
$SEACR ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.001 norm relaxed ${projPath}/seacr/${clusterName}_SEACRpeaks_relaxedPeaks.bed
done
