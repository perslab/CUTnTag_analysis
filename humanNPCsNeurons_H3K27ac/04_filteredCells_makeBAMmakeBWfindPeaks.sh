#!/bin/bash

#04 Filtered cells: Sinto make bam, deeptools make bigwig, find peaks
outsPath="/projects/oliver/SCOaaP/outs/SCOP_2021_0139/scATACseq/03_PipelineOut/210218_kirkeby_cutNtag/cellranger/EndoPool_3/outs"
hashPath="/projects/oliver/SCOaaP/outs/SCOP_2021_0139/scATACseq/03_PipelineOut/210218_kirkeby_cutNtag/kite/HashPool_3/featurecounts"
projPath="/nfsdata/projects/jph712/cutntag_analysis/humanNPCsNeurons_H3K27ac"
ExptName="humanNPCsNeurons_H3K27ac"

nproc=30

##mouse
#mappable_genome=2652783500 #mappable genome size: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
#blacklist="/nfsdata/projects/jph712/cutntag_analysis/mm10/mm10-blacklist.v2.bed"

##human:
mappable_genome=2805636331 #mappable genome size: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
blacklist="/nfsdata/projects/jph712/cutntag_analysis/hg38/hg38-blacklist.v2.bed"





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

cd $projPath/bam; sinto filterbarcodes -b ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam -p $nproc -c $projPath/bam/FilteredCellBarcodes.tsv ; cd $projPath
samtools index $projPath/bam/${ExptName}_FilteredCells.bam

#######################
#04b deeptools bigwig
#######################
mkdir $projPath/bigwig -p

#all cellranger fragments:
bamCoverage --bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam --outFileName ${projPath}/bigwig/${ExptName}_AllCellrangerFragments_nucFree.bigWig --outFileFormat bigwig --numberOfProcessors $nproc   --minFragmentLength 20 --maxFragmentLength 121 --blackListFileName $blacklist -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $mappable_genome
bamCoverage --bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam --outFileName ${projPath}/bigwig/${ExptName}_AllCellrangerFragments_Nucleosomal.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 120 --maxFragmentLength 800 --blackListFileName $blacklist -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $mappable_genome
bamCoverage --bam ${projPath}/bam/${ExptName}_AllCellrangerFragments.bam --outFileName ${projPath}/bigwig/${ExptName}_AllCellrangerFragments_extendReads.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 20 --maxFragmentLength 800 --blackListFileName $blacklist -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $mappable_genome

#Filtered cells:
bamCoverage --bam $projPath/bam/${ExptName}_FilteredCells.bam --outFileName ${projPath}/bigwig/${ExptName}_filteredCells_nucFree.bigWig --outFileFormat bigwig --numberOfProcessors $nproc   --minFragmentLength 20 --maxFragmentLength 121 --blackListFileName $blacklist -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $mappable_genome
bamCoverage --bam $projPath/bam/${ExptName}_FilteredCells.bam --outFileName ${projPath}/bigwig/${ExptName}_filteredCells_Nucleosomal.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 120 --maxFragmentLength 800 --blackListFileName $blacklist -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $mappable_genome
bamCoverage --bam $projPath/bam/${ExptName}_FilteredCells.bam --outFileName ${projPath}/bigwig/${ExptName}_filteredCells_extendReads.bigWig --outFileFormat bigwig --numberOfProcessors $nproc  --minFragmentLength 20 --maxFragmentLength 800 --blackListFileName $blacklist -bs 1 --ignoreDuplicates --ignoreForNormalization chrX chrM --extendReads --normalizeUsing RPGC --effectiveGenomeSize $mappable_genome



#######################
#04c MACS2
#######################
mkdir $projPath/macs2/BinClusters/narrow -p
mkdir $projPath/macs2/BinClusters/broad -p

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
SEACR_1.3.sh ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.1 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.1peaks.bed
SEACR_1.3.sh ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.01 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.01peaks.bed
SEACR_1.3.sh ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.001 norm stringent ${projPath}/seacr/${clusterName}_SEACRpeaks_stringentTop0.001peaks.bed
SEACR_1.3.sh ${projPath}/seacr/bed/${clusterName}.fragments.bedgraph 0.001 norm relaxed ${projPath}/seacr/${clusterName}_SEACRpeaks_relaxedPeaks.bed
done
