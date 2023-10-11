#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=16G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=intersectTEsWindowsSV

#SBATCH --output=intersectTEsWindowsSV.output

#SBATCH --error=intersectTEsWindowsSV.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10
source activate /homes/users/mcoronado/.conda/envs/5GP

assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
tissues="head gut ovary"
histones="H3K9me3 H3K27ac"
inDIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics/DATA/"
outDIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics/ANALYSES/"

> ${outDIR}/intersectTesHistoneMarks/intersectTEsWindowsSignalValueTrack.bed 
for assembly in $assemblies
do
for tissue in $tissues
do
for histone in $histones
do
mkdir -p ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/
bedtools intersect -a ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindows.bed -b $inDIR/ChIPseq/${assembly}/${tissue}/${histone}/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean.bed -wa -wb > ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsWindowsSignalValueTrack.bed
sort -k1,1 -k2,2n ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsWindowsSignalValueTrack.bed > ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsWindowsSignalValueTrack.sort.bed
awk  -v tissue="$tissue" -v strain="$assembly" -v histone="$histone" '{print $0 "\t" tissue "\t" strain "\t" histone }' ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsWindowsSignalValueTrack.sort.bed >> ${outDIR}/intersectTesHistoneMarks/intersectTEsWindowsSignalValueTrack.bed 

done
done
done

sort -k1,1 -k2,2n ${outDIR}/intersectTesHistoneMarks/intersectTEsWindowsSignalValueTrack.bed > ${outDIR}/intersectTesHistoneMarks/intersectTEsWindowsSignalValueTrack.sort.bed


> ${outDIR}/intersectTesHistoneMarks/intersectTEsFilterWindowsSignalValueTrack.bed 
for assembly in $assemblies
do
for tissue in $tissues
do
for histone in $histones
do
mkdir -p ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/
bedtools intersect -a ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindowsFilter.bed -b $inDIR/ChIPseq/${assembly}/${tissue}/${histone}/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean.bed -wa -wb > ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsFilterWindowsSignalValueTrack.bed
sort -k1,1 -k2,2n ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsFilterWindowsSignalValueTrack.bed > ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsFilterWindowsSignalValueTrack.sort.bed
awk  -v tissue="$tissue" -v strain="$assembly" -v histone="$histone" '{print $0 "\t" tissue "\t" strain "\t" histone }' ${outDIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${assembly}/${tissue}/${histone}/intersectTEsFilterWindowsSignalValueTrack.sort.bed >> ${outDIR}/intersectTesHistoneMarks/intersectTEsFilterWindowsSignalValueTrack.bed 

done
done
done

sort -k1,1 -k2,2n ${outDIR}/intersectTesHistoneMarks/intersectTEsFilterWindowsSignalValueTrack.bed  > ${outDIR}/intersectTesHistoneMarks/intersectTEsFilterWindowsSignalValueTrack.sort.bed 

grep -f <(cat ../breakpointsTEs/TEsTSD50bp.lst |sed  "s/$/\t/") ${outDIR}/intersectTesHistoneMarks/intersectTEsFilterWindowsSignalValueTrack.sort.bed  > ${outDIR}/intersectTesHistoneMarks/intersectTEsFilterBreakpointsWindowsSignalValueTrack.sort.bed 