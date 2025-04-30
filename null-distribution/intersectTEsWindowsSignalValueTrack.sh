#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=64G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=intersect_windows_TE

#SBATCH --output=logs/intersect_windows_TE_%a.output

#SBATCH --error=logs/intersect_windows_TE_%a.error


module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

histones="H3K9me3 H3K27ac"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment"

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/null_distribution/R01-GB-null_distribution_v2/"
DATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/"

random_region_n=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $DIR/number_random_regions.tab)
tissue=$(echo "$random_region_n" | cut -f1)
strain=$(echo "$random_region_n" | cut -f2)
position2=$(echo "$random_region_n" | cut -f3)
n=$(echo "$random_region_n" | cut -f4)
histone=$(echo "$random_region_n" | cut -f5) # manually set H3K27ac 
effect=$(echo "$random_region_n" | cut -f6)

> ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed 


mkdir -p ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/
bedtools intersect -a ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindows.bed -b $DATA/ChIPseq/${strain}/${tissue}/${histone}/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean.bed.gz -wa -wb > ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed 
sort -k1,1 -k2,2n ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed  > ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/${histone}_${effect}_${strain}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed
awk  -v tissue="$tissue" -v strain="$strain" -v histone="$histone" '{print $0 "\t" tissue "\t" strain "\t" histone }' ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/${histone}_${effect}_${strain}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed >> ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed 


sort -k1,1 -k2,2n ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed > ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed

gzip -f ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed

rm -rf ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed
rm -rf ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/${histone}_${effect}_${strain}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed
rm -rf ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone}/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed

# version H3K27ac
histone2=H3K27ac

> ${DIR}/intersectTesHistoneMarks/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed 


mkdir -p ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/
bedtools intersect -a ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindows.bed -b $DATA/ChIPseq/${strain}/${tissue}/${histone2}/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean.bed.gz -wa -wb > ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed 
sort -k1,1 -k2,2n ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed  > ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/${histone2}_${effect}_${strain}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed
awk  -v tissue="$tissue" -v strain="$strain" -v histone2="$histone2" '{print $0 "\t" tissue "\t" strain "\t" histone2 }' ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/${histone2}_${effect}_${strain}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed >> ${DIR}/intersectTesHistoneMarks/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed 


sort -k1,1 -k2,2n ${DIR}/intersectTesHistoneMarks/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed > ${DIR}/intersectTesHistoneMarks/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed

gzip -f ${DIR}/intersectTesHistoneMarks/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed


rm -rf ${DIR}/intersectTesHistoneMarks/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed
rm -rf ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/${histone2}_${effect}_${strain}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed
rm -rf ${DIR}/intersectTesHistoneMarks/intersectSignalValueTrack/TE/${strain}/${tissue}/${histone2}/${histone2}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.bed
