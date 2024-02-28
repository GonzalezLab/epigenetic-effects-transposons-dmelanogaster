#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=18G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=plot

#SBATCH --output=logs/plot_%a.output

#SBATCH --error=logs/plot_%a.error

module load Miniconda3/4.9.2
source activate epigenetics

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/R01-enrichment"
cd $DIR
histones="H3K27ac H3K9me3"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

for assembly in $assemblies
do
#bedtools makewindows -b $DIR/DATA/chrSizes/$assembly.chrom.sizes -w 10000 > $DIR/DATA/$assembly.chrom.window.10kb.bed
for tissue in $tissues
do
for histone in $histones
do

gunzip $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.bed.gz

bedmap --echo  --mean $DIR/DATA/$assembly.chrom.window.10kb.bed $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.bed  | tr '|' '\t' | awk -v OFS='\t' '$3 = $3 OFS "."' | sed 's/NAN/0.000000/' > $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean10kb.bed
#gzip $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.bed
#gzip $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean10kb.bed
done
done
done