#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=18G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=heatmap

#SBATCH --output=logs/heatmap_%a.output

#SBATCH --error=logs/heatmap_%a.error

module load Miniconda3/4.9.2
source activate epigenetics

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/R01-heatmap"
cd $DIR
histones="H3K27ac H3K9me3"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

for assembly in $assemblies
do
#sed -i -r 's/\S+\t//2' $DIR/DATA/chrSizes/$assembly.chrom.sizes
for tissue in $tissues
do
for histone in $histones
do
mkdir -p $DIR/RESULT/$assembly/$tissue/$histone
call-macs2_signal_track/$assembly/$tissue/$histone/${assembly}_${histone}_${tissue}.bedGraph 
cut -f1-3,5 $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean.bed > $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/${assembly}_${histone}_${tissue}.bedGraph
./bedGraphToBigWig $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/${assembly}_${histone}_${tissue}.bedGraph $DIR/DATA/chrSizes/$assembly.chrom.sizes $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/${assembly}_${histone}_${tissue}.bw
#gzip $DIR/DATA/call-macs2_signal_track/$assembly/$tissue/$histone/allReps.clean_permseq.filter_sorted.nodup_x_ctl.pooled.fc.signal.mean.sort.mean.bed
#gzip $DIR/DATA/
computeMatrix reference-point -S $DIR/DATA/call-macs2_signal_track/$assembly/${tissue}/${histone}/${assembly}_${histone}_$tissue.bw -R $DIR/DATA/TE-Clean-Nested/$assembly.bed -b 20000 -a 20000 --outFileName $DIR/RESULT/$assembly/$tissue/$histone/matrix.TE.dat --outFileNameMatrix $DIR/RESULT/$assembly/$tissue/$histone/matrix.TE.txt --referencePoint center --numberOfProcessors 4
plotHeatmap --matrixFile $DIR/RESULT/$assembly/$tissue/$histone/matrix.TE.dat --outFileName $DIR/RESULT/$assembly/$tissue/$histone/TE.pdf --sortRegions descend --sortUsing median --plotFileFormat pdf --xAxisLabel TE  --regionsLabel "TEs" --samplesLabel "$histone $tissue"

done
done

computeMatrix reference-point -S $DIR/DATA/call-macs2_signal_track/$assembly/head/H3K9me3/${assembly}_H3K9me3_head.bw $DIR/DATA/call-macs2_signal_track/$assembly/head/H3K27ac/${assembly}_H3K27ac_head.bw $DIR/DATA/call-macs2_signal_track/$assembly/gut/H3K9me3/${assembly}_H3K9me3_gut.bw $DIR/DATA/call-macs2_signal_track/$assembly/gut/H3K27ac/${assembly}_H3K27ac_gut.bw $DIR/DATA/call-macs2_signal_track/$assembly/ovary/H3K9me3/${assembly}_H3K9me3_ovary.bw $DIR/DATA/call-macs2_signal_track/$assembly/ovary/H3K27ac/${assembly}_H3K27ac_ovary.bw -R $DIR/DATA/TE-Clean-Nested/$assembly.bed -b 20000 -a 20000 --outFileName $DIR/RESULT/$assembly/matrix.TE.dat --outFileNameMatrix $DIR/RESULT/$assembly/matrix.TE.txt --referencePoint center --numberOfProcessors 4

plotHeatmap --matrixFile $DIR/RESULT/$assembly/matrix.TE.dat --outFileName $DIR/RESULT/$assembly/TE.pdf --sortRegions descend --sortUsing median --plotFileFormat pdf --xAxisLabel TE  --regionsLabel "TEs" --samplesLabel "H3K9me3 head" "H3K27ac head" "H3K9me3 gut" "H3K27ac gut" "H3K9me3 ovary" "H3K27ac ovary"

done

