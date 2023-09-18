#!/bin/bash

module load Miniconda3/4.7.10
source activate /homes/users/mcoronado/.conda/envs/5GP

assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
tissues="head gut ovary"
histones="H3K9me3 H3K27ac"
inDIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics/DATA/"
outDIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics/ANALYSES/"

for strain in $assemblies
do
> ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed
done

while read TEinfo
do
TE=$(echo "$TEinfo" | cut -f1 -d' ' )
strain=$(echo "$TEinfo" | cut -f3 -d' ' )
chr=$(cat $outDIR/breakpointsTEs/breakpoint/$TE/$strain/$TE.breakpoint.$strain.bed  | cut -f1 )
breakpoint=$(cat $outDIR/breakpointsTEs/breakpoint/$TE/$strain/$TE.breakpoint.$strain.bed  | cut -f2 )
breakpointStart=$(( $breakpoint - 20000 ))
strand=$(cat $outDIR/breakpointsTEs/breakpoint/$TE/$strain/$TE.breakpoint.$strain.bed  | cut -f5 )
for s in `seq 20 -1 1`
do
echo -ne "$chr\t$breakpointStart"  >> ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed
breakpointStart=$(($breakpointStart + 1000))
echo -e "\t$breakpointStart\t${TE}\t1000\t$strand\t-${s}kb\tbreakpoint" >> ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed
done
breakpointS=$(( $breakpoint - 1 ))
echo -e "$chr\t$breakpointS\t$breakpoint\t${TE}\t1\t$strand\tbreakpoint\tbreakpoint" >> ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed

for s in `seq 1 20`
do
echo -ne "$chr\t$breakpoint" >> ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed
breakpoint=$(($breakpoint + 1000))
echo -e "\t$breakpoint\t${TE}\t1000\t$strand\t${s}kb\tbreakpoint" >> ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed
done

done < ${outDIR}/breakpointsTEs/statusBreakpointAnalyzable50bp.log

#cd slidingWindowsBed
#cut -f4 */*_nonOverlappingWindowsBreakpoints.bed | sort -u |wc -l