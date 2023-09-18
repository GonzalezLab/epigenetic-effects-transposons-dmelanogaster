#!/bin/bash

module load Miniconda3/4.7.10
source activate /homes/users/mcoronado/.conda/envs/5GP

assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
tissues="head gut ovary"
histones="H3K9me3 H3K27ac"
inDIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics/DATA/"
outDIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics/ANALYSES/"

# Keep all TEs analyzable (4823)
for assembly in $assemblies
do
mkdir -p ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}

while read TEbed
do
chr=$(echo "$TEbed" | cut -f1 )
TE=$(echo "$TEbed" | cut -f4 )
length=$(echo "$TEbed" | cut -f5 )
strand=$(echo "$TEbed" | cut -f6 )
family=$(echo "$TEbed" | cut -f7 )
startTE=$(echo "$TEbed" | cut -f 2)
start=$(($startTE - 20000))
endTE=$(echo "$TEbed" | cut -f 3)
end=$(echo "$TEbed" | cut -f 3)

for s in `seq 20 -1 1`
do
echo -ne "$chr\t$start"
start=$(($start + 1000))
echo -e "\t$start\t$TE\t1000\t$strand\t$family\t-${s}kb\tTE"
done

echo -e "$chr\t$startTE\t$endTE\t$TE\t$length\t$strand\t$family\tTE\tTE"

for s in `seq 1 20`
do
echo -ne "$chr\t$end"
end=$(($end + 1000))
echo -e "\t$end\t$TE\t1000\t$strand\t$family\t${s}kb\tTE"
done

done < ${outDIR}/TEsDescriptiveAnalyses/TE-Clean-Nested/${assembly}.bed  > ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindows.bed

done

#cd slidingWindowsBed
#cut -f4 */*_nonOverlappingWindows.bed | sort -u |wc -l

# Only keep the ones with breakpoints + fixed (2425+945=3370 in total)
for assembly in $assemblies
do
grep -f <(cat ../breakpointsTEs/TEsTSD50bp.lst |sed  "s/$/\t/") ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindows.bed > ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindowsFilter.tmp.bed
grep -f <(grep "5$" ${outDIR}/TEsDescriptiveAnalyses/TE-Frequency/allStrains.tsv | cut -f5 |sort -u | sed  "s/$/\t/") ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindows.bed  > ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindowsFilter.tmp2.bed
cat  ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindowsFilter.tmp.bed ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindowsFilter.tmp2.bed > ${outDIR}/intersectTesHistoneMarks/slidingWindowsBed/${assembly}/${assembly}_nonOverlappingWindowsFilter.bed
done

#cd slidingWindowsBed
#cut -f4 */*_nonOverlappingWindowsFilter.bed | sort -u |wc -l