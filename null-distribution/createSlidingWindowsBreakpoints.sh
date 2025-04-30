#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=3G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=createSlidingWindowsBP

#SBATCH --output=logs/createSlidingWindowsBP_%a.output

#SBATCH --error=logs/createSlidingWindowsBP_%a.error

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
histone=$(echo "$random_region_n" | cut -f5)
effect=$(echo "$random_region_n" | cut -f6)

#> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_nonOverlappingWindowsBreakpoints.bed

absentStrains=$(grep -v $strain $DIR/genomeAssemblies/genomeAssembly.txt)

for absentStrain in $absentStrains
do
> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/${strain}_${absentStrain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindowsBreakpoints.bed

done


while read TEinfo
do
TE=$(echo "$TEinfo" | cut -f1 )
chr=$(echo "$TEinfo" | cut -f7 )
start=$(echo "$TEinfo" | cut -f8 )
end=$(echo "$TEinfo" | cut -f9  )
position=$(echo "$TEinfo" | cut -f6,10-13)

while read absentStrain
do
(
mkdir "${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/"


# this now needs to be a grep in a big file
#chr2=$(cut -f1 $DIR/breakpoint/${TE}.${chr}.${start}.${end}/${absentStrain}/${TE}.breakpoint.${absentStrain}.${chr}.${start}.${end}.bed)
chr2=$(cat $DIR/breakpoint/TEs_epigenetics_${histone}_${effect}_${absentStrain}_${tissue}.location_${position2}.random.tab | grep -P "\t${TE}:${chr}.${start}.${end}::" | cut -f1)
if  [[ -z "$chr2" ]]; then
echo "$TE $start $end not found" >> problem.txt
fi
#breakpoint=$(cut -f2 $DIR/breakpoint/${TE}.${chr}.${start}.${end}/${absentStrain}/${TE}.breakpoint.${absentStrain}.${chr}.${start}.${end}.bed)
breakpoint=$(cat $DIR/breakpoint/TEs_epigenetics_${histone}_${effect}_${absentStrain}_${tissue}.location_${position2}.random.tab | grep -P "\t${TE}:${chr}.${start}.${end}::" | cut -f2)

breakpointStart=$(( $breakpoint - 20000 ))
for s in `seq 20 -1 1`
do
echo -ne "$chr2\t$breakpointStart"  >> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/${strain}_${absentStrain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindowsBreakpoints.bed
breakpointStart=$(($breakpointStart + 1000))
echo -e "\t$breakpointStart\t${TE}\t1000\t-${s}kb\tbreakpoint\t$chr\t$start\t$end\t$position" >> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/${strain}_${absentStrain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindowsBreakpoints.bed
done
breakpointS=$(( $breakpoint - 1 ))
echo -e "$chr2\t$breakpointS\t$breakpoint\t${TE}\t1\tbreakpoint\tbreakpoint\t$chr\t$start\t$end\t$position" >> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/${strain}_${absentStrain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindowsBreakpoints.bed

for s in `seq 1 20`
do
echo -ne "$chr2\t$breakpoint" >> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/${strain}_${absentStrain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindowsBreakpoints.bed
breakpoint=$(($breakpoint + 1000))
echo -e "\t$breakpoint\t${TE}\t1000\t${s}kb\tbreakpoint\t$chr\t$start\t$end\t$position" >> ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${absentStrain}/${strain}_${absentStrain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindowsBreakpoints.bed
done
) &
done <<< "$absentStrains"
wait

done < $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.final.clean.tab
