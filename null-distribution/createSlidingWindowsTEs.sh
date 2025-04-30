#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=8G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=createSlidingWindowsTEs

#SBATCH --output=logs/createSlidingWindowsTEs_%a.output

#SBATCH --error=logs/createSlidingWindowsTEs_%a.error


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

mkdir -p ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}

grep -P "\t1$"  $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.breakpointDetected.tab  >  $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.final.tab

absentStrains=$(grep -v $strain $DIR/genomeAssemblies/genomeAssembly.txt)

for absentStrain in $absentStrains
do
> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.problem.txt
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
chr2=$(cat $DIR/breakpoint/TEs_epigenetics_${histone}_${effect}_${absentStrain}_${tissue}.location_${position2}.random.tab | grep -P "\t${TE}:${chr}.${start}.${end}::" | cut -f1)
if  [[ -z "$chr2" ]]; then
echo -e "${TE}\t${chr}\t${start}\t${end}" >> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.problem.txt
fi
) &
done <<< "$absentStrains"
wait

done < $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.final.tab

awk 'FNR==NR {key[$1"\t"$2"\t"$3"\t"$4]; next} !($1"\t"$7"\t"$8"\t"$9 in key)' \
"$DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.problem.txt" \
"$DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.final.tab" > $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.final.clean.tab


while read TEbed
do
chr=$(echo "$TEbed" | cut -f7 )
TE=$(echo "$TEbed" | cut -f1 )
startTE=$(echo "$TEbed" | cut -f 8)
start=$(($startTE - 20000))
endTE=$(echo "$TEbed" | cut -f 9)
end=$(echo "$TEbed" | cut -f 9)
position=$(echo "$TEbed" | cut -f6,10-13)

for s in `seq 20 -1 1`
do
echo -ne "$chr\t$start"
start=$(($start + 1000))
echo -e "\t$start\t$TE\t1000\t-${s}kb\tTE\t$chr\t$startTE\t$endTE\t$position"
done

echo -e "$chr\t$startTE\t$endTE\t$TE\t1000\tTE\tTE\t$chr\t$startTE\t$endTE\t$position"

for s in `seq 1 20`
do
echo -ne "$chr\t$end"
end=$(($end + 1000))
echo -e "\t$end\t$TE\t1000\t${s}kb\tTE\t$chr\t$startTE\t$endTE\t$position"
done

done <  $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position2}.random.final.clean.tab   > ${DIR}/intersectTesHistoneMarks/slidingWindowsBed/${strain}/${strain}_${tissue}_${histone}_${effect}_location_${position2}_nonOverlappingWindows.bed



