#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=20

#SBATCH --mem=32G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=enrichment

#SBATCH --output=logs/enrichment_H3K27ac_%a.output

#SBATCH --error=logs/enrichment_H3K27ac_%a.error

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
#histone=$(echo "$random_region_n" | cut -f5) # manually set H3K27ac 
effect=$(echo "$random_region_n" | cut -f6)

> $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.H3K27ac.txt

MAX_PARALLEL=20
COUNT=0

#cont=$(grep -v $'\t'  $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.txt | awk 'NF' | tail -n 1 | cut -f1 -d' ' | grep -f - -n  $DIR/DATA/TEs_epigenetics_H3K9me3_${effect}_${strain}_${tissue}.location_${position2}.random.final.clean.tab  | cut -f1 -d':')
#echo "${strain} ${tissue} ${position2} $cont"  > $DIR/intersectTesHistoneMarks/log_pos_H3K9me3_${effect}_${strain}_${tissue}.location_${position2}.txt

cont=$(cat $DIR/intersectTesHistoneMarks/log_pos_H3K9me3_${effect}_${strain}_${tissue}.location_${position2}.txt | cut -f 4 -d' ')


cat $DIR/DATA/TEs_epigenetics_H3K9me3_${effect}_${strain}_${tissue}.location_${position2}.random.final.clean.tab | head -n ${cont} | while read TEbed; do
(
  chr=$(echo "${TEbed}" | cut -f7)
  start=$(echo "${TEbed}" | cut -f8)
  end=$(echo "${TEbed}" | cut -f9)
  TE=$(echo "${TEbed}" | cut -f1)
  tissue=$(echo "${TEbed}" | cut -f3)
  gene=$(echo "${TEbed}" | cut -f10)
  position=$(echo "${TEbed}" | cut -f6)
  zsc=$(echo "${TEbed}" | cut -f11)
  pval=$(echo "${TEbed}" | cut -f12)
  log=$(echo "${TEbed}" | cut -f13)
  strainTE=$(echo "${TEbed}" | cut -f2)
  strainsBreakpoint=$(grep -v "$strainTE" $DIR/genomeAssemblies/genomeAssembly.txt)
    #for histone in H3K9me3 H3K27ac; do

  for histone in H3K27ac; do
    mkdir -p ${DIR}/intersectTesHistoneMarks/tmp/${TE}.${chr}.${start}.${end}/${tissue}/${histone}/
    zcat ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${tissue}_location_${position2}_intersectTEsWindowsSignalValueTrack.sort.bed.gz \
    | awk -v TE="$TE" -v start="$start" -v end="$end" -v chr="$chr" '$4 == TE && $8 == chr && $9 == start && $10 == end' \
    > ${DIR}/intersectTesHistoneMarks/tmp/${TE}.${chr}.${start}.${end}/${tissue}/${histone}/TE_${TE}_${tissue}_${strainTE}_${histone}.bed

    while read strainBreakpoint; do
    (
      zcat ${DIR}/intersectTesHistoneMarks/${histone}_${effect}_${strain}_${strainBreakpoint}_${tissue}_location_${position2}_intersectBreakpointsWindowsSignalValueTrack.sort.bed.gz \
      | awk -v TE="$TE" -v start="$start" -v end="$end" -v chr="$chr" '$4 == TE && $8 == chr && $9 == start && $10 == end' \
      > ${DIR}/intersectTesHistoneMarks/tmp/${TE}.${chr}.${start}.${end}/${tissue}/${histone}/breakpoint_${TE}_${tissue}_${strainBreakpoint}_${histone}.bed

      Rscript $DIR/enrichment3.R "$TE" "$chr" "$start" "$end" "$tissue" "$histone" "$strainTE" "$strainBreakpoint" "$gene" "$position" "$zsc" "$pval" "$log"
    ) &
    done < <(echo "$strainsBreakpoint")

    wait  # wait for all strainBreakpoint jobs for this histone
  done

  # Cleanup after both histones processed
  rm -rf ${DIR}/intersectTesHistoneMarks/tmp/${TE}.${chr}.${start}.${end}

) &

((COUNT++))
if [[ $COUNT -ge $MAX_PARALLEL ]]; then
  wait
  COUNT=0
fi

done
wait  # final wait for all outer jobs

# cat $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.txt | tr ' ' '\t' > $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab

# IMPORTANT: FIX THE \N PROBLEM WITH SOME LINES -- now it does a cat with \t and \n so I htink lines hould be ok. i should only run a sort

# do the sort

#awk 'NF' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.txt | sort -k1.8n -k6 > $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab
