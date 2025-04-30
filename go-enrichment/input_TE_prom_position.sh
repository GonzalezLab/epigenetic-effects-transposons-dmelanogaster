#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=2G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=get_breakpoints

#SBATCH --output=logs/input_prom_position.output

#SBATCH --error=logs/input_prom_position.error


module load Miniconda3/4.9.2
source activate epigenetics

DIR=/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/R01-GB-GOWINDA


tissues=("head" "gut" "ovary")
effects=("All" "Positive" "Negative")

for t in "${tissues[@]}"; do
for z in "${effects[@]}"; do
while read -r gene_TE; do
TE=$(echo "$gene_TE" | cut -f1 -d' ')
strain=$(echo "$gene_TE" | cut -f2 -d' ')

chr=$(awk -v TE="$TE" '$4 == TE' TE_reference_genome/${strain}_Ref_Coord.bed | cut -f 1)
start=$(awk -v TE="$TE" '$4 == TE' TE_reference_genome/${strain}_Ref_Coord.bed | cut -f 2)
end=$(awk -v TE="$TE" '$4 == TE' TE_reference_genome/${strain}_Ref_Coord.bed | cut -f3)
mid=$(( (start + end) / 2 ))

echo -e "$chr\t$mid\t$mid\t$TE\t$strain" 
done < "INPUT/${t}.${z}.TE_genes.lst" > "INPUT/${t}.${z}.TE_genes_position.bed"
done
done

cat INPUT/*.TE_genes.lst | cut -f3 -d' ' | sort -u > INPUT/all_genes_DE.lst

grep -f INPUT/all_genes_DE.lst $DIR//dmel-all-r6.31.gtf > $DIR/dmel-filter.gtf

cat TE_reference_genome/*_Ref_Coord.bed | cut -f1,2,3,4 | sort -u > TE_reference_genome/all_insertions.tmp.bed 

while read -r info; do
TE=$(echo "$info" | cut -f4)
start=$(echo "$info" | cut -f2)
end=$(echo "$info" | cut -f3)
chr=$(echo "$info" | cut -f1)

mid=$(( (start + end) / 2 ))

echo -e "$chr\t$mid\t$mid\t$TE" 

done < TE_reference_genome/all_insertions.tmp.bed  > TE_reference_genome/all_insertions.bed 
