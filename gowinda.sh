#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=40

#SBATCH --mem=8G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=get_breakpoints

#SBATCH --output=logs/gowinda.output

#SBATCH --error=logs/gowinda.error


module load Miniconda3/4.9.2
source activate epigenetics

DIR=/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/R01-GB-GOWINDA



tissues=("head" "gut" "ovary")
effects=("All" "Positive" "Negative")

for t in "${tissues[@]}"; do
for z in "${effects[@]}"; do

java -jar $DIR//Gowinda-1.12.jar \
--snp-file $DIR/TE_reference_genome/all_insertions.bed \
--candidate-snp-file "INPUT/${t}.${z}.TE_genes_position.bed" \
--gene-set-file $DIR//funcassociate_go_associations.txt \
--annotation-file $DIR/dmel-filter.gtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition gene \
--threads 8 \
--output-file $DIR/OUTPUT/${t}.${z}.TE_genes_position.gene.txt \
--mode gene \
--min-genes 2 &

java -jar $DIR//Gowinda-1.12.jar \
--snp-file $DIR/TE_reference_genome/all_insertions.bed \
--candidate-snp-file "INPUT/${t}.${z}.TE_genes_position.bed" \
--gene-set-file $DIR//funcassociate_go_associations.txt \
--annotation-file $DIR/dmel-filter.gtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition updownstream1000 \
--threads 8 \
--output-file $DIR/OUTPUT/${t}.${z}.TE_genes_position.updownstream1000.txt \
--mode gene \
--min-genes 2 &

java -jar $DIR//Gowinda-1.12.jar \
--snp-file $DIR/TE_reference_genome/all_insertions.bed \
--candidate-snp-file "INPUT/${t}.${z}.TE_genes_position.bed" \
--gene-set-file $DIR//funcassociate_go_associations.txt \
--annotation-file $DIR/dmel-filter.gtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition updownstream2000 \
--threads 8 \
--output-file $DIR/OUTPUT/${t}.${z}.TE_genes_position.updownstream2000.txt \
--mode gene \
--min-genes 2 &

java -jar $DIR//Gowinda-1.12.jar \
--snp-file $DIR/TE_reference_genome/all_insertions.bed \
--candidate-snp-file "INPUT/${t}.${z}.TE_genes_position.bed" \
--gene-set-file $DIR//funcassociate_go_associations.txt \
--annotation-file $DIR/dmel-filter.gtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition updownstream5000 \
--threads 8 \
--output-file $DIR/OUTPUT/${t}.${z}.TE_genes_position.updownstream5000.txt \
--mode gene \
--min-genes 2 &

java -jar $DIR//Gowinda-1.12.jar \
--snp-file $DIR/TE_reference_genome/all_insertions.bed \
--candidate-snp-file "INPUT/${t}.${z}.TE_genes_position.bed" \
--gene-set-file $DIR//funcassociate_go_associations.txt \
--annotation-file $DIR/dmel-filter.gtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition updownstream20000 \
--threads 8 \
--output-file $DIR/OUTPUT/${t}.${z}.TE_genes_position.updownstream20000.txt \
--mode gene \
--min-genes 2 &


wait

done
done

analysis="gene updownstream1000 updownstream2000 updownstream5000 updownstream20000"

for t in "${tissues[@]}"; do
for z in "${effects[@]}"; do
for a in $analysis; do
awk -F'\t' '$5 < 0.05 { print "'$t'\t'$z'\t'$a'\t"$1"\t"$5"\t"$9"\t"$10 }' "$DIR/OUTPUT/${t}.${z}.TE_genes_position.${a}.txt" 
done
done
done > processed_gowinda.fdr0.05.txt


analysis="updownstream20000"

for t in "${tissues[@]}"; do
for z in "${effects[@]}"; do
for a in $analysis; do
awk -F'\t' ' { print "'$t'\t'$z'\t"$1"\t"$5"\t"$9"\t"$10 }' "$DIR/OUTPUT/${t}.${z}.TE_genes_position.${a}.txt" 
done
done
done > processed_gowinda.$analysis.txt
