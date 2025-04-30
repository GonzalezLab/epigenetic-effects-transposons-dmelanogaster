#!/bin/bash

# set the partition where the job will run


# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=16G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=intersect_%j

#SBATCH --output=logs/intersect_%a.output

#SBATCH --error=logs/intersect_%a.error

module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

#histones="H3K9me3 H3K27ac"
# tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
#effects="increment"

#histones="H3K9me3"
#tissues="ovary"
#effects="increment"
#assemblies="TOM-007"

wDir="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/R01-GB-expression_analysis/"
cd $wDir

task_file="${wDir}/tasks.txt"

# rm -f "$task_file"  # Remove existing task file if any
# for assembly in $assemblies; do
#     for tissue in $tissues; do
#         echo "$assembly $tissue" >> "$task_file"
#     done
# done

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$task_file")
read assembly tissue <<< "$line"


# > $wDir/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${assembly}.bed

# bedtools window -w "1000" -a $wDir/INPUT/${assembly}_TE_unique.bed -b $wDir/INPUT/geneAnnotationsLiftoff/${assembly}/${assembly}_genes.gtf |  cut -f 4,18 | cut -f -1,2 -d' ' | sed "s/gene_id \"//g" | sed "s/\";//g"  >> $wDir/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${assembly}.bed

# cat $wDir/INPUT/geneAnnotationsLiftoff/total_genes.lst | grep -f - $wDir/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${assembly}.bed > $wDir/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${assembly}_filtered.bed


> $wDir/INPUT/intersect/intersect2kb_genes_TE_${tissue}_${assembly}.bed

bedtools window -w "2000" -a $wDir/INPUT/${assembly}_TE_unique.bed -b $wDir/INPUT/geneAnnotationsLiftoff/${assembly}/${assembly}_genes.gtf |  cut -f 4,18 | cut -f -1,2 -d' ' | sed "s/gene_id \"//g" | sed "s/\";//g"  >> $wDir/INPUT/intersect/intersect2kb_genes_TE_${tissue}_${assembly}.bed

cat $wDir/INPUT/geneAnnotationsLiftoff/total_genes.lst | grep -f - $wDir/INPUT/intersect/intersect2kb_genes_TE_${tissue}_${assembly}.bed > $wDir/INPUT/intersect/intersect2kb_genes_TE_${tissue}_${assembly}_filtered.bed


# only TEs in 1