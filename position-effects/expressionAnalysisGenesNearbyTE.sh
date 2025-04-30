#!/bin/bash

# set the partition where the job will run


# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=4G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=expressionAnalysis_%j

#SBATCH --output=logs/expressionAnalysis_%a.output

#SBATCH --error=logs/expressionAnalysis_%a.error

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
#  cut -f9 */*_gene* | cut -f 2 -d' ' | tr -d '"' | tr -d ';' | sort | uniq -c |  awk '$1 == 5 {print $2}'>total_genes.lst
task_file="${wDir}/tasks.txt"

# rm -f "$task_file"  # Remove existing task file if any
# for assembly in $assemblies; do
#     for tissue in $tissues; do
#         echo "$assembly $tissue" >> "$task_file"
#     done
# done

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$task_file")
read assembly tissue <<< "$line"

cd $wDir

# for assembly in $assemblies
# do
# for tissue in $tissues
# do
#for histone in $histones
#do
#for effect in $effects
#do
#mkdir -p $wDir/zscoreTMM_geneTE/$assembly/$tissue/geneTE1kb
mkdir -p $wDir/zscoreTMM_geneTE/$assembly/$tissue/geneTE2kb

mkdir -p tmp2/assembly/tissue

# > zscoreTMM_geneTE/${assembly}/${tissue}//geneTE1kb/zscore_${tissue}-${assembly}.tab
> zscoreTMM_geneTE/${assembly}/${tissue}//geneTE2kb/zscore_${tissue}-${assembly}.tab

while read gene
do
#TE=$(echo "$gene" | cut -f1)
#gene=$(echo "$gene" |cut -f2)
position=NA
distance=NA
mkdir -p tmp2/${assembly}/$tissue/
#> tmp/${assembly}/$tissue//$gene.FPKM.matrix
#> tmp/${assembly}/$tissue//$gene.TPM.matrix
> tmp2/${assembly}/$tissue//$gene.TMM.matrix

strainTEpresent=$assembly
strainTEabsent=$(echo "$assemblies" | tr ' ' '\n' | grep -v "^$assembly$" | tr '\n' ' ' )

#awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${assembly}/FPKM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "TE"}' >> tmp/${assembly}/$tissue//$gene.FPKM.matrix
#awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${assembly}/TPM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "TE"}' >> tmp/${assembly}/$tissue//$gene.TPM.matrix
awk -v gene="$gene" ' $1 == gene ' ../expression_analysis/RNA-seq/${assembly}/TMM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "TE"}' >> tmp2/${assembly}/$tissue//$gene.TMM.matrix


for strainNoTE in $strainTEabsent
do
#awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${strainNoTE}/FPKM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "noTE"}' >> tmp/${assembly}/$tissue//$gene.FPKM.matrix
#awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${strainNoTE}/TPM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "noTE"}' >> tmp/${assembly}/$tissue//$gene.TPM.matrix
awk -v gene="$gene" ' $1 == gene ' ../expression_analysis/RNA-seq/${strainNoTE}/TMM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "noTE"}' >> tmp2/${assembly}/$tissue//$gene.TMM.matrix
done #< <(echo "$strainTEabsent")

# R

Rscript zscore_geneTE2.R tmp2/${assembly}/$tissue//$gene.TMM.matrix $assembly $tissue $gene $position $distance

#done < <(cut -f2 $wDir/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${assembly}_filtered.bed | sort -u)
done < <(cut -f2 $wDir/INPUT/intersect/intersect2kb_genes_TE_${tissue}_${assembly}_filtered.bed | sort -u)


# > zscore_all_pos.tab
# for histone in $histones
# do
# for effect in $effects
# do
# cat zscoreTMM/*/*/${histone}/${effect}/* > zscore_${histone}_${effect}.tab
# cat zscoreTMM/*/*/${histone}/${effect}/* >> zscore_all_pos.tab
# done
# done

# #cat zscoreTMM/*/*/*/*/* > zscore_all.tab
