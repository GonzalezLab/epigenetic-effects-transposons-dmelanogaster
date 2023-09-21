#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=16G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=expressionAnalysis

#SBATCH --output=expressionAnalysis.output

#SBATCH --error=expressionAnalysis.error

module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

histones="H3K9me3 H3K27ac bivalent"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment depletion"

#histones="H3K9me3"
#tissues="ovary"
#effects="increment"
#assemblies="TOM-007"

wDir="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/expression_analysis/"

cd $wDir

for assembly in $assemblies
do
for tissue in $tissues
do
for histone in $histones
do
for effect in $effects
do
mkdir -p zscoreTMM/$assembly/$tissue/$histone/$effect

> zscoreTMM/${assembly}/${tissue}/$histone/$effect/zscore_${histone}-${tissue}-${effect}-${assembly}.tab
while read TEgene
do
TE=$(echo "$TEgene" | cut -f1)
gene=$(echo "$TEgene" |cut -f2)
position=$(echo "$TEgene" | cut -f3 )
distance=$(echo "$TEgene" | cut -f4 )
mkdir -p tmp/${assembly}/$tissue/$histone/$effect
> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.FPKM.matrix
> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TPM.matrix
> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TMM.matrix

strainTEpresent=$assembly
strainTEabsent=$(awk -v TE="$TE" ' $1 == TE ' ${wDir}/listEffects/${assembly}/${tissue}/list_TEs_${effect}_${histone}.lst | cut -f2 )

awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${assembly}/FPKM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "TE"}' >> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.FPKM.matrix
awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${assembly}/TPM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "TE"}' >> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TPM.matrix
awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${assembly}/TMM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "TE"}' >> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TMM.matrix


while read strainNoTE
do
awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${strainNoTE}/FPKM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "noTE"}' >> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.FPKM.matrix
awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${strainNoTE}/TPM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "noTE"}' >> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TPM.matrix
awk -v gene="$gene" ' $1 == gene ' ${wDir}/RNA-seq/${strainNoTE}/TMM_matrix_def.tab | awk -v tissue="$tissue" ' $5 == tissue ' | awk '{print $0 "\t" "noTE"}' >> tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TMM.matrix
done < <(echo "$strainTEabsent")


# R

Rscript zscore.R tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TPM.matrix tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.FPKM.matrix tmp/${assembly}/$tissue/$histone/$effect/$TE.$gene.TMM.matrix $assembly $tissue $histone $effect $TE $gene $position $distance

done < ${wDir}/intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}_effect_position.lst
done
done
done
done


> zscore_all_pos.tab
for histone in $histones
do
for effect in $effects
do
cat zscoreTMM/*/*/${histone}/${effect}/* > zscore_${histone}_${effect}.tab
cat zscoreTMM/*/*/${histone}/${effect}/* >> zscore_all_pos.tab
done
done

#cat zscoreTMM/*/*/*/*/* > zscore_all.tab
