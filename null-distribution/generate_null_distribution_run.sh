#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=8G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=expressionAnalysis

#SBATCH --output=logs/expressionAnalysis.output

#SBATCH --error=logs/expressionAnalysis.error

module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

histones="H3K9me3 H3K27ac bivalent"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment depletion"


DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/null_distribution/R01-GB-null_distribution_v2/"
DATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/"

# ALREADY DONE

# #### DEFINE HETEROCHROMATIC REGIONS

# samtools faidx $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta > $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta.fai
# cut -f1,2 $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta.fai | grep -Ew "2L|2R|3L|3R|X"> $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.sizes.genome

# bedtools complement -i $DIR/DATA/eucromatic_regions_reference.bed -g $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.sizes.genome > $DIR/DATA/heterochromatin/heterochromatin_regions_reference.bed

# bedtools getfasta -fi $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta \
# -bed $DIR/DATA/heterochromatin/heterochromatin_regions_reference.bed \
# -fo $DIR/DATA/heterochromatin/heterochromatin_regions_reference.fasta

# module load GCCcore/11.2.0
# module load minimap2/2.22

# for strain in $assemblies
# do

# minimap2 -ax asm5 \
# --secondary=no \
# ~/scratch/5GenomesProject/DATA/genomestrain/$strain/$strain.$pos.fasta \
# $DIR/DATA/heterochromatin/heterochromatin_regions_reference.fasta > $DIR/DATA/heterochromatin/heterochromMapping.$strain.sam

# minimap2 -x asm5 \
# --secondary=no \
# ~/scratch/5GenomesProject/DATA/genomestrain/$strain/$strain.$pos.fasta \
# $DIR/DATA/heterochromatin/heterochromatin_regions_reference.fasta > $DIR/DATA/heterochromatin/heterochromMapping.$strain.paf

# Rscript heterochromatin.R $strain

# For each strain manually extract the coordinates of the heterochromatin
# grep "^2L" heterochromMapping.$strain.paf | awk '$6 == "2L" ' | grep ":0-" | cut -f1-10 | sort -k8h
# grep "^2L" heterochromMapping.$strain.paf | awk '$6 == "2L" ' | grep -v ":0-" | cut -f1-10 | sort -k8h
# grep "^2R" heterochromMapping.$strain.paf | awk '$6 == "2R" ' | grep  ":0-" | cut -f1-10 | sort -k8h
# grep "^2R" heterochromMapping.$strain.paf | awk '$6 == "2R" ' | grep -v ":0-" | cut -f1-10 | sort -k8h
# grep "^3L" heterochromMapping.$strain.paf | awk '$6 == "3L" ' | grep  ":0-" | cut -f1-10 | sort -k8h
# grep "^3L" heterochromMapping.$strain.paf | awk '$6 == "3L" ' | grep -v ":0-" | cut -f1-10 | sort -k8h
# grep "^X" heterochromMapping.$strain.paf | awk '$6 == "X" ' | grep  ":0-" | cut -f1-10 | sort -k8h
# grep "^X" heterochromMapping.$strain.paf | awk '$6 == "X" ' | grep -v ":0-" | cut -f1-10 | sort -k8h
# done

#### GENOMIC AREAS AVAILABLE 

# Intergenic: extract genes from gff, -1000 and +1000 in the start and end position, bedtools complementary, remove regions overlapping with TEs (+-1000 start and end TEs), the remaining regions are intergenic regions that are +1kb away from genes and TEs.

# CDS, intron, 3'UTR and 5'UTR: extract these features and remove the ones that overlap a TE +- 1kb

# Upstream and downstream: extract 2kb upstream and downstream for genes, remove TEs (+1 kb)

mkdir "$DIR/DATA/genomic_areas"

# Intergenic

for strain in $assemblies
do
mkdir "$DIR/DATA/genomic_areas/$strain"
mkdir -p "$DIR/tmp/TEannotations/$strain/"

awk -v OFS='\t' -v bp=1000 '$3 == "gene" { $4 = $4 - bp; $5 = $5 + bp; print; }' "$DATA/GeneTransference/${strain}/${strain}.gff" > $DIR/DATA/genomic_areas/$strain/genes_1kb.gff
bedtools complement -i $DIR/DATA/genomic_areas/$strain/genes_1kb.gff -g /lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/${strain}/${strain}.chrom.sizes >  $DIR/DATA/genomic_areas/$strain/intergenic_regions.bed
awk -v OFS='\t' -v bp=2000 ' { $2 = $2 - bp; $3 = $3 + bp; print; }' "/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/TEannotations/${strain}.bed" > $DIR/tmp/TEannotations/$strain/$strain.2kb.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/intergenic_regions.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/intergenic_regions.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/intergenic_regions.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '$3 - $2 > 2000 { print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/intergenic_regions.clean.bed

done


# Introns

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "intron" ' "$DATA/GeneTransference/${strain}/${strain}v2.gff" | awk ' $5 > $4 ' | cut -f 1,4,5,9 | awk -F'\t' '{
  match($4, /gene_id=FBgn[0-9]+/, m);
  if (m[0] != "") {
    split(m[0], a, "=");
    print $1, $2, $3, a[2];
  }
}' OFS='\t' | sort -u > $DIR/DATA/genomic_areas/$strain/introns.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/introns.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/introns.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/introns.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2, $4 }' >  $DIR/DATA/genomic_areas/$strain/introns.clean.bed

done

tissues="head gut ovary"
for strain in $assemblies
do
for tissue in $tissues
do
## add gene expr
# Loop through each line in the file
while read intron
do
# Define the corresponding z-score file
zscore_file="../../R01-GB-expression_analysis/zscoreTMM/${strain}/${tissue}/zscore_${tissue}-${strain}.tab"
gene_id=$(echo "$intron" | cut -f 5)
# Check if the file exists
if [[ -f "$zscore_file" ]]; then
# Find the gene's expression level
expression_level=$(awk -v gene="$gene_id" '$1 == gene {print $0}' "$zscore_file" | cut -f4- -d' ' | tr ' ' '\t')

# If found, print the results
if [[ -n "$expression_level" ]]; then
echo -e "$intron\t$expression_level"
else
echo -e "$intron\tNA\tNA\tNA"
fi
else
echo -e "$intron\tFile not found: $zscore_file"
fi
done < $DIR/DATA/genomic_areas/$strain/introns.clean.bed | grep -v NA | grep -v NaN | grep -v Inf > "$DIR/DATA/genomic_areas/$strain/introns.clean.expr.$tissue.bed" 
done
done

# CDS

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "CDS" ' "$DATA/GeneTransference/${strain}/${strain}.gff" | awk ' $5 > $4 ' | cut -f 1,4,5,9  | awk -F'\t' '{
  match($4, /gene_id=FBgn[0-9]+/, m);
  if (m[0] != "") {
    split(m[0], a, "=");
    print $1, $2, $3, a[2];
  }
}' OFS='\t' | sort -u  > $DIR/DATA/genomic_areas/$strain/CDS.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/CDS.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/CDS.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/CDS.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2, $4 }' >  $DIR/DATA/genomic_areas/$strain/CDS.clean.bed

done


tissues="head gut ovary"
for strain in $assemblies
do
for tissue in $tissues
do
## add gene expr
# Loop through each line in the file
while read CDS
do
# Define the corresponding z-score file
zscore_file="../../R01-GB-expression_analysis/zscoreTMM/${strain}/${tissue}/zscore_${tissue}-${strain}.tab"
gene_id=$(echo "$CDS" | cut -f 5)
# Check if the file exists
if [[ -f "$zscore_file" ]]; then
# Find the gene's expression level
expression_level=$(awk -v gene="$gene_id" '$1 == gene {print $0}' "$zscore_file" | cut -f4- -d' ' | tr ' ' '\t')

# If found, print the results
if [[ -n "$expression_level" ]]; then
echo -e "$CDS\t$expression_level"
else
echo -e "$CDS\tNA\tNA\tNA"
fi
else
echo -e "$CDS\tFile not found: $zscore_file"
fi
done < $DIR/DATA/genomic_areas/$strain/CDS.clean.bed | grep -v NA | grep -v NaN | grep -v Inf > "$DIR/DATA/genomic_areas/$strain/CDS.clean.expr.$tissue.bed" 
done
done

# 5'UTR

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "5UTR" ' "$DATA/GeneTransference/${strain}/${strain}.gff" | awk ' $5 > $4 ' | cut -f 1,4,5,9  | awk -F'\t' '{
  match($4, /gene_id=FBgn[0-9]+/, m);
  if (m[0] != "") {
    split(m[0], a, "=");
    print $1, $2, $3, a[2];
  }
}' OFS='\t' | sort -u  > $DIR/DATA/genomic_areas/$strain/5UTR.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/5UTR.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/5UTR.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/5UTR.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }' | awk -v OFS='\t' ' {
  if ($5 in max && $4 > max[$5]) {
    max[$5] = $4
    line[$5] = $0
  } else if (!($5 in max)) {
    max[$5] = $4
    line[$5] = $0
  }
}
END {
  for (g in line) print line[g]
}
' | sort -u >  $DIR/DATA/genomic_areas/$strain/5UTR.clean.bed # I keep the longest 5UTR

done



tissues="head gut ovary"
for strain in $assemblies
do
for tissue in $tissues
do
## add gene expr
# Loop through each line in the file
while read UTR5
do
# Define the corresponding z-score file
zscore_file="../../R01-GB-expression_analysis/zscoreTMM/${strain}/${tissue}/zscore_${tissue}-${strain}.tab"
gene_id=$(echo "$UTR5" | cut -f 5)
# Check if the file exists
if [[ -f "$zscore_file" ]]; then
# Find the gene's expression level
expression_level=$(awk -v gene="$gene_id" '$1 == gene {print $0}' "$zscore_file" | cut -f4- -d' ' | tr ' ' '\t')

# If found, print the results
if [[ -n "$expression_level" ]]; then
echo -e "$UTR5\t$expression_level"
else
echo -e "$UTR5\tNA\tNA\tNA"
fi
else
echo -e "$UTR5\tFile not found: $zscore_file"
fi
done < $DIR/DATA/genomic_areas/$strain/5UTR.clean.bed | grep -v NA | grep -v NaN | grep -v Inf > "$DIR/DATA/genomic_areas/$strain/5UTR.clean.expr.$tissue.bed" 
done
done

# 3'UTR

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "3UTR" ' "$DATA/GeneTransference/${strain}/${strain}.gff" | awk ' $5 > $4 ' | cut -f 1,4,5,9  | awk -F'\t' '{
  match($4, /gene_id=FBgn[0-9]+/, m);
  if (m[0] != "") {
    split(m[0], a, "=");
    print $1, $2, $3, a[2];
  }
}' OFS='\t' | sort -u  > $DIR/DATA/genomic_areas/$strain/3UTR.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/3UTR.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/3UTR.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/3UTR.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }'  | awk -v OFS='\t' ' {
  if ($5 in max && $4 > max[$5]) {
    max[$5] = $4
    line[$5] = $0
  } else if (!($5 in max)) {
    max[$5] = $4
    line[$5] = $0
  }
}
END {
  for (g in line) print line[g]
}
' | sort -u >  $DIR/DATA/genomic_areas/$strain/3UTR.clean.bed

done



tissues="head gut ovary"
for strain in $assemblies
do
for tissue in $tissues
do
## add gene expr
# Loop through each line in the file
while read UTR3
do
# Define the corresponding z-score file
zscore_file="../../R01-GB-expression_analysis/zscoreTMM/${strain}/${tissue}/zscore_${tissue}-${strain}.tab"
gene_id=$(echo "$UTR3" | cut -f 5)
# Check if the file exists
if [[ -f "$zscore_file" ]]; then
# Find the gene's expression level
expression_level=$(awk -v gene="$gene_id" '$1 == gene {print $0}' "$zscore_file" | cut -f4- -d' ' | tr ' ' '\t')

# If found, print the results
if [[ -n "$expression_level" ]]; then
echo -e "$UTR3\t$expression_level"
else
echo -e "$UTR3\tNA\tNA\tNA"
fi
else
echo -e "$UTR3\tFile not found: $zscore_file"
fi
done < $DIR/DATA/genomic_areas/$strain/3UTR.clean.bed | grep -v NA | grep -v NaN | grep -v Inf > "$DIR/DATA/genomic_areas/$strain/3UTR.clean.expr.$tissue.bed" 
done
done


# Upstream

# Intergenic

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "gene" ' "$DATA/GeneTransference/${strain}/${strain}.gff" > $DIR/DATA/genomic_areas/$strain/genes.gff

# Define the GFF file path
gff_file="$DIR/DATA/genomic_areas/$strain/genes.gff"
> $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed
> $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed
# Loop through the GFF file and extract genes based on strand
while IFS=$'\t' read -r chrom source feature start end no strand rest attr; do
# Calculate upstream and downstream coordinates based on strand
#echo $chrom $start $end $strand
if [[ "$strand" == *"-"* ]]; then
# Gene is on the minus strand
upstream_start=$end
upstream_end=$((end + 1000))
downstream_start=$((start - 1000))
downstream_end=$start
else
# Gene is on the plus strand
upstream_start=$((start - 1000))
upstream_end=$start
downstream_start=$end
downstream_end=$((end + 1000))
fi
gene=$(echo "$attr" | cut -f1 -d';' | cut -f2 -d'=' )
# Output genes to the respective files
echo -e "$chrom\t$upstream_start\t$upstream_end\t$gene" >> $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed
echo -e "$chrom\t$downstream_start\t$downstream_end\t$gene" >> $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed
done < "$gff_file"
done

for strain in $assemblies
do

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed -b $DIR/DATA/genomic_areas/$strain/genes.gff >  $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp2.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp2.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }' >  $DIR/DATA/genomic_areas/$strain/upstream.clean.bed

done



tissues="head gut ovary"
for strain in $assemblies
do
for tissue in $tissues
do
## add gene expr
# Loop through each line in the file
while read upstream
do
# Define the corresponding z-score file
zscore_file="../../R01-GB-expression_analysis/zscoreTMM/${strain}/${tissue}/zscore_${tissue}-${strain}.tab"
gene_id=$(echo "$upstream" | cut -f 5)
# Check if the file exists
if [[ -f "$zscore_file" ]]; then
# Find the gene's expression level
expression_level=$(awk -v gene="$gene_id" '$1 == gene {print $0}' "$zscore_file" | cut -f4- -d' ' | tr ' ' '\t')

# If found, print the results
if [[ -n "$expression_level" ]]; then
echo -e "$upstream\t$expression_level"
else
echo -e "$upstream\tNA\tNA\tNA"
fi
else
echo -e "$upstream\tFile not found: $zscore_file"
fi
done < $DIR/DATA/genomic_areas/$strain/upstream.clean.bed | grep -v NA | grep -v NaN | grep -v Inf > "$DIR/DATA/genomic_areas/$strain/upstream.clean.expr.$tissue.bed" 
done
done


for strain in $assemblies
do

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed -b $DIR/DATA/genomic_areas/$strain/genes.gff >  $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp2.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp2.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }' >  $DIR/DATA/genomic_areas/$strain/downstream.clean.bed

done



tissues="head gut ovary"
for strain in $assemblies
do
for tissue in $tissues
do
## add gene expr
# Loop through each line in the file
while read downstream
do
# Define the corresponding z-score file
zscore_file="../../R01-GB-expression_analysis/zscoreTMM/${strain}/${tissue}/zscore_${tissue}-${strain}.tab"
gene_id=$(echo "$downstream" | cut -f 5)
# Check if the file exists
if [[ -f "$zscore_file" ]]; then
# Find the gene's expression level
expression_level=$(awk -v gene="$gene_id" '$1 == gene {print $0}' "$zscore_file" | cut -f4- -d' ' | tr ' ' '\t')

# If found, print the results
if [[ -n "$expression_level" ]]; then
echo -e "$downstream\t$expression_level"
else
echo -e "$downstream\tNA\tNA\tNA"
fi
else
echo -e "$downstream\tFile not found: $zscore_file"
fi
done < $DIR/DATA/genomic_areas/$strain/downstream.clean.bed | grep -v NA | grep -v NaN | grep -v Inf > "$DIR/DATA/genomic_areas/$strain/downstream.clean.expr.$tissue.bed" 
done
done 


### console create number_random_regions

tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
positions="3UTR 5UTR CDS downstream introns upstream"
for strain in $assemblies
do
for tissue in $tissues
do
for position in $positions
do
lines=$(cat "$DIR/DATA/genomic_areas/$strain/$position.clean.expr.$tissue.bed"  | cut -f5 | sort -u | wc -l)
echo -e "$tissue\t$strain\t$position\t$lines\tH3K9me3\tincrement"
done
done
done > $DIR/number_random_regions.tab

tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
positions="3UTR 5UTR CDS downstream introns upstream"
for strain in $assemblies
do
for tissue in $tissues
do
for position in $positions
do
lines=$(cat "$DIR/DATA/genomic_areas/$strain/$position.clean.expr.$tissue.bed"  | cut -f5 | sort -u | wc -l)
echo -e "$tissue\t$strain\t$position\t$lines\tH3K27ac\tincrement"
done
done
done > $DIR/number_random_regions_H3K27ac.tab
