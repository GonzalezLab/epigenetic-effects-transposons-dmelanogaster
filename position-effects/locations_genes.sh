DIR=/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/R01-GB-expression_analysis/
DATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/"

module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

histones="H3K9me3 H3K27ac bivalent"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment depletion"


mkdir -p "$DIR/DATA/genomic_areas"

# Introns

for strain in $assemblies
do
mkdir -p "$DIR/DATA/genomic_areas/$strain"

awk -v OFS='\t' '$3 == "intron" ' "$DATA/GeneTransference/${strain}/${strain}v2.gff" | awk ' $5 > $4 ' | cut -f 1,4,5,9 | awk -F'\t' '{
  match($4, /gene_id=FBgn[0-9]+/, m);
  if (m[0] != "") {
    split(m[0], a, "=");
    print $1, $2, $3, a[2];
  }
}' OFS='\t' | sort -u > $DIR/DATA/genomic_areas/$strain/introns.bed
#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/introns.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/introns.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/introns.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2, $4 }' >  $DIR/DATA/genomic_areas/$strain/introns.clean.bed

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
#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/CDS.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/CDS.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/CDS.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2, $4 }' >  $DIR/DATA/genomic_areas/$strain/CDS.clean.bed

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
#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/5UTR.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/5UTR.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/5UTR.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }' | awk -v OFS='\t' ' {
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
#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/3UTR.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/3UTR.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/3UTR.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }'  | awk -v OFS='\t' ' {
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


# down and upstream 

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

#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed -b $DIR/DATA/genomic_areas/$strain/genes.gff >  $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp.bed

#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp2.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }' >  $DIR/DATA/genomic_areas/$strain/upstream.clean.bed
 
done



for strain in $assemblies
do

#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed -b $DIR/DATA/genomic_areas/$strain/genes.gff >  $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp.bed

#bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp2.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2,$4 }' >  $DIR/DATA/genomic_areas/$strain/downstream.clean.bed

done

for strain in $assemblies
do
cd $DIR/DATA/genomic_areas/$strain/
> all_positions.clean.bed
for file in 3UTR.clean.bed 5UTR.clean.bed CDS.clean.bed downstream.clean.bed  introns.clean.bed upstream.clean.bed; do
echo $file
feature=$(basename "$file" .clean.bed)
awk -v f="$feature" '{print $0 "\t" f}' "$file" >> all_positions.clean.bed
done
done

cd $DIR
for tissue in $tissues
do
for strain in $assemblies
do
  echo $tissue $strain
echo $tissue $strain > log_analysis
> $DIR/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${strain}_locations.tmp.bed
while read gene_TE
do
gene=$(echo "$gene_TE" | cut -f2)
TE=$(echo "$gene_TE" | cut -f1)

result=$(bedtools intersect -a <(grep -w "$TE" /lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/TEannotations/TE-Clean-Nested/$strain.bed) \
-b $DIR/DATA/genomic_areas/$strain/all_positions.clean.bed -wa -wb)

if [ -z "$result" ]; then
echo "$TE $gene not found" >> $DIR/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${strain}_locations.tmp.bed
else
echo "$result" | cut -f 4,12,13 >> $DIR/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${strain}_locations.tmp.bed
fi


done < INPUT/intersect/intersect1kb_genes_TE_${tissue}_${strain}_filtered.bed

done
done

cd $DIR
for tissue in $tissues
do
for strain in $assemblies
do
cat $DIR/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${strain}_locations.tmp.bed | sort -u > $DIR/INPUT/intersect/intersect1kb_genes_TE_${tissue}_${strain}_locations.bed
done
done
# zscoreTMM_geneTE/AKA-017/head/geneTE1kb/zscore_head-AKA-017.tab

# 2L_1133260_1133260_1360 FBgn0031304

# TE=2L_1133260_1133260_1360
# gene=FBgn0031304

# 2L_1133260_1133260_1360

# bedtools intersect -a <(grep 2L_1133260_1133260_1360 /lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/TEannotations/TE-Clean-Nested/$strain.bed) -b $DIR/DATA/genomic_areas/$strain/all_positions.clean.bed -wa -wb
