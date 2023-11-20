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


DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/null_distribution/"
DATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/"

#### DEFINE HETEROCHROMATIC REGIONS

samtools faidx $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta > $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta.fai
cut -f1,2 $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta.fai | grep -Ew "2L|2R|3L|3R|X"> $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.sizes.genome

bedtools complement -i $DIR/DATA/eucromatic_regions_reference.bed -g $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.sizes.genome > $DIR/DATA/heterochromatin/heterochromatin_regions_reference.bed

bedtools getfasta -fi $DIR/DATA/heterochromatin/dmel-all-chromosome-r6.31.fasta \
-bed $DIR/DATA/heterochromatin/heterochromatin_regions_reference.bed \
-fo $DIR/DATA/heterochromatin/heterochromatin_regions_reference.fasta

module load GCCcore/11.2.0
module load minimap2/2.22

for strain in $assemblies
do

minimap2 -ax asm5 \
--secondary=no \
~/scratch/5GenomesProject/DATA/genomestrain/$strain/$strain.$pos.fasta \
$DIR/DATA/heterochromatin/heterochromatin_regions_reference.fasta > $DIR/DATA/heterochromatin/heterochromMapping.$strain.sam

minimap2 -x asm5 \
--secondary=no \
~/scratch/5GenomesProject/DATA/genomestrain/$strain/$strain.$pos.fasta \
$DIR/DATA/heterochromatin/heterochromatin_regions_reference.fasta > $DIR/DATA/heterochromatin/heterochromMapping.$strain.paf

Rscript heterochromatin.R $strain

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
bedtools complement -i $DIR/DATA/genomic_areas/$strain/genes_1kb.gff -g ~/scratch/5GenomesProject/DATA/genomestrain/$strain/$strain.chrom.sizes >  $DIR/DATA/genomic_areas/$strain/intergenic_regions.bed
awk -v OFS='\t' -v bp=2000 ' { $2 = $2 - bp; $3 = $3 + bp; print; }' "$DATA/TEannotations/$strain.$pos.bed" > $DIR/tmp/TEannotations/$strain/$strain.2kb.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/intergenic_regions.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/intergenic_regions.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/intergenic_regions.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '$3 - $2 > 2000 { print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/intergenic_regions.clean.bed

done


# Introns

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "intron" ' "$DATA/GeneTransference/${strain}/${strain}v2.gff" | awk ' $5 > $4 ' | cut -f 1,4,5  > $DIR/DATA/genomic_areas/$strain/introns.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/introns.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/introns.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/introns.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/introns.clean.bed

done

# CDS

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "CDS" ' "$DATA/GeneTransference/${strain}/${strain}.gff" | awk ' $5 > $4 ' | cut -f 1,4,5  > $DIR/DATA/genomic_areas/$strain/CDS.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/CDS.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/CDS.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/CDS.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/CDS.clean.bed

done

# 5'UTR

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "5UTR" ' "$DATA/GeneTransference/${strain}/${strain}.gff" | awk ' $5 > $4 ' | cut -f 1,4,5  > $DIR/DATA/genomic_areas/$strain/5UTR.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/5UTR.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/5UTR.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/5UTR.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/5UTR.clean.bed

done

# 3'UTR

for strain in $assemblies
do

awk -v OFS='\t' '$3 == "3UTR" ' "$DATA/GeneTransference/${strain}/${strain}.gff" | awk ' $5 > $4 ' | cut -f 1,4,5  > $DIR/DATA/genomic_areas/$strain/3UTR.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/3UTR.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/3UTR.clean.tmp.bed
bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/3UTR.clean.tmp.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/3UTR.clean.bed

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
while IFS=$'\t' read -r chrom source feature start end no strand rest; do
# Calculate upstream and downstream coordinates based on strand
echo $chrom $start $end $strand
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

# Output genes to the respective files
echo -e "$chrom\t$upstream_start\t$upstream_end" >> $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed
echo -e "$chrom\t$downstream_start\t$downstream_end" >> $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed
done < "$gff_file"
done

for strain in $assemblies
do

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.tmp.bed -b $DIR/DATA/genomic_areas/$strain/genes.gff >  $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp2.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/upstream.clean.tmp2.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/upstream.clean.bed

done

for strain in $assemblies
do

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.tmp.bed -b $DIR/DATA/genomic_areas/$strain/genes.gff >  $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp.bed -b $DIR/tmp/TEannotations/$strain/$strain.2kb.bed >  $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp2.bed

bedtools intersect -v -a $DIR/DATA/genomic_areas/$strain/downstream.clean.tmp2.bed -b $DIR/DATA/heterochromatin/heterochromatin_$strain.$pos.bed | awk -v OFS='\t' '{ print $1, $2, $3, $3 - $2 }' >  $DIR/DATA/genomic_areas/$strain/downstream.clean.bed

done

# After generating the regions that will be used for the null distribution, we will:
# generate list for TEs and effects per body part
# select each TE of the list, strain present, and strain absent list
# select a random position according to the TE position
# get 500pb upstream and downstream that position, and use the getBreakpoints.sh code to find the same position in the "absent strains"
# with this we will obtain the null-distribution positions

histones="H3K9me3 H3K27ac bivalent"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment depletion noEnrichment"

for histone in $histones
do
for tissue in $tissues
do
for effect in $effects
do
if [[ $effect == "noEnrichment" ]]; then
cat $DIR/DATA/TEs_epigenetics_effects.tab | grep $effect | cut -f1-3 | sort -u > $DIR/DATA/TEs_epigenetics_$effect.tab
else
grep $histone $DIR/DATA/TEs_epigenetics_effects.tab | grep $effect | cut -f1-3 | sort -u > $DIR/DATA/TEs_epigenetics_${histone}_$effect.tab
> $DIR/DATA/TEs_epigenetics_${histone}_$effect.locations.tab
while IFS=$'\t' read -r TE strain tissue; do


result=$( grep $TE $DIR/DATA/TEs_epigenetics_effects_genes.tab | grep $tissue | grep $effect | grep $strain)

if [ -n "$result" ]; then
if [[ $histone == "bivalent" ]]; then
histone2=Bivalent
else
histone2="$histone"
fi
awk -v TE="$TE" -v strain="$strain" -v tissue="$tissue" -v effect="$effect" -v histone="$histone2" -F'\t' '{
if ($1 == TE && $2 == strain && $3 == tissue && $4 == effect && $5 == histone) {
print $1 "\t" $2 "\t" $3 "\t" $4  "\t" $5  "\t" $7
} 

}' $DIR/DATA/TEs_epigenetics_effects_genes.tab >> $DIR/DATA/TEs_epigenetics_${histone}_$effect.locations.tab
else
echo -e "$TE\t$strain\t$tissue\t$effect\t$histone\tintergenic" >> $DIR/DATA/TEs_epigenetics_${histone}_$effect.locations.tab
fi

# Match TE information with location file and append position information
 # >> $DIR/DATA/TEs_epigenetics_$effect.locations.tab

done < $DIR/DATA/TEs_epigenetics_${histone}_$effect.tab
fi
done
done
done

cd $DIR

for strain in $assemblies
do
mkdir -p $DIR/tmp/$strain
mkdir $DIR/tmp/$strain/intersect

> $DIR/tmp/$strain/intersect/intersect_genes_TE_noEnrichment_${strain}.bed

while read TE
do
grep -w  "$TE" "$DATA/TEannotations/$strain.$pos.bed" > $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed
bedtools window -w "1000" -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b geneAnnotationsLiftoff/${strain}/${strain}_genes.gtf >> $DIR/tmp/$strain/intersect/intersect_genes_TE_noEnrichment_${strain}.bed
done < <(cat $DIR/DATA/TEs_epigenetics_noEnrichment.tab | grep $strain | cut -f 1 | sort -u)
done


for strain in $assemblies
do
> $DIR/tmp/$strain/intersect//list_genes_TE_noEnrichment_${strain}.lst
while read result
do
TE=$(echo "$result" | cut -f4)
gene=$(echo "$result" | cut -f 16 | cut -f 2 -d' ' |tr -d '"' |tr -d ';')
echo -e "$TE\t$gene" >> $DIR/tmp/$strain/intersect//list_genes_TE_noEnrichment_${strain}.lst
done < $DIR/tmp/$strain/intersect/intersect_genes_TE_noEnrichment_${strain}.bed
done

# position
for strain in $assemblies
do
mkdir $DIR/tmp/${strain}/gene

> $DIR/tmp/$strain/intersect/list_genes_TE_noEnrichment_${strain}_effect_position.lst
while read pair
do
TE=$(echo "$pair" | cut -f1)
gene=$(echo "$pair" | cut -f2)
intersect=$(awk -v TE="$TE" ' $4 == TE ' $DIR/tmp/$strain/intersect/intersect_genes_TE_noEnrichment_${strain}.bed | grep -w $gene )
TEstart=$(echo "$intersect" | cut -f 2)
TEend=$(echo "$intersect" | cut -f 3)
genestart=$(echo "$intersect" | cut -f 11)
geneend=$(echo "$intersect" | cut -f 12)
geneDir=$(echo "$intersect" | cut -f 14)
echo "$strain $tissue $histone $effect $gene $TE"
if [ "$TEend" -le "$genestart" -a "$geneDir" == "+" ]; then
dist=$(( $genestart - $TEend + 1 ))
echo -e "$TE\t$gene\tupstream\t$dist"

elif [ "$TEend" -gt "$geneend" -a "$TEstart" -lt "$genestart"  ];then
grep -w "$gene" geneAnnotationsLiftoff/${strain}/${strain}.gtf > $DIR/tmp/${strain}/gene/${gene}.gtf
pos=$(bedtools intersect -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b  $DIR/tmp/${strain}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
echo -e "$TE\t$gene\tCDS\t0"

elif [ "$TEend" -gt "$genestart" -a "$TEstart" -lt "$genestart" -a "$geneDir" == "+" ];then
grep -w "$gene" geneAnnotationsLiftoff/${strain}/${strain}.gtf > $DIR/tmp/${strain}/gene/${gene}.gtf
pos=$(bedtools intersect -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b  $DIR/tmp/${strain}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="5UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

elif [ "$TEstart" -ge "$geneend" -a "$geneDir" == "-" ]; then
dist=$(( $TEstart - $geneend + 1 ))
echo -e "$TE\t$gene\tupstream\t$dist"

elif [ "$TEstart" -lt "$geneend" -a "$TEend" -gt "$geneend" -a "$geneDir" == "-" ];then
grep -w "$gene" geneAnnotationsLiftoff/${strain}/${strain}.gtf > $DIR/tmp/${strain}/gene/${gene}.gtf
pos=$(bedtools intersect -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b  $DIR/tmp/${strain}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" |  grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="5UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

elif [ "$TEend" -le "$geneend" -a "$TEstart" -ge "$genestart" ];then
grep -w "$gene" geneAnnotationsLiftoff/${strain}/${strain}.gtf > $DIR/tmp/${strain}/gene/${gene}.gtf
pos=$(bedtools intersect -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b  $DIR/tmp/${strain}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" |grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u)
if [ -z $pos ]; then
pos=intron
fi
echo -e "$TE\t$gene\t$pos\t0"

elif [ "$TEend" -le "$genestart" -a "$geneDir" == "-" ]; then

dist=$(( $genestart - $TEend + 1 ))
echo -e "$TE\t$gene\tdownstream\t$dist"

elif [ "$TEend" -gt "$genestart" -a "$TEstart" -lt "$genestart" -a "$geneDir" == "-" ];then
grep -w "$gene" geneAnnotationsLiftoff/${strain}/${strain}.gtf > $DIR/tmp/${strain}/gene/${gene}.gtf

pos=$(bedtools intersect -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b  $DIR/tmp/${strain}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="3UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

elif [ "$TEstart" -ge "$geneend" -a "$geneDir" == "+" ]; then
dist=$(( $TEstart - $geneend + 1 ))
echo -e "$TE\t$gene\tdownstream\t$dist"

elif [ "$TEstart" -lt "$geneend" -a "$TEend" -gt "$geneend" -a "$geneDir" == "+" ];then
grep -w "$gene" geneAnnotationsLiftoff/${strain}/${strain}.gtf > $DIR/tmp/${strain}/gene/${gene}.gtf

pos=$(bedtools intersect -a $DIR/tmp/${strain}/TE_noEnrichment_${strain}_${TE}.bed -b  $DIR/tmp/${strain}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="3UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

else
echo -e "$TE\t$gene\tproblem"
fi >> $DIR/tmp/$strain/intersect/list_genes_TE_noEnrichment_${strain}_effect_position.lst


done < $DIR/tmp/$strain/intersect//list_genes_TE_noEnrichment_${strain}.lst
done


assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment depletion noEnrichment"
> $DIR/DATA/TEs_epigenetics_noEnrichment.locations.tab

for strain in $assemblies
do


while IFS=$'\t' read -r TE strain tissue; do


result=$( grep $TE $DIR/tmp/$strain/intersect/list_genes_TE_noEnrichment_${strain}_effect_position.lst )

if [ -n "$result" ]; then

awk -v TE="$TE" -v strain="$strain" -v tissue="$tissue" -F'\t' '{
if ($1 == TE ) {
print $1 "\t" strain "\t" tissue "\tnoEnrichment\tNA\t" $3
} 

}' $DIR/tmp/$strain/intersect/list_genes_TE_noEnrichment_${strain}_effect_position.lst  >> $DIR/DATA/TEs_epigenetics_noEnrichment.locations.tab
else
echo -e "$TE\t$strain\t$tissue\tnoEnrichment\tNA\tintergenic" >> $DIR/DATA/TEs_epigenetics_noEnrichment.locations.tab
fi

# Match TE information with location file and append position information
 # >> $DIR/DATA/TEs_epigenetics_$effect.locations.tab

done < <(grep $strain $DIR/DATA/TEs_epigenetics_noEnrichment.tab)

done

# manually changed 2L_16607264_16607264_roo for just 3UTR, and not 3TR,CDS
histones="H3K9me3 H3K27ac bivalent"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment depletion"


for effect in $effects
do
for histone in $histones
do
echo $effect $histone
> $DIR/DATA/TEs_epigenetics_${histone}_${effect}.locations.random.tab
ls $DIR/DATA/TEs_epigenetics_${histone}_${effect}.locations.tab

while IFS=$'\t' read -r TE strain tissue effect2 histone2 position; do

if [[ $position == "intergenic" ]]; then
position=intergenic_regions
elif [[ $position == "Intron" ]]; then
position=introns
elif [[ $position == "5'UTR" ]]; then
position=5UTR
elif [[ $position == "3'UTR" ]]; then
position=3UTR
elif [[ $position == "Downstream" ]]; then
position=downstream
elif [[ $position == "Upstream" ]]; then
position=upstream
fi

# Count the number of lines in the file (number of intergenic regions)
num_lines=$(wc -l < $DIR/DATA/genomic_areas/$strain/$position.clean.bed)

# Generate a random line number within the range of lines
random_line_number=$((1 + $RANDOM % num_lines))

# Get the intergenic region at the random line
random_region=$(sed -n "${random_line_number}p" $DIR/DATA/genomic_areas/$strain/$position.clean.bed)

# Extract the chromosome, start, and end coordinates
chromosome=$(echo "$random_region" | awk '{print $1}')
start=$(echo "$random_region" | awk '{print $2}')
end=$(echo "$random_region" | awk '{print $3}')

# Generate a random 1bp coordinate within the selected intergenic region
random_coordinate=$((start + $RANDOM % (end - start + 1)))

# Print the random 1bp coordinate
echo -e "$TE\t$strain\t$tissue\t$effect2\t$histone2\t$position\t$chromosome\t$random_coordinate\t$((random_coordinate + 1))" >> $DIR/DATA/TEs_epigenetics_${histone}_${effect}.locations.random.tab

done < "$DIR/DATA/TEs_epigenetics_${histone}_${effect}.locations.tab"

done
done

effect=noEnrichment
> $DIR/DATA/TEs_epigenetics_$effect.locations.random.tab

while IFS=$'\t' read -r TE strain tissue effect histone position; do

if [[ $position == "intergenic" ]]; then
position=intergenic_regions
elif [[ $position == "intron" ]]; then
position=introns
fi

# Count the number of lines in the file (number of intergenic regions)
num_lines=$(wc -l < $DIR/DATA/genomic_areas/$strain/$position.clean.bed)

# Generate a random line number within the range of lines
random_line_number=$((1 + $RANDOM % num_lines))

# Get the intergenic region at the random line
random_region=$(sed -n "${random_line_number}p" $DIR/DATA/genomic_areas/$strain/$position.clean.bed)

# Extract the chromosome, start, and end coordinates
chromosome=$(echo "$random_region" | awk '{print $1}')
start=$(echo "$random_region" | awk '{print $2}')
end=$(echo "$random_region" | awk '{print $3}')

# Generate a random 1bp coordinate within the selected intergenic region
random_coordinate=$((start + $RANDOM % (end - start + 1)))

# Print the random 1bp coordinate
echo -e "$TE\t$strain\t$tissue\t$effect\t$histone\t$position\t$chromosome\t$random_coordinate\t$((random_coordinate + 1))" >> $DIR/DATA/TEs_epigenetics_$effect.locations.random.tab

done < $DIR/DATA/TEs_epigenetics_$effect.locations.tab


for strain in $assemblies
do
grep "$strain" $DIR/DATA/TEs_epigenetics_noEnrichment.locations.random.tab | cut -f2 -d':' > $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab 
grep "$strain" $DIR/DATA/TEs_epigenetics_H*_*.locations.random.tab | cut -f2 -d':' >> $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab 

grep "$strain" $DIR/DATA/TEs_epigenetics_bivalent_*.locations.random.tab | cut -f2 -d':' >> $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab 

done

# Remove noEnrichment if there are increment or bivalent


for strain in $assemblies
do
grep -E "increment|depletion" $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab | cut -f1-3 > $DIR/DATA/TEs_epigenetics_${strain}.increment_depletion.locations.random.tab
mv $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab  $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tmp.tab
grep "mix" $DIR/DATA/TEs_epigenetics_effects.tab | awk -v strain="$strain" '$2 == strain ' | sort -u | cut -f 1,3 > $DIR/DATA/TEs_${strain}_mix_remove.lst
Rscript $DIR/filter_mix_increment_depletion.R $strain

done

# Get breakpoints

# Run script: get_breakpoints.sh

