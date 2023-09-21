#!/bin/bash

module load Miniconda3/4.9.2
source activate epigenetics

tissues="gut head ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="depletion increment"
histones="bivalent H3K27ac H3K9me3"

wDir="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/expression_analysis/"

cd $wDir

for assembly in $assemblies
do
mkdir -p tmp/$assembly
mkdir intersect
awk ' $3 == "gene" ' geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > geneAnnotationsLiftoff/${assembly}/${assembly}_genes.gtf

for tissue in $tissues
do
for effect in $effects
do
for histone in $histones
do
> intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed
while read TE
do
spread=$(grep "$tissue" ../epigenetic_effects/epigeneticEffectsTesParsed_v5.tab | awk -v histone="$histone" ' $13 == histone' | awk -v effect="$effect" ' $4 == effect'	 | grep -w "$TE" | cut -f8 | sort -h | head -n1 )
echo -e "$TE\t$spread"
grep -w  "$TE" TEannotations/${assembly}/TE_${assembly}.bed > tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed
bedtools window -w "${spread}000" -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b geneAnnotationsLiftoff/${assembly}/${assembly}_genes.gtf >> intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed
done < <(grep "$tissue" ../epigenetic_effects/epigeneticEffectsTesParsed_v5.tab | awk -v histone="$histone" ' $13 == histone' | awk -v effect="$effect" ' $4 == effect'	| awk -v assembly="$assembly" ' $3 == assembly ' | cut -f 1 | sort -u)

nTEs=$(cat intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed | cut -f4 |sort -u |wc -l)
nGenes=$(cat intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed | cut -f16 |sort -u |wc -l)
echo -e "$assembly\t$tissue\t$histone\t$effect\t$nlistTEs\t$nTEs ($nGenes)"
> intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}.lst

while read result
do
TE=$(echo "$result" | cut -f4)
spread=$(grep "$tissue" ../epigenetic_effects/epigeneticEffectsTesParsed_v5.tab | awk -v histone="$histone" ' $13 == histone' | awk -v effect="$effect" ' $4 == effect'	 | grep -w "$TE" | cut -f8 | sort -h | head -n1 )
gene=$(echo "$result" | cut -f 16 | cut -f 2 -d' ' |tr -d '"' |tr -d ';')
echo -e "$TE\t$gene\t$spread" >> intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}.lst
done < intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed

done
done
done
done


for assembly in $assemblies
do
for tissue in $tissues
do
for effect in $effects
do
for histone in $histones
do
> intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}_effect.lst
while read pair
do
TE=$(echo "$pair" | cut -f1)
gene=$(echo "$pair" | cut -f2)
intersect=$(awk -v TE="$TE" ' $4 == TE ' intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed | grep -w $gene )
TEstart=$(echo "$intersect" | cut -f 2)
TEend=$(echo "$intersect" | cut -f 3)
genestart=$(echo "$intersect" | cut -f 11)
geneend=$(echo "$intersect" | cut -f 12)
geneDir=$(echo "$intersect" | cut -f 14)

if [ "$TEend" -le "$genestart" -a "$geneDir" == "+" ]; then
echo -e "$TE\t$gene\tupstream"
elif [ "$TEend" -gt "$genestart" -a "$TEstart" -lt "$genestart" -a "$geneDir" == "+" ];then
echo -e "$TE\t$gene\tupstream"
elif [ "$TEstart" -ge "$geneend" -a "$geneDir" == "-" ]; then
echo -e "$TE\t$gene\tupstream"
elif [ "$TEstart" -lt "$geneend" -a "$TEend" -gt "$geneend" -a "$geneDir" == "-" ];then
echo -e "$TE\t$gene\tupstream"
elif [ "$TEend" -le "$geneend" -a "$TEstart" -ge "$genestart" ];then
echo -e "$TE\t$gene\twithin"
elif [ "$TEend" -le "$genestart" -a "$geneDir" == "-" ]; then
echo -e "$TE\t$gene\tdownstream"
elif [ "$TEend" -gt "$genestart" -a "$TEstart" -lt "$genestart" -a "$geneDir" == "-" ];then
echo -e "$TE\t$gene\tdownstream"
elif [ "$TEstart" -ge "$geneend" -a "$geneDir" == "+" ]; then
echo -e "$TE\t$gene\tdownstream"
elif [ "$TEstart" -lt "$geneend" -a "$TEend" -gt "$geneend" -a "$geneDir" == "+" ];then
echo -e "$TE\t$gene\tdownstream"
else
echo -e "$TE\t$gene\tproblem"
fi >> intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}_effect.lst  


done < intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}.lst
done
done
done
done

# position
for assembly in $assemblies
do
mkdir tmp/${assembly}/gene
for tissue in $tissues
do
for effect in $effects
do
for histone in $histones
do
> intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}_effect_position.lst
while read pair
do
TE=$(echo "$pair" | cut -f1)
gene=$(echo "$pair" | cut -f2)
intersect=$(awk -v TE="$TE" ' $4 == TE ' intersect/intersect_genes_TE_${tissue}_${histone}_${effect}_${assembly}.bed | grep -w $gene )
TEstart=$(echo "$intersect" | cut -f 2)
TEend=$(echo "$intersect" | cut -f 3)
genestart=$(echo "$intersect" | cut -f 11)
geneend=$(echo "$intersect" | cut -f 12)
geneDir=$(echo "$intersect" | cut -f 14)
echo "$assembly $tissue $histone $effect $gene $TE"
if [ "$TEend" -le "$genestart" -a "$geneDir" == "+" ]; then
dist=$(( $genestart - $TEend + 1 ))
echo -e "$TE\t$gene\tupstream\t$dist"

elif [ "$TEend" -gt "$geneend" -a "$TEstart" -lt "$genestart"  ];then
grep -w "$gene" geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > tmp/${assembly}/gene/${gene}.gtf
pos=$(bedtools intersect -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b  tmp/${assembly}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
echo -e "$TE\t$gene\tCDS\t0"

elif [ "$TEend" -gt "$genestart" -a "$TEstart" -lt "$genestart" -a "$geneDir" == "+" ];then
grep -w "$gene" geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > tmp/${assembly}/gene/${gene}.gtf
pos=$(bedtools intersect -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b  tmp/${assembly}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="5UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

elif [ "$TEstart" -ge "$geneend" -a "$geneDir" == "-" ]; then
dist=$(( $TEstart - $geneend + 1 ))
echo -e "$TE\t$gene\tupstream\t$dist"

elif [ "$TEstart" -lt "$geneend" -a "$TEend" -gt "$geneend" -a "$geneDir" == "-" ];then
grep -w "$gene" geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > tmp/${assembly}/gene/${gene}.gtf
pos=$(bedtools intersect -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b  tmp/${assembly}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" |  grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="5UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

elif [ "$TEend" -le "$geneend" -a "$TEstart" -ge "$genestart" ];then
grep -w "$gene" geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > tmp/${assembly}/gene/${gene}.gtf
pos=$(bedtools intersect -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b  tmp/${assembly}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" |grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u)
if [ -z $pos ]; then
pos=intron
fi
echo -e "$TE\t$gene\t$pos\t0"

elif [ "$TEend" -le "$genestart" -a "$geneDir" == "-" ]; then

dist=$(( $genestart - $TEend + 1 ))
echo -e "$TE\t$gene\tdownstream\t$dist"

elif [ "$TEend" -gt "$genestart" -a "$TEstart" -lt "$genestart" -a "$geneDir" == "-" ];then
grep -w "$gene" geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > tmp/${assembly}/gene/${gene}.gtf

pos=$(bedtools intersect -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b  tmp/${assembly}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="3UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

elif [ "$TEstart" -ge "$geneend" -a "$geneDir" == "+" ]; then
dist=$(( $TEstart - $geneend + 1 ))
echo -e "$TE\t$gene\tdownstream\t$dist"

elif [ "$TEstart" -lt "$geneend" -a "$TEend" -gt "$geneend" -a "$geneDir" == "+" ];then
grep -w "$gene" geneAnnotationsLiftoff/${assembly}/${assembly}.gtf > tmp/${assembly}/gene/${gene}.gtf

pos=$(bedtools intersect -a tmp/${assembly}/TE_${tissue}_${histone}_${effect}_${assembly}_${TE}.bed -b  tmp/${assembly}/gene/${gene}.gtf -wb -wa |grep -vP "\tgene\t" |grep -vP "\tmRNA\t" | grep -vP "\tRNA\t" | grep -vP "\tncRNA\t" | grep -vP "\texon\t" | grep -vP "\tpseudogene\t" | grep -vP "\tstop_codon\t" | cut -f 10 |sort -u |tr '\n' ';' )
if [ -z $pos ]; then
pos="3UTR;"
fi
echo -e "$TE\t$gene\t${pos::-1}\t0"

else
echo -e "$TE\t$gene\tproblem"
fi >> intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}_effect_position.lst  


done < intersect/list_genes_TE_${tissue}_${histone}_${effect}_${assembly}.lst
done
done
done
done
