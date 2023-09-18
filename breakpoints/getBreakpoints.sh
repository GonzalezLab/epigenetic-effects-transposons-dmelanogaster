#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=16G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=getBreakpoints

#SBATCH --output=getBreakpoints.output

#SBATCH --error=getBreakpoints.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10
source activate /homes/users/mcoronado/.conda/envs/5GP

assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

DIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics"

# Get list of TEs (unique and pol)
tail -n +2 ${DIR}/ANALYSES/TEsDescriptiveAnalyses/TE-Frequency/allStrains.tsv | awk ' $12 != "5" ' | cut -f 5 | sort -u > ${DIR}/ANALYSES/breakpointsTEs/listUniquePolTEs.lst
# Copy TE libraries
mkdir TEannotations
cp ${DIR}/ANALYSES/TEsDescriptiveAnalyses/TE-Frequency/*.tsv TEannotations 
# Copy genome assemblies
mkdir genomeAssemblies
for assembly in $assemblies
do
cp /homes/users/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/$assembly/$assembly.fasta genomeAssemblies
done
> status.log
while read TE
do
echo "Processing $TE"
mkdir -p tmp/$TE/bed
mkdir -p tmp/$TE/sam
mkdir -p tmp/$TE/fasta
mkdir -p tmp/$TE/log
mkdir -p breakpoint/$TE

nTE=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | wc -l )
if [ $nTE -eq 1 ]
then
echo "Unique TE $TE  $nTE"
# get info strain 
strain=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f 11 )
# get TE from the strain and create a tmp bed file
grep $TE TEannotations/$strain.tsv | grep $strain | awk -v OFS='\t'  '{print $2,$3,$4,$5,$6,$7,$1}' > tmp/$TE/bed/$TE.$strain.bed
grep $TE TEannotations/$strain.tsv | grep $strain | awk -v OFS='\t'  '{print $2,$3,$3,$5"-f1",$6,$7,$1}' > tmp/$TE/bed/$TE.f1.tmp.$strain.bed
grep $TE TEannotations/$strain.tsv | grep $strain | awk -v OFS='\t'  '{print $2,$4,$4,$5"-f2",$6,$7,$1}' > tmp/$TE/bed/$TE.f2.tmp.$strain.bed
# bedops to get 500 upstream and 500 downstream
bedops --range -500:500 -u tmp/$TE/bed/$TE.$strain.bed > tmp/$TE/bed/$TE.f1-TE-f2.$strain.bed
bedops --range -500:0 -u tmp/$TE/bed/$TE.f1.tmp.$strain.bed > tmp/$TE/bed/$TE.f1.$strain.bed
bedops --range 0:500 -u tmp/$TE/bed/$TE.f2.tmp.$strain.bed > tmp/$TE/bed/$TE.f2.$strain.bed
# bedtools: get TE flanking regions
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.f1.$strain.bed -name > tmp/$TE/fasta/$TE.f1.$strain.fasta
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.f2.$strain.bed -name > tmp/$TE/fasta/$TE.f2.$strain.fasta
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.f1-TE-f2.$strain.bed -name > tmp/$TE/fasta/$TE.f1-TE-f2.$strain.fasta
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.$strain.bed -name > tmp/$TE/fasta/$TE.$strain.fasta
cat tmp/$TE/fasta/$TE.f1.$strain.fasta > tmp/$TE/fasta/$TE.f1-f2.$strain.fasta
tail -n+2 tmp/$TE/fasta/$TE.f2.$strain.fasta >> tmp/$TE/fasta/$TE.f1-f2.$strain.fasta
# minimap2: map the TE flanking regions in the strains without the TE
strainsNoTE=$(grep -v $strain /homes/users/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/genomeAssembly.txt)
echo "$TE present $strain 1" >> status.log
while read -r strainNoTE
do
mkdir tmp/$TE/sam/$strainNoTE/
mkdir tmp/$TE/log/$strainNoTE/
mkdir -p tmp/$TE/bed/minimap2/$strainNoTE/
mkdir -p breakpoint/$TE/$strainNoTE

minimap2 -ax splice -G100k -B2 -t 4 -C5 genomeAssemblies/${strainNoTE}.fasta tmp/$TE/fasta/$TE.f1.$strain.fasta > tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam
minimap2 -ax splice -G100k -B2 -t 4 -C5 genomeAssemblies/${strainNoTE}.fasta tmp/$TE/fasta/$TE.f2.$strain.fasta > tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam
minimap2 -ax splice -G100k -B2 -t 4 -C5 genomeAssemblies/${strainNoTE}.fasta tmp/$TE/fasta/$TE.f1-f2.$strain.fasta > tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam
# find the coordinates
bedtools bamtobed -i tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam > tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed
bedtools bamtobed -i tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam > tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed
bedtools bamtobed -i tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam > tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.bed


# some quality checks (log)
echo "Analysis of TE: $TE in $strainNoTE" > tmp/$TE/log/$strainNoTE/info.log

nf1=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam)
nf2=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam)
nf12=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam)
if [ $nf1 -eq 1 -a $nf2 -eq 1 -a $nf12 -eq 1 ]
then
echo "Perfect breakpoint detection"

# parse breakpoint: we consider the breakpoint in the end of the coordinate of the F1
cut -f1,3,4,5,6 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed > breakpoint/$TE/$strainNoTE/$TE.breakpoint.$strainNoTE.bed
cut -f1,3,4,5,6 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed > tmp/$TE/bed/minimap2/$strainNoTE/$TE.breakpoint.$strainNoTE.bed

# get stats
bpf1=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
bpf2=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
bpf12=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.bed)
endf1=$(cut -f3 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
startf2=$(cut -f2  tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
TSD=$((endf1-startf2))

echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
 
echo "Size match f1: $bpf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f2: $bpf2" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f1-f2: $bpf12" >>  tmp/$TE/log/$strainNoTE/info.log

echo "Overlap  between f1 and f2 (TSD): $TSD" >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 1" >> status.log
elif [ $nf1 -gt 1 -o $nf2 -gt 1 -o $nf12 -gt 1 ]
then
echo "Problem in the unique detection" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 0" >> status.log
elif [ $nf1 -eq 0 -o $nf2 -eq 0 -o $nf12 -eq 0 ]
then
echo "No detection"
echo "$TE absent $strainNoTE 0" >> status.log
fi

done < <(echo "$strainsNoTE")

else
echo "Pol TE $TE $nTE"

# get info strain: if there are more than one, we pick the first one
strainAll=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f 11 )
strain=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f 11 | head -n 1)
# get TE from the strain and create a tmp bed file
grep $TE TEannotations/$strain.tsv | grep $strain | awk -v OFS='\t'  '{print $2,$3,$4,$5,$6,$7,$1}' > tmp/$TE/bed/$TE.$strain.bed
grep $TE TEannotations/$strain.tsv | grep $strain | awk -v OFS='\t'  '{print $2,$3,$3,$5"-f1",$6,$7,$1}' > tmp/$TE/bed/$TE.f1.tmp.$strain.bed
grep $TE TEannotations/$strain.tsv | grep $strain | awk -v  OFS='\t'  '{print $2,$4,$4,$5"-f2",$6,$7,$1}' > tmp/$TE/bed/$TE.f2.tmp.$strain.bed
# bedops to get 500 upstream and 500 downstream
bedops --range -500:500 -u tmp/$TE/bed/$TE.$strain.bed > tmp/$TE/bed/$TE.f1-TE-f2.$strain.bed
bedops --range -500:0 -u tmp/$TE/bed/$TE.f1.tmp.$strain.bed > tmp/$TE/bed/$TE.f1.$strain.bed
bedops --range 0:500 -u tmp/$TE/bed/$TE.f2.tmp.$strain.bed > tmp/$TE/bed/$TE.f2.$strain.bed
# bedtools: get TE flanking regions
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.f1.$strain.bed -name > tmp/$TE/fasta/$TE.f1.$strain.fasta
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.f2.$strain.bed -name > tmp/$TE/fasta/$TE.f2.$strain.fasta
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.f1-TE-f2.$strain.bed -name > tmp/$TE/fasta/$TE.f1-TE-f2.$strain.fasta
bedtools getfasta -fi genomeAssemblies/${strain}.fasta -bed tmp/$TE/bed/$TE.$strain.bed -name > tmp/$TE/fasta/$TE.$strain.fasta
cat tmp/$TE/fasta/$TE.f1.$strain.fasta > tmp/$TE/fasta/$TE.f1-f2.$strain.fasta
tail -n+2 tmp/$TE/fasta/$TE.f2.$strain.fasta >> tmp/$TE/fasta/$TE.f1-f2.$strain.fasta
# minimap2: map the TE flanking regions in the strains without the TE
strainsNoTE=$(grep -v "$strainAll" /homes/users/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/genomeAssembly.txt)
for strainTE in $strainAll
do
echo "$TE present $strainTE 1" >> status.log
done

while read -r strainNoTE
do
mkdir tmp/$TE/sam/$strainNoTE/
mkdir tmp/$TE/log/$strainNoTE/
mkdir -p tmp/$TE/bed/minimap2/$strainNoTE/
mkdir -p breakpoint/$TE/$strainNoTE

minimap2 -ax splice -G100k -B2 -t 4 -C5 genomeAssemblies/${strainNoTE}.fasta tmp/$TE/fasta/$TE.f1.$strain.fasta > tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam
minimap2 -ax splice -G100k -B2 -t 4 -C5 genomeAssemblies/${strainNoTE}.fasta tmp/$TE/fasta/$TE.f2.$strain.fasta > tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam
minimap2 -ax splice -G100k -B2 -t 4 -C5 genomeAssemblies/${strainNoTE}.fasta tmp/$TE/fasta/$TE.f1-f2.$strain.fasta > tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam
# find the coordinates
bedtools bamtobed -i tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam > tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed
bedtools bamtobed -i tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam > tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed
bedtools bamtobed -i tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam > tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.bed


# some quality checks (log)
echo "Analysis of TE: $TE in $strainNoTE" > tmp/$TE/log/$strainNoTE/info.log

nf1=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam)
nf2=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam)
nf12=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam)
if [ $nf1 -eq 1 -a $nf2 -eq 1 -a $nf12 -eq 1 ]
then
echo "Perfect breakpoint detection"

# parse breakpoint: we consider the breakpoint in the end of the coordinate of the F1
cut -f1,3,4,5,6 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed > breakpoint/$TE/$strainNoTE/$TE.breakpoint.$strainNoTE.bed
cut -f1,3,4,5,6 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed > tmp/$TE/bed/minimap2/$strainNoTE/$TE.breakpoint.$strainNoTE.bed

# get stats
bpf1=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
bpf2=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
bpf12=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.bed)
endf1=$(cut -f3 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
startf2=$(cut -f2  tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
TSD=$((endf1-startf2))

echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
 
echo "Size match f1: $bpf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f2: $bpf2" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f1-f2: $bpf12" >>  tmp/$TE/log/$strainNoTE/info.log

echo "Overlap  between f1 and f2 (TSD): $TSD" >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 1" >> status.log
elif [ $nf1 -gt 1 -o $nf2 -gt 1 -o $nf12 -gt 1 ]
then
echo "Problem in the unique detection" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 0" >> status.log
elif [ $nf1 -eq 0 -o $nf2 -eq 0 -o $nf12 -eq 0 ]
then
echo "No detection"
echo "$TE absent $strainNoTE 0" >> status.log
fi

done < <(echo "$strainsNoTE")

fi

done < <(cat ${DIR}/ANALYSES/breakpointsTEs/listUniquePolTEs.lst)




