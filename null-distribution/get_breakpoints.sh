#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=8G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=get_breakpoints

#SBATCH --output=logs/get_breakpoints_%a.output

#SBATCH --error=logs/get_breakpoints_%a.error

module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

histones="H3K9me3 H3K27ac"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
effects="increment"

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/null_distribution/R01-GB-null_distribution_v2/"
DATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/"


#> $DIR/DATA/TEs_epigenetics_${histone}_${effect}.locations.random.tab

# # head
# tissue=head
# position=upstream
# n=201
# strain=AKA-017
# histone=H3K9me3
# effect=increment 

random_region_n=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $DIR/number_random_regions.tab)
tissue=$(echo "$random_region_n" | cut -f1)
strain=$(echo "$random_region_n" | cut -f2)
position=$(echo "$random_region_n" | cut -f3)
n=$(echo "$random_region_n" | cut -f4)
histone=$(echo "$random_region_n" | cut -f5)
effect=$(echo "$random_region_n" | cut -f6)

#cp $DIR/DATA/TEs_epigenetics_${strain}.locations.random.original.tab $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab

#cp $DIR/DATA/TEs_epigenetics_${strain}.locations.random.tab $DIR/DATA/TEs_epigenetics_${strain}.locations.random.original.tab

> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.breakpointDetected.tab 

> $DIR/status_${histone}_${effect}_${strain}_${tissue}.location_${position}.log

while read TEinfo
do
TE=$(echo "$TEinfo" | cut -f1)
echo "Processing region $TE $strain"

# get TE random new position
pos=$(echo "$TEinfo" | cut -f 7-9 | tr '\t' '.')
chr=$(echo "$TEinfo" | cut -f 7)
start=$(echo "$TEinfo" | cut -f 8)
end=$(echo "$TEinfo" | cut -f 9)
effect=$(echo "$TEinfo" | cut -f 4)
position=$(echo "$TEinfo" | cut -f 6)
gene=$(echo "$TEinfo" | cut -f 10)
expr=$(echo "$TEinfo" | cut -f 11-)

mkdir -p $DIR/tmp/$TE.$pos.$strain.$tissue/bed
mkdir -p $DIR/tmp/$TE.$pos.$strain.$tissue/sam
mkdir -p $DIR/tmp/$TE.$pos.$strain.$tissue/fasta
mkdir -p $DIR/tmp/$TE.$pos.$strain.$tissue/log
mkdir -p $DIR/breakpoint/$TE.$pos


echo "$TEinfo" | cut -f 7-9 |  awk -v TE="$TE" -v pos="$pos" '{print $0 "\t" TE ":" pos}' > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.$strain.$pos.bed
echo "$TEinfo" | cut -f 7-9 |  awk -v TE="$TE" -v pos="$pos" '{print $0 "\t" TE ":" pos}' > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.$strain.$pos.f1.bed
echo "$TEinfo" | cut -f 7-9 |  awk -v TE="$TE" -v pos="$pos" '{print $0 "\t" TE ":" pos}' > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.$strain.$pos.f2.bed

# bedops to get 500 upstream and 500 downstream
bedops --range -500:499 -u $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.$strain.$pos.bed > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.f1-TE-f2.$strain.$pos.bed
bedops --range -500:-1 -u $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.$strain.$pos.f1.bed > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.f1.$strain.$pos.bed
bedops --range 0:499 -u $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.$strain.$pos.f2.bed > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.f2.$strain.$pos.bed
# bedtools: get TE flanking regions
bedtools getfasta -fi $DIR/genomeAssemblies/${strain}.fasta -bed $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.f1.$strain.$pos.bed -name > $DIR/tmp/$TE.$pos.$strain.$tissue/fasta/$TE.f1.$strain.$pos.fasta
bedtools getfasta -fi $DIR/genomeAssemblies/${strain}.fasta -bed $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.f2.$strain.$pos.bed -name > $DIR/tmp/$TE.$pos.$strain.$tissue/fasta/$TE.f2.$strain.$pos.fasta
bedtools getfasta -fi $DIR/genomeAssemblies/${strain}.fasta -bed $DIR/tmp/$TE.$pos.$strain.$tissue/bed/$TE.f1-TE-f2.$strain.$pos.bed -name > $DIR/tmp/$TE.$pos.$strain.$tissue/fasta/$TE.f1-f2.$strain.$pos.fasta
# minimap2: map the TE flanking regions in the strains without the TE
strainsNoTE=$(grep -v $strain $DIR/genomeAssemblies/genomeAssembly.txt)
echo "$TE present $strain $pos 1 $effect $histone" >> $DIR/status_${histone}_${effect}_${strain}_${tissue}.location_${position}.log

while read -r strainNoTE
do
mkdir $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/
mkdir $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/
mkdir -p $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/
mkdir -p breakpoint/$TE.$pos/$strainNoTE

minimap2 -ax splice -G100k -B2 -t 2 -C5 $DIR/genomeAssemblies/${strainNoTE}.fasta $DIR/tmp/$TE.$pos.$strain.$tissue/fasta/$TE.f1.$strain.$pos.fasta > $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f1.$strainNoTE.$pos.sam &
minimap2 -ax splice -G100k -B2 -t 2 -C5 $DIR/genomeAssemblies/${strainNoTE}.fasta $DIR/tmp/$TE.$pos.$strain.$tissue/fasta/$TE.f2.$strain.$pos.fasta > $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f2.$strainNoTE.$pos.sam &
minimap2 -ax splice -G100k -B2 -t 2 -C5 $DIR/genomeAssemblies/${strainNoTE}.fasta $DIR/tmp/$TE.$pos.$strain.$tissue/fasta/$TE.f1-f2.$strain.$pos.fasta > $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.$pos.sam &
wait

# find the coordinates
samtools view -b $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f1.$strainNoTE.$pos.sam | bedtools bamtobed -i stdin > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.$pos.bed &
samtools view -b  $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f2.$strainNoTE.$pos.sam | bedtools bamtobed -i stdin  > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.$pos.bed &
samtools view -b  $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.$pos.sam | bedtools bamtobed -i stdin  > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.$pos.bed &
wait

# some quality checks (log)
echo "Analysis of TE: $TE in $strainNoTE ($pos)" > $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log

nf1=$(samtools view -F 0x04 -c $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f1.$strainNoTE.$pos.sam)
nf2=$(samtools view -F 0x04 -c $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f2.$strainNoTE.$pos.sam)
nf12=$(samtools view -F 0x04 -c $DIR/tmp/$TE.$pos.$strain.$tissue/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.$pos.sam)
if [ $nf1 -eq 1 -a $nf2 -eq 1 -a $nf12 -eq 1 ]
then
tag=1
echo "Perfect breakpoint detection"
# parse breakpoint: we consider the breakpoint in the end of the coordinate of the F1
#cut -f1,3,4,5,6 $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.$pos.bed > breakpoint/$TE.$pos/$strainNoTE/$TE.breakpoint.$strainNoTE.$pos.bed
cut -f1,3,4,5,6 $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.$pos.bed >> breakpoint/TEs_epigenetics_${histone}_${effect}_${strainNoTE}_${tissue}.location_${position}.random.tab
#cut -f1,3,4,5,6 $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.$pos.bed > $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.breakpoint.$strainNoTE.$pos.bed

# get stats
bpf1=$(awk '{print($3-$2)}' $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.$pos.bed)
bpf2=$(awk '{print($3-$2)}' $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.$pos.bed)
bpf12=$(awk '{print($3-$2)}' $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.$pos.bed)
endf1=$(cut -f3 $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.$pos.bed)
startf2=$(cut -f2  $DIR/tmp/$TE.$pos.$strain.$tissue/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.$pos.bed)
TSD=$((endf1-startf2))

echo "Number of alignments f1: $nf1"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Number of alignments f2: $nf2"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Number of alignments f1-f2: $nf12"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
 
echo "Size match f1: $bpf1"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Size match f2: $bpf2" >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Size match f1-f2: $bpf12" >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log

echo "Overlap  between f1 and f2 (TSD): $TSD" >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "$TE absent $strainNoTE $pos 1" >> $DIR/status_${histone}_${effect}_${strain}_${tissue}.location_${position}.log
elif [ $nf1 -gt 1 -o $nf2 -gt 1 -o $nf12 -gt 1 ]
then
tag=0
echo "Problem in the unique detection" >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Number of alignments f1: $nf1"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Number of alignments f2: $nf2"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "Number of alignments f1-f2: $nf12"  >>  $DIR/tmp/$TE.$pos.$strain.$tissue/log/$strainNoTE/info.$pos.log
echo "$TE absent $strainNoTE $pos 0" >> $DIR/status_${histone}_${effect}_${strain}_${tissue}.location_${position}.log
elif [ $nf1 -eq 0 -o $nf2 -eq 0 -o $nf12 -eq 0 ]
then
tag=0
echo "No detection"
echo "$TE absent $strainNoTE $pos 0" >> $DIR/status_${histone}_${effect}_${strain}_${tissue}.location_${position}.log
fi

#rm -rf breakpoint/$TE.$pos/$strainNoTE
done < <(echo "$strainsNoTE")

#sed -i "/$TEinfo/ s/$/\t$tag/" $DIR/DATA/TEs_epigenetics_${histone}_${effect}_$strain.locations.random.tab 
echo -e "$TEinfo\t$tag" >> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.breakpointDetected.tab 

rm -rf $DIR/tmp/$TE.$pos.$strain.$tissue
rm -rf breakpoint/$TE.$pos


done < <(cat "$DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.tab" )


# create the fina file per tissue and position
histone=H3K9me3
effect=increment
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
positions="3UTR 5UTR CDS downstream introns upstream"
# # # Get final list (use breakpointDetected)
# # # console
# for strain in $assemblies
# do
# for tissue in $tissues
# do
# for position in $positions
#do
grep -P "\t1$"  $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.breakpointDetected.tab  >  $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.final.tab
# done
# done
# done
# done

# cat $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_*.random.final.tab > "$DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.locations.random.final.tab"