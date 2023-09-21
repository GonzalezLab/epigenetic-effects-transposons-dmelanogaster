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

> status.log

while read TE
do
echo "Processing $TE"
nTE=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | wc -l )
if [ $nTE -eq 1 ]
then
echo "Unique TE $TE $nTE"
# get info strain 
strain=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f 11 )
strainsNoTE=$(grep -v $strain /homes/users/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/genomeAssembly.txt)
statusTE=$(awk -v TE="$TE" ' $4 == TE ' /gpfs42/projects/lab_jgonzalez/gonzalez_lab/santi/Drosophila/TEAnnot/common_name_problematic_flag/${strain}_TE_Annotation_problematic_flag.bed |cut -f 8)
echo "$TE present $strain 1 TE $statusTE" >> status.log

while read -r strainNoTE
do
# some quality checks (log)
echo "Analysis of TE: $TE in $strainNoTE" > tmp/$TE/log/$strainNoTE/info.log

nf1=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam)
nf2=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam)
nf12=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam)

if [ $nf1 -eq 1 -a $nf2 -eq 1 -a $nf12 -eq 1 ]
then
echo "Perfect breakpoint detection"

# get stats
bpf1=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
bpf2=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
bpf12=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.bed)
endf1=$(cut -f3 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
startf2=$(cut -f2  tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
TSD=$((startf2-endf1))


echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
 
echo "Size match f1: $bpf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f2: $bpf2" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f1-f2: $bpf12" >>  tmp/$TE/log/$strainNoTE/info.log

echo "Overlap  between f1 and f2 (TSD): $TSD" >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 1 $TSD $statusTE" >> status.log

elif [ $nf1 -gt 1 -o $nf2 -gt 1 -o $nf12 -gt 1 ]
then
echo "Problem in the unique detection" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 0 CC $statusTE" >> status.log

elif [ $nf1 -eq 0 -o $nf2 -eq 0 -o $nf12 -eq 0 ]
then
echo "No detection"
echo "$TE absent $strainNoTE 0 NA $statusTE" >> status.log
fi

done < <(echo "$strainsNoTE")

else
echo "Pol TE $TE $nTE"

# get info strain: if there are more than one, we pick the first one
strainAll=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f 11 )
strain=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f 11 | head -n 1) # get one of them
strainsNoTE=$(grep -v "$strainAll" /homes/users/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/genomeAssembly.txt)


for strainTE in $strainAll
do
statusTE=$(awk -v TE="$TE" ' $4 == TE ' /gpfs42/projects/lab_jgonzalez/gonzalez_lab/santi/Drosophila/TEAnnot/common_name_problematic_flag/${strainTE}_TE_Annotation_problematic_flag.bed |cut -f 8)
echo "$TE present $strainTE 1 TE $statusTE" >> status.log
done

statusTE=$(awk -v TE="$TE" ' $4 == TE ' /gpfs42/projects/lab_jgonzalez/gonzalez_lab/santi/Drosophila/TEAnnot/common_name_problematic_flag/${strain}_TE_Annotation_problematic_flag.bed |cut -f 8)


while read -r strainNoTE
do

# some quality checks (log)
echo "Analysis of TE: $TE in $strainNoTE" > tmp/$TE/log/$strainNoTE/info.log

nf1=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1.$strainNoTE.sam)
nf2=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f2.$strainNoTE.sam)
nf12=$(samtools view -F 0x04 -c tmp/$TE/sam/$strainNoTE/$TE.f1-f2.$strainNoTE.sam)

if [ $nf1 -eq 1 -a $nf2 -eq 1 -a $nf12 -eq 1 ]
then
echo "Perfect breakpoint detection"

# get stats
bpf1=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
bpf2=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
bpf12=$(awk '{print($3-$2)}' tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1-f2.$strainNoTE.bed)
endf1=$(cut -f3 tmp/$TE/bed/minimap2/$strainNoTE/$TE.f1.$strainNoTE.bed)
startf2=$(cut -f2  tmp/$TE/bed/minimap2/$strainNoTE/$TE.f2.$strainNoTE.bed)
TSD=$((startf2-endf1))

echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
 
echo "Size match f1: $bpf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f2: $bpf2" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Size match f1-f2: $bpf12" >>  tmp/$TE/log/$strainNoTE/info.log

echo "Overlap  between f1 and f2 (TSD): $TSD" >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 1 $TSD $statusTE" >> status.log

elif [ $nf1 -gt 1 -o $nf2 -gt 1 -o $nf12 -gt 1 ]
then
echo "Problem in the unique detection" >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1: $nf1"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f2: $nf2"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "Number of alignments f1-f2: $nf12"  >>  tmp/$TE/log/$strainNoTE/info.log
echo "$TE absent $strainNoTE 0 CC $statusTE" >> status.log

elif [ $nf1 -eq 0 -o $nf2 -eq 0 -o $nf12 -eq 0 ]
then
echo "No detection"
echo "$TE absent $strainNoTE 0 NA $statusTE" >> status.log
fi

done < <(echo "$strainsNoTE")

fi

done < <(cat ${DIR}/ANALYSES/breakpointsTEs/listUniquePolTEs.lst)


cut -f1 -d' ' status.log | sort -u | wc -l # 3878
cat  status.log | grep absent | awk '$4 == 1'  | cut -f1 -d' ' | sort -u | wc -l # 3036
cut -f1,4 -d' ' status.log | sort -u | cut -f1 -d' ' | uniq -c | grep "1 " | wc -l # 2791
cat  status.log | grep absent | awk '$4 == 1' | awk '$5 >= -50' | awk '$5 <= 50 '  | cut -f1 -d' ' | sort -u | wc -l


[mcoronado@mr-login2 breakpointsTEs]$ cat <(cat status.log | grep absent | awk '$4 == 1' | awk '$5 >= -10 && $5 <= 10 ' | cut -f1 -d' ' | sort -u ) <(cat status.log | grep absent | awk '$4 == 0' | cut -f1 -d' ' |sort -u ) | sort -u |wc -l
2972
[mcoronado@mr-login2 breakpointsTEs]$ comm -12 <(cat status.log | grep absent | awk '$4 == 1' | awk '$5 >= -10 && $5 <= 10 ' | cut -f1 -d' ' | sort -u ) <(cat status.log | grep absent | awk '$4 == 0' | cut -f1 -d' ' |sort -u ) | sort -u |wc -l
136
[mcoronado@mr-login2 breakpointsTEs]$ cat status.log | grep absent | awk '$4 == 1' | awk '$5 < -10 || $5 > 10 ' | cut -f1 -d' ' |sort -u | wc -l
1143
[mcoronado@mr-login2 breakpointsTEs]$ cat status.log | grep absent | awk '$4 == 0' | cut -f1 -d' ' |sort -u | wc -l         1087
[mcoronado@mr-login2 breakpointsTEs]$ cat status.log | grep absent | awk '$4 == 1' | awk '$5 >= -10 && $5 <= 10 ' | cut -f1 -d' ' | sort -u | wc -l
2021

 cat status.log | grep absent | awk '$4 == 1' | awk '$5 >= -50' | awk '$5 <= 50 ' > statusBreakpointAnalyzable50bp.log





