#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=4

#SBATCH --mem=16G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=recoverBreakpoints

#SBATCH --output=recoverBreakpoints.output

#SBATCH --error=recoverBreakpoints.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10
source activate /homes/users/mcoronado/.conda/envs/5GP

assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

DIR="/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/roleTEsEpigenetics"

# After running getBreakpoints.sh, I need to check how many of the TEs without breakpoints 

grep " 0$" status.log > TEproblems.lst

while read TEproblem
do
TE=$(echo "$TEproblem" | cut -f1 -d' ' )
strain=$(echo "$TEproblem" | cut -f3 -d' ' )
nF1=$(cat tmp/$TE/log/$strain/info.log | grep "f1:" | cut -f 2 -d':')
nF2=$(cat tmp/$TE/log/$strain/info.log | grep " f2:" | cut -f 2 -d':')
nF12=$(cat tmp/$TE/log/$strain/info.log | grep "f1-f2:" | cut -f 2 -d':')
chrReal=$(awk -v TE="$TE" ' $5 == TE ' TEannotations/allStrains.tsv | cut -f2 | sort -u)
chr_F1=noOk
chr_F2=noOk
chr_F12=noOk

if [ -z $nF1 ]
then
nF1=0
elif [ $nF1 -eq 1 ]
then
chr_F1=$(cut -f1 tmp/$TE/bed/minimap2/$strain/$TE.f1.$strain.bed)
if [ $chr_F1 == $chrReal ]
then
chr_F1="ok"
else
chr_F1="noOk"
fi
fi

if [ -z $nF2 ]
then
nF2=0
elif [ $nF2 -eq 1 ]
then
chr_F2=$(cut -f1 tmp/$TE/bed/minimap2/$strain/$TE.f2.$strain.bed)
if [ $chr_F2 == $chrReal ]
then
chr_F2="ok"
else
chr_F2="noOk"
fi
fi

if [ -z $nF12 ]
then
nF12=0
elif [ $nF12 -eq 1 ]
then
chr_F12=$(cut -f1 tmp/$TE/bed/minimap2/$strain/$TE.f1-f2.$strain.bed)
if [ $chr_F12 == $chrReal ]
then
chr_F12="ok"
else
chr_F12="noOk"
fi
fi


echo $TE $strain $nF1-$chr_F1 $nF2-$chr_F2 $nF12-$chr_F12 
done < TEproblems.lst   > TEproblemsCheckN.log

# añadir TSD a status.log
cp status.log statusTSD.log
while read TEinfo
do
TE=$(echo "$TEinfo" | cut -f1 -d' ' )
strain=$(echo "$TEinfo" | cut -f3 -d' ' )
TSD=$(cat tmp/$TE/log/$strain/info.log | grep "TSD" | cut -f 2 -d':')
sed -i "s/$TE absent $strain 1/$TE absent $strain 1$TSD/g" statusTSD.log
done < <(grep absent status.log |grep 1$ )

# info total
cut -f1 -d' ' statusTSD.log | sort -u | wc -l # 3878
cat  statusTSD.log | grep absent | awk '$4 == 1'  | cut -f1 -d' ' | sort -u | wc -l # 3036
cut -f1,4 -d' ' statusTSD.log | sort -u | cut -f1 -d' ' | uniq -c | grep "1 " | wc -l # 2791

# info TSD smaller than 10 bp
cat  statusTSD.log | grep absent | awk '$4 == 1' | awk '$5 < 10' | cut -f1 -d' ' | sort -u | wc -l # 2615


# hacer copia de status.log y reempalzar 0 con 1 cuando se cumpla
cp statusTSD.log statusAfterRecoverBreakpoints.log
while read TEproblem
do
TE=$(echo "$TEproblem" | cut -f1 -d' ' )
strain=$(echo "$TEproblem" | cut -f2 -d' ' )
statusF1=$(echo "$TEproblem" | cut -f3 -d' ' )
statusF2=$(echo "$TEproblem" | cut -f4 -d' ' )
sed -i "s/$TE absent $strain 0/$TE absent $strain 1 NA/g" statusAfterRecoverBreakpoints.log
if [[ $statusF1 == "1-ok" ]]; then
cut -f1,3,4,5,6 tmp/$TE/bed/minimap2/$strain/$TE.f1.$strain.bed > breakpoint/$TE/$strain/$TE.breakpoint.$strain.bed

elif [[ $statusF2 == "1-ok" ]]; then
cut -f1,2,4,5,6 tmp/$TE/bed/minimap2/$strain/$TE.f2.$strain.bed > breakpoint/$TE/$strain/$TE.breakpoint.$strain.bed
fi
done < <(awk '$3 == "1-ok" || $4 == "1-ok" ' TEproblemsCheckN.log)

# info total
cut -f1 -d' ' statusTSD.log | sort -u | wc -l # 3878
cat  statusAfterRecoverBreakpoints.log | grep absent | awk '$4 == 1'  | cut -f1 -d' ' | sort -u | wc -l # 3036
cut -f1,4 -d' ' statusAfterRecoverBreakpoints.log | sort -u | cut -f1 -d' ' | uniq -c | grep "1 " | wc -l # 3207



# The TEs that we can analyze are:
cat  statusAfterRecoverBreakpoints.log | grep absent | grep "1$" | cut -f1 -d' ' | sort -u | wc -l
# So we create a new bed file with the subset of analyzable TEs

for assembly in $assemblies
do
TEs=$(grep "present $assembly" statusAfterRecoverBreakpoints.log  | cut -f1 -d' ')
grep -P -f "\t$TEs\t" $DIR/ANALYSES/TEsDescriptiveAnalyses/TE-Clean-Nested/$assembly.bed
done