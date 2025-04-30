#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=20

#SBATCH --mem=16G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=epigeneticsEffectsTesFilter

#SBATCH --output=logs/epigeneticsEffectsTesFilter_%a.output

#SBATCH --error=logs/epigeneticsEffectsTesFilter_%a.error

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

#SLURM_ARRAY_TASK_ID=$1

random_region_n=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $DIR/number_random_regions.tab)
tissue=$(echo "$random_region_n" | cut -f1)
strain=$(echo "$random_region_n" | cut -f2)
position2=$(echo "$random_region_n" | cut -f3)
n=$(echo "$random_region_n" | cut -f4)
histone=$(echo "$random_region_n" | cut -f5) # manually set H3K27ac 
effect=$(echo "$random_region_n" | cut -f6)

cat $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.txt $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.H3K27ac.txt | tr ' ' '\t' | sort -u > $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.join.txt

awk 'NF' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.tmp.join.txt | sort -k1.8n -k6 | sed 's/\(region_[0-9][0-9]*\)/\n\1/g' | awk 'NF' > $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab

#> epigeneticEffectsTEs.clean.consistent.tab
#> epigeneticEffectsTEs.clean.consistent.noEnrichment.tab
#> epigeneticEffectsTEs.clean.consistent.polymorphicDiscordances.tab

# > epigeneticEffectsTEs.clean.consistent_v3.tab

# while read TE
# do
# for histone in $histones
# do
# for tissue in $tissues
# do
# strains=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | cut -f2 | sort -u)
# nStatus=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | cut -f7 | sort -u | wc -l)
# status=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | cut -f7 | sort -u )
# if [[ $nStatus == 1 && $status != "noEnrichment" ]]; then
# awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue >> epigeneticEffectsTEs.clean.consistent_v3.tab
# fi
# #if [[ $nStatus == 1 && $status == "noEnrichment" ]]; then
# #awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue >> epigeneticEffectsTEs.clean.consistent.noEnrichment.tab
# #fi
# done
# done
# #done < <(tail -n+2 $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep -vP "1\tN" | cut -f1 | sort -u )
# done < <(tail -n+2 $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| cut -f1 | sort -u )


# while read TE
# do
# for tissue in $tissues
# do
# nHistone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_v3.tab  | grep $tissue | cut -f 6 | sort -u| wc -l)
# histone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_v3.tab  | grep $tissue | cut -f 6 | sort -u )
# status=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_v3.tab | grep $tissue  | cut -f 6,7 | sort -u | tr '\n' ';' | tr '\t' '_')
# freq=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_v3.tab | grep $tissue  | cut -f 10 | sort -u)

# if [[ nHistone -eq 2 ]]; then
# echo $TE $tissue $status $freq
# elif [[ nHistone -eq 1 ]]; then
# echo $TE $tissue $status $freq
# fi

# done

# done < <(cat epigeneticEffectsTEs.clean.consistent_v3.tab | cut -f1 | sort -u)



#cd $DIR/intersectTesHistoneMarks

#cat $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.*.*.*.clean.tab > $DIR/intersectTesHistoneMarks/$DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab 


> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab

MAX_PARALLEL=20
COUNT=0

while read TEinfo
do
(
TE=$(echo "$TEinfo" | cut -f1)
chr=$(echo "$TEinfo" | cut -f3)
start=$(echo "$TEinfo" | cut -f4)
end=$(echo "$TEinfo" | cut -f5)
histone=$(echo "$TEinfo" | cut -f2)

strains=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | cut -f2 | sort -u)
nStrains=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | cut -f2 | sort -u |wc -l )
nStatus=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | cut -f7  | sort -u | wc -l)
status=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end |  cut -f7  | sort -u )
# Si sólo hay 1 status, y el status no es ni enrichment ni NA
if [[ $nStatus == 1 && $status != "noEnrichment" && $status != "NA" ]]; then
n=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -vP "\tNA\t" | wc -l)
# Si hay por lo menos datos en 3 comparaciones
if [[ $n -ge 3 ]]; then
awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue  | grep -w $chr | grep -w $start | grep -w $end | awk '{print $0  "\tconsistent"}' >> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab
else
echo "No llega a 3"
fi
elif [[ $nStatus == 1 && $status == "noEnrichment" ]]; then
awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue  | grep -w $chr | grep -w $start | grep -w $end | awk '{print $0  "\tconsistent_noEnrichment"}' >> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab
fi

# Si hay 1 cepa sólo pero 2 status
if [[ $nStrains == 1 && $nStatus == 2 ]];then
# Cuántas veces está el min
min=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | cut -f 7 |sort | uniq -c | sort |  sed 's/^\s*//' |head -n1 |cut -f1 -d' ' )
# Cuántas veces está el max
max=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f1 -d' ' )
# Cuál es el status asociado a max
maxValue=$(awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f2 -d' ' )
# Si el min está 1 vez y el máximo 3, y el maximo no es ni no enrichment ni NA
if [[ $min == 1 && $max == 3 && $maxValue != "noEnrichment"  && $maxValue != "NA" ]]; then
awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue  | grep -w $chr | grep -w $start | grep -w $end| grep $maxValue | awk '{print $0  "\tthreshold=1"}' >> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab
#awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1"}'>> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab
elif [[ $min == 1 && $max == 3 && $maxValue == "noEnrichment"  && $maxValue != "NA" ]]; then
awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | grep -w $chr | grep -w $start | grep -w $end | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}' >> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab
#awk -v TE="$TE" ' $1 == TE ' $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab
fi
fi

#done < <(tail -n+2 $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab| grep -vP "1\tN" | cut -f1 | sort -u )
  ) &

  ((COUNT++))
  if [[ $COUNT -ge $MAX_PARALLEL ]]; then
    wait
    COUNT=0
  fi
done < <(cat "$DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab" | cut -f1,10,11,12,6 | sort -u)

wait
#done < <(cat $DIR/intersectTesHistoneMarks/epigeneticEffectsTEs.$strain.$tissue.$position2.clean.tab | cut -f1,10,11,12,6 | sort -u )

grep -v noEnrichment RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.tab > RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab

while read TEinfo
do
(
TE=$(echo "$TEinfo" | cut -f1)
chr=$(echo "$TEinfo" | cut -f2)
start=$(echo "$TEinfo" | cut -f3)
end=$(echo "$TEinfo" | cut -f4)

nhistone=$(awk -v TE="$TE" ' $1 == TE ' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f6 | sort -u |wc -l)
histone=$(awk -v TE="$TE" ' $1 == TE ' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f6 | sort -u)
neffects=$(awk -v TE="$TE" ' $1 == TE ' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f7 | sort -u |wc -l)
effects=$(awk -v TE="$TE" ' $1 == TE ' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f7 | sort -u )
if [[ $nhistone == 1 ]]; then

if [[ $effects == "increment" ]];then
result=$(awk -v TE="$TE" '$1 == TE' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 < 0) print "mix" }' | grep mix | sort -u | wc -l)
if [[ $result -eq 1 ]];then
echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
else
echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
fi
elif [[ $effects == "depletion" ]];then
result=$(awk -v TE="$TE" '$1 == TE' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 > 0) print "mix" }' | grep mix | sort -u | wc -l )
if [[ $result -eq 1 ]];then
echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
else
echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
fi
fi


elif [[ $nhistone == 2 && $neffects == 1 ]]; then
histone="bivalent"
if [[ $effects == "increment" ]];then
result=$(awk -v TE="$TE" '$1 == TE' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 < 0) print "mix" }' | grep mix | sort -u | wc -l)
if [[ $result -eq 1 ]];then
echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
else
echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
fi
elif [[ $effects == "depletion" ]];then
result=$(awk -v TE="$TE" '$1 == TE' RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 > 0) print "mix" }' | grep mix | sort -u | wc -l)
if [[ $result -eq 1 ]];then
echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
else
echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
fi
fi
	
elif [[ $nhistone == 2 && $neffects == 2  ]]; then
echo -e "$TE\t$tissue\tmix\tmix\t$chr\t$start\t$end"	
fi
  ) &

  ((COUNT++))
  if [[ $COUNT -ge $MAX_PARALLEL ]]; then
    wait
    COUNT=0
  fi
done < <(cut -f1,10,11,12 RESULTS/epigeneticEffectsTEs.clean.consistent.threshold=1_v3.$strain.$tissue.$position2.effects.tab | sort -u) > RESULTS/status.threshold.$strain.$tissue.$position2.tab
wait
# I think this is not needed?

# grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab| grep consistent > epigeneticEffectsTEs.clean.consistent_effects_v3.tab

# grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab > epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab



# while read TEinfo
# do
# TE=$(echo "$TEinfo" | cut -f1)
# chr=$(echo "$TEinfo" | cut -f2)
# start=$(echo "$TEinfo" | cut -f3)
# end=$(echo "$TEinfo" | cut -f4)

# for tissue in $tissues
# do
# nhistone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f6 | sort -u |wc -l)
# histone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f6 | sort -u)
# neffects=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f7 | sort -u |wc -l)
# effects=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | grep -w $chr | grep -w $start | grep -w $end| cut -f7 | sort -u )
# if [[ $nhistone == 1 ]]; then

# if [[ $effects == "increment" ]];then
# result=$(awk -v TE="$TE" '$1 == TE' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 < 0) print "mix" }' | grep mix | sort -u | wc -l)
# if [[ $result -eq 1 ]];then
# echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
# else
# echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
# fi
# elif [[ $effects == "depletion" ]];then
# result=$(awk -v TE="$TE" '$1 == TE' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 > 0) print "mix" }' | grep mix | sort -u | wc -l )
# if [[ $result -eq 1 ]];then
# echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
# else
# echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
# fi
# fi


# elif [[ $nhistone == 2 && $neffects == 1 ]]; then
# histone="bivalent"
# if [[ $effects == "increment" ]];then
# result=$(awk -v TE="$TE" '$1 == TE' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 < 0) print "mix" }' | grep mix | sort -u | wc -l)
# if [[ $result -eq 1 ]];then
# echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
# else
# echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
# fi
# elif [[ $effects == "depletion" ]];then
# result=$(awk -v TE="$TE" '$1 == TE' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab | grep "$tissue" | grep -w "$chr" | grep -w "$start" | grep -w "$end" | cut -f9 | awk '{if($1 > 0) print "mix" }' | grep mix | sort -u | wc -l)
# if [[ $result -eq 1 ]];then
# echo -e "$TE\t$tissue\t$histone\tmix\t$chr\t$start\t$end"
# else
# echo -e "$TE\t$tissue\t$histone\t$effects\t$chr\t$start\t$end"
# fi
# fi
	
# elif [[ $nhistone == 2 && $neffects == 2  ]]; then
# echo -e "$TE\t$tissue\tmix\tmix\t$chr\t$start\t$end"	
# fi

# done	
# done < <(cut -f1,10,11,12 epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab | sort -u) > status.threshold
