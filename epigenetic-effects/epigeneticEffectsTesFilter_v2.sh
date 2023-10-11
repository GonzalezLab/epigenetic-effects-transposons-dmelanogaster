#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=2

#SBATCH --mem=8G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=enrichment2

#SBATCH --output=enrichment2.output

#SBATCH --error=enrichment2.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10
source activate /homes/users/mcoronado/.conda/envs/5GP


histones="H3K27ac H3K9me3"
tissues="head gut ovary"

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
# strains=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f2 | sort -u)
# nStatus=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f7 | sort -u | wc -l)
# status=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f7 | sort -u )
# if [[ $nStatus == 1 && $status != "noEnrichment" ]]; then
# awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue >> epigeneticEffectsTEs.clean.consistent_v3.tab
# fi
# #if [[ $nStatus == 1 && $status == "noEnrichment" ]]; then
# #awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue >> epigeneticEffectsTEs.clean.consistent.noEnrichment.tab
# #fi
# done
# done
# #done < <(tail -n+2 epigeneticEffectsTEs.clean.tab | grep -vP "1\tN" | cut -f1 | sort -u )
# done < <(tail -n+2 epigeneticEffectsTEs.clean.tab | cut -f1 | sort -u )


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



> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab

while read TE
do
for histone in $histones
do
for tissue in $tissues
do
strains=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f2 | sort -u)
nStrains=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f2 | sort -u |wc -l )
nStatus=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f7  | sort -u | wc -l)
status=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f7  | sort -u )
# Si sólo hay 1 status, y el status no es ni enrichment ni NA
if [[ $nStatus == 1 && $status != "noEnrichment" && $status != "NA" ]]; then
n=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep -vP "\tNA\t" | wc -l)
# Si hay por lo menos datos en 3 comparaciones
if [[ $n -ge 3 ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tconsistent"}' >> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
else
echo "No llega a 3"
fi
elif [[ $nStatus == 1 && $status == "noEnrichment" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tconsistent_noEnrichment"}' >> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
fi

# Si hay 1 cepa sólo pero 2 status
if [[ $nStrains == 1 && $nStatus == 2 ]];then
# Cuántas veces está el min
min=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort | uniq -c | sort |  sed 's/^\s*//' |head -n1 |cut -f1 -d' ' )
# Cuántas veces está el max
max=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f1 -d' ' )
# Cuál es el status asociado a max
maxValue=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f2 -d' ' )
# Si el min está 1 vez y el máximo 3, y el maximo no es ni no enrichment ni NA
if [[ $min == 1 && $max == 3 && $maxValue != "noEnrichment"  && $maxValue != "NA" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
#awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
elif [[ $min == 1 && $max == 3 && $maxValue == "noEnrichment"  && $maxValue != "NA" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
#awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
fi

# Si hay 2 cepas y dos status
#elif [[ $nStrains == 2 && $nStatus -ge 2 ]]; then
elif [[ $nStrains == 2 && $nStatus == 2 ]]; then
# maxValues=""
# tags=""
# tag="exit"
# # Para cada cepa
# while read strain
# do
# # Número de status
# nStatus=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $strain | grep $tissue | cut -f 7 |sort -u |wc -l)
# # Veces que está el min
# min=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $strain | grep $tissue | cut -f 7 |sort | uniq -c | sort |  sed 's/^\s*//' |head -n1 |cut -f1 -d' ' )
# # Veces que está el max
# max=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $strain | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f1 -d' ' )
# # Max value
# maxValue=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone  | grep $strain | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f2 -d' ' )
# # Vector con los max values
# maxValues=$(echo $maxValues $maxValue)
# # Si hay más de 1 status, el min es 1 y hay mas de 1 max, el max value es no enrichment ni NA
# if [[ nStatus -gt 1 && $min == 1 && $max > 1 && $maxValue != "noEnrichment"  && $maxValue != "NA" ]]; then
# tags=$(echo ok $tags)
# # Si solo hay 1 status, max value es mayor q 1 y el max value es no enrichment ni NA
# elif [[ $nStatus == 1 && $max > 1  && $maxValue != "noEnrichment"  && $maxValue != "NA" ]]; then
# tags=$(echo ok $tags)
# elif [[ $nStatus == 1 && $max > 1  && $maxValue == "noEnrichment"  && $maxValue != "NA" ]]; then
# tags=$(echo ok $tags)
# elif [[ $nStatus -gt 1 && $max > 1  && $maxValue == "noEnrichment"  && $maxValue != "NA" ]]; then
# tags=$(echo ok $tags)
# else
# tags=$(echo exit $tags)
# fi

# done < <(echo "$strains")

min=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| head -n1 | cut -f1 -d' ')
max=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| tail -n1 | cut -f1 -d' ')
minValue=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| head -n1 | cut -f2 -d' ')
maxValue=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| tail -n1 | cut -f2 -d' ')
# maxValues=$(echo "$maxValues" | tr ' ' '\n' |sort -u | wc -l)
# tag=$(echo "$tags" | tr ' ' '\n' |sort -u )

if [[ $min == 1 && $max == 5 && $maxValue != "noEnrichment" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
elif [[ $min == 1 && $max == 5 && $maxValue == "noEnrichment" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
fi
# if [[ $maxValues == 1 && $tag == "ok" && $maxValue != "noEnrichment" ]]; then
# #maxValue=$(echo "$maxValue" | tr ' ' '\n' |sort -u | wc -l)
# #awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# elif [[ $maxValues == 1 && $tag == "ok" && $maxValue == "noEnrichment" ]]; then
# #maxValue=$(echo "$maxValue" | tr ' ' '\n' |sort -u | wc -l)
# awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# #awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# fi

# Si hay 4 cepas y 2 status
elif [[ $nStrains == 4 && $nStatus == 2 ]]; then
# Cuántas veces está el min
min=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone   | grep $tissue | cut -f 7 |sort | uniq -c | sort |  sed 's/^\s*//' |head -n1 |cut -f1 -d' ' )
# Cuántas veces está el max
max=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone   | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f1 -d' ' )
# Cuántas veces está el max value
maxValue=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone   | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f2 -d' ' )
# Si el min está una vez, el máx 3, y el max value no es ni no enrichment ni NA
if [[ $min == 1 && $max == 3 && $maxValue != "noEnrichment"  && $maxValue != "NA" ]]; then
#awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue  | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
elif [[ $min == 1 && $max == 3 && $maxValue == "noEnrichment"  && $maxValue != "NA" ]]; then
#awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue  | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
fi
# Si hay 3 cepas y más de 2 status
#elif [[ $nStrains == 3 && $nStatus -ge 2 ]]; then
elif [[ $nStrains == 3 && $nStatus == 2 ]]; then

min=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| head -n1 | cut -f1 -d' ')
max=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| tail -n1 | cut -f1 -d' ')
minValue=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| head -n1 | cut -f2 -d' ')
maxValue=$( awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | cut -f 7 |sort |uniq -c | sed 's/^\s*//' |sort| tail -n1 | cut -f2 -d' ')
# maxValues=$(echo "$maxValues" | tr ' ' '\n' |sort -u | wc -l)
# tag=$(echo "$tags" | tr ' ' '\n' |sort -u )

if [[ $min == 1 && $max == 5 && $maxValue != "noEnrichment" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
elif [[ $min == 1 && $max == 5 && $maxValue == "noEnrichment" ]]; then
awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
fi
# tags=""
# tag="exit"
# # Para cada cepa
# while read strain
# do
# # Número de status
# nStatus=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $strain | grep $tissue | cut -f 7 |sort -u |wc -l)
# # Veces que está el min
# min=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $strain | grep $tissue | cut -f 7 |sort | uniq -c | sort |  sed 's/^\s*//' |head -n1 |cut -f1 -d' ' )
# # Veces que está el max
# max=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $strain | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f1 -d' ' )
# # Max value
# maxValue=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone  | grep $strain | grep $tissue | cut -f 7 |sort | uniq -c | sort | sed 's/^\s*//' | tail -n1 |cut -f2 -d' ' )
# # Si hay más de 1 status, el min es 1 y hay mas de 1 max, el max value es no enrichment ni NA
# if [[ $nStatus == 1 && $max == 2 && $maxValue != "noEnrichment"  && $maxValue != "NA" ]]; then
# tags=$(echo ok:$maxValue $tags)
# # Si solo hay 1 status, max value es mayor q 1 y el max value es no enrichment ni NA
# elif [[ $nStatus == 1 && $max == 2 && $maxValue == "noEnrichment"  && $maxValue != "NA" ]]; then
# tags=$(echo ok:$maxValue $tags)
# else
# tags=$(echo exit $tags)
# fi

# done < <(echo "$strains")

# tag=$(echo "$tags" | tr ' ' '\n' |sort | uniq -c | sed 's/^\s*//' | grep "2 ok" |wc -l )
# maxValue=$( echo "$tags" | tr ' ' '\n' |sort | uniq -c | sed 's/^\s*//' | grep "2 ok" | cut -f2 -d':' )
# if [[ $tag == 1 && $maxValue != "noEnrichment" ]]; then
# #maxValue=$(echo "$maxValue" | tr ' ' '\n' |sort -u | wc -l)
# awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab

# #awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# elif [[ $tag == 1 && $maxValue == "noEnrichment" ]]; then
# #maxValue=$(echo "$maxValue" | tr ' ' '\n' |sort -u | wc -l)
# awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | grep $maxValue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# #awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue | awk '{print $0  "\tthreshold=1_noEnrichment"}'>> epigeneticEffectsTEs.clean.consistent.threshold=1_v3.tab
# fi


fi
#if [[ $nStatus == 1 && $status == "noEnrichment" ]]; then
#awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.tab | grep $histone | grep $tissue >> epigeneticEffectsTEs.clean.consistent.noEnrichment.tab
#fi


done
done
#done < <(tail -n+2 epigeneticEffectsTEs.clean.tab | grep -vP "1\tN" | cut -f1 | sort -u )
done < <(tail -n+2 epigeneticEffectsTEs.clean.tab  | cut -f1 | sort -u )

# Sacar los bivalentes: grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab| cut -f1,5,6 | sort -u | cut -f1,2| uniq -c 
#grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab| cut -f1,5,6 | sort -u | cut -f1,2| uniq -c | sed 's/^\s*//' | tr ' ' '\t' > bivalent.status
#grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab| grep consistent| cut -f1,5,6 | sort -u | cut -f1,2| uniq -c | sed 's/^\s*//' | tr ' ' '\t' > bivalent.status.consistent

grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab| grep consistent > epigeneticEffectsTEs.clean.consistent_effects_v3.tab

grep -v noEnrichment epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab > epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab


while read TE
do
for tissue in $tissues
do
nhistone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_effects_v3.tab |grep $tissue | cut -f6 | sort -u |wc -l)
histone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_effects_v3.tab |grep $tissue | cut -f6 | sort -u)
neffects=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_effects_v3.tab |grep $tissue | cut -f7 | sort -u |wc -l)
effects=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent_effects_v3.tab |grep $tissue | cut -f7 | sort -u )
if [[ $nhistone == 1 ]]; then
echo -e "$TE\t$tissue\t$histone\t$effects"
elif [[ $nhistone == 2 && $neffects == 1 ]]; then
echo -e "$TE\t$tissue\tbivalent\t$effects"	
elif [[ $nhistone == 2 && $neffects == 2  ]]; then
echo -e "$TE\t$tissue\tmix\tmix"	
fi

done	
done < <(cut -f1 epigeneticEffectsTEs.clean.consistent_effects_v3.tab | sort -u) > status.consistent


while read TE
do
for tissue in $tissues
do
nhistone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | cut -f6 | sort -u |wc -l)
histone=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | cut -f6 | sort -u)
neffects=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | cut -f7 | sort -u |wc -l)
effects=$(awk -v TE="$TE" ' $1 == TE ' epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab  |grep $tissue | cut -f7 | sort -u )
if [[ $nhistone == 1 ]]; then
echo -e "$TE\t$tissue\t$histone\t$effects"
elif [[ $nhistone == 2 && $neffects == 1 ]]; then
echo -e "$TE\t$tissue\tbivalent\t$effects"	
elif [[ $nhistone == 2 && $neffects == 2  ]]; then
echo -e "$TE\t$tissue\tmix\tmix"	
fi

done	
done < <(cut -f1 epigeneticEffectsTEs.clean.consistent.threshold\=1_effects_v3.tab | sort -u) > status.threshold

# no enrichment
grep "_noEnrichment" epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab | cut -f1,5,6   |  sort -u | cut -f1,2,3 | sort -u |cut -f1 | uniq -c | sort -r | grep " 6 " | wc -l
77


grep "_noEnrichment" epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab | cut -f1,5,6  | grep H3K27ac |  sort -u | cut -f1,2,3 | sort -u |cut -f1 | uniq -c | sort -r | grep " 3 " | head

[mcoronado@mr-login2 intersectTesHistoneMarks]$ comm -12 <(grep "_noEnrichment" epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab | cut -f1,5,6  | grep H3K9me3 |  sort -u | cut -f1,2,3 | sort -u |cut -f1 | uniq -c | sort -r | grep " 3 "  | cut -f 8 -d ' ' |sort )  <(grep "_noEnrichment" epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab | cut -f1,5,6  | grep H3K27ac |  sort -u | cut -f1,2,3 | sort -u |cut -f1 | uniq -c | sort -r | grep " 3 "  | cut -f 8 -d ' ' |sort ) |wc -l
77
[mcoronado@mr-login2 intersectTesHistoneMarks]$ grep "_noEnrichment" epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab | cut -f1,5,6  | grep H3K9me3 |  sort -u | cut -f1,2,3 | sort -u |cut -f1 | uniq -c | sort -r | grep " 3 "  | cut -f 8 -d ' ' |sort |wc -l
226
[mcoronado@mr-login2 intersectTesHistoneMarks]$ grep "_noEnrichment" epigeneticEffectsTEs.clean.consistent.threshold\=1_v3.tab | cut -f1,5,6  | grep H3K27ac |  sort -u | cut -f1,2,3 | sort -u |cut -f1 | uniq -c | sort -r | grep " 3 "  | cut -f 8 -d ' ' |sort | wc -l
414
