#!/bin/bash

# set the partition where the job will run

# set the number of nodes

#SBATCH --cpus-per-task=1

#SBATCH --mem=8G

# mail alert at start, end and abortion of execution

#SBATCH --job-name=randomRegions

#SBATCH --output=logs/randomRegions_%a.output

#SBATCH --error=logs/randomRegions_%a.error

module load Miniconda3/4.9.2
source activate epigenetics
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

# histones="H3K9me3 H3K27ac"
# tissues="head gut ovary"
# assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
# effects="increment"

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/null_distribution/R01-GB-null_distribution_v2/"
DATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/epigenetics/DATA/"



# # head
# tissue=head
# position=upstream
# n=201
# strain=AKA-017
# histone=H3K9me3
# effect=increment 

#random_region_n=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $DIR/number_random_regions.tab) # for H3K9me3 already run
random_region_n=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $DIR/number_random_regions_H3K27ac.tab)

tissue=$(echo "$random_region_n" | cut -f1)
strain=$(echo "$random_region_n" | cut -f2)
position=$(echo "$random_region_n" | cut -f3)
n=$(echo "$random_region_n" | cut -f4)
histone=$(echo "$random_region_n" | cut -f5)
effect=$(echo "$random_region_n" | cut -f6)

> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.tab
s=1
# new version
while read gene
do
random_region=$(grep -w "$gene" "$DIR/DATA/genomic_areas/$strain/$position.clean.expr.$tissue.bed" | shuf -n 1)

# Extract fields
chromosome=$(echo "$random_region" | awk '{print $1}')
start=$(echo "$random_region" | awk '{print $2}')
end=$(echo "$random_region" | awk '{print $3}')
gene=$(echo "$random_region" | cut -f 5)
expr=$(echo "$random_region" | cut -f 6-)

# Generate a random 1bp coordinate
random_coordinate=$((start + RANDOM % (end - start + 1)))

# Print to file
echo -e "region_$s\t$strain\t$tissue\tincrement\t$histone\t$position\t$chromosome\t$random_coordinate\t$((random_coordinate + 1))\t$gene\t$expr" >> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.tab
s=$(( $s + 1 ))

done < <(cut -f5 "$DIR/DATA/genomic_areas/$strain/$position.clean.expr.$tissue.bed" | sort -u)


# # Initialize an array to store used gene IDs
# used_genes=()

# for s in $(seq 1 $n)
# do
#     while true; do
#         # Count the number of lines in the file
#         num_lines=$(wc -l < "$DIR/DATA/genomic_areas/$strain/$position.clean.expr.$tissue.bed")

#         # Generate a random line number
#         random_line_number=$((1 + RANDOM % num_lines))

#         # Get the intergenic region at the random line
#         random_region=$(sed -n "${random_line_number}p" "$DIR/DATA/genomic_areas/$strain/$position.clean.expr.$tissue.bed")

#         # Extract fields
#         chromosome=$(echo "$random_region" | awk '{print $1}')
#         start=$(echo "$random_region" | awk '{print $2}')
#         end=$(echo "$random_region" | awk '{print $3}')
#         gene=$(echo "$random_region" | cut -f 5)
#         expr=$(echo "$random_region" | cut -f 6-)

#         # Check if gene has already been used
#         if [[ " ${used_genes[*]} " == *" $gene "* ]]; then
#             continue  # try another region
#         fi

#         # If it's a new gene, add it to the list and proceed
#         used_genes+=("$gene")

#         # Generate a random 1bp coordinate
#         random_coordinate=$((start + RANDOM % (end - start + 1)))

#         # Print to file
#         echo -e "region_$s\t$strain\t$tissue\tincrement\tH3K9me3\t$position\t$chromosome\t$random_coordinate\t$((random_coordinate + 1))\t$gene\t$expr" >> $DIR/DATA/TEs_epigenetics_${histone}_${effect}_${strain}_${tissue}.location_${position}.random.tab
#         break  # move to next s
#     done
# done
# Get breakpoints

# Run script: get_breakpoints.sh