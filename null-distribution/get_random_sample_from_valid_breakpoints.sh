# randomly sample breakpoint positions

random_region_n=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $DIR/number_random_regions.tab)


declare -A target_counts=(
  ["upstream"]=283
  ["downstream"]=216
  ["5UTR"]=906
  ["CDS"]=0
  ["introns"]=5418
  ["3UTR"]=0
)



used_genes_file="$DIR/DATA/random_sampled_areas/used_genes_gut_AKA-017_H3K9me3_increment.txt"
> "$used_genes_file"  # Empty the file


for pos in "${!target_counts[@]}"; do
  count=${target_counts[$pos]}
  input_file="TEs_epigenetics_H3K9me3_increment_AKA-017_gut.location_${pos}.random.tab"

  cat $input_file | grep -v -f $used_genes_file | shuf -n $count > $DIR/DATA/random_sampled_areas/random_sampled_H3K9me3_increment_AKA-017_gut_$pos.tab

  cut -f 10 $DIR/DATA/random_sampled_areas/random_sampled_H3K9me3_increment_AKA-017_gut_$pos.tab >> $used_genes_file

done
