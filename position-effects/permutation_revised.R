
library(data.table)
library(ggplot2)

run_permutation_test_tissue <- function(strain, histone_mark,  n_simulations = 1000, tissueT="head", dist="1kb",seed=123) {
  set.seed(seed)  # Ensure reproducibility
  
  # Define file paths for tissues
  base_path <- "zscoreTMM_geneTE"
  #tissues <- c("head", "gut", "ovary")
  tissues<-tissueT
  
  df_positions <- fread(paste0("INPUT/intersect/intersect1kb_genes_TE_",tissueT,"_",strain,"_locations.bed"), header = F, sep = "\t")
  colnames(df_positions) <- c("TE", "Gene", "Position")
  # Define custom priority (lower is higher priority)
  priority <- c("CDS" = 1, "5UTR" = 2, "3UTR" = 2, "introns" = 3, "upstream" = 4, "downstream" = 4)
  
  # Assign rank based on Position
  df_positions <- df_positions %>%
    mutate(Priority = priority[Position]) %>%
    group_by(TE, Gene) %>%
    slice_min(Priority, with_ties = FALSE) %>%
    ungroup()
  
TE_gene_pos<-  df_positions %>%
    distinct(TE, Gene, Position) %>%  # ensure unique TE-Gene-Position combinations
    group_by(Position) %>%
    summarise(Num_TE_Gene_Pairs = n()) %>%
    arrange(Num_TE_Gene_Pairs)
  
  # Read and combine data for all tissues
  data_list <- lapply(tissues, function(tissue) {
    file_path <- sprintf("%s/%s/%s/geneTE%s/zscore_%s-%s.tab", base_path, strain, tissue, dist, tissue, strain)
    if (file.exists(file_path)) {
      df <- fread(file_path, header = FALSE)
      na.omit(df)
    } else {
      message(sprintf("File not found: %s", file_path))
      NULL
    }
  })
  
  
  # Remove NULL values (files that don't exist)
  data_list <- Filter(Negate(is.null), data_list)
  
  # Combine into one dataframe
  combined_data <- rbindlist(data_list)
  
  df_TE<-fread(paste0("INPUT/intersect/intersect1kb_genes_TE_",tissueT,"_",strain,"_filtered.bed"), header=F)
  combined_data<-merge(combined_data,
                       df_TE,
                       by.x="V1",
                       by.y="V2")
  combined_data_pos <- unique(merge(combined_data, df_positions, by.x=c("V1","V1.y"), by.y=c("Gene", "TE")))
  
  zscore <- if (histone_mark == "H3K27ac") "Positive" else "Negative"
  
  # Prepare input files
  count_file <- paste0("../null_distribution/R01-GB-null_distribution_v2/RESULTS/counts_stable/zscoreTMM_TEs_effects_", effect, "_", histone_mark, "_", strain, "_count.tab")
  observed_file <- paste0("../null_distribution/R01-GB-null_distribution_v2/RESULTS/counts_stable/zscoreTMM_TEs_effects_", effect, "_", histone_mark, "_", strain, "_",zscore,"_count.tab")

  # Read counts and clean names
  sample_size <- fread(count_file)[, .(position, count = get(tissueT))]
  observed_value <- fread(observed_file)[, .(position, count = get(tissueT))]
  
  # Normalize position names
  for (dt in list(sample_size, observed_value)) {
    dt[position == "Intron", position := "introns"]
    dt[position == "5'UTR", position := "5UTR"]
    dt[position == "3'UTR", position := "3UTR"]
    dt[position == "Upstream", position := "upstream"]
    dt[position == "Downstream", position := "downstream"]
  }
  
  colnames(TE_gene_pos) <- c("position","count")
  TE_gene_pos$type<-"TE_gene_all"
  sample_size$type<-"Epigenetic effects"
  observed_value$type<-"Expression effects"
  stable<-rbind(TE_gene_pos, sample_size, observed_value)
  library(tidyr)
  
  stable_wide <- stable %>%
    pivot_wider(
      names_from = type,
      values_from = count,
      values_fill = 0  # fill missing combinations with 0
    )
  # Run the simulations
  set.seed(123)  # For reproducibility
  counts <- replicate(n_simulations, {
    sampled_subset <- rbindlist(lapply(1:nrow(sample_size), function(i) {
      pos <- sample_size$position[i]
      n <- sample_size$count[i]
      candidates <- combined_data_pos[Position == pos]
      if (nrow(candidates) >= n) {
        candidates[sample(.N, n)]
      } else {
        candidates[sample(.N, n, replace = TRUE)]
      }
    }))
    if (histone_mark == "H3K27ac") {
      sum(sampled_subset$V5 < 0.05 & sampled_subset$V4 > 0)
    } else {
      sum(sampled_subset$V5 < 0.05 & sampled_subset$V4 < 0)
    }
  })
  
  # Convert to dataframe
  counts_df <- data.frame(counts = counts)
  observed_total <- sum(observed_value$count)
  
  # Plot
  p<-ggplot(counts_df, aes(x = counts)) +
    geom_density(fill = "gray90", color = "gray90") +
    geom_vline(xintercept = observed_total, color = "firebrick4", linetype = "dashed") +
    theme_minimal(base_size = 16) +
    labs(
      title = paste(strain, "-", tissueT, "- Genes", ifelse(histone_mark == "H3K27ac", "upregulated", "downregulated")),
      x = "Frequency",
      y = "Density"
    ) +
    theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))
  
  # Empirical p-value
  p_value <- mean(counts >= observed_total)
  cat("Tissue", tissueT, "Strain", strain, "Histone mark", histone_mark, "\n")
  
  
  return(list(p_value = p_value,plot=p, stable_wide=stable_wide))
  
  
}

run_permutation_test_all <- function(strain, histone_mark, sample_size = c(head=32,gut=29,ovary=25), n_simulations = 1000, observed_value = c(head=17,gut=9,ovary=15), dist="1kb",seed=123) {
  set.seed(seed)  # Ensure reproducibility
  
  # Define file paths for tissues
  base_path <- "zscoreTMM_geneTE"
  tissues <- c("head", "gut", "ovary")
  # Read and combine data for all tissues
  data_list <- lapply(tissues, function(tissue) {
    file_path <- sprintf("%s/%s/%s/geneTE%s/zscore_%s-%s.tab", base_path, strain, tissue, dist, tissue, strain)
    if (file.exists(file_path)) {
      df <- fread(file_path, header = FALSE)
      na.omit(df)
    } else {
      message(sprintf("File not found: %s", file_path))
      NULL
    }
  })
  
  # Remove NULL values (files that don't exist)
  data_list <- Filter(Negate(is.null), data_list)
  
  # Combine into one dataframe
  combined_data <- rbindlist(data_list)
  
  # Read and combine data for all tissues
  data_list_TE <- lapply(tissues, function(tissue) {
    file_path <- sprintf("INPUT/intersect/intersect%s_genes_TE_%s_%s_filtered.bed", dist, tissue, strain)
    if (file.exists(file_path)) {
      df <- fread(file_path, header = FALSE)
      na.omit(df)
    } else {
      message(sprintf("File not found: %s", file_path))
      NULL
    }
  })
  
  # Remove NULL values (files that don't exist)
  data_list_TE <- Filter(Negate(is.null), data_list_TE)
  
  # Combine into one dataframe
  combined_data_TE <- rbindlist(data_list_TE)
  combined_data_TE<-unique(combined_data_TE)
  
  # Define condition based on histone mark
  if (histone_mark == "H3K9me3") {
    condition <- combined_data$V5 < 0.05 & combined_data$V4 < 0
  } else if (histone_mark == "H3K27ac") {
    condition <- combined_data$V5 < 0.05 & combined_data$V4 > 0
  } else {
    stop("Invalid histone mark! Choose 'H3K9me3' or 'H3K27ac'.")
  }
  
  combined_data<-merge(combined_data,
                       combined_data_TE,
                       by.x="V1",
                       by.y="V2")
  
  # Run simulations
  counts <- replicate(n_simulations, {
    sample_head_TE <- sample(unique(combined_data[combined_data$V2=="head",]$V1.y),sample_size["head"])
    combined_data_head_sampled <- combined_data[combined_data$V1.y%in%sample_head_TE & combined_data$V2=="head",]
    #sampled_data <- combined_data[sample(.N, sample_size)]
    
    sample_ovary_TE <- sample(unique(combined_data[combined_data$V2=="ovary",]$V1.y),sample_size["ovary"])
    combined_data_ovary_sampled <- combined_data[combined_data$V1.y%in%sample_ovary_TE & combined_data$V2=="ovary",]
    
    sample_gut_TE <- sample(unique(combined_data[combined_data$V2=="gut",]$V1.y),sample_size["gut"])
    combined_data_gut_sampled <- combined_data[combined_data$V1.y%in%sample_gut_TE & combined_data$V2=="gut",]
    
    sampled_data <- rbind(combined_data_gut_sampled,combined_data_head_sampled,combined_data_ovary_sampled)
    
    if (histone_mark == "H3K9me3") {
      length(unique(sampled_data[sampled_data$V5 < 0.05 & sampled_data$V4 < 0,]$V1.y))
    } else if (histone_mark == "H3K27ac") {
      length(unique(sampled_data[sampled_data$V5 < 0.05 & sampled_data$V4 > 0,]$V1.y))    }
    
  })
  
  # Convert to data frame for plotting
  counts_df <- data.frame(counts = counts/sum(sample_size))
  
  # Plot results
  if (histone_mark == "H3K9me3") {
    plot <- ggplot(counts_df, aes(x = counts)) +
      geom_density(fill = "gray90", color = "gray90") +
      geom_vline(xintercept = sum(observed_value)/sum(sample_size), color = "firebrick4", linetype = "dashed") +
      theme_minimal(base_size = 16) +
      labs(
        title = sprintf("%s - %s: Genes nearby TEs (%s) that are downregulated", strain, histone_mark, dist),
        x = "Frequency",
        y = "Density"
      ) +
      theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))
  } else if (histone_mark == "H3K27ac") {
    plot <- ggplot(counts_df, aes(x = counts)) +
      geom_density(fill = "gray90", color = "gray90") +
      geom_vline(xintercept = sum(observed_value)/sum(sample_size), color = "firebrick4", linetype = "dashed") +
      theme_minimal(base_size = 16) +
      labs(
        title = sprintf("%s - %s: Genes nearby TEs (%s) that are upregulated", strain, histone_mark, dist),
        x = "Frequency",
        y = "Density"
      ) +
      theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))
    
  }
  print(plot)
  
  # Compute p-value
  p_value <- mean(counts >= sum(observed_value))
  message(sprintf("P-value for %s (%s): %.5f", strain, histone_mark, p_value))
  
  return(list(plot = plot, p_value = p_value))
}

tissues<-c("head","gut","ovary")

for (histonem in histones) {
for (s in strains) {
  for (tissue in tissues) {
   
  
  print(run_permutation_test_tissue(strain = s, histone_mark = histonem, n_simulations = 1000, dist="1kb", tissueT = tissue,seed=123))
  }
}
}

run_permutation_test_tissue("AKA-017", "H3K9me3", n_simulations = 1000, dist="1kb", tissueT = "head",seed=123)




run_permutation_test_tissue("AKA-017", "H3K9me3",sample_size = 25, observed_value = 15, n_simulations = 1000, dist="1kb", tissueT = "ovary",seed=123)
run_permutation_test_tissue("AKA-017", "H3K9me3",sample_size = 29, observed_value = 9, n_simulations = 1000, dist="1kb", tissueT = "gut",seed=123)

run_permutation_test_all("AKA-017", "H3K9me3",sample_size = c(head=32,gut=29,ovary=25), n_simulations = 1000, observed_value = c(head=17,gut=9,ovary=15), dist="1kb",seed=123)
