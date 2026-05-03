# Load required packages for genomic analysis of assembled genomes as long as they possess .GFF and GFT files
   # Specific packages for exploratory genomic analysis e.g.
   # Basic genome statistics
   # Visualize genome comparisons
   # Gene content analysis
   # Synteny analysis preparation <<<<<>>>>> Generally this kind of analysis requires multiple contigs or chromossomes and blast results
   # Phylogenetics analysis 
   # Virulence factors search as well as other gene products (self criteria)

                      
library(Biostrings)
library(rtracklayer)
library("genoPlotR")
library("ape")
library(ggplot2)                        #Check out manual for an array of plots
library(dplyr)
                             # In case they aren't installed

install.packages("devtools")
install.packages("genoPlotR")
install.packages("ape")
install.packages("ggplot2")
install.packages("dplyr")
                            # If the package is not available download form github (devtools required)

devtools::install_github("username/repository")

                            # Use different aproaches suited to your files 

                            #library(rtracklayer)

              # Read GFF/GTF
              #gff_data <- import("annotation.gff")
              #gtf_data <- import("annotation.gtf")

              # Read GBFF (requires different approach)
              #library(genbankr)
              #gbff_data <- readGenBank("sequence.gbff")

              # Or read GBFF as text and parse
              #gb_lines <- readLines("sequence.gbff")



                            # This will return the package path if installed, or error if not
tryCatch({
  find.package("ggplot2")
  print("Package is installed")
}, error = function(e) {
  print("Package is not installed")
})


                             # Create project directory structure
dir.create("Sporothrix_spp")  # your own
setwd("")                     # your own
dir.create("results")
dir.create("plots")
                              #List files in wd

list.files("C:/Users/kelto/Downloads/Sporothrix spp/sporothrixspp") #file inside files to be listed and readed



                             # Read genome sequences

sporothrix_genomes <- list()
genome_files <- list.files(pattern = "\\.fna$|\\.fasta$")

for(file in genome_files){
  species_name <- gsub("\\.fna$|\\.fasta$", "", file)
  sporothrix_genomes[[species_name]] <- readDNAStringSet(file)
  cat("Loaded", species_name, "-", length(sporothrix_genomes[[species_name]]), "sequences\n")
}

                            # Read annotation files
sporothrix_annotations <- list()
gff_files <- list.files(pattern = "\\.gff$|\\.gff3$")

for(file in gff_files){
  species_name <- gsub("\\.gff$|\\.gff3$", "", file)
  sporothrix_annotations[[species_name]] <- import(file)      #error here, fixed with reentering first line
  cat("Loaded", species_name, "annotations -", 
      length(sporothrix_annotations[[species_name]]), "features\n")
}



                            #Basic genome statistics 
                            # Calculate basic genome statistics

genome_stats <- data.frame()

for(species in names(sporothrix_genomes)){
  genome <- sporothrix_genomes[[species]]
  annotation <- sporothrix_annotations[[species]]
  
  total_length <- sum(width(genome))
  n_contigs <- length(genome)
  gc_content <- sum(letterFrequency(genome, "GC")) / total_length
  n_genes <- sum(annotation$type == "gene", na.rm = TRUE)
  n_cds <- sum(annotation$type == "CDS", na.rm = TRUE)
  
  stats <- data.frame(
    Species = species,
    Total_length = total_length,
    N_contigs = n_contigs,
    GC_content = round(gc_content, 4),
    N_genes = n_genes,
    N_CDS = n_cds
  )
  
  genome_stats <- rbind(genome_stats, stats)
}

print(genome_stats)


                          #Vizualise genome comparisons 
                          # Plot genome statistics

library(ggplot2)

# Genome size comparison

p1 <- ggplot(genome_stats, aes(x = Species, y = Total_length/1e6, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(title = "Genome Size Comparison", y = "Genome Size (Mb)") +
  theme_minimal()

# GC content comparison

p2 <- ggplot(genome_stats, aes(x = Species, y = GC_content, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(title = "GC Content Comparison", y = "GC Content") +
  theme_minimal()

# Save plots
ggsave("plots/genome_sizes.png", p1, width = 8, height = 6)
ggsave("plots/gc_content.png", p2, width = 8, height = 6)


                    #Gene content analysis
                    # Compare gene features across species

gene_features <- data.frame()

for(species in names(sporothrix_annotations)){
  annotation <- sporothrix_annotations[[species]]
  genes <- annotation[annotation$type == "gene"]
  cds <- annotation[annotation$type == "CDS"]
  
  # Gene lengths
  gene_lengths <- width(genes)
  cds_lengths <- width(cds)
  
  features <- data.frame(
    Species = species,
    Feature = "Genes",
    Mean_length = mean(gene_lengths, na.rm = TRUE),
    Median_length = median(gene_lengths, na.rm = TRUE),
    Count = length(genes)
  )
  
  gene_features <- rbind(gene_features, features)
}

print(gene_features)


      # Prepare data for synteny analysis (if you have multiple contigs/chromosomes)
      # This requires BLAST results - we'll simulate or you can generate real BLAST output

# For now, let's create a simple comparison plot
comparison_plot <- function(stats_df, variable, title) {
  ggplot(stats_df, aes_string(x = "Species", y = variable, fill = "Species")) +
    geom_bar(stat = "identity") +
    labs(title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create multiple comparison plots
plots <- list()
plots[[1]] <- comparison_plot(genome_stats, "Total_length", "Total Genome Length")
plots[[2]] <- comparison_plot(genome_stats, "GC_content", "GC Content")
plots[[3]] <- comparison_plot(genome_stats, "N_genes", "Number of Genes")

# Arrange plots
library(gridExtra)
grid.arrange(grobs = plots, ncol = 2)






      # Extract conserved genes for phylogenetic analysis
      # For Sporothrix, we might look at housekeeping genes

      # Simple approach: use whole genome alignment or conserved single-copy genes
      # This requires external tools like OrthoFinder or BUSCO

      # Alternative: Calculate genomic distances based on k-mers
library(ape)

calculate_genomic_distance <- function(genome_list) {
      # Simple k-mer based distance (for demonstration)
      # In practice, use proper alignment methods
  kmer_matrix <- matrix(0, nrow = length(genome_list), ncol = length(genome_list))
  rownames(kmer_matrix) <- names(genome_list)
  colnames(kmer_matrix) <- names(genome_list)
  
  for(i in 1:length(genome_list)) {
    for(j in 1:length(genome_list)) {
      if(i != j) {
        # Simple comparison - in real analysis use proper distance metrics
        kmer_matrix[i,j] <- 1 - (min(length(genome_list[[i]]), 
                                     length(genome_list[[j]])) / 
                                   max(length(genome_list[[i]]), 
                                       length(genome_list[[j]])))
      }
    }
  }
  return(as.dist(kmer_matrix))
}

# Build simple tree (if you have multiple genomes)
if(length(sporothrix_genomes) > 2) {
  genomic_dist <- calculate_genomic_distance(sporothrix_genomes)
  sporothrix_tree <- nj(genomic_dist)
  plot(sporothrix_tree, main = "Sporothrix Genomic Relationships")
}



# Write statistics to file
write.csv(genome_stats, "results/genome_statistics.csv", row.names = FALSE)
write.csv(gene_features, "results/gene_features.csv", row.names = FALSE)

# Save session info for reproducibility
sink("results/session_info.txt")
sessionInfo()
sink()




              # Search for known Sporothrix virulence factors in annotations
virulence_genes <- c("proteinase", "lipase", "melanin", "heat shock protein")

find_virulence_genes <- function(annotation, keywords) {
  hits <- list()
  for(keyword in keywords) {
    gene_hits <- annotation[grep(keyword, annotation$product, ignore.case = TRUE)]
    if(length(gene_hits) > 0) {
      hits[[keyword]] <- gene_hits
    }
  }
  return(hits)
}

# Search in each species
virulence_results <- list()
for(species in names(sporothrix_annotations)) {
  virulence_results[[species]] <- find_virulence_genes(
    sporothrix_annotations[[species]], virulence_genes)
  cat("Found", length(unlist(virulence_results[[species]])), 
      "potential virulence genes in", species, "\n")
}





          # Create a summary data frame of virulence gene counts
virulence_summary <- data.frame()

for(species in names(virulence_results)) {
  species_results <- virulence_results[[species]]
  
  for(keyword in names(species_results)) {
    count <- length(species_results[[keyword]])
    if(count > 0) {
      summary_row <- data.frame(
        Species = species,
        Virulence_Factor = keyword,
        Count = count,
        Genes = paste(species_results[[keyword]]$ID, collapse = "; ")
      )
      virulence_summary <- rbind(virulence_summary, summary_row)
    }
  }
}

# Print summary
print(virulence_summary)






library(ggplot2)
library(dplyr)

# Plot 1: Total virulence factors by species
p1 <- ggplot(virulence_summary, aes(x = Species, y = Count, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Virulence Factors by Sporothrix Species",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("plots/virulence_factors_total.png", p1, width = 8, height = 6)











# Plot 2: Stacked by virulence factor type
p2 <- ggplot(virulence_summary, aes(x = Species, y = Count, fill = Virulence_Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Virulence Factor Composition by Species",
       y = "Number of Genes",
       fill = "Virulence Factor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
ggsave("plots/virulence_factors_stacked.png", p2, width = 10, height = 6)













# Plot 3: Heatmap style
p3 <- ggplot(virulence_summary, aes(x = Virulence_Factor, y = Species, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), color = "white", size = 4) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Virulence Factor Distribution Heatmap",
       x = "Virulence Factor",
       y = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)
ggsave("plots/virulence_factors_heatmap.png", p3, width = 8, height = 6)













# Plot 4: Faceted by virulence factor
p4 <- ggplot(virulence_summary, aes(x = Species, y = Count, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Virulence_Factor, scales = "free_y") +
  labs(title = "Virulence Factor Distribution Across Species",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

print(p4)
ggsave("plots/virulence_factors_faceted.png", p4, width = 12, height = 8)














# Export summary table
write.csv(virulence_summary, 
          "results/virulence_factors_summary.csv", 
          row.names = FALSE)

# Export detailed gene information
detailed_virulence <- data.frame()

for(species in names(virulence_results)) {
  species_results <- virulence_results[[species]]
  
  for(keyword in names(species_results)) {
    genes <- species_results[[keyword]]
    if(length(genes) > 0) {
      for(i in 1:length(genes)) {
        gene_info <- data.frame(
          Species = species,
          Virulence_Factor = keyword,
          Gene_ID = ifelse(!is.null(genes[i]$ID), genes[i]$ID, "NA"),
          Product = ifelse(!is.null(genes[i]$product), genes[i]$product, "NA"),
          Start = start(genes[i]),
          End = end(genes[i]),
          Strand = strand(genes[i])
        )
        detailed_virulence <- rbind(detailed_virulence, gene_info)
      }
    }
  }
}

# Export detailed data
write.csv(detailed_virulence, 
          "results/virulence_factors_detailed.csv", 
          row.names = FALSE)

cat("Exported results to CSV files:\n")
cat("- virulence_factors_summary.csv\n")
cat("- virulence_factors_detailed.csv\n")
