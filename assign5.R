# Install and Load Required Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, dplyr, tidyr, 
  Biostrings, DECIPHER, seqinr, rentrez, msa, 
  ape, phangorn, treeio, picante, 
  rgbif, sf, 
  vegan, cluster, factoextra, 
  future, furrr, parallel
)

# Set Seeds for Reproducibility
set.seed(123)  # For community analysis
set.seed(42)   # For phylogenetic analysis

# 1. GBIF Occurrence Data Retrieval
occ_search_results <- occ_search(scientificName = "Anolis", hasCoordinate = TRUE, limit = 500)
occurrence_data <- occ_search_results$data %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(scientificName)) %>%
  distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  dplyr::select(species, decimalLatitude, decimalLongitude, country, stateProvince, iucnRedListCategory)


# View the result
head(occurrence_data)

# 2. Sequence Data Retrieval
#Anolis_16S_search <- entrez_search(db = "nuccore", term = "Anolis AND 16S", retmax = 300)
#Anolis_16S_search1 <- entrez_search(db = "nuccore", term = "Anolis AND 16S AND 400:650[SLEN]", retmax = 300)

# Fetch sequences
#Anolis_16S_sequences <- entrez_fetch(db = "nuccore", id = #Anolis_16S_search1$ids, rettype = "fasta", retmode = "text")

# Print sequences to console
#cat(Anolis_16S_sequences, sep = "\n")

# Save sequences to FASTA file
#writeLines(Anolis_16S_sequences, "Anolis_16S_sequences.fasta")

# Load sequences as DNAStringSet
library(Biostrings)
Anolis_16S_summ <- readDNAStringSet("./data/Anolis_16S_sequences.fasta")
# Quality Control and Filtering of Sequences

# a. Remove sequences with 'N' (missing bases)
Anolis_16S_summ_cleaned <- Anolis_16S_summ[!grepl("N", Anolis_16S_summ)]

# b. Remove sequences with too many gaps (>50%)
max_gap_percentage <- 0.5
Anolis_16S_summ_cleaned <- Anolis_16S_summ_cleaned[sapply(Anolis_16S_summ_cleaned, function(seq) {
  sum(strsplit(as.character(seq), NULL)[[1]] == "-") / length(seq) < max_gap_percentage
})]

# c. Remove duplicate sequences
Anolis_16S_summ_cleaned <- unique(Anolis_16S_summ_cleaned)

# d. Remove sequences shorter than 100 bases
min_length <- 100
Anolis_16S_summ_cleaned <- Anolis_16S_summ_cleaned[nchar(Anolis_16S_summ_cleaned) > min_length]

# Create a data frame with cleaned sequences
dfAnolis_16S <- data.frame(Anolis16S_Title = names(Anolis_16S_summ_cleaned), 
                           Anolis_16S_Sequence = paste(Anolis_16S_summ_cleaned))
dfAnolis_16S$Species_Name <- word(dfAnolis_16S$Anolis16S_Title, 2L, 3L)

dfAnolis_16S <- dfAnolis_16S %>%
  mutate(Species_Name = gsub("^A\\.\\s*", "Anolis ", Species_Name)) %>%
  mutate(Species_Name = word(Species_Name, 1, 2)) %>%
  filter(grepl("^Anolis\\s\\w+", Species_Name)) %>%
  select(Anolis16S_Title, Species_Name, Anolis_16S_Sequence)



# Check the resulting dataframe
head(dfAnolis_16S)

dfAnolis_16S <- dfAnolis_16S[, c("Anolis16S_Title", "Species_Name", "Anolis_16S_Sequence")]

# Create the DNAStringSet from unique species
dna_sequences <- DNAStringSet(dfAnolis_16S$Anolis_16S_Sequence)
names(dna_sequences) <- dfAnolis_16S$Species_Name

# Remove duplicate species names
dfAnolis_16S_unique <- dfAnolis_16S[!duplicated(dfAnolis_16S$Species_Name), ]
dna_sequences_unique <- DNAStringSet(dfAnolis_16S_unique$Anolis_16S_Sequence)
names(dna_sequences_unique) <- dfAnolis_16S_unique$Species_Name

# Perform multiple sequence alignment
library(muscle)
aligned_sequences_unique <- muscle::muscle(dna_sequences_unique)

# Convert to phyDat format
library(phangorn)
aligned_phyDat_unique <- as.phyDat(aligned_sequences_unique)  

# Compute the distance matrix for phylogenetic tree construction
dist_matrix <- dist.ml(aligned_phyDat_unique, model = "JC69")

# Create a neighbor-joining tree
nj_tree <- nj(dist_matrix)

# Plot the NJ tree
plot(nj_tree, main = "Neighbor-Joining Phylogenetic Tree of Anolis species", 
     cex = 0.4, edge.width = 0.3, no.margin = TRUE)


# Maximum Likelihood Optimization
ml_tree <- pml(nj_tree, aligned_phyDat_unique)
ml_tree_optimized <- optim.pml(ml_tree, model = "GTR", optInv = TRUE, optGamma = TRUE)

# Print the summary of the optimized tree
summary(ml_tree_optimized)

# Plot the tree to visualize the phylogenetic relationships
#plot(ml_tree_optimized$tree)

# Set graphical parameters before plotting
par(mar = c(1, 1, 1, 1))  # Adjust margins for better spacing
plot(ml_tree_optimized$tree, cex = 0.4)  # Adjust 'cex' to scale tree label size

# Extract species names from both dataframes
species_occurrence <- unique(occurrence_data$species)  # Species from occurrence data
species_sequences <- unique(dfAnolis_16S_unique$Species_Name)  # Species from sequence data

# Find common species
common_species <- intersect(species_occurrence, species_sequences)

# Filter the dataframes to keep only rows with common species
occurrence_common <- occurrence_data %>% 
  filter(species %in% common_species)

sequence_common <- dfAnolis_16S_unique %>% 
  filter(Species_Name %in% common_species)

# Check the results
print(common_species)
head(occurrence_common)
head(sequence_common)

# Load necessary libraries
library(plotly)
library(RColorBrewer)

# Assuming 'occurrence_common' is your dataset with the relevant columns
# Check unique countries in your dataset
unique_countries <- unique(occurrence_common$country)

# Generate a color palette with rainbow colors to cover all countries
colors <- rainbow(length(unique_countries))

# Create a named vector for mapping country names to colors
country_colors <- setNames(colors, unique_countries)

# Assign a color to each row based on its country
occurrence_common$color <- country_colors[occurrence_common$country]

# Create the Plotly map
plot_ly(occurrence_common, 
        type = 'scattermapbox',  
        lat = ~decimalLatitude,  
        lon = ~decimalLongitude,  
        mode = 'markers',  
        marker = list(color = ~color, size = 8, opacity = 0.7),  # Use the 'color' column for marker color
        text = ~paste0("Species: ", species, "<br>",  
                       "Country: ", country, "<br>",  
                       "State/Province: ", stateProvince, "<br>",  
                       "Latitude: ", decimalLatitude, "<br>",  
                       "Longitude: ", decimalLongitude),  
        hoverinfo = 'text') %>%  
  layout(  
    mapbox = list(
      style = 'open-street-map', 
      zoom = 3,  # Adjust zoom level to focus on the Americas
      center = list(lat = 0, lon = -90),  # Center the map over the Americas
      bearing = 0,  # North up orientation
      pitch = 0  # No tilt
    ),  
    margin = list(l = 0, r = 0, t = 0, b = 0)  # Remove margins
  )


library(dplyr)
library(tidyr)
library(tibble)

# Combine occurrence data with phylogenetic information
community_matrix <- occurrence_common %>%
  group_by(country, species) %>%
  summarise(n = n(), .groups = "drop") %>%  # Drop grouping after summarizing
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "country")

# Check the resulting community matrix
head(community_matrix)


# Prepare data for phylogenetic community structure
# Combine occurrence data with phylogenetic information
library(picante)

# Ensure tree tip labels match community matrix species names
pruned_tree <- keep.tip(ml_tree_optimized$tree, colnames(community_matrix))

# Calculate phylogenetic community structure metrics
phylo_community_structure <- ses.mpd(community_matrix, cophenetic(pruned_tree))

# Remove rows with missing or invalid 'mpd.obs' values
phylo_community_structure_clean <- phylo_community_structure %>%
  filter(!is.na(mpd.obs) & mpd.obs >= 0)

# Visualize the phylogenetic community structure with the clean dat

ggplot(phylo_community_structure_clean, aes(x = reorder(row.names(phylo_community_structure_clean), mpd.obs), y = mpd.obs)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Mean Pairwise Phylogenetic Distance by Country",
       x = "Country", y = "MPD Observed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(mpd.obs, 2)), vjust = -0.3, size = 3)


library(ape)
library(picante)

# Method 1: Midpoint rooting
rooted_tree <- midpoint(pruned_tree)

# Verify the rooted tree
is.rooted(rooted_tree)

# Calculate Phylogenetic Diversity
pd_results <- pd(community_matrix, rooted_tree)

# Print results
print(pd_results)

# Visualization
library(ggplot2)

# Prepare data for plotting
pd_plot_data <- data.frame(
  Country = rownames(pd_results),
  PD = pd_results$PD,
  SR = pd_results$SR
)

# Create the plot
ggplot(pd_plot_data, aes(x = Country, y = PD)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  labs(title = "Phylogenetic Diversity by Country",
       x = "Country", 
       y = "Phylogenetic Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(PD, 2)), 
            vjust = -0.3, 
            size = 3)


# Make sure ggplot2 is loaded
library(ggplot2)
library(ggtree)
library(tidytree)
library(RColorBrewer)

# Convert tree to treedata object
tree_data <- as_tibble(rooted_tree)

# Create community annotation data
community_annotations <- data.frame(
  species = colnames(community_matrix),
  community = apply(community_matrix, 2, function(x) names(which.max(x)))
)

# Create color palette for communities
n_communities <- length(unique(community_annotations$community))
community_colors <- setNames(
  brewer.pal(n = max(3, min(n_communities, 9)), "Set3"),
  unique(community_annotations$community)
)

# Plot enhanced phylogenetic tree
p <- ggtree(rooted_tree) %<+% community_annotations +
  geom_tippoint(aes(color = community), size = 3) +
  scale_color_manual(values = community_colors) +
  theme_tree2() +
  ggplot2::labs(title = "Phylogenetic Tree with Community Distribution",
                color = "Community")

# Add scale bar
p <- p + geom_treescale()
plot(p)

# Compute ses.mpd once
ses_mpd_results <- ses.mpd(community_matrix, 
                           cophenetic(rooted_tree),
                           null.model = "taxa.labels",
                           abundance.weighted = FALSE,
                           runs = 999)

# Clean up and filter results
df_clean <- na.omit(ses_mpd_results)

# Create data frame for analysis
mpd_clean <- data.frame(
  Community = rownames(df_clean),
  MPD_Observed = df_clean$mpd.obs,
  MPD_Z_Score = df_clean$mpd.obs.z,
  MPD_P_Value = df_clean$mpd.obs.p
)

# Interpretation based on Z-scores
mpd_clean$Interpretation <- case_when(
  mpd_clean$MPD_P_Value < 0.05 & mpd_clean$MPD_Z_Score < -1.96 ~ "Phylogenetic Clustering",
  mpd_clean$MPD_P_Value < 0.05 & mpd_clean$MPD_Z_Score > 1.96 ~ "Phylogenetic Overdispersion",
  TRUE ~ "Random Phylogenetic Structure"
)

# Plotting the Z-scores
ggplot(mpd_clean, aes(x = reorder(Community, MPD_Z_Score), y = MPD_Z_Score, fill = Interpretation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("blue", "gray", "gray")) +
  labs(title = "Phylogenetic Community Structure",
       x = "Community", 
       y = "MPD Standardized Effect Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted", color = "red")
# Additional insights: Calculate summary statistics grouped by interpretation
summary_stats <- mpd_clean %>%
  group_by(Interpretation) %>%
  summarise(
    Count = n(),
    Mean_Z_Score = mean(MPD_Z_Score),
    Mean_P_Value = mean(MPD_P_Value)
  )

# Print summary stats
cat("Summary of Phylogenetic Community Structure:\n")
print(summary_stats)

