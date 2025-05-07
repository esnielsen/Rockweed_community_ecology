## Get RW-bot spp overlap from layer MARINe data (Jaccard similarity in presences as proxy for distributional overlap)

library(readxl)
library("xlsx")
library(dplyr)

# read table
photo.lay <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/LTM_layer_data/photolayerdata_20231218.xlsx", 
                        sheet = "photolayerraw_download")

# filter to remove non-living species
photo.all <- filter(photo.lay, !top_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))
photo.all <- filter(photo.all, !bottom_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET", "NONE"))

# Combine 'top_layer' and 'bottom_layer' into a single column of species per site
df <- photo.all %>%
  rowwise() %>%
  mutate(all_species = paste(top_layer, bottom_layer, sep = ",")) %>%
  ungroup()

# Split species into separate rows
df_long <- df %>%
  separate_rows(all_species, sep = ",") %>%
  distinct(ltm_latitude, ltm_longitude, all_species) # Remove duplicate species per site

## Create a wide-format dataframe with one column per species, showing pres/abs 
df_wide <- df_long %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = all_species, values_from = presence, values_fill = list(presence = 0))

# Define a function to calculate Jaccard similarity
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)  # Intersection: sites where both are present
  union <- sum(x | y)         # Union: sites where at least one is present
  if (union == 0) {
    return(NA)  # Avoid division by zero
  } else {
    return(intersection / union)
  }
}

# Columns to calculate Jaccard similarity with
target_columns <- c("SILCOM", "FUCGAR", "PELLIM")

# Initialize an empty list to store results
jaccard_list <- list()

# Loop through each target column
for (target in target_columns) {
  # Calculate Jaccard similarity for the current target column
  jaccard_results <- df_wide %>%
    select(-ltm_latitude, -ltm_longitude) %>%
    summarise(across(everything(), ~ jaccard_similarity(.x, !!sym(target)))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "species",
      values_to = "jaccard_value"
    ) %>%
    mutate(target_species = target)  # Add a column to identify the target species
  
  # Append to the list
  jaccard_list[[target]] <- jaccard_results
}

# Combine all results into a single long-format dataframe
jaccard_long_combined <- bind_rows(jaccard_list)
RW_spp_pres_cors <- as.data.frame(jaccard_long_combined)

# Convert the long-format dataframe to wide format
pres_cors_wide <- RW_spp_pres_cors %>%
  pivot_wider(
    names_from = target_species,    # Use 'target_species' values as new column names
    values_from = jaccard_value     # Fill the wide-format columns with Jaccard similarity values
  )

# Import bot spp DF - from biodiv.indices.R
bot_spp_cov_wide.table <- read.csv("~/bot_spp_cov_wide.table.csv")

# Combine DFs by spp code
merged_df <- pres_cors_wide %>%
  inner_join(bot_spp_cov_wide.table, by = c("species" = "species_code"))

# export table (to create final Table 1 for manuscript)
write.xlsx(merged_df, file = "bot_spp_pres_cov_corrs.xlsx")
