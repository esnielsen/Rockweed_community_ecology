### Code to analyze and plotsignficant RW ~ understory spp correlations

library(readxl)
library(ggplot2)
library(dplyr) 
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(maps)
library(viridis)
library(pglm)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(Polychrome)
library(pals)
library(patchwork)

# read raw data
photo.lay <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/LTM_layer_data/photolayerdata_20231218.xlsx", 
                        sheet = "photolayerraw_download")

# filter to only living organisms in both top and bottom layer
photo.all <- filter(photo.lay, !top_layer %in% c("UNIDEN", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))
photo.all <- filter(photo.all, !bottom_layer %in% c("UNIDEN", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET", "NONE", "ROCK", "SAND"))

# filter to RW spp
photo.all <-filter(photo.all, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis"))

# We have to swap top_layer record with bottom_layer 'NONE' values so we don't lose understory data
photo.all$bottom_layer[photo.all$bottom_layer == 'NONE'] <- photo.all$top_layer[photo.all$bottom_layer == 'NONE']

# edit again to say if top = bottom, call top layer as "NONE"
photo.all <- photo.all  %>% mutate(top_layer = if_else(bottom_layer==top_layer, 'NONE', top_layer))

### swap bot and top layers for those with RW bot layer (so we only have RW as top layer)

# Step 1: Identify rows to be edited before modification
rows_to_edit_before <- photo.all %>%
  filter(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL"))

# Print rows before editing
cat("Rows before editing:\n")
print(rows_to_edit_before)

# Step 2: Apply the modifications
photo.all <- photo.all %>%
  mutate(
    # Swap bottom_layer and top_layer where conditions are met
    new_bottom_layer = ifelse(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL"),
                              top_layer, bottom_layer),
    new_top_layer = ifelse(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL"),
                           bottom_layer, top_layer)
  ) %>%
  select(-bottom_layer, -top_layer) %>%
  rename(bottom_layer = new_bottom_layer, top_layer = new_top_layer)

# Step 3: Identify rows that were edited after modification
rows_to_edit_after <- photo.all %>%
  filter(top_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL") & bottom_layer == "NONE")

# calc % cov per top spp
top_per_cov<-photo.all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(top_layer)
top_per_cov$RW_cov <- as.numeric(top_per_cov$n ) / 100 

# calc % cov per bot spp
bot_per_cov<-photo.all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(bottom_layer)
bot_per_cov$bot_cov <- as.numeric(bot_per_cov$n ) / 100

# add identifier colums to combine DFS
top_per_cov$site_quad_sp <- paste(top_per_cov$marine_site_name, "_", top_per_cov$quadrat_code, "_", top_per_cov$marine_common_year)
bot_per_cov$site_quad_sp <- paste(bot_per_cov$marine_site_name, "_", bot_per_cov$quadrat_code, "_", bot_per_cov$marine_common_year)

#combine into single DF (this now changed to include spp as to_layer column - use this in Loop!)
plot_df_list <- list(bot_per_cov, top_per_cov)
plot_spp_cov <- plot_df_list %>% reduce(full_join, by='site_quad_sp')

# remove rows with non-RW top_layer
plot_spp_cov <-filter(plot_spp_cov, top_layer %in%  c("SILCOM", "PELLIM", "FUCGAR"))


##### create plot of each understory spp with signif interaction (% cov regressed to % cov of RW layer)

# Set colors
mypalette <- createPalette(75, c("#ff0000", "#010101", "#00ff00", "#0000ff"))
names(mypalette) <- sort(unique(plot_spp_cov$bottom_layer))

# Define function to fit linear model and tidy it
fit_and_tidy_model <- function(data) {
  model <- glm(botcov2 ~ RW_cov, family = beta_family(link = "logit"), data = data) 
  tidy_model <- broom::tidy(model)
  return(tidy_model)
}

# Function to create plot for each species
create_species_plot <- function(df, species_name) {
  df <- df %>%
    mutate(botcov2 = replace(bot_cov, bot_cov == 1, 0.99),
           marine_common_year.x = as.factor(marine_common_year.x))
  
  # Fit models and filter significant ones
  fitted_models <- df %>%
    group_by(bottom_layer) %>%
    do(model = fit_and_tidy_model(.))
  
  sig_models <- fitted_models %>%
    filter(model$p.value[2] < 0.05)
  
  filt_data <- df %>%
    filter(bottom_layer %in% sig_models$bottom_layer)
  
  # Create ggplot
  plot <- ggplot(filt_data, aes(x = RW_cov, y = botcov2, color = bottom_layer)) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Rockweed cover", y = "Understory cover", title = species_name) +
    scale_color_manual(values = mypalette) +
    coord_cartesian(ylim = c(0, 0.35)) +  # Set fixed y-axis limits
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "italic")
    )
  
  return(plot)
}

# Combine data for shared legend
sil_data <- plot_spp_cov %>% filter(top_layer == 'SILCOM') %>% mutate(species = "S. compressa")
pel_data <- plot_spp_cov %>% filter(top_layer == 'PELLIM') %>% mutate(species = "P. limitata")
fuc_data <- plot_spp_cov %>% filter(top_layer == 'FUCGAR') %>% mutate(species = "F. distichus")

# Combine all species data into one dataframe
combined_data <- bind_rows(sil_data, pel_data, fuc_data)

# Create individual plots
sil_plot <- create_species_plot(sil_data, "S. compressa")
pel_plot <- create_species_plot(pel_data, "P. limitata")
fuc_plot <- create_species_plot(fuc_data, "F. distichus")

# Combine plots into a single figure with a shared legend
combined_plot <- (sil_plot | fuc_plot | pel_plot) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8), # Reduce legend text size
    legend.key.size = unit(0.5, "lines")  # Reduce legend key size
  )

# Adjust legend to have multiple rows
combined_plot <- combined_plot & guides(color = guide_legend(nrow = 3))

# Add a title to the combined plot
combined_plot <- combined_plot + 
  plot_annotation(title = "Species-specific Understory Cover by Rockweed Cover")

# Save the combined plot to a PDF
ggsave("combined_species_intxn_Yaxes_F.dist.pdf", combined_plot, width = 12, height = 8)

