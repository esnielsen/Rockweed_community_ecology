##### Code to run PERMANOVA and SIMPER analyses on rockweed layer data


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
library(reshape2)
library(vegan)
library(RColorBrewer)
library(data.table)


photo.lay <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/LTM_layer_data/photolayerdata_20231218.xlsx", 
                        sheet = "photolayerraw_download")


# We have to swap top_layer record with bottom_layer 'NONE' values so we don't lose understory data
photo.lay$bottom_layer[photo.lay$bottom_layer == 'NONE'] <- photo.lay$top_layer[photo.lay$bottom_layer == 'NONE']

# edit again to say if top = bottom, call top layer as "NONE"
photo.lay <- photo.lay  %>% mutate(top_layer = if_else(bottom_layer==top_layer, 'NONE', top_layer))

#filter to remove non-living species in top layer
photo.lay <- filter(photo.lay, !top_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))

# remove non-living spp in bottom layer
photo.lay <- filter(photo.lay, !bottom_layer %in% c("UNIDEN",  "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))

# Filter again to only include target assemblages = RW
photo.lay <-filter(photo.lay, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis"))


# calculate per cov per quadrat (with counts that we divide by 100)
top_per_cov<-photo.lay %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(top_layer)
top_per_cov$n <- as.numeric(top_per_cov$n ) / 100

# do again for bottom layer spp
bot_per_cov<-photo.lay %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code, top_layer) %>%
  dplyr::count(bottom_layer)
bot_per_cov$n <- as.numeric(bot_per_cov$n ) / 100

bot_per_cov <- bot_per_cov %>% 
  rename("bot_cov" = "n")

# add top_layer % cov to the DF
top_bot_plots <- merge(top_per_cov, bot_per_cov)
top_bot_plots <- top_bot_plots %>% 
  rename("top_cov" = "n")

#filter to RW top cov
RW_top_cov_plot <-filter(top_bot_plots, top_layer %in%  c("SILCOM", "PELLIM", "FUCGAR"))


# filter to each RW spp
sil_diff_plots <-filter(RW_top_cov_plot, top_layer %in%  "SILCOM")
sil_diff_plots <-filter(sil_diff_plots, target_assemblage %in% "silvetia")

pel_diff_plots <-filter(RW_top_cov_plot, top_layer %in%  "PELLIM")
pel_diff_plots <-filter(pel_diff_plots, target_assemblage %in% "pelvetiopsis")

fuc_diff_plots <-filter(RW_top_cov_plot, top_layer %in% "FUCGAR")
fuc_diff_plots <-filter(fuc_diff_plots, target_assemblage %in%  "fucus")


#make wide
sil_diff_plot_w <- sil_diff_plots %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
#change NAs to zeros
sil_diff_plot_w <- sil_diff_plot_w %>% replace(is.na(.), 0)

# do for other spp

pel_diff_plot_w <- pel_diff_plots %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
pel_diff_plot_w <- pel_diff_plot_w %>% replace(is.na(.), 0)

fuc_diff_plot_w <- fuc_diff_plots %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
fuc_diff_plot_w <- fuc_diff_plot_w %>% replace(is.na(.), 0)


## permanova with region and only high/low cov
sil.manova <- adonis2(sil_diff_plot_w[, 11:70] ~ top_cov * georegion + marine_site_name, 
                      data = sil_diff_plot_w, 
                      permutations = 999, 
                      strata = sil_diff_plot_w$quadrat_code,
                      by = "terms")
sil.manova

pel.manova <- adonis2(pel_diff_plot_w[ , 11:44] ~ top_cov * georegion + marine_site_name, data = pel_diff_plot_w, permutations = 999, strata = pel_diff_plot_w$quadrat_code, by = "terms")
pel.manova

fuc.manova <- adonis2(fuc_diff_plot_w[ , 11:45] ~ top_cov * georegion + marine_site_name, data = fuc_diff_plot_w, permutations = 999, strata = fuc_diff_plot_w$quadrat_code, by = "terms")
fuc.manova


# filter to spp with > 5% cover avg per RW spp 

bot_spp_cov_wide.table <- read.csv("C:/Users/erica.nielsen/Desktop/Synz/MARINe/PCI_MARINe_synth/Fucoids/layer_DSR_analyses/bot_spp_cov_wide.table.csv")

# sil filt
sil_filt <- bot_spp_cov_wide.table %>%
  filter(S..compressa > 4) %>%
  select(1)  # Selects the first column
# filter wide DF by filt DF
sil_diff_filt <- sil_diff_plots %>%
  filter(bottom_layer %in% sil_filt$species_code)

# fuc filt
fuc_filt <- bot_spp_cov_wide.table %>%
  filter(F..garderni > 4) %>%
  select(1)  # Selects the first column
fuc_diff_filt <- fuc_diff_plots %>%
  filter(bottom_layer %in% fuc_filt$species_code)

# pel filt
pel_filt <- bot_spp_cov_wide.table %>%
  filter(P..limitata > 4) %>%
  select(1)  # Selects the first column
pel_diff_filt <- pel_diff_plots %>%
  filter(bottom_layer %in% pel_filt$species_code)


# make wide again
sil_diff_filt_w <- sil_diff_filt %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
sil_diff_filt_w <- sil_diff_filt_w %>% replace(is.na(.), 0)

pel_diff_filt_w <- pel_diff_filt %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
pel_diff_filt_w <- pel_diff_filt_w %>% replace(is.na(.), 0)

fuc_diff_filt_w <- fuc_diff_filt %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
fuc_diff_filt_w <- fuc_diff_filt_w %>% replace(is.na(.), 0)

# run permanova (THESE NOT FINAL RESULTS - used above!!)
sil.f.manova <- adonis2(sil_diff_filt_w[ , 11:21] ~ top_cov * georegion, data = sil_diff_filt_w, permutations = 999, by = "terms")
sil.f.manova

pel.f.manova <- adonis2(pel_diff_filt_w[ , 11:19] ~ top_cov * georegion, data = pel_diff_filt_w, permutations = 999, by = "terms")
pel.f.manova

fuc.f.manova <- adonis2(fuc_diff_filt_w[ , 11:23] ~ top_cov * georegion, data = fuc_diff_filt_w, permutations = 999, by = "terms")
fuc.f.manova


#### perform simper
#######################################################################################################################

## create groups for categorical predictor

#make new column to group RW cov into high/low (above/below 50% cov)
sil_diff_plot_w <- sil_diff_plot_w %>%
mutate(cov_group = case_when(top_cov>=0.67~"high_cover", between(top_cov, 0.33, 0.67)~"med_cover",
top_cov<=0.32~"low_cover"))

fuc_diff_plot_w <- fuc_diff_plot_w %>%
  mutate(cov_group = case_when(top_cov>=0.67~"high_cover", between(top_cov, 0.33, 0.67)~"med_cover",
                               top_cov<=0.32~"low_cover"))

pel_diff_plot_w <- pel_diff_plot_w %>%
  mutate(cov_group = case_when(top_cov>=0.67~"high_cover", between(top_cov, 0.33, 0.67)~"med_cover",
                               top_cov<=0.32~"low_cover"))


# turn spp df into matrix
mat_sil_spp <- as.matrix(sil_diff_plot_w[ , 12:70])

# run simper, with RW cov group as pred var
simp_sil <- simper(mat_sil_spp, sil_diff_plot_w$cov_group)

# get spp contrib to low vs high cov
spp_cont_sil <- as.data.frame(simp_sil$high_cover_low_cover)

# create abund diff column
spp_cont_sil <- spp_cont_sil %>%
  mutate(abund.diff = (ava - avb) * 100)

## PEL SIMPER
# 1st turn spp df into matrix
mat_pel_spp <- as.matrix(pel_diff_plot_w[ , 11:44])

# run simper, with RW cov group as pred var
simp_pel <- simper(mat_pel_spp, pel_diff_plot_w$cov_group)

# get spp contrib to low vs high cov
spp_cont_pel <- as.data.frame(simp_pel$low_cover_high_cover)

# create abund diff column
spp_cont_pel <- spp_cont_pel %>%
  mutate(abund.diff = (avb - ava) * 100)


## fuc SIMPER
# 1st turn spp df into matrix
mat_fuc_spp <- as.matrix(fuc_diff_plot_w[ , 11:45])

# run simper, with RW cov group as pred var
simp_fuc <- simper(mat_fuc_spp, fuc_diff_plot_w$cov_group)

# get spp contrib to low vs high cov
spp_cont_fuc <- as.data.frame(simp_fuc$low_cover_high_cover)

# create abund diff column
spp_cont_fuc <- spp_cont_fuc %>%
  mutate(abund.diff = (avb - ava) * 100)

all_spp_simp <- bind_rows(lst(spp_cont_sil, spp_cont_pel, spp_cont_fuc), .id = 'id')

write.csv(all_spp_simp, "all_spp_simp.NEW.csv")


## plot simper results
# bar graph

# set RW spp order
filt_spp_SIMP <- filt_spp_SIMP %>%
  mutate(rockweed = factor(rockweed, levels = c("S. compressa", "F. distichus", "P. limitata")))

library(cowplot)

# Find the max y-value to ensure consistent y-axis limits
max_y <- max(filt_spp_SIMP$abund.diff)

# Create separate plots while ensuring the same y scale and removing extra legends
plots_list <- filt_spp_SIMP %>%
  split(.$rockweed) %>%
  lapply(function(df) {
    ggplot(df, aes(x = reorder(species, -average), y = abund.diff, fill = species)) +
      geom_bar(stat = "identity") +
      labs(x = "Understory species",
           y = "Difference in % cover",
           title = unique(df$rockweed)) +
      theme_minimal() +
      scale_fill_manual(values = mypalette) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ylim(min(filt_spp_SIMP$abund.diff), max(filt_spp_SIMP$abund.diff)) + 
      guides(fill = "none") # Remove individual legends
  })

# Combine plots and attach the legend
final_plot <- plot_grid(plotlist = plots_list, ncol = 3)

# Save the final plot
save_plot("simper.plot.pdf", final_plot, base_width = 15, base_height = 7)

