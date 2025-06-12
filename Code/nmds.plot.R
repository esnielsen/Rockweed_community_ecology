### code to run NMDS plots assessing understory change with RW cover change

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
library(tibble)
library(ggrepel)


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
 
# filter to spp with > 5% cover avg per RW (bc then we are looking at how main players in community change w/RW)

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

#make wide
sil_diff_plot_w <- sil_diff_filt %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
#change NAs to zeros
sil_diff_plot_w <- sil_diff_plot_w %>% replace(is.na(.), 0)

# do for other spp

pel_diff_plot_w <- pel_diff_filt %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
pel_diff_plot_w <- pel_diff_plot_w %>% replace(is.na(.), 0)

fuc_diff_plot_w <- fuc_diff_filt %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
fuc_diff_plot_w <- fuc_diff_plot_w %>% replace(is.na(.), 0)



#### Nmds with all points

# Select only species cover columns (columns 11-21)
species_data <- sil_diff_filt_w[, 11:21]

# Perform NMDS
nmds_sil <- metaMDS(comm = sil_diff_filt_w[ , 11:21],  # Define the community data 
                    distance = "bray",       # Specify a bray-curtis distance
                    try = 100)               # Number of iterations

#get 'stress', lower is better (<0.05 means good data representation)
nmds_sil$stress 

# Get NMDS coordinates
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
nmds_scores = as.data.frame(scores(nmds_sil)$sites)

#add columns to data frame 
nmds_scores$RW_cover = sil_diff_filt_w$top_cov
nmds_scores$Region = sil_diff_filt_w$georegion
nmds_scores$Site = sil_diff_filt_w$marine_site_name

# Fit species vectors
species_fit <- envfit(nmds_sil, species_data)

species_scores <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores$species <- rownames(species_scores)

# Create NMDS plot
sil.nmds.plot<-ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, size = RW_cover, color = Region, alpha = 0.5)) +
  geom_segment(data = species_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), fontface = "bold") +
  mdthemes::md_theme_classic() + theme(plot.margin = unit(c(6,0,6,0), "pt"))+
  labs(title = "*S. compressa*", x = "NMDS1", y = "NMDS2") +
  scale_color_manual(values = c("CA Central" = "darkorange", "CA South" = "magenta"), name = "Region") +
  guides(size = guide_legend(title = "Rockweed Cover"), color = guide_legend(title = "Georegion"))+guides(alpha = "none")


pdf("sil.nmds.NEW.unfilt.pdf")
sil.nmds.plot
dev.off()

## do for Fuc 
##

# Select only species cover columns (columns 11-21)
species_data <- fuc_diff_filt_w[, 11:24]

# Perform NMDS
nmds_fuc <- metaMDS(comm = fuc_diff_filt_w[ , 11:24],  # Define the community data 
                    distance = "bray",       # Specify a bray-curtis distance
                    try = 100)               # Number of iterations

#get 'stress', lower is better (<0.05 means good data representation)
nmds_fuc$stress 

# Get NMDS coordinates
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
nmds_scores = as.data.frame(scores(nmds_fuc)$sites)

#add columns to data frame 
nmds_scores$RW_cover = fuc_diff_filt_w$top_cov
nmds_scores$Region = fuc_diff_filt_w$georegion
nmds_scores$Site = fuc_diff_filt_w$marine_site_name

# Fit species vectors
species_fit <- envfit(nmds_fuc, species_data)

species_scores <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores$species <- rownames(species_scores)

# Create NMDS plot
fuc.nmds.plot<-ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, size = RW_cover, color = Region, alpha = 0.5)) +
  geom_segment(data = species_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), fontface = "bold") +
  mdthemes::md_theme_classic() + theme(plot.margin = unit(c(6,0,6,0), "pt"))+
  labs(title = "*F. distichus*", x = "NMDS1", y = "NMDS2") +
  scale_color_manual(values = c("CA Central" = "darkorange", "CA North" = "lightseagreen"), name = "Region") +
  guides(size = guide_legend(title = "Rockweed Cover"), color = guide_legend(title = "Georegion"))+guides(alpha = "none")


pdf("fuc.nmds.NEW.unfilt.pdf")
fuc.nmds.plot
dev.off()


### do for pellim
##


# Select only species cover columns (columns 11-21)
species_data <- pel_diff_filt_w[, 11:19]

# Perform NMDS
nmds_pel <- metaMDS(comm = pel_diff_filt_w[ , 11:19],  # Define the community data 
                    distance = "bray",       # Specify a bray-curtis distance
                    try = 100)               # Number of iterations

#get 'stress', lower is better (<0.05 means good data representation)
nmds_pel$stress 

# Get NMDS coordinates
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
nmds_scores = as.data.frame(scores(nmds_pel)$sites)

#add columns to data frame 
nmds_scores$RW_cover = pel_diff_filt_w$top_cov
nmds_scores$Region = pel_diff_filt_w$georegion
nmds_scores$Site = pel_diff_filt_w$marine_site_name

# Fit species vectors
species_fit <- envfit(nmds_pel, species_data)

species_scores <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores$species <- rownames(species_scores)

# Create NMDS plot
pel.nmds.plot<-ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, size = RW_cover, color = Region, alpha = 0.5)) +
  geom_segment(data = species_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), fontface = "bold") +
  mdthemes::md_theme_classic() + theme(plot.margin = unit(c(6,0,6,0), "pt"))+
  labs(title = "*P. limitata*", x = "NMDS1", y = "NMDS2") +
  scale_color_manual(values = c("CA Central" = "darkorange", "CA North" = "lightseagreen"), name = "Region") +
  guides(size = guide_legend(title = "Rockweed Cover"), color = guide_legend(title = "Georegion"))+guides(alpha = "none")


pdf("pel.nmds.NEW.unfilt.pdf")
pel.nmds.plot
dev.off()


## combine into single plot

library(cowplot)

# Extract legends
legend1 <- get_legend(sil.nmds.plot + theme(legend.position="right"))
legend2 <- get_legend(fuc.nmds.plot + theme(legend.position="right"))

# Combine legends into one grid
combined_legend <- plot_grid(legend1, legend2, ncol = 1)

# arrange the three plots in a single row
prow <- plot_grid( sil.nmds.plot + theme(legend.position="none"),
                   fuc.nmds.plot + theme(legend.position="none"),
                   pel.nmds.plot + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("a)", "b)", "c)"),
                   hjust = -1,
                   nrow = 1
)

final_plot <- plot_grid(prow, combined_legend, ncol = 1, rel_heights = c(3, 1))
final_plot

save_plot("combined.nmds.legend.pdf", final_plot, base_width = 15, base_height = 7)


