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

# read dat
photo.lay <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/LTM_layer_data/photolayerdata_20231218.xlsx", 
    sheet = "photolayerraw_download")

# filter to only include sites with at least 7 years sampled

#get number years sampled
phot_count <- photo.lay %>%                           
 group_by(marine_site_name) %>%
 dplyr::summarise(count = n_distinct(marine_common_year))

photo.yr<-merge(photo.lay,  phot_count, by.x = "marine_site_name", by.y = "marine_site_name")

# filtering to include at least 7 years of sampling, calling 'NorCal' because only those sites pass filters
photo.NorCal <- photo.yr %>% filter(count > 7)

# We have to swap top_layer record with bottom_layer 'NONE' values so we don't lose understory data
photo.NorCal$bottom_layer[photo.NorCal$bottom_layer == 'NONE'] <- photo.NorCal$top_layer[photo.NorCal$bottom_layer == 'NONE']

# edit again to say if top = bottom, call top layer as "NONE"
photo.NorCal <- photo.NorCal  %>% mutate(top_layer = if_else(bottom_layer==top_layer, 'NONE', top_layer))

#filter to remove non-living species in top layer
photo.NorCal <- filter(photo.NorCal, !top_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))

# just filter unidentified spp in bottom layer
photo.NorCal <- filter(photo.NorCal, !bottom_layer %in% c("UNIDEN"))

# Filter again to only include target assemblages = RW
photo.NorCal <-filter(photo.NorCal, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis", "hesperophycus"))


# calculate per cov per quadrat (with counts that we divide by 100)
top_per_cov<-photo.NorCal %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(top_layer)
top_per_cov$n <- as.numeric(top_per_cov$n ) / 100

# do again for bottom layer spp
bot_per_cov<-photo.NorCal %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, season_name, target_assemblage, quadrat_code, top_layer) %>%
  dplyr::count(bottom_layer)
bot_per_cov$n <- as.numeric(bot_per_cov$n ) / 100

## filter to plots with > 70% change over the years
# filt top_cov to only include RW of interest
RW_top_cov_plot <-filter(top_per_cov, top_layer %in%  c("SILHES", "SILCOM", "PELLIM", "PELHYB", "FUCGAR", "HESCAL"))

# get diffs btwn the highest and lowest % covs per site/plot
RW_cov_diff <- RW_top_cov_plot %>% 
    group_by(marine_site_name,quadrat_code, target_assemblage, top_layer) %>%
    summarise(Diff = n[which.max(n)] - n[which.min(n)])

# sort to get sites/plot per spp w highest diffs
d <- data.table(RW_cov_diff, key="Diff")
d2 <- d[, head(.SD, 3), by=Diff]
View(d2)

## READ:
# I then looked at the bottom of the 'd2' dataframe and took the 5 sites per spp with highest difference in % cov over time and copied those into an excel sheet titled 'top_cov_diff_plots

top_cov_diff_plots <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/PCI_MARINe_synth/Fucoids/layer_DSR_analyses/top_cov_diff_plots.xlsx")

# remove white space of site names for join later
bot_per_cov$marine_site_name <- gsub("\\s+", "_", bot_per_cov$marine_site_name)
top_per_cov$marine_site_name <- gsub("\\s+", "_", top_per_cov$marine_site_name)

# filter bot_cov by plots of interest
bot_diff_plots <- merge(bot_per_cov, top_cov_diff_plots)
bot_diff_plots <- bot_diff_plots %>% 
       rename("bot_cov" = "n")

# add top_layer % cov to the DF
top_diff_plots <- merge(top_per_cov, top_cov_diff_plots)
top_diff_plots <- top_diff_plots %>% 
       rename("top_cov" = "n")
diff_plots_DF <- inner_join(bot_diff_plots, top_diff_plots, relationship = "many-to-many")

#make new column to group RW cov into high/low (above/below 50% cov)
diff_plots_DF <- diff_plots_DF %>%
      mutate(cov_group = case_when(top_cov>=0.5~"high_cover",
                                     top_cov<0.5~"low_cover"))

# filter to each RW spp
sil_diff_plots <-filter(diff_plots_DF, top_layer %in%  "SILCOM")
sil_diff_plots <-filter(sil_diff_plots, target_assemblage %in% "silvetia")
                                                                  
pel_diff_plots <-filter(diff_plots_DF, top_layer %in%  "PELLIM")
pel_diff_plots <-filter(pel_diff_plots, target_assemblage %in% "pelvetiopsis")

fuc_diff_plots <-filter(diff_plots_DF, top_layer %in% "FUCGAR")
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

# 
nmds_sil <- metaMDS(comm = sil_diff_plot_w[ , 12:51],  # Define the community data 
                        distance = "bray",       # Specify a bray-curtis distance
                        try = 100)               # Number of iterations

#get 'stress', lower is better (<0.05 means good data representation)
nmds_sil$stress 

#get coordinates
nmds_sil$points %>% head()

#save coords
nmds_sil.pts <- data.frame(nmds_sil$points)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
sil.scores = as.data.frame(scores(nmds_sil)$sites)

#add columns to data frame 
sil.scores$RW_cover = sil_diff_plot_w$top_cov
sil.scores$Year = sil_diff_plot_w$marine_common_year
sil.scores$Site = sil_diff_plot_w$marine_site_name


### Get metadata ready
sil.env <- sil.scores[, c("RW_cover", "Site")]

sil.envfit <- envfit(nmds_sil, sil.env, permutations = 999)

site.scrs <- as.data.frame(scores(nmds_sil, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Site = sil.scores$Site) #add site names as variable if you want to display on plot

head(site.scrs)

### Get spp data ready
sil<-sil_diff_plot_w[ , 12:51]
sil.spp.fit <- envfit(nmds_sil, sil, permutations = 999) # this fits species vectors

  
spp.scrs <- as.data.frame(scores(sil.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = sil.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

head(spp.scrs)

### make env dataset
env.scores <- as.data.frame(scores(sil.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names

env.scores <- cbind(env.scores, pval = sil.envfit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

head(sig.env.scrs)

### Get spp data ready
sil<-sil_diff_plot_w[ , 12:51]
sil.spp.fit <- envfit(nmds_sil, sil, permutations = 999) # this fits species vectors


### make env dataset
env.scores <- as.data.frame(scores(sil.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names

env.scores <- cbind(env.scores, pval = sil.envfit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

### make spp dataset
ord <- metaMDS(sil)
ord

species.scores <- as.data.frame(scores(ord, "species"))
species.scores$Species <- rownames(species.scores)

# Get spp taxonomy 
bot_spp_groups <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/PCI_MARINe_synth/Fucoids/layer_DSR_analyses/bot_spp_groups.xlsx")
bot_spp_groups <- na.omit(bot_spp_groups)

# Merge species scores with taxa information
botgroup <-bot_spp_groups %>% remove_rownames %>% column_to_rownames(var="species")
spps_score_group <- merge(species.scores, botgroup, by = "row.names")  # Ensure correct merge key

# Extract NMDS site and species scores
site.scores <- as.data.frame(scores(ord, "sites"))
site.scores$Site <- rownames(site.scores)

species.scores <- as.data.frame(scores(ord, "species"))
species.scores$Species <- rownames(species.scores)

# Merge species scores with taxa information
spps_score_group <- merge(species.scores, botgroup, by = "row.names") 

### Create plot

# Define custom taxa colors
okabe <- c("tan4", "mediumseagreen","purple4", "darkolivegreen", "red3", "grey30" )
names(okabe) <- sort(unique(bot_spp_groups$taxa))  # Assign colors to taxa names

# Add RW_cover to site scores
site.scores$RW_cover <- sil.env$RW_cover

# Generate contour data for RW_cover using kde2d
kde <- with(site.scores, MASS::kde2d(NMDS1, NMDS2, h = c(0.5, 0.5), n = 100))  # Adjust 'h' for smoothing
contour_data <- data.frame(expand.grid(x = kde$x, y = kde$y), z = as.vector(kde$z))

# Generate plot
# Plot NMDS with ggplot
sil.nmds.plot <- ggplot() +
  # Contour lines for RW_cover
  geom_contour_filled(data = contour_data, aes(x = x, y = y, z = z), alpha = 0.5) +
  scale_fill_viridis_d(name = "RW Cover", option = "plasma", direction = -1) + # Plasma color scale for contours
  
  # Add NMDS site points (grey, with transparency)
  geom_point(data = site.scores, aes(x = NMDS1, y = NMDS2), color = "grey", alpha = 0.6, size = 3) +
  
  # New color scale for taxa
  ggnewscale::new_scale_color() +
  
  # Add species labels (bold & colored by taxa)
  geom_text_repel(data = spps_score_group, 
                  aes(x = NMDS1, y = NMDS2, label = Species, color = taxa), 
                  size = 4, fontface = "bold", box.padding = 0.3, max.overlaps = 20)+
  scale_color_manual(name = "Taxa", values = okabe) +  # Apply custom colors
  
  # Equal aspect ratio for NMDS
  coord_fixed() +
  
 # Theme and labels
  mdthemes::md_theme_classic() + theme(plot.margin = unit(c(6,0,6,0), "pt"))+
  labs(title = "*S. compressa*", x = "NMDS1", y = "NMDS2")

#pdf("sil.nmds.NEW.pdf")
sil.nmds.plot
#dev.off()

## Do for Fucus
# run nmds
nmds_fuc <- metaMDS(comm = fuc_diff_plot_w[ , 12:38],  # Define the community data 
                    distance = "bray",       # Specify a bray-curtis distance
                    try = 100)               # Number of iterations

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
fuc.scores = as.data.frame(scores(nmds_fuc)$sites)

#add columns to data frame 
fuc.scores$RW_cover = fuc_diff_plot_w$top_cov
fuc.scores$Year = fuc_diff_plot_w$marine_common_year
fuc.scores$Site = fuc_diff_plot_w$marine_site_name

### Do spp and env fits
fuc <-fuc_diff_plot_w[ , 12:38]
fuc.env <- fuc.scores[, c("RW_cover", "Site")]

fuc.spp.fit <- envfit(nmds_fuc, fuc, permutations = 999) # this fits species vectors
fuc.envfit <- envfit(nmds_fuc, fuc.env, permutations = 999)

### make env dataset
env.scores <- as.data.frame(scores(fuc.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names

env.scores <- cbind(env.scores, pval = fuc.envfit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

### make spp dataset
ord <- metaMDS(fuc)
ord

species.scores <- as.data.frame(scores(ord, "species"))
species.scores$Species <- rownames(species.scores)

# Merge species scores with taxa information
botgroup <-bot_spp_groups %>% remove_rownames %>% column_to_rownames(var="species")
spps_score_group <- merge(species.scores, botgroup, by = "row.names")  # Ensure correct merge key

# Extract NMDS site and species scores
site.scores <- as.data.frame(scores(ord, "sites"))
site.scores$Site <- rownames(site.scores)

species.scores <- as.data.frame(scores(ord, "species"))
species.scores$Species <- rownames(species.scores)

# Merge species scores with taxa information
spps_score_group <- merge(species.scores, botgroup, by = "row.names") 

# Define custom taxa colors
okabe <- c("tan4", "mediumseagreen","purple4", "darkolivegreen", "red3", "grey30" )
names(okabe) <- sort(unique(bot_spp_groups$taxa))  # Assign colors to taxa names

# Add RW_cover to site scores
site.scores$RW_cover <- fuc.env$RW_cover

# Generate contour data for RW_cover using kde2d
kde <- with(site.scores, MASS::kde2d(NMDS1, NMDS2, h = c(0.5, 0.5), n = 100))  # Adjust 'h' for smoothing
contour_data <- data.frame(expand.grid(x = kde$x, y = kde$y), z = as.vector(kde$z))

# Generate plot
# Plot NMDS with ggplot
fuc.nmds.plot <- ggplot() +
  # Contour lines for RW_cover
  geom_contour_filled(data = contour_data, aes(x = x, y = y, z = z), alpha = 0.5) +
  scale_fill_viridis_d(name = "RW Cover", option = "plasma", direction = -1) + # Plasma color scale for contours
  
  # Add NMDS site points (grey, with transparency)
  geom_point(data = site.scores, aes(x = NMDS1, y = NMDS2), color = "grey", alpha = 0.6, size = 3) +
  
  # New color scale for taxa
  ggnewscale::new_scale_color() +
  
  # Add species labels (bold & colored by taxa)
  geom_text_repel(data = spps_score_group, 
                  aes(x = NMDS1, y = NMDS2, label = Species, color = taxa), 
                  size = 4, fontface = "bold", box.padding = 0.3, max.overlaps = 20)+
  scale_color_manual(name = "Taxa", values = okabe) +  # Apply custom colors
  
  # Equal aspect ratio for NMDS
  coord_fixed() +
  
  # Theme and labels
  mdthemes::md_theme_classic() + theme(plot.margin = unit(c(6,0,6,0), "pt"))+
  labs(title = "*F. gardneri*", x = "NMDS1", y = "NMDS2")

#pdf("fuc.nmds.NEW.pdf")
fuc.nmds.plot
#dev.off()


### Do for Pellim


# run nmds
nmds_pel <- metaMDS(comm = pel_diff_plot_w[ , 12:28],  # Define the community data 
                    distance = "bray",       # Specify a bray-curtis distance
                    try = 100)               # Number of iterations

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
pel.scores = as.data.frame(scores(nmds_pel)$sites)

#add columns to data frame 
pel.scores$RW_cover = pel_diff_plot_w$top_cov
pel.scores$Year = pel_diff_plot_w$marine_common_year
pel.scores$Site = pel_diff_plot_w$marine_site_name

### Do spp and env fits
pel <-pel_diff_plot_w[ , 12:28]
pel.env <- pel.scores[, c("RW_cover", "Site")]

pel.spp.fit <- envfit(nmds_pel, pel, permutations = 999) # this fits species vectors
pel.envfit <- envfit(nmds_pel, pel.env, permutations = 999)

### make env dataset
env.scores <- as.data.frame(scores(pel.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names

env.scores <- cbind(env.scores, pval = fuc.envfit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

### make spp dataset
ord <- metaMDS(pel)
ord

species.scores <- as.data.frame(scores(ord, "species"))
species.scores$Species <- rownames(species.scores)

# Merge species scores with taxa information
botgroup <-bot_spp_groups %>% remove_rownames %>% column_to_rownames(var="species")
spps_score_group <- merge(species.scores, botgroup, by = "row.names")  # Ensure correct merge key

# Extract NMDS site and species scores
site.scores <- as.data.frame(scores(ord, "sites"))
site.scores$Site <- rownames(site.scores)

species.scores <- as.data.frame(scores(ord, "species"))
species.scores$Species <- rownames(species.scores)

# Merge species scores with taxa information
spps_score_group <- merge(species.scores, botgroup, by = "row.names") 

# Define custom taxa colors
okabe <- c("tan4", "mediumseagreen","purple4", "darkolivegreen", "red3", "grey30" )
names(okabe) <- sort(unique(bot_spp_groups$taxa))  # Assign colors to taxa names

# Add RW_cover to site scores
site.scores$RW_cover <- pel.env$RW_cover

# Generate contour data for RW_cover using kde2d
kde <- with(site.scores, MASS::kde2d(NMDS1, NMDS2, h = c(0.5, 0.5), n = 100))  # Adjust 'h' for smoothing
contour_data <- data.frame(expand.grid(x = kde$x, y = kde$y), z = as.vector(kde$z))

# Generate plot
# Plot NMDS with ggplot
pel.nmds.plot <- ggplot() +
  # Contour lines for RW_cover
  geom_contour_filled(data = contour_data, aes(x = x, y = y, z = z), alpha = 0.5) +
  scale_fill_viridis_d(name = "RW Cover", option = "plasma", direction = -1) + # Plasma color scale for contours
  
  # Add NMDS site points (grey, with transparency)
  geom_point(data = site.scores, aes(x = NMDS1, y = NMDS2), color = "grey", alpha = 0.6, size = 3) +
  
  # New color scale for taxa
  ggnewscale::new_scale_color() +
  
  # Add species labels (bold & colored by taxa)
  geom_text_repel(data = spps_score_group, 
                  aes(x = NMDS1, y = NMDS2, label = Species, color = taxa), 
                  size = 4, fontface = "bold", box.padding = 0.3, max.overlaps = 20)+
  scale_color_manual(name = "Taxa", values = okabe) +  # Apply custom colors
  
  # Equal aspect ratio for NMDS
  coord_fixed() +
  
  # Theme and labels
  mdthemes::md_theme_classic() + theme(plot.margin = unit(c(6,0,6,0), "pt"))+
  labs(title = "*P. limitata*", x = "NMDS1", y = "NMDS2")

#pdf("pel.nmds.NEW.pdf")
pel.nmds.plot
#dev.off()

## Combine into a single plot for all 3 spp

library(cowplot)

# arrange the three plots in a single row
prow <- plot_grid( sil.nmds.plot + theme(legend.position="none"),
           fuc.nmds.plot + theme(legend.position="none"),
           pel.nmds.plot + theme(legend.position="none"),
           align = 'vh',
           labels = c("a)", "b)", "c)"),
           hjust = -1,
           nrow = 1
           )

legend <- get_legend(sil.nmds.plot + theme(legend.position="right"))

p <- plot_grid( prow, ncol = 1)
p

save_plot("combined.nmds.pdf", p, base_width = 15, base_height = 7)

l <- plot_grid(legend, ncol = 1)

save_plot("comb.nmds.legend.pdf", l, base_height = 7)
