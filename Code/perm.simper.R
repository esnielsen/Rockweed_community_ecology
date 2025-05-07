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

# read raw data
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
photo.NorCal <-filter(photo.NorCal, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis"))


# calculate per cov per quadrat (with counts that we divide by 100)
top_per_cov<-photo.NorCal %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(top_layer)
top_per_cov$n <- as.numeric(top_per_cov$n ) / 100

# do again for bottom layer spp
bot_per_cov<-photo.NorCal %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code, top_layer) %>%
  dplyr::count(bottom_layer)
bot_per_cov$n <- as.numeric(bot_per_cov$n ) / 100


# running simper & permanova on all sites

bot_per_cov <- bot_per_cov %>% 
  rename("bot_cov" = "n")

# add top_layer % cov to the DF
top_bot_plots <- merge(top_per_cov, bot_per_cov)
top_bot_plots <- top_bot_plots %>% 
  rename("top_cov" = "n")

#filter to RW top cov
RW_top_cov_plot <-filter(top_bot_plots, top_layer %in%  c("SILCOM", "PELLIM", "FUCGAR"))

# make new column of high medium low
#make new column to group RW cov into high/low (above/below 50% cov)
RW_top_cov_plot <- RW_top_cov_plot %>%
  mutate(cov_group = case_when(top_cov>=0.67~"high_cover", between(top_cov, 0.33, 0.67)~"med_cover",
                               top_cov<=0.32~"low_cover"))


# filter to each RW spp
sil_diff_plots <-filter(RW_top_cov_plot, top_layer %in%  "SILCOM")
sil_diff_plots <-filter(sil_diff_plots, target_assemblage %in% "silvetia")

pel_diff_plots <-filter(RW_top_cov_plot, top_layer %in%  "PELLIM")
pel_diff_plots <-filter(pel_diff_plots, target_assemblage %in% "pelvetiopsis")

fuc_diff_plots <-filter(RW_top_cov_plot, top_layer %in% "FUCGAR")
fuc_diff_plots <-filter(fuc_diff_plots, target_assemblage %in%  "fucus")


#make wide per RW spp
sil_diff_plot_w <- sil_diff_plots %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
#change NAs to zeros
sil_diff_plot_w <- sil_diff_plot_w %>% replace(is.na(.), 0)

pel_diff_plot_w <- pel_diff_plots %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
pel_diff_plot_w <- pel_diff_plot_w %>% replace(is.na(.), 0)

fuc_diff_plot_w <- fuc_diff_plots %>%
  pivot_wider(names_from = bottom_layer, values_from = bot_cov)
fuc_diff_plot_w <- fuc_diff_plot_w %>% replace(is.na(.), 0)

# filter out medium cover rows
sil_h.l_cov <-filter(sil_diff_plot_w, !cov_group %in% "med_cover")
pel_h.l_cov <-filter(pel_diff_plot_w, !cov_group %in% "med_cover")
fuc_h.l_cov <-filter(fuc_diff_plot_w, !cov_group %in% "med_cover")


### R script for later
pairwise_permanova <- function(sp_matrix, group_var, dist = "bray", adj = "fdr", perm = 10000) {
  
  require(vegan)
  
  ## list contrasts
  group_var <- as.character(group_var)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
  
  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA
  )
  
  for (i in seq(nrow(contrasts))) {
    sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i] 
    contrast_matrix <- sp_matrix[sp_subset,]
    
    ## fit contrast using adonis
    fit <- vegan::adonis2(
      contrast_matrix ~ group_var[sp_subset],
      method = dist, 
      perm = perm
    )
    
    contrasts$R2[i] <- round(fit$R2[1], digits = 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
    contrasts$df1[i] <- fit$Df[1]
    contrasts$df2[i] <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
  }
  
  ## adjust p-values for multiple comparisons
  contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)
  
  return(list(
    contrasts = contrasts, 
    "p-value adjustment" = adj, 
    permutations = perm
  ))
}



# permanova w/o region & three coverage groups
sil.manova <- adonis2(sil_diff_plot_w[ , 12:71] ~ sil_diff_plot_w$cov_group, method = "euclidean", data= sil_diff_plot_w)
sil.manova

pairwise_permanova(sil_diff_plot_w[ , 12:71], sil_diff_plot_w$cov_group)

pel.manova <- adonis2(pel_diff_plot_w[ , 12:50] ~ pel_diff_plot_w$cov_group, method = "euclidean", data= pel_diff_plot_w)
pel.manova

pairwise_permanova(pel_diff_plot_w[ , 12:50], pel_diff_plot_w$cov_group)

fuc.manova <- adonis2(fuc_diff_plot_w[ , 12:53] ~ fuc_diff_plot_w$cov_group, method = "euclidean", data= fuc_diff_plot_w)
fuc.manova

pairwise_permanova(fuc_diff_plot_w[ , 12:53], fuc_diff_plot_w$cov_group)


## permanova with region and only high/low cov
sil.2.manova <- adonis2(sil_h.l_cov[ , 12:71] ~ cov_group * georegion, data = sil_h.l_cov, permutations = 999)
sil.2.manova

pel.2.manova <- adonis2(pel_h.l_cov[ , 12:50] ~ cov_group * georegion, data = pel_h.l_cov, permutations = 999)
pel.2.manova

fuc.2.manova <- adonis2(fuc_h.l_cov[ , 12:53] ~ cov_group * georegion, data = fuc_h.l_cov, permutations = 999)
fuc.2.manova


#### SIL simper

# 1st turn spp df into matrix
mat_sil_spp <- as.matrix(sil_diff_plot_w[ , 12:70])

# run simper, with RW cov group as pred var
simp_sil <- simper(mat_sil_spp, sil_diff_plot_w$cov_group)

# Filter for species contributing more than 3%
spp_cont <- as.data.frame(simp_sil$low_cover_high_cover)

sil_simp_spp <- spp_cont[spp_cont$average > 0.03, ]


## PEL SIMPER
# 1st turn spp df into matrix
mat_pel_spp <- as.matrix(pel_diff_plot_w[ , 11:49])

# run simper, with RW cov group as pred var
simp_pel <- simper(mat_pel_spp, pel_diff_plot_w$cov_group)

# Filter for species contributing more than 3%
spp_cont <- as.data.frame(simp_pel$low_cover_high_cover)

pel_simp_spp <- spp_cont[spp_cont$average > 0.03, ]


## FUC SIMPER
# 1st turn spp df into matrix
mat_fuc_spp <- as.matrix(fuc_diff_plot_w[ , 11:52])

# run simper, with RW cov group as pred var
simp_fuc <- simper(mat_fuc_spp, fuc_diff_plot_w$cov_group)

# Filter for species contributing more than 3%
spp_cont <- as.data.frame(simp_fuc$high_cover_low_cover)

fuc_simp_spp <- spp_cont[spp_cont$average > 0.03, ]

all_spp_simp <- bind_rows(lst(sil_simp_spp, pel_simp_spp, fuc_simp_spp), .id = 'id')

write.csv(all_spp_simp, "all_spp_simp.csv")
