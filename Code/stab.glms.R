# Script to run models for stability

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
library(codyn)

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

#filter to remove non-living species
photo.NorCal <- filter(photo.NorCal, !top_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))
photo.NorCal <- filter(photo.NorCal, !bottom_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))

# Filter again to only include target assemblages = RW
photo.NorCal <-filter(photo.NorCal, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis"))

# We have to swap top_layer record with bottom_layer 'NONE' values so we don't lose understory data
photo.NorCal$bottom_layer[photo.NorCal$bottom_layer == 'NONE'] <- photo.NorCal$top_layer[photo.NorCal$bottom_layer == 'NONE']

# edit again to say if top = bottom, call top layer as "NONE"
photo.NorCal <- photo.NorCal  %>% mutate(top_layer = if_else(bottom_layer==top_layer, 'NONE', top_layer))

## edit again to make occurrences with RW bottom and NONE top to be swapped (RW only as top layer)
# Step 1: Identify rows to be edited before modification
rows_to_edit_before <- photo.NorCal %>%
  filter(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR"))

# Step 2: Apply the modifications
photo.NorCal <- photo.NorCal %>%
  mutate(
    # Swap bottom_layer and top_layer where conditions are met
    new_bottom_layer = ifelse(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL"),
                              top_layer, bottom_layer),
    new_top_layer = ifelse(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL"),
                           bottom_layer, top_layer)
  ) %>%
  select(-bottom_layer, -top_layer) %>%
  rename(bottom_layer = new_bottom_layer, top_layer = new_top_layer)


# calculate per cov per quadrat (with counts that we divide by 100)
top_per_cov<-photo.NorCal %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, season_name, georegion, target_assemblage, quadrat_code) %>%
  dplyr::count(top_layer)
top_per_cov$n <- as.numeric(top_per_cov$n ) / 100

RW_top_cov_all <-filter(top_per_cov, top_layer %in%  c("SILCOM", "PELLIM", "FUCGAR"))

# average percent cover for each site/year (avg over quadrats)
top_per_cov_avg<-top_per_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, georegion, top_layer) %>%
  dplyr::summarize(Mean = mean(n, na.rm=TRUE))

# do again for bottom layer spp
bot_per_cov<-photo.NorCal %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, season_name, georegion, target_assemblage, quadrat_code) %>%
  dplyr::count(bottom_layer)
bot_per_cov$n <- as.numeric(bot_per_cov$n ) / 100

bot_quad_cov <- bot_per_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, season_name, georegion, marine_common_year,quadrat_code, target_assemblage) %>%
  dplyr::summarize(total_U_cov = sum(n, na.rm=TRUE))

RW_quad_cov <- RW_top_cov_all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, season_name, georegion, marine_common_year, quadrat_code, top_layer, target_assemblage) %>%
  dplyr::summarize(total_RW_cov = sum(n, na.rm=TRUE))

bot_plot_cov <- bot_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, georegion, quadrat_code, target_assemblage) %>%
  dplyr::summarize(bot_cov = mean(total_U_cov, na.rm=TRUE))

RW_plot_cov <- RW_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, marine_common_year, georegion, quadrat_code, top_layer) %>%
  dplyr::summarize(RW_cov = mean(total_RW_cov, na.rm=TRUE))

# need to combine site, year, assemblg into single factor to use as reps in codyn script

# bot community cover stability
bot_plot_cov$site_quad_sp <- paste(bot_plot_cov$marine_site_name, "_", bot_plot_cov$quadrat_code)

bot_plot_stab <- community_stability(bot_plot_cov, 
                                     time.var = "marine_common_year", 
                                     abundance.var = "bot_cov", 
                                     replicate.var = "site_quad_sp")
names(bot_plot_stab)[2] <- "all_bot_stab"


# RW stability

RW_plot_cov_spp <- RW_plot_cov %>%  pivot_wider(names_from = "top_layer", 
                                                   values_from = "RW_cov")

RW_plot_cov_spp$site_quad_sp <- paste(RW_plot_cov_spp$marine_site_name, "_", RW_plot_cov_spp$quadrat_code)

sil_stab <- community_stability(RW_plot_cov_spp, 
                                    time.var = "marine_common_year", 
                                    abundance.var = "SILCOM", 
                                    replicate.var = "site_quad_sp")
names(sil_stab)[2] <- "sil_RW_stab"

fuc_stab <- community_stability(RW_plot_cov_spp, 
                                time.var = "marine_common_year", 
                                abundance.var = "FUCGAR", 
                                replicate.var = "site_quad_sp")
names(fuc_stab)[2] <- "fuc_RW_stab"

pel_stab <- community_stability(RW_plot_cov_spp, 
                                time.var = "marine_common_year", 
                                abundance.var = "PELLIM", 
                                replicate.var = "site_quad_sp")
names(pel_stab)[2] <- "pel_RW_stab"

#combine spp stab DFs
plot_df_list <- list(sil_stab, fuc_stab, pel_stab)
RW_spp_stab <- plot_df_list %>% reduce(full_join, by='site_quad_sp')


# create RW avg cov 
# calc sum of RW cov per site/year - first need to filter top_per_cov to only RW
RW_top_cov_all <-filter(top_per_cov, top_layer %in%  c( "SILCOM", "PELLIM", "FUCGAR"))

RW_quad_cov <- RW_top_cov_all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, season_name, georegion, marine_common_year, quadrat_code, top_layer) %>%
  dplyr::summarize(total_RW_cov = sum(n, na.rm=TRUE))

#average % cover across seasons
RW_tot_avg_cov <- RW_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, quadrat_code, top_layer) %>%
  dplyr::summarize(avg_tot_RW_cov = mean(total_RW_cov, na.rm=TRUE))

#reshape to get a column per RW cov
RW_quad_cov_spp <- RW_tot_avg_cov %>%  pivot_wider(names_from = "top_layer", names_prefix = "RW_cov_",
                                                values_from = "avg_tot_RW_cov")

RW_quad_cov_spp$site_quad_sp <- paste(RW_quad_cov_spp$marine_site_name, "_", RW_quad_cov_spp$quadrat_code)


## combine into single DF
plot_df_list <- list(bot_plot_stab, bot_avg_plot_SR, RW_spp_stab, RW_quad_cov_spp)
plot_DF_all <- plot_df_list %>% reduce(full_join, by='site_quad_sp')

# make longer with spp as group for RW stab
plot_df_long <- pivot_longer(
  plot_DF_all, 
  cols = ends_with("RW_stab"), 
  names_to = "rw_sp1", 
  values_to = "RW.stab"
)

# make longer with spp as group for RW cov
plot_df_long2 <- pivot_longer(
  plot_DF_all, 
  cols = starts_with("RW_cov_"), 
  names_to = "rw_sp2", 
  values_to = "avg_RW_cov"
)

# combine the DFs & remove the duplicated cols
plot_DF_all2 <- cbind(plot_df_long, plot_df_long2)
plot_DF_all2 <-  plot_DF_all2[!duplicated(colnames(plot_DF_all2))]

plot_DF_all2_filt <- plot_DF_all2 %>% filter_all(all_vars(!is.infinite(.)))

####  Run understory stab models per spp
# do w/o georegion, including mean bot SR as predictor
groups <- unique(plot_DF_all2_filt$rw_sp1)

stab_models <- list() 
stab_AICc_tables <- list()

for (group in groups) {
  spp_group <- subset(plot_DF_all2_filt, rw_sp1 == group)
  
  spp_group <- spp_group %>% drop_na(RW.stab)
  
  m1.mean.SR <- glmmTMB(all_bot_stab ~  Mean.SR + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m1a.mean.SR <- glmmTMB(all_bot_stab ~  Mean.SR*georegion.x + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m2.RW.stab <- glmmTMB(all_bot_stab ~  RW.stab + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m2a.RW.stab <- glmmTMB(all_bot_stab ~  RW.stab*georegion.x + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m3.mean.RW.cov <- glmmTMB(all_bot_stab ~ avg_RW_cov + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m3a.mean.RW.cov <- glmmTMB(all_bot_stab ~ avg_RW_cov*georegion.x + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m4.mSR.RWstab <- glmmTMB(all_bot_stab ~ RW.stab + Mean.SR + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m5.mSR.mRWcov <- glmmTMB(all_bot_stab ~ avg_RW_cov + Mean.SR + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m6.RWstab.mRWcov <- glmmTMB(all_bot_stab ~ avg_RW_cov + RW.stab + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  m7.mSR.RWstab.mRWcov <- glmmTMB(all_bot_stab ~ avg_RW_cov + RW.stab + Mean.SR + (1 | marine_site_name.x),  data = spp_group, family=Gamma(link = "log"))
  
  stab_models[[as.character(group)]] <- list(m1.mean.SR = m1.mean.SR, m1a.mean.SR=m1a.mean.SR, m2.RW.stab=m2.RW.stab, m2a.RW.stab= m2a.RW.stab, m3.mean.RW.cov=m3.mean.RW.cov, m3a.mean.RW.cov=m3a.mean.RW.cov, m4.mSR.RWstab=m4.mSR.RWstab, m5.mSR.mRWcov=m5.mSR.mRWcov, m6.RWstab.mRWcov=m6.RWstab.mRWcov, m7.mSR.RWstab.mRWcov=m7.mSR.RWstab.mRWcov)
  stab_AICc_table <- bbmle::AICctab(m1.mean.SR, m1a.mean.SR, m2.RW.stab, m2a.RW.stab, m3.mean.RW.cov, m3a.mean.RW.cov, m4.mSR.RWstab, m5.mSR.mRWcov, m6.RWstab.mRWcov, m7.mSR.RWstab.mRWcov)
  stab_AICc_tables[[as.character(group)]] <- stab_AICc_table }

print(stab_AICc_tables[[as.character(groups[1])]]) #change the # after groups to get the AIC table per RW spp -> choose best model to include in code below

# run summay for each RW spp best model
summary(stab_models$sil_RW_stab$m3.mean.RW.cov)

sjPlot::tab_model(stab_models$sil_RW_stab$m3.mean.RW.cov,  stab_models$fuc_RW_stab$m1.mean.SR, stab_models$pel_RW_stab$m3a.mean.RW.cov, show.icc = F, show.aic = T, show.zeroinf = T, show.ngroups = F)

sjPlot::plot_model(stab_models$sil_RW_stab$m3.mean.RW.cov, type = "pred")
sjPlot::plot_model(stab_models$fuc_RW_stab$m1.mean.SR, type = "pred")
sjPlot::plot_model(stab_models$pel_RW_stab$m3a.mean.RW.cov, type = "pred")


