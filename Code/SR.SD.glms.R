## R code to run models assessing RW cover on understory cover, richness and diversity

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

# read raw data
photo.lay <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/LTM_layer_data/photolayerdata_20231218.xlsx", 
                        sheet = "photolayerraw_download")

# filter to remove non-living species
photo.all <- filter(photo.lay, !top_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))
photo.all <- filter(photo.all, !bottom_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET", "NONE"))

# filter to RW spp
photo.all <-filter(photo.all, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis"))

# We have to swap top_layer record with bottom_layer 'NONE' values so we don't lose understory data
photo.all$bottom_layer[photo.all$bottom_layer == 'NONE'] <- photo.all$top_layer[photo.all$bottom_layer == 'NONE']

# edit again to say if top = bottom, call top layer as "NONE"
photo.all <- photo.all  %>% mutate(top_layer = if_else(bottom_layer==top_layer, 'NONE', top_layer))

### swap bot and top layers for those with RW bot layer (so we only have RW as top layer)

# Step 1: Identify rows to be edited before modification
rows_to_edit_before <- photo.all %>%
  filter(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR"))

# Print rows before editing
cat("Rows before editing:\n")
print(rows_to_edit_before)

# Step 2: Apply the modifications
photo.all <- photo.all %>%
  mutate(
    # Swap bottom_layer and top_layer where conditions are met
    new_bottom_layer = ifelse(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR"),
                              top_layer, bottom_layer),
    new_top_layer = ifelse(top_layer == "NONE" & bottom_layer %in% c("SILCOM", "PELLIM", "FUCGAR"),
                           bottom_layer, top_layer)
  ) %>%
  select(-bottom_layer, -top_layer) %>%
  rename(bottom_layer = new_bottom_layer, top_layer = new_top_layer)

# Step 3: Identify rows that were edited after modification
rows_to_edit_after <- photo.all %>%
  filter(top_layer %in% c("SILCOM", "PELLIM", "FUCGAR") & bottom_layer == "NONE")


# filter out rows with NONE bottom layer
photo.all <- filter(photo.all, !bottom_layer %in% c("NONE"))

# calc % cov per top spp
top_per_cov<-photo.all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(top_layer)
top_per_cov$n <- as.numeric(top_per_cov$n ) / 100 

# calc % cov per bot spp
bot_per_cov<-photo.all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code) %>%
  dplyr::count(bottom_layer)
bot_per_cov$n <- as.numeric(bot_per_cov$n ) / 100

# calc sum of RW cov per site/year - first need to filter top_per_cov to only RW
RW_top_cov_all <-filter(top_per_cov, top_layer %in%  c( "SILCOM", "PELLIM", "FUCGAR"))

# get total RW cover (sum across spp)
RW_quad_cov <- RW_top_cov_all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, season_name, marine_common_year, quadrat_code, target_assemblage, top_layer) %>%
  dplyr::summarize(total_RW_cov = sum(n, na.rm=TRUE))

# filter out non-living bottom records
bot_per_cov <- filter(bot_per_cov, !bottom_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET", "NONE"))

# get total bottom cover (sum across spp) - this 
bot_quad_cov <- bot_per_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, season_name, marine_common_year, quadrat_code, target_assemblage) %>%
  dplyr::summarize(total_U_cov = sum(n, na.rm=TRUE))

# get bot cov per quadrat/year/assemblage (avg across seasons included in previous line)
bot_plot_cov <- bot_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage) %>%
  dplyr::summarize(bot_cov = mean(total_U_cov, na.rm=TRUE))

# add identifier column
bot_plot_cov$site_quad_sp <- paste(bot_plot_cov$marine_site_name, "_", bot_plot_cov$quadrat_code, "_", bot_plot_cov$marine_common_year)

## calcualte RW cover
RW_plot_cov <- RW_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage, top_layer) %>%
  dplyr::summarize(RW_cov = mean(total_RW_cov, na.rm=TRUE))

# add identifier column
RW_plot_cov$site_quad_sp <- paste(RW_plot_cov$marine_site_name, "_", RW_plot_cov$quadrat_code, "_", RW_plot_cov$marine_common_year)

## calcualte species richness
plot_bot_SR <- bot_per_cov %>% group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage)%>%
  dplyr::summarise(richness=n_distinct(bottom_layer))

# add identifier column
plot_bot_SR$site_quad_sp <- paste(plot_bot_SR$marine_site_name, "_", plot_bot_SR$quadrat_code, "_", plot_bot_SR$marine_common_year)

## calcualte shannon div
#first take avg across seasons
bot_spp_cov <- bot_per_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage, bottom_layer) %>%
  dplyr::summarize(avg_cov = mean(n, na.rm=TRUE))

# calc SD
s.divs <-bot_spp_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage) %>%
  filter(avg_cov>0) %>%
  summarise(N=sum(avg_cov),
            shannon.di=-sum((avg_cov/sum(avg_cov))*log(avg_cov/sum(avg_cov))),
            exp.shannon.di=exp(-sum((avg_cov/sum(avg_cov))*log(avg_cov/sum(avg_cov)))),
            simpson.di=1-sum((avg_cov/sum(avg_cov))^2),
            inv.simpson.di=1/sum((avg_cov/sum(avg_cov))^2)) %>%
  arrange(-shannon.di)

# add to SR DF
div.all <- left_join(s.divs, plot_bot_SR)

# combine into single DF 
plot_df_list <- list(bot_plot_cov, RW_plot_cov, div.all)

plot_div_cov <- purrr::reduce(.x = plot_df_list, merge, by = c("site_quad_sp", "target_assemblage"), all = F)

# plot histogram of shannon/simpson div
ggplot(plot_div_cov, aes(x=exp.shannon.di)) + geom_histogram()

ggplot(plot_div_cov.f, aes(x=inv.simpson.di)) + geom_histogram()

## plot histograms, facet by RW spp
# give spp labels for plot
spp_names <- c(
  `FUCGAR` = "F. gardneri",
  `PELLIM` = "P. limitata",
  `SILCOM` = "S. compressa")

SD.histos <- ggplot(plot_div_cov, aes(x=exp.shannon.di)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 10) +
  facet_grid(top_layer ~ ., labeller = as_labeller(spp_names))+ theme_minimal() + labs(x = "Shannon diversity exponential")+
  theme(plot.title = element_text(face = "italic"))

SR.histos <- ggplot(plot_div_cov, aes(x=richness)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 10) +
  facet_grid(top_layer ~ ., labeller = as_labeller(spp_names))+ theme_minimal() + labs(x = "Shannon diversity exponential")+
  theme(plot.title = element_text(face = "italic"))

cov.histos <- ggplot(plot_div_cov, aes(x=bot_cov)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 10) +
  facet_grid(top_layer ~ ., labeller = as_labeller(spp_names))+ theme_minimal() + labs(x = "Shannon diversity exponential")+
  theme(plot.title = element_text(face = "italic"))


## run models PER SPP

# Initialize lists to store model summaries and models
model_summaries <- list()
models <- list()

# Loop through each group
groups <- unique(div_cov_all_DF$top_layer)

for (group in groups) {
  
  # Subset data for the current group
  spp_group <- div_cov_all_DF %>%
    filter(top_layer == group)
  
  # Check if the subsetted data has rows
  if (nrow(spp_group) > 0) {
    
    # Define the GLMM formula
    formula <- richness ~ RW_cov + (1 | marine_site_name/quadrat_code)
    
    # Fit the GLMM using lmer from the lme4 package
    bSR.RWcov <- glmmTMB(formula, data = spp_group)
    
    # Store the summary in the list
    model_summary <- summary(bSR.RWcov)
    model_summaries[[group]] <- model_summary
    
    # Store the model in the list
    models[[group]] <- bSR.RWcov
    
  } else {
    warning(paste("No data available for group:", group))
  }
}


# Print summaries for each group

for (i in seq_along(groups)) {
  
  cat("Summary for group:", groups[i], "\n")
  
  print(model_summaries[[groups[i]]])
  
  cat("\n")
}

diagnose(models$PELLIM)

simulationOutput <- simulateResiduals(fittedModel = models$SILCOM, plot = T)
simulationOutput <- simulateResiduals(fittedModel = models$FUCGAR, plot = T)
simulationOutput <- simulateResiduals(fittedModel = models$PELLIM, plot = T)

sjPlot::tab_model(models$SILCOM, models$FUCGAR, models$PELLIM, show.icc = F, show.aic = T, show.zeroinf = T, show.ngroups = F)


sjPlot::plot_model(models$SILCOM, type = "pred", title = "Silcom")
sjPlot::plot_model(models$FUCGAR, type = "pred", title = "Fucus")
sjPlot::plot_model(models$PELLIM, type = "pred", title = "Pellim")

# diagnostic for temporal correlation
sil.df <- div_cov_all_DF %>% filter(top_layer == 'SILCOM')
fuc.df <- div_cov_all_DF %>% filter(top_layer == 'FUCGAR')
pel.df <- div_cov_all_DF %>% filter(top_layer == 'PELLIM')

pdf("botSR.resid.plots.year.pdf")
boxplot(residuals(models$"SILCOM") ~ sil.df$marine_common_year)
boxplot(residuals(models$"FUCGAR") ~ fuc.df$marine_common_year)
boxplot(residuals(models$"PELLIM") ~ pel.df$marine_common_year)
dev.off()

