### Script to create and plot biodivserity indices (SD = shannon diversity, SR = species richness) from the RW layer data

library(readxl)
library(ggplot2)
library(dplyr) 
library(tidyr)
library(purrr)
library(broom)
library(maps)
library(viridis)

# Read raw data
photo.lay <- read_excel("~/photolayerdata_20231218.xlsx", 
                        sheet = "photolayerraw_download")

#filter to remove non-living species
photo.all <- filter(photo.lay, !top_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "OTHSUB", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET"))
photo.all <- filter(photo.all, !bottom_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET", "NONE"))

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
  filter(top_layer %in% c("SILCOM", "PELLIM", "FUCGAR", "HESCAL") & bottom_layer == "NONE")

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

###
# get list of top spp with avg cov
top_all_avg_cov <- top_per_cov %>%
  group_by(target_assemblage, top_layer) %>%
  dplyr::summarize(top_per_cov = mean(n, na.rm=TRUE))

colnames(top_all_avg_cov)[2] <- "species_code"

photolayerdata_20231218 <- read_excel("C:/Users/erica.nielsen/Desktop/Synz/MARINe/LTM_layer_data/photolayerdata_20231218.xlsx", 
                                      sheet = "marine_lumping_codes")

top_spp_tab <- merge(top_all_avg_cov, photolayerdata_20231218, by="species_code")

#export table
#writexl::write_xlsx(top_spp_tab, "avg.cov.top.spp.xlsx")

###
# get list of bot spp with avg cover (Table 1)
bot_all_avg_cov <- bot_per_cov %>%
  group_by(target_assemblage, bottom_layer) %>%
  dplyr::summarize(bot_per_cov = mean(n, na.rm=TRUE))

colnames(bot_all_avg_cov)[2] <- "species_code"

bot_spp_tab <- merge(bot_all_avg_cov, photolayerdata_20231218, by="species_code")

#export table
#writexl::write_xlsx(bot_spp_tab, "avg.cov.bot.spp.xlsx")


##### Calc species richness
#####

# calc sum of RW cov per site/year - first need to filter top_per_cov to only RW
RW_top_cov_all <-filter(top_per_cov, top_layer %in%  c("SILHES", "SILCOM", "PELLIM", "PELHYB", "FUCGAR"))

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

# get RW cover per quadrat/year/assemblage
RW_plot_cov <- RW_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage, top_layer) %>%
  dplyr::summarize(RW_cov = mean(total_RW_cov, na.rm=TRUE))

RW_plot_cov$site_quad_sp <- paste(RW_plot_cov$marine_site_name, "_", RW_plot_cov$quadrat_code, "_", RW_plot_cov$marine_common_year)

# calc bot richness per quadrat/year/assemblage
plot_bot_SR <- bot_per_cov %>% group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage)%>%
  dplyr::summarise(richness=n_distinct(bottom_layer))

plot_bot_SR$site_quad_sp <- paste(plot_bot_SR$marine_site_name, "_", plot_bot_SR$quadrat_code, "_", plot_bot_SR$marine_common_year)

# calc top richness per quadrat/year/assemblage
# photo.all.top <- filter(photo.all, !top_layer %in% "NONE")

#plot_top_SR <-photo.all.top %>% group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage)%>%
 # dplyr::summarise(top_SR=n_distinct(top_layer))


####### calcualte shannon div

#first take avg across seasons
bot_spp_cov <- bot_per_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage, bottom_layer) %>%
  dplyr::summarize(avg_cov = mean(n, na.rm=TRUE))

s.divs <-bot_spp_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage) %>%
  filter(avg_cov>0) %>%
  summarise(N=sum(avg_cov),
            shannon.di=-sum((avg_cov/sum(avg_cov))*log(avg_cov/sum(avg_cov))),
            exp.shannon.di=exp(-sum((avg_cov/sum(avg_cov))*log(avg_cov/sum(avg_cov)))),
            simpson.di=1-sum((avg_cov/sum(avg_cov))^2),
            inv.simpson.di=1/sum((avg_cov/sum(avg_cov))^2)) %>%
  arrange(-shannon.di)

div.all <- left_join(s.divs, plot_bot_SR)

# biplot of the diff diversity measures
library(ggfortify)
fit <- princomp(div.all[,c(9:13)], cor=TRUE)

autoplot( fit, data=div.all, colour="target_assemblage", loadings=TRUE, loadings.label = TRUE, loadings.colour = "black",  loadings.label.colour = "black")
 

#combine into single DF (bottom % cov, RW top % cov, SR, SD)
plot_df_list <- list(bot_plot_cov, RW_plot_cov, div.all)

plot_div_cov <- purrr::reduce(.x = plot_df_list, merge, by = c("site_quad_sp", "target_assemblage"), all = F)

#### make boxplots of div/richness per RW

setwd("C:/Users/erica.nielsen/Desktop/Synz/MARINe/PCI_MARINe_synth/Fucoids/layer_DSR_analyses")

div.all$georegion = factor(div.all$georegion, levels=c('CA North','CA Central','CA South'))
div.all$target_assemblage = factor(div.all$target_assemblage, levels=c('silvetia','fucus','pelvetiopsis'))


pdf("SR.boxplots.pdf")
ggplot(div.all, aes(x=factor(target_assemblage), y=richness, fill = georegion))+
  geom_boxplot()+  scale_fill_manual(values = c("CA North" ='lightseagreen', "CA Central" ='darkorange', 'CA South' = 'magenta')) +
  theme( legend.position = "none", axis.title.x=element_blank() ) +facet_wrap(target_assemblage~., scales = "free_x")+theme_minimal()
dev.off()

pdf("SD.boxplots.pdf")
ggplot(div.all, aes(x=factor(target_assemblage), y=exp.shannon.di, fill = georegion))+
  geom_boxplot()+ scale_fill_manual(values = c("CA North" ='lightseagreen', "CA Central" ='darkorange', 'CA South' = 'magenta')) +
  theme( legend.position = "none" )+facet_wrap(target_assemblage~., scales = "free_x")+theme_minimal()
dev.off()

### Run ANOVA to assess differences in SR and SD between rockweed species

# SR
SR.mod <- aov(richness ~ target_assemblage, data = div.all)

summary(SR.mod)

# Use Tukey's HSD test for post-hoc comparisons
tukey_test <- TukeyHSD(SR.mod)


# Print the results of the Tukey's HSD test
print(tukey_test)

# SD
SD.mod <- aov(exp.shannon.di ~ target_assemblage, data = div.all)

summary(SD.mod)

# Use Tukey's HSD test for post-hoc comparisons
tukey_test <- TukeyHSD(SD.mod)


# Print the results of the Tukey's HSD test
print(tukey_test)



ggplot(div_cov_all_DF, aes(x=target_assemblage, y=exp.shannon.di, color = target_assemblage)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90))
ggplot(div_cov_all_DF, aes(x=target_assemblage, y=richness, color = target_assemblage)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90))

ggplot(div_cov_all_DF, aes(x=georegion, y=exp.shannon.di, color = georegion)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90))
ggplot(div_cov_all_DF, aes(x=georegion, y=richness, color = georegion)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90))

