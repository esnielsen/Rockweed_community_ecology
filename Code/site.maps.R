#### R script to create map of sample sites for RW layer data analyses

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
library(ggplot2)
library(maps)
library(dplyr)

# read raw data
photo.lay <- read_excel("~/photolayerdata_20231218.xlsx", 
                        sheet = "photolayerraw_download")

#get number years sampled
phot_count <- photo.lay %>%                           
  group_by(marine_site_name) %>%
  dplyr::summarise(count = n_distinct(marine_common_year))

#new DF with yrs sampled column
photo.yr <- merge(photo.lay,  phot_count, by.x = "marine_site_name", by.y = "marine_site_name")

#filter to RW species
photo.all <-filter(photo.yr, target_assemblage %in%  c("silvetia", "fucus", "pelvetiopsis"))

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
#photo.all <- filter(photo.all, !bottom_layer %in% c("NONE"))

# calc % cov per top spp
top_per_cov<-photo.all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code, count) %>%
  dplyr::count(top_layer)
top_per_cov$n <- as.numeric(top_per_cov$n ) / 100 

# calc % cov per bot spp
bot_per_cov<-photo.all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, season_name, target_assemblage, quadrat_code, count) %>%
  dplyr::count(bottom_layer)
bot_per_cov$n <- as.numeric(bot_per_cov$n ) / 100

# calc sum of RW cov per site/year - first need to filter top_per_cov to only RW
RW_top_cov_all <-filter(top_per_cov, top_layer %in%  c( "SILCOM", "PELLIM", "FUCGAR"))

# calc top_SR of other algal spp, so we filter out RWs and other non-canopy spp (like barnacles & anemones)
#non_RW_top_cov <- filter(top_per_cov, !top_layer %in% c("SILHES", "SILCOM", "PELLIM", "PELHYB", "FUCGAR", "HESCAL", "MYTTRO", "MOPSPP", "ANTXAN", "NONCRU", "BALGLA", "CRUCOR", "POLPOL", "MYTCAL", "LIMPET", "CHTDAL", "ARTCOR", "ANTELE", "NONE", "TETRUB"))

# get total RW cover (sum across spp)
RW_quad_cov <- RW_top_cov_all %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, season_name, marine_common_year, quadrat_code, target_assemblage, top_layer, count) %>%
  dplyr::summarize(total_RW_cov = sum(n, na.rm=TRUE))

# filter out non-living bottom records
bot_per_cov <- filter(bot_per_cov, !bottom_layer %in% c("UNIDEN", "ROCK", "TAR", "SAND", "DEABAL", "DEACHT", "DEACRA", "DEADCB", "DEAINV", "DEAMCA", "DEAMTR", "DEASBB", "DEASEM", "DEATET", "NONE"))

# get total bottom cover (sum across spp) -
bot_quad_cov <- bot_per_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, season_name, marine_common_year, quadrat_code, target_assemblage, count) %>%
  dplyr::summarize(total_U_cov = sum(n, na.rm=TRUE))

# get bot cov per quadrat/year/assemblage (avg across seasons included in previous line)
bot_plot_cov <- bot_quad_cov %>%
  group_by(ltm_latitude, ltm_longitude, marine_site_name, georegion, marine_common_year, quadrat_code, target_assemblage, count) %>%
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

# calcualte shannon div
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

#combine into single DF (this now changed to include spp as to_layer column - use this in Loop!)
plot_df_list <- list(bot_plot_cov, RW_plot_cov, div.all)

plot_div_cov <- purrr::reduce(.x = plot_df_list, merge, by = c("site_quad_sp", "target_assemblage"), all = F)

#make final DF
div_cov_all_DF <- plot_div_cov %>% drop_na() %>% filter_all(all_vars(!is.infinite(.)))

#### make maps
setwd("C:/Users/erica.nielsen/Desktop/Synz/MARINe/PCI_MARINe_synth/Fucoids/layer_DSR_analyses/site_maps")

#get state shapefile
states <- map_data("state")
wc <- states %>%
  filter(region %in% c("california"))

### Plot for Pellim
photo.s <-filter(div_cov_all_DF, top_layer %in%  c("PELLIM"))

pel.map <- ggplot() +
  geom_polygon(data = wc, aes(x=long, y = lat, group = group)) + coord_fixed(1.3)  +
  geom_point(data=photo.s, mapping = aes(x = ltm_longitude, y = ltm_latitude, color=georegion.x), cex=4)+
  scale_color_manual(values = c("CA Central" ='darkorange', "CA North" ='lightseagreen'))+
  theme(plot.background=element_blank(), panel.background=element_blank(), axis.text=element_text(colour='black',size=14), axis.title.x=element_blank(), axis.title.y=element_blank())


### Plot for Fucus
photo.s <-filter(div_cov_all_DF, top_layer %in%  c("FUCGAR"))

fuc.map <- ggplot() +
  geom_polygon(data = wc, aes(x=long, y = lat, group = group)) + coord_fixed(1.3)  +
  geom_point(data=photo.s, mapping = aes(x = ltm_longitude, y = ltm_latitude, color=georegion.x), cex=4)+
  scale_color_manual(values = c("CA Central" ='darkorange', "CA North" ='lightseagreen'))+
  theme(plot.background=element_blank(), panel.background=element_blank(), axis.text=element_text(colour='black',size=14), axis.title.x=element_blank(), axis.title.y=element_blank())

### Plot for Silcom
photo.s <-filter(div_cov_all_DF, top_layer %in%  c("SILCOM"))

sil.map <- ggplot() +
  geom_polygon(data = wc, aes(x=long, y = lat, group = group)) + coord_fixed(1.3)  +
  geom_point(data=photo.s, mapping = aes(x = ltm_longitude, y = ltm_latitude, color=georegion.x), cex=4)+
  scale_color_manual(values = c("CA Central" ='darkorange', 'CA South' = 'magenta'))+
  theme(plot.background=element_blank(), panel.background=element_blank(), axis.text=element_text(colour='black',size=14), axis.title.x=element_blank(), axis.title.y=element_blank())

pdf("pel.site.maps.pdf")
pel.map
dev.off()

pdf("fuc.site.maps.pdf")
fuc.map
dev.off()

pdf("sil.site.maps.pdf")
sil.map
dev.off()
