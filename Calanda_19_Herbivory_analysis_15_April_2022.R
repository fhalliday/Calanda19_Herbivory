
# packages----
require(tidyverse)
require(piecewiseSEM)
require(nlme)
require(interactions)
require(ggbeeswarm)
require(lme4)
require(V.PhyloMaker)
require(car)
require(ggeffects)
require(MuMIn)
require(cowplot)
require(emmeans)

# ggplot theme----
theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) 
theme_set(theme_figs)

update_geom_defaults('point', c(size = 1.5))
update_geom_defaults('errorbar', c(size = 0.5))

pd <- position_dodge(0.3)



# data----
calanda19.com <- read.csv("Calanda_19_data/Community_disease_load.csv") %>% 
  mutate(plot = str_sub(SubplotID, 1, 4)) %>% 
  
  # add plot level data
  left_join(
    read.csv("Calanda_19_data/Plot_elevation.csv") %>% 
      select(-Notes)) %>%
  
  # add plant diversity data
  left_join(
    read.csv("Calanda_19_data/Calanda19_diversity.csv") %>% 
      select(SubplotID = PlotID, plant.richness:hill.simpson)
  ) %>% 
  
  # add phylogenetic diversity
  left_join(
    read.csv("Calanda_19_data/Calanda_19_mpd_no_trees.csv") %>% 
      select(SubplotID = PlotID, mpd.obs, mpd.obs.z)
  ) %>% 
  
  # add traits
  left_join(
    read.csv("Calanda_19_data/Calanda_2019_community_traits_no_trees.csv") %>% 
      select(SubplotID = PlotID, Chlorophyll:RaoQ)
  ) %>% 
  
  # add data-loggers
  left_join(
    read.csv("Calanda_19_data/Calanda_summarized_tms_2019.csv") %>% 
      select(site = Site, mean.min.soil.t:mean.moisture)
  ) %>% 
  
  # some transformations
  mutate(
    l.richness = car::logit(plant.richness), 
    l.simpson = car::logit(hill.simpson),
    Chlo.z = scale(Chlorophyll),
    CN.z = scale(CN),
    f1.z = scale(f1),
    H.z = scale(Height),
    Lon.z = scale(Leaf_lifespan),
    N.z = scale(Leaf_N),
    P.z = scale(Leaf_P),
    Amax.z = scale(Photosynthetic_rate),
    SLA.z = scale(SLA),
    SM.z = scale(Seed_mass),
    # rescale variables
    soil_surf_t_lowsd_scaled = mean.soil.surface.t-(mean(mean.soil.surface.t)-sd(mean.soil.surface.t)),
    soil_surf_t_highsd_scaled = mean.soil.surface.t-(mean(mean.soil.surface.t)+sd(mean.soil.surface.t)),
    soil_surf_t_centered = scale(mean.soil.surface.t, scale = F),
    soil_t_centered = scale(mean.soil.t, scale = F),
    air_t_centered = scale(mean.air.t, scale = F),
    air_t_lowsd_scaled = mean.air.t-(mean(mean.air.t)-sd(mean.air.t)),
    air_t_highsd_scaled = mean.air.t-(mean(mean.air.t)+sd(mean.air.t)),
    moisture_centered = scale(mean.moisture, scale = F),
    elevation_centered = scale(elevation, scale = F),
    sqrt.Disease.l = sqrt(Disease.l),
    f1_scaled = scale(f1, scale = F),
    richness_centered = scale(plant.richness, scale = F),
    sqrt.Herbivory.l = sqrt(Herbivory.l),
    sqrt.chew.l = sqrt(Chew.l),
    sqrt.scrape.l = sqrt(mScrape.l + iScrape.l),
    sqrt.mine.l = sqrt(Mine.l),
    sqrt.window.l = sqrt(Window.l),
    sqrt.thrips.l = sqrt(thrips.l),
    sqrt.skeleton.l = sqrt(mScrape.l + iScrape.l + Window.l),
    sqrt.gall.l = sqrt(Gall.l),
    FDiv_centered = scale(FDiv, scale = F)) %>% 
  
  # some additional variables
  mutate(
    Meadow = factor(meadow),
    Site = site,
    PlotID = plot
  )

summary(calanda19.com)
names(calanda19.com)

# correlations ----
# look at data
calanda19.com %>% 
  dplyr::select(elevation, mean.soil.t, mean.soil.surface.t, mean.air.t, mean.moisture, Disease.l, plant.richness, hill.shannon, f1, Seed_mass, Height, mpd.obs.z, FRic) %>% 
  pairs()

calanda19.com %>% 
  dplyr::select(elevation, mean.soil.t, mean.soil.surface.t, mean.air.t, mean.moisture, Disease.l, plant.richness, hill.shannon, f1, Seed_mass, Height, mpd.obs.z, FRic) %>% 
  cor()

# patterns of herbivory across sites
# Herbivory.l
# Chew.l
# mScrape.l
# iScape.l needs to be combined with thrips.l
# Mine.l
# Gall.l
# Window.l

calanda19.herb <- calanda19.com %>% 
  mutate(Thrips.l = iScrape.l + thrips.l,
         Total_herbivory = Herbivory.l, Chewing = Chew.l, Molluscs = mScrape.l, Thrips = Thrips.l, Mining = Mine.l, Galling = Gall.l, Windowpaning = Window.l) %>% 
  gather(key = "Herbivore_type",
         value = "Community_load", 
         Total_herbivory, Chewing, Molluscs, Thrips, Mining, Galling, Windowpaning) %>% 
  mutate(sqrt_load = sqrt(Community_load))
summary(calanda19.herb)

calanda19.herb %>% 
  ggplot(aes(x = fct_reorder(Site,elevation), y = Community_load, color = fct_reorder(Herbivore_type, Community_load,.desc = TRUE), group = Herbivore_type)) +
  geom_point(position = position_quasirandom(dodge.width = .5)) +
  scale_y_sqrt() +
  theme(legend.position = "bottom")

calanda19.herb %>% 
  ggplot(aes(x = fct_reorder(Site,elevation), y = Community_load)) +
  facet_wrap(~ fct_reorder(Herbivore_type, Community_load,.desc = TRUE), scales = "free") +
  geom_point(position = position_quasirandom(), alpha = .35) +
  stat_summary() +
  scale_y_sqrt() +
  theme(legend.position = "bottom") +
  labs(x = "Site ordered by elevation", y = "Community load")

calanda19.com %>% 
  ggplot(aes(x = fct_reorder(Site, elevation), y = Herbivory.l)) +
  geom_point(position = position_quasirandom(), alpha = .35) +
  stat_summary() +
  scale_y_sqrt() +
  theme(legend.position = "bottom") +
  labs(x = "Site ordered by elevation", y = "Herbivore community load")

calanda19.com %>% 
  ggplot(aes(x = fct_reorder(Site, elevation), y = Chew.l)) +
  geom_point(position = position_quasirandom(), alpha = .35) +
  stat_summary() +
  scale_y_sqrt() +
  theme(legend.position = "bottom") +
  labs(x = "Site ordered by elevation", y = "Chewing community load",
       title = "Chewing insects")

calanda19.com %>% 
  mutate(Thrips.l = iScrape.l + thrips.l) %>% 
  ggplot(aes(x = fct_reorder(Site, elevation), y = Thrips.l)) +
  geom_point(position = position_quasirandom(), alpha = .35) +
  stat_summary() +
  scale_y_sqrt() +
  theme(legend.position = "bottom") +
  labs(x = "Site ordered by elevation", y = "Thrips community load",
       title = "Thrips")

calanda19.herb %>% 
  ggplot(aes(x = mean.soil.surface.t, y = Community_load)) +
  facet_wrap(~ fct_reorder(Herbivore_type, Community_load,.desc = TRUE), scales = "free") +
  geom_point(alpha = .35) +
  geom_smooth(method = "lm") +
  scale_y_sqrt() +
  theme(legend.position = "bottom") +
  labs(x = "Soil surface temperature", y = "Community load")

# are functional and phylo div correlated?
calanda19.com %>% 
  ggplot(aes(x = mpd.obs.z, y = FDiv)) +
  geom_point(alpha = 0.35) +
  labs(x = "phylogenetic diversity", y = "functional diversity")

cor(calanda19.com$mpd.obs.z, calanda19.com$FDiv, method = "pearson")

calanda19.com %>% 
  ggplot(aes(x = plant.richness, y = FDiv)) +
  geom_point(alpha = 0.35) +
  labs(x = "taxonomic diversity", y = "functional diversity")

cor(calanda19.com$plant.richness, calanda19.com$FDiv, method = "pearson")

calanda19.com %>% 
  ggplot(aes(x = plant.richness, y = mpd.obs.z)) +
  geom_point(alpha = 0.35) +
  labs(x = "taxonomic diversity", y = "phylogenetic diversity")

cor(calanda19.com$plant.richness, calanda19.com$mpd.obs.z, method = "pearson")


# check model assumptions ----
# check that component models fit assumptions
com.her.ri <- lme(Herbivory.l ~ plant.richness * elevation + FDiv*elevation + mpd.obs.z*elevation, random = ~1|Meadow/Site/PlotID,  data = calanda19.com)
# check residuals
resids.fig <- function(mod, df) {
  residdf <- dplyr::mutate(df, resids = residuals(mod, type = 'normalized'),
                           fits = fitted(mod))
  fig2 <-ggplot(residdf, aes(x = fits, y = resids)) + geom_point() +
    labs(x = 'Fitted values', y = '')
  
  fig3 <- ggplot(residdf) + stat_qq(aes(sample = resids)) +
    labs(x = 'Theoretical Quantiles', y = 'Sample Quantiles')
  
  # qqline plot = FALSE, according to James should work
  
  fig4 <- ggplot(residdf, aes(x = resids)) + geom_histogram(aes(y=..density..), colour = 'grey50') +
    labs(x = 'Residuals', y = 'Frequency') + scale_y_continuous(expand = c(0, 0)) +
    stat_function(fun = dnorm, color = "red", args = list(mean = mean(residdf$resids),
                                                          sd = sd(residdf$resids)))
  grid::grid.draw(rbind(ggplotGrob(fig2), ggplotGrob(fig3), ggplotGrob(fig4), size = 'first'))
  
  return(summary(mod))
}
variables.fig <- function(df, mod, variable){
  df %>% 
    mutate(resids = residuals(mod, type = 'normalized')) %>% 
    ggplot(aes_string(x=variable, y="resids", color=variable)) +
    xlab(variable) +
    ylab("residuals") +
    # geom_boxplot() +
    geom_point(position=position_jitter(h=0, w=0.4))
}

resids.fig(com.her.ri, calanda19.com)
# appears to violate assumptions of heteroscedasticity and normality

variables.fig(calanda19.com, com.her.ri, "elevation")
variables.fig(calanda19.com, com.her.ri, "plant.richness")
variables.fig(calanda19.com, com.her.ri, "mpd.obs.z")
variables.fig(calanda19.com, com.her.ri, "FDiv")
variables.fig(calanda19.com, com.her.ri, "Site")

# square root transform the response
com.her.ri2 <- lme(sqrt.Herbivory.l ~ plant.richness * elevation + 
                     FDiv*elevation + 
                     mpd.obs.z*elevation, 
                   random = ~1|Meadow/Site/PlotID,  data = calanda19.com)
# check residuals
resids.fig(com.her.ri2, calanda19.com)
# that normalized the residuals somewhat, but still appears to violate assumptions of heterscedasticity
# still appears to violate assumptions of heteroscedasticity

variables.fig(calanda19.com, com.her.ri2, "plant.richness")
variables.fig(calanda19.com, com.her.ri2, "elevation")
variables.fig(calanda19.com, com.her.ri2, "mpd.obs.z")
variables.fig(calanda19.com, com.her.ri2, "FDiv")
variables.fig(calanda19.com, com.her.ri2, "Site")

# model variance separately among sites
com.her.ri3 <- lme(sqrt.Herbivory.l ~ plant.richness * elevation + 
                     FDiv.z*elevation +
                     mpd.obs.z*elevation +
                     FDiv * elevation,
                   random = ~1|Meadow/Site/PlotID,
                   weights = varIdent(form = ~ 1 | Site),
                   data = calanda19.com,
                   control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
resids.fig(com.her.ri3, calanda19.com)
# a bit better
variables.fig(calanda19.com, com.her.ri3, "plant.richness")
variables.fig(calanda19.com, com.her.ri3, "elevation")
variables.fig(calanda19.com, com.her.ri3, "mpd.obs.z") # this one still isn't perfect -- possibly due to outliers?
variables.fig(calanda19.com, com.her.ri3, "FDiv")
variables.fig(calanda19.com, com.her.ri3, "Site")

# plant richness model
pr.ri <- lme(plant.richness ~ elevation, 
             random = ~1|Meadow/Site/PlotID,
             weights = varIdent(form = ~ 1 | Site),
             data = calanda19.com,
             control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

resids.fig(pr.ri, calanda19.com)
# looks ok...

variables.fig(calanda19.com, pr.ri, "elevation")
# looks pretty good
variables.fig(calanda19.com, pr.ri, "mean.moisture")

# functional diversity model
FDiv.ri <- lme(FDiv ~ elevation, 
             random = ~1|Meadow/Site/PlotID,
             weights = varIdent(form = ~ 1 | Site),
             data = calanda19.com,
             control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

resids.fig(FDiv.ri, calanda19.com)
# seems fine

# phylogenetic diversity model
mpd.ri <- lme(mpd.obs.z ~ elevation, 
              random = ~1|Meadow/Site/PlotID,
              weights = varIdent(form = ~ 1 | Site),
              data = calanda19.com,
              control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

resids.fig(mpd.ri, calanda19.com)
# this one seems a bit odd, possibly due to the presence of outliers?

# 1. does air temp modify community structure?----
calanda19.mpdmod <- lme(mpd.obs.z ~ mean.air.t, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.mpdmod)
# significant effect of elevation on mean-pairwise phylogenetic diversity
MuMIn::r.squaredGLMM(calanda19.mpdmod)

calanda19.rmod <- lme(plant.richness ~ mean.air.t, 
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.rmod)
# significant but very weak effect of elevation on richness
MuMIn::r.squaredGLMM(calanda19.rmod)

calanda19.fdivmod <- lme(FDiv ~ mean.air.t, 
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.fdivmod)
# marginal and weak
MuMIn::r.squaredGLMM(calanda19.fdivmod)


r.t <- ggeffects::ggpredict(calanda19.rmod2, "mean.air.t") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = mean.air.t, y = plant.richness), position = position_jitter(), shape = 1)+
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "\nPlant richness\n", x = "Air temperature", title = "")

# pdf("Figures/richness_temp.pdf", height = 4, width = 6)
r.t
# dev.off()


mpd.t <- ggeffects::ggpredict(calanda19.mpdmod, "mean.air.t") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = mean.air.t, y = mpd.obs.z), position = position_jitter(), shape = 1) + 
  geom_hline(yintercept = 0, lty = 2, color = "grey50") +
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Plant phylogenetic \ndiversity (mpd.obs.z)", x = "Air temperature", title = "")

# pdf("Figures/mpdz_temp.pdf", height = 4, width = 6)
mpd.t
# dev.off()


fd.t <- ggeffects::ggpredict(calanda19.fdivmod, "mean.air.t") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = mean.air.t, y = FDiv), position = position_jitter(), shape = 1) + 
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Plant community \nfunctional diversity", x = "Air temperature", title = "")

# pdf("Figures/fdiv_elevation.pdf", height = 4, width = 6)
fd.t
# dev.off()


# pdf("Figures/moderators__air_temp.pdf", height = 7, width = 4)
plot_grid(r.t, mpd.t, fd.t, ncol = 1)
# dev.off()

# add herbivory to this figure
calanda19.herbmod <- lme(sqrt.Herbivory.l ~ mean.air.t, 
                         random = ~1|Meadow/Site/PlotID,
                         weights = varIdent(form = ~ 1 | Site),
                         data = calanda19.com,
                         control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

h.t <- ggeffects::ggpredict(
  calanda19.herbmod, "mean.air.t") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = mean.air.t, y = sqrt.Herbivory.l), position = position_jitter(), shape = 1) + 
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "CWM herbivory \n (square-root transformed)", x = "Air temperature", title = "")

# pdf("Figures/Figure S1.pdf", height = 10, width = 5)
plot_grid(r.t, fd.t, mpd.t, h.t, ncol=1, labels = "auto")
# dev.off()

# as jpg
# jpeg('Figures/Fig S1.jpg',
#      height = 10, width = 5 , res = 250, units = "in")
# plot_grid(r.t, fd.t, mpd.t, h.t, ncol=1, labels = "auto")
# dev.off()


# 2. How does temp combine with community structure to affect herbivory? ----

# note that this model excludes any residual correlation among FDiv, mpd.obs.z, and plant.richness, because we a priori selected variables that are independent from one another
model.c19h <- psem(
  
  lme(sqrt.Herbivory.l ~ air_t_centered  + mpd.obs.z + FDiv + plant.richness +
        plant.richness * air_t_centered + 
        FDiv * air_t_centered +
        mpd.obs.z * air_t_centered,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(mpd.obs.z ~ air_t_centered, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ air_t_centered, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(FDiv ~ air_t_centered, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(air_t_centered ~ elevation, 
      data = calanda19.com,
      random = ~1|Meadow,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)))


model.c19h.fit <- summary(model.c19h, .progressBar = T)
model.c19h.fit

# assess model by rescaling temp
model.c19h_low <- psem(
  
  lme(sqrt.Herbivory.l ~ air_t_lowsd_scaled  + mpd.obs.z + FDiv + plant.richness +
        plant.richness * air_t_lowsd_scaled + 
        FDiv * air_t_lowsd_scaled +
        mpd.obs.z * air_t_lowsd_scaled,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(mpd.obs.z ~ air_t_lowsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ air_t_lowsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(FDiv ~ air_t_lowsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(air_t_lowsd_scaled ~ elevation, 
      data = calanda19.com,
      random = ~1|Meadow,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)))
# note that warning about NAs refers to data that we won't ever use, so i think it's ok.

model.c19h_low.fit <- summary(model.c19h_low, .progressBar = T)
model.c19h_low.fit


model.c19h_high <- psem(
  
  lme(sqrt.Herbivory.l ~ air_t_highsd_scaled  + mpd.obs.z + FDiv + plant.richness +
        plant.richness * air_t_highsd_scaled + 
        FDiv * air_t_highsd_scaled +
        mpd.obs.z * air_t_highsd_scaled,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(mpd.obs.z ~ air_t_highsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ air_t_highsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(FDiv ~ air_t_highsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(air_t_highsd_scaled ~ elevation, 
      data = calanda19.com,
      random = ~1|Meadow,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  # residual correlations among variables
  # mpd.obs.z is independent of species richness, so no need to include that here
  FDiv %~~% mpd.obs.z)
# note that warning about NAs refers to data that we won't ever use, so i think it's ok.

model.c19h_high.fit <- summary(model.c19h_high, .progressBar = T)
model.c19h_high.fit

# Plot the results of the herbivory paths in the SEM ----

calanda19.hmod <- lme(sqrt.Herbivory.l ~ mean.air.t + plant.richness + mpd.obs.z + FDiv +
                         plant.richness*mean.air.t + mpd.obs.z*mean.air.t + FDiv*mean.air.t,
                       random = ~1|Meadow/Site/PlotID,
                       weights = varIdent(form = ~ 1 | Site),
                       data = calanda19.com,
                       control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.hmod)
# r squared
MuMIn::r.squaredGLMM(calanda19.hmod)

fd_modplot <- emmeans::emtrends(calanda19.hmod, ~mean.air.t, var = "FDiv", at = list(mean.air.t = seq(min(calanda19.com$mean.air.t),max(calanda19.com$mean.air.t), .1))) %>% 
  data.frame() %>% 
  ggplot(aes(x = mean.air.t, y = FDiv.trend)) +
  geom_ribbon(fill = "grey", aes(ymin = lower.CL, ymax = upper.CL)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.com, aes(y = FDiv), sides = "b") +
  labs(y = "Effect of functional \ndiversity on herbivory\n", x = "Mean air temperature")

# pdf("Figures/FDiv_moderation.pdf", height = 4, width = 5)
fd_modplot
# dev.off()

mpd_modplot <- emmeans::emtrends(calanda19.hmod, ~mean.air.t, var = "mpd.obs.z", at = list(mean.air.t = seq(min(calanda19.com$mean.air.t),max(calanda19.com$mean.air.t), .1))) %>% 
  data.frame() %>% 
  ggplot(aes(x = mean.air.t, y = mpd.obs.z.trend)) +
  geom_ribbon(fill = "grey", aes(ymin = lower.CL, ymax = upper.CL)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.com, aes(y = 0), sides = "b") +
  labs(y = "Effect of phylogenetic \ndiversity on herbivory", x = "Mean air temperature")

# pdf("Figures/mpd_moderation.pdf", height = 4, width = 5)
mpd_modplot
# dev.off()# 

# A different way to look at the results:

seq(min(calanda19.com$mpd.obs.z),max(calanda19.com$mpd.obs.z), .1)

# phylogenetic diversity (lines represent mean air temp, one sd below and one sd above)
mpd_rawplot <- ggeffects::ggpredict(calanda19.hmod, terms = c("mpd.obs.z[-4:4 by=.1]", "mean.air.t")) %>% 
  data.frame() %>% 
  filter(
    # get rid of values of mpd in group 13.84 that are less than the minimum value that was observed in all observations at or below that temperature
    !(group %in% 13.84 & x < filter(calanda19.com, mean.air.t<=13.84) %>% summarize(min(mpd.obs.z)) %>% as.numeric()),
    # get rid of values of mpd in group 13.84 that are more than the maximum value that was observed in all observations at or below that temperature
    !(group %in% 13.84 & x > filter(calanda19.com, mean.air.t<=13.84) %>% summarize(max(mpd.obs.z)) %>% as.numeric()),
    
    # remove predictions outside the range of observations for the group from 13.84 - 17.42 (called group 15.63 here)
    !(group %in% 15.63 & x < filter(calanda19.com, mean.air.t>13.84 & mean.air.t < 17.42) %>% summarize(min(mpd.obs.z)) %>% as.numeric()),
    !(group %in% 15.63 & x > filter(calanda19.com, mean.air.t>13.84 & mean.air.t < 17.42) %>% summarize(max(mpd.obs.z)) %>% as.numeric()),
    
    # remove predictions outside the range of observations for the group 17.42
    !(group %in% 17.42 & x < filter(calanda19.com, mean.air.t >= 17.42) %>% summarize(min(mpd.obs.z)) %>% as.numeric()),
    !(group %in% 17.42 & x > filter(calanda19.com, mean.air.t >= 17.42) %>% summarize(max(mpd.obs.z)) %>% as.numeric()),
  ) %>% 
  # rename 'group' to mean.air.t
  mutate(mean.air.t = as.numeric(as.character(group))) %>% 
  ggplot(aes(x = x, y = predicted)) +
  geom_point(data = calanda19.com, aes(x = mpd.obs.z, y = sqrt.Herbivory.l, color = mean.air.t), alpha = 0.6) + 
  # geom_ribbon(alpha = .2, aes(ymin = conf.low, ymax = conf.high, group = mean.air.t)) +
  geom_line(aes(color = mean.air.t, group = group), size = 1) +
  scale_color_viridis_c() +
  theme(legend.position = c(.9,.8),
        legend.title = element_text(hjust = 0),
        legend.direction="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0)) +
  labs(y = "CWM Herbivory \n(square-root transformed)", x = "Phylogenetic diversity", color = "Mean air \ntemperature")

mpd_rawplot


# functional diversity
fd_rawplot <- ggeffects::ggpredict(calanda19.hmod, terms = c("FDiv[0.4:0.97 by=.01]", "mean.air.t")) %>% 
  data.frame() %>% 
  filter(
    # get rid of values of mpd in group 13.84 that are less than the minimum value that was observed in all observations at or below that temperature
    !(group %in% 13.84 & x < filter(calanda19.com, mean.air.t<=13.84) %>% summarize(min(FDiv)) %>% as.numeric()),
    # get rid of values of mpd in group 13.84 that are more than the maximum value that was observed in all observations at or below that temperature
    !(group %in% 13.84 & x > filter(calanda19.com, mean.air.t<=13.84) %>% summarize(max(FDiv)) %>% as.numeric()),
    
    # remove predictions outside the range of observations for the group from 13.84 - 17.42 (called group 15.63 here)
    !(group %in% 15.63 & x < filter(calanda19.com, mean.air.t>13.84 & mean.air.t < 17.42) %>% summarize(min(FDiv)) %>% as.numeric()),
    !(group %in% 15.63 & x > filter(calanda19.com, mean.air.t>13.84 & mean.air.t < 17.42) %>% summarize(max(FDiv)) %>% as.numeric()),
    
    # remove predictions outside the range of observations for the group 17.42
    !(group %in% 17.42 & x < filter(calanda19.com, mean.air.t >= 17.42) %>% summarize(min(FDiv)) %>% as.numeric()),
    !(group %in% 17.42 & x > filter(calanda19.com, mean.air.t >= 17.42) %>% summarize(max(FDiv)) %>% as.numeric()),
  ) %>% 
  # rename 'group' to mean.air.t
  mutate(mean.air.t = as.numeric(as.character(group))) %>% 
  ggplot(aes(x = x, y = predicted)) +
  geom_point(data = calanda19.com, aes(x = FDiv, y = sqrt.Herbivory.l, color = mean.air.t), alpha = 0.6) + 
  # geom_ribbon(alpha = .2, aes(ymin = conf.low, ymax = conf.high, group = mean.air.t)) +
  geom_line(aes(color = mean.air.t, group = group), size = 1) +
  scale_color_viridis_c() +
  theme(legend.position = c(.4,.9),
        legend.title = element_text(hjust = 0),
        legend.direction="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0)) +
  labs(y = "CWM Herbivory \n(square-root transformed)", x = "Functional diversity", color = "Mean air \ntemperature")

fd_rawplot

# pdf("Figures/Fig 3.pdf", height = 6, width = 8)
#  plot_grid(fd_rawplot, fd_modplot, mpd_rawplot, mpd_modplot, labels = "auto")
# dev.off()

# as jpg
# jpeg('Figures/Fig 3.jpg',
#      height = 6, width = 8 , res = 250, units = "in")
# plot_grid(fd_rawplot, fd_modplot, mpd_rawplot, mpd_modplot, labels = "auto")
# dev.off()


# 3. consider different herbivore guilds using a multi-response regression ----
calanda19.long <- calanda19.com %>% 
  dplyr::select(SubplotID, PlotID, Site, Meadow, elevation, mean.soil.t, mean.soil.surface.t, mean.air.t, plant.richness, mpd.obs.z, FDiv, sqrt.chew.l:sqrt.gall.l) %>% 
  # Create scaled variables (scaled by max) -- omit galling insects because there aren't observations in every site
  mutate(Chew = sqrt.chew.l * 100 / max(sqrt.chew.l), 
         Thrips =  sqrt.thrips.l * 100 / max(sqrt.thrips.l),
         Skeleton = sqrt.skeleton.l * 100 / max(sqrt.skeleton.l),
         Mine = sqrt.mine.l * 100 / max(sqrt.mine.l)) %>%  
  # Format data frame for multivariate analysis
  gather(key = Var, value = Y, Chew, Thrips, Skeleton, Mine) %>% 
  mutate(Var = factor(Var),
         PlotID = factor(PlotID),
         VarNum = case_when(
           Var == "Chew" ~ 1,
           Var == "Thrips" ~ 2,
           Var == "Skeleton" ~ 3,
           Var == "Mine" ~ 4))

# Multivariate response models
c19.m.multi <- lme(Y ~ Var + 
                     # main effects
                     mean.air.t:Var + plant.richness:Var + mpd.obs.z:Var + FDiv:Var +
                     # interactions
                     mean.air.t:plant.richness:Var + mean.air.t:mpd.obs.z:Var + mean.air.t:FDiv:Var -1, 
                   data = calanda19.long, 
                   na.action = na.omit,
                   random = ~1 | Meadow/Site/PlotID/SubplotID,
                   correlation = corSymm(form = ~ -1| Meadow/Site/PlotID/SubplotID),
                   weights = varComb(varIdent(form = ~ 1 | Site), varIdent(form = ~ -1 | VarNum)),
                   control = list(maxIter=5000000, msMaxIter=500000), method = "ML")

# is that extra correlation necessarry?
c19.m.multi2 <- lme(Y ~ Var + 
                      # main effects
                      mean.air.t:Var + plant.richness:Var + mpd.obs.z:Var + FDiv:Var +
                      # interactions
                      mean.air.t:plant.richness:Var + mean.air.t:mpd.obs.z:Var + mean.air.t:FDiv:Var -1, 
                    data = calanda19.long, 
                    na.action = na.omit,
                    random = ~1 | Meadow/Site/PlotID/SubplotID,
                    correlation = corSymm(form = ~ -1| Meadow/Site/PlotID/SubplotID),
                    weights = varIdent(form = ~ -1 | VarNum),
                    control = list(maxIter=500000000, msMaxIter=50000000), method = "ML")

# Multivariate response ANOVAs ----
car::Anova(c19.m.multi)
car::Anova(c19.m.multi2)
AIC(c19.m.multi,
    c19.m.multi2)
# still an improvement to include the variance structure
# Extract pseudo R-squared
r.squaredGLMM(c19.m.multi)

# Construct multivariate response figure----
c19.ha <- emmeans::emtrends(c19.m.multi, ~mean.air.t | Var, var = "FDiv", at = list(mean.air.t = seq(min(calanda19.long$mean.air.t),max(calanda19.long$mean.air.t), .1))) %>% 
  data.frame() %>% 
  ggplot(aes(x = mean.air.t, y = FDiv.trend)) +
  facet_grid(~Var) + 
  geom_ribbon(fill = "grey", color = "white", aes(ymin = lower.CL, ymax = upper.CL), alpha = .7) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.long, aes(y = FDiv), sides = "b") +
  labs(y = "Functional diversity effect", x = "Mean air temperature")

c19.hb <- emmeans::emtrends(c19.m.multi, ~mean.air.t | Var, var = "mpd.obs.z", at = list(mean.air.t = seq(min(calanda19.long$mean.air.t),max(calanda19.long$mean.air.t), .1))) %>% 
  data.frame() %>% 
  ggplot(aes(x = mean.air.t, y = mpd.obs.z.trend)) +
  facet_grid(~Var) + 
  geom_ribbon(fill = "grey", color = "white", aes(ymin = lower.CL, ymax = upper.CL), alpha = .7) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.long, aes(y = FDiv), sides = "b") +
  labs(y = "Phylogenetic diversity effect\n", x = "Mean air temperature")

c19.hc <- emmeans::emtrends(c19.m.multi, ~mean.air.t | Var, var = "plant.richness", at = list(mean.air.t = seq(min(calanda19.long$mean.air.t),max(calanda19.long$mean.air.t), .1))) %>% 
  data.frame() %>% 
  ggplot(aes(x = mean.air.t, y = plant.richness.trend)) +
  facet_grid(~Var) + 
  geom_ribbon(fill = "grey", color = "white", aes(ymin = lower.CL, ymax = upper.CL), alpha = .7) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.long, aes(y = FDiv), sides = "b") +
  labs(y = "Species richness effect\n", x = "Mean air temperature")

# pdf("Figures/Fig 4.pdf", height = 7, width = 8)
# plot_grid(c19.hc, c19.ha, c19.hb, ncol = 1, labels = "auto")
# dev.off()

# as jpg
# jpeg('Figures/Fig 4.jpg',
#      height = 7, width = 8 , res = 250, units = "in")
# plot_grid(c19.hc, c19.ha, c19.hb, ncol = 1, labels = "auto")
# dev.off()

# some additional plots ----
# Plot net effect of elevation
calanda19.com %>% 
  # create "skeleton" damage type
  mutate(Skeleton = mScrape.l + iScrape.l + Window.l,
         Chew = Chew.l,
         Thrips = thrips.l,
         Mine = Mine.l) %>% 
  # long-form for easy plotting
  gather(key = Var, value = Y, Chew, Thrips, Skeleton, Mine) %>% 
  mutate(Var = factor(Var),
         PlotID = factor(PlotID),
         VarNum = case_when(
           Var == "Chew" ~ 1,
           Var == "Thrips" ~ 2,
           Var == "Skeleton" ~ 3,
           Var == "Mine" ~ 4)) %>% 
    mutate(meadow2 = case_when(
    Meadow == "Im Bofel" ~ "Im Bofel",
    Meadow == "Arella" ~ "Arella",
    Meadow == "Nesselboden" ~ "Nesselboden",
    Meadow == "Oberberg" ~ "Oberberg / Underalp",
    Meadow == "Under Alp" ~ "Oberberg / Underalp") %>% 
      fct_relevel(., "Im Bofel", "Arella", "Nesselboden")
    ) %>% 
  ggplot(aes(x = meadow2, y = Y, color = meadow2)) +
  facet_wrap(~Var, nrow = 1, scales = "free_y") + 
  geom_boxplot(position = pd, width = .6) +
  geom_point(shape = 1, alpha = .8, position = position_quasirandom(dodge.width = .3)) +
  scale_y_sqrt()

calanda19.long %>% 
  mutate(meadow2 = case_when(
    Meadow == "Im Bofel" ~ "Im Bofel",
    Meadow == "Arella" ~ "Arella",
    Meadow == "Nesselboden" ~ "Nesselboden",
    Meadow == "Oberberg" ~ "Oberberg / Underalp",
    Meadow == "Under Alp" ~ "Oberberg / Underalp") %>% 
      fct_relevel(., "Im Bofel", "Arella", "Nesselboden")
  ) %>% 
  ggplot(aes(x = Var, y = Y, color = meadow2)) +
  stat_summary(position = pd, size = .7) +
  geom_point(shape = 1, alpha = .8, position = position_quasirandom(dodge.width = .3))

calanda19.long %>% 
  ggplot(aes(x = elevation, y = Y, color = Var)) + 
  geom_point(shape = 1, alpha = .8) +
  geom_smooth(se = F, span = 2) + 
  facet_wrap(~Var, nrow = 1)

fig5 <- calanda19.com %>% 
  # create "skeleton" damage type
  mutate(Skeleton = mScrape.l + iScrape.l + Window.l,
         Chew = Chew.l,
         Thrips = thrips.l,
         Mine = Mine.l) %>% 
  # long-form for easy plotting
  gather(key = Var, value = Y, Chew, Thrips, Skeleton, Mine) %>% 
  mutate(Var = factor(Var),
         PlotID = factor(PlotID),
         VarNum = case_when(
           Var == "Chew" ~ 1,
           Var == "Thrips" ~ 2,
           Var == "Skeleton" ~ 3,
           Var == "Mine" ~ 4)) %>% 
  ggplot(aes(x = elevation, y = Y)) +
  facet_wrap(~Var, nrow = 1, scales = "free_y") + 
  geom_point(shape = 1, alpha = .8) +
  geom_smooth(se = F, method = "lm", color = "black") +
  scale_y_sqrt() + 
  labs(y = "CWM Herbivory (%)", x = "Elevation (m.a.s.l)")

# pdf("Figures/Fig 5.pdf", height = 3, width = 8)
# fig5
# dev.off()

# as jpg
# jpeg('Figures/Fig 5.jpg',
#      height = 3, width = 8 , res = 250, units = "in")
# print(fig5)
# dev.off()

# Look at patterns in the raw herbivory data, but note that these distributions
# could be a bit misleading since sampling of species is not perfectly random
# (more abundant species are more likely to be measured in our sampling design),
# which is why we analyze at community scale.

cd <- read.csv("Calanda_19_data/Calanda FULL Community Disease survey CLEAN no blanks plotID and speciesID fixed 5.11.2019.csv") %>% 
  mutate(skeleton = mScrape + iScrape + Window,
         meadow = substr(SubplotID, 1, 1)) %>% 
  mutate(Meadow = case_when(
    meadow == "I" ~ "Im Bofel",
    meadow == "A" ~ "Arella",
    meadow == "N" ~ "Nesselboden",
    meadow == "O" ~ "Oberberg/Underalp",
    meadow == "U" ~ "Oberberg/Underalp"
  )) %>% 
  gather(key = herbivore_guild, value = damage, Chew, Mine, Window, thrips, skeleton, Gall) %>% 
  mutate(damage2 = ifelse(damage>0,1,0)) %>% 
  group_by(Meadow, herbivore_guild) %>% 
  summarize(damage2 = sum(damage2),
            nleaves = n()) %>% 
  mutate(frequency_of_damage = damage2/nleaves)

cd %>% 
  ggplot(aes(x = fct_relevel(Meadow, "Im Bofel", "Arella"), y = frequency_of_damage)) + 
  facet_wrap(~herbivore_guild,ncol = 1, scales = "free_y") +
  geom_col() +
  labs(y = "frequency of damage (leaves)", x = "meadow arranged by elevation")

# plant level incicdence

cd_pl <- read.csv("Calanda_19_data/Calanda FULL Community Disease survey CLEAN no blanks plotID and speciesID fixed 5.11.2019.csv") %>% 
  mutate(Skeleton = mScrape + iScrape + Window,
         Thrips = thrips,
         meadow = substr(SubplotID, 1, 1)) %>% 
  mutate(Meadow = case_when(
    meadow == "I" ~ "Im Bofel",
    meadow == "A" ~ "Arella",
    meadow == "N" ~ "Nesselboden",
    meadow == "O" ~ "Oberberg",
    meadow == "U" ~ "Underalp"
  )) %>% 
  gather(key = herbivore_guild, value = damage, Chew, Mine, Thrips, Skeleton) %>% 
  mutate(Unique_plant = paste(SubplotID, SPP, PlantID)) %>% 
  group_by(Meadow, herbivore_guild, Unique_plant) %>% 
  summarize(damage = mean(damage)) %>% 
  ungroup() %>% 
  mutate(damage2 = ifelse(damage>0,1,0)) %>% 
  group_by(Meadow, herbivore_guild) %>% 
  summarize(damage2 = sum(damage2),
            nplants = n()) %>% 
  mutate(frequency_of_damage = damage2/nplants)

cd_pl %>% 
  ggplot(aes(x = fct_relevel(Meadow, "Im Bofel", "Arella"), y = frequency_of_damage)) + 
  facet_wrap(~herbivore_guild,ncol = 1, scales = "free_y") +
  geom_col() +
  labs(y = "frequency of damage (plants)", x = "meadow arranged by elevation")

figs2 <- cd_pl %>% 
  ggplot(aes(x = fct_relevel(Meadow, "Im Bofel", "Arella"), y = frequency_of_damage, fill = fct_reorder(herbivore_guild, -frequency_of_damage))) +
  geom_col(position = position_dodge(width = .9)) +
  labs(y = "frequency of damage (plants)", x = "meadow arranged by elevation") +
  theme(legend.position = c(.2,.9)) +
  colorblindr::scale_fill_OkabeIto()

# pdf("Figures/Fig s2.pdf", height = 4, width = 6)
# figs2
# dev.off()

# as jpg
# jpeg('Figures/Fig S2.jpg',
#      height = 4, width = 6 , res = 250, units = "in")
# print(figs2)
# dev.off()

cd2 <- read.csv("Calanda_19_data/Calanda FULL Community Disease survey CLEAN no blanks plotID and speciesID fixed 5.11.2019.csv") %>% 
  mutate(skeleton = mScrape + iScrape + Window,
         meadow = substr(SubplotID, 1, 1)) %>% 
  mutate(Meadow = case_when(
    meadow == "I" ~ "Im Bofel",
    meadow == "A" ~ "Arella",
    meadow == "N" ~ "Nesselboden",
    meadow == "O" ~ "Oberberg/Underalp",
    meadow == "U" ~ "Oberberg/Underalp"
  ))

cd2 %>% 
  ggplot(aes(x = Chew)) +
  facet_wrap(~fct_relevel(Meadow, "Im Bofel", "Arella")) +
  geom_histogram()

cd2 %>% 
  ggplot(aes(x = skeleton)) +
  facet_wrap(~fct_relevel(Meadow, "Im Bofel", "Arella")) +
  geom_histogram()

cd2 %>% 
  ggplot(aes(x = thrips)) +
  facet_wrap(~fct_relevel(Meadow, "Im Bofel", "Arella")) +
  geom_histogram()

