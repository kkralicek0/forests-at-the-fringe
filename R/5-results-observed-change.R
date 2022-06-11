###############################################################################
## Project - Forests at the fringe                                         ----
## Script  - 5-results-observed-change.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Forests at the Fringe: comparing observed
## change to projected climate change impacts for five species in the Pacific 
## Northwest' by Karin Kralicek, Tara Barrett, Jay Ver Hoef, and Temesgen
## Hailemariam.
## 
## About - this script:
## - Summarize and visualize observed change estimates within PECA-domains to
##   compare with projected climate change impacts and assess if data suggests
##   that projected impacts are already being experienced by species... 
## - Produce figures/tables for manuscript and appendix.
##   (note - 'projected impacts' results are analyzed in the previous script)
##
###############################################################################
# data manipulation
library(dplyr)      # for ... E V E R Y T H I N G ...
library(tidyr)
library(purrr)      # o! beloved map...
library(sf)         # spatial data for plotting
library(ggplot2)    # plots
library(patchwork)  # plots
library(parallel)

# also using lemon to place legend within facets for a few figures, but doing
# it as lemon::reposition_legend() so not loading entire library here..
# library(lemon)

# == Paths =============================================================== ####
# path to Rdata
path_top <- paste0("C:/Users/karinkralicek/Documents/", 
                   "projects/active/SI - forest at the fringe/code/")

# paths to intermediate data
path_output <- paste0(path_top, "x - output/")
path_data  <- paste0(path_top, "x - data/")
path_output_ch2 <- paste0("C:/Users/karinkralicek/Documents/CUI/",
                          "ch2 output sample/") # for spatial data..

# where to save figures (at correct dpi, etc.)
path_images <- paste0(path_top, "images/sort/")

# == load data =========================================================== ####
# - load observed change estimates --------------------------------------- ####
# domain estimates from script 3
load(paste0(path_output, "FF-3-domain-estimates.Rdata"))

est_overall <- est_overall %>% map(~.x %>% ungroup)
est_elevlat <- est_elevlat %>% map_depth(2, ~.x %>% ungroup)
est_peca_cna <- est_peca_cna %>% 
  map_depth(2, function(sp.x) {
    sp.x %>% 
      map(~.x %>% bind_rows(.id = "peca_domain")) %>% 
      bind_rows(.id = "mcmc_n") %>%
      ungroup
  }) %>%
  transpose %>%
  # pretty scenario labels for plotting
  map(function(sp.x) {
    sp.x %>%
      bind_rows(.id = "scenario") %>% 
      mutate(scenario = case_when(
        scenario == "CCSM4_rcp45" ~ "RCP 4.5 CCSM4",
        scenario == "HadGEM2_ES_rcp45" ~ "RCP 4.5 HadGEM2-ES",
        scenario == "CCSM4_rcp85" ~ "RCP 8.5 CCSM4",
        scenario == "HadGEM2_ES_rcp85" ~ "RCP 8.5 HadGEM2-ES")) %>%
      mutate(scenario = factor(scenario, 
                               levels = c("RCP 4.5 CCSM4",
                                          "RCP 4.5 HadGEM2-ES",
                                          "RCP 8.5 CCSM4",
                                          "RCP 8.5 HadGEM2-ES")))
  })

# - load elev-lat groups ------------------------------------------------- ####
# Load these for some specific summaries about prevalence within naive-div
elevlat <- readRDS(paste0(path_data, "FF-3-elev-lat-groups.rds"))

# - load spatial data ---------------------------------------------------- ####
# For mapping elev/lat (naive-) and PECA-divisions
plot_polygons_public.sf <- readRDS(
  paste0(path_output_ch2, "ch2-8-plot_polygons_public_sf.rds"))
study_area.sf <- readRDS(
  paste0(path_output_ch2, "ch2-8-study_area_sf.rds"))

# == Color palettes & plotting fns ======================================= ####
clim_ds_palette <- c("#999999", "#af8dc3", "#7fbf7b", "#762a83", "#1b7837") %>%
  setNames(c("current", "CCSM4_rcp45", "HadGEM2_ES_rcp45", 
             "CCSM4_rcp85","HadGEM2_ES_rcp85"))

peca_domain_palette <- c('#377eb8', "#009E73", "#D55E00", "#E69F00") %>%
  setNames(c("E", "P", "C", "A"))

# set general text-size & figure width arguments
# - these based on the min allowable for Frontiers' in F&GC in their author's
#   guidelines
FF_font_size <- 8
FF_fig_width_2col <- 180 #mm
FF_fig_width_1col <- 85  #mm
FF_width_unit <- "mm"
FF_min_dpi <- 300

# == Tidy some `NaN`s in PECA data ======================================= ####
# About: sometimes we were unable to make an estimate for a PECA-domain 
# that was projected to exist under an MCMC sample; such cases produced 
# `NaN` for one of the two reasons outlined below:
# - Case 1: an estimate (and SEM) will have a value of `NaN` when the est
#   was calculated as 0/0 - which will happen for stat_id's 2|4|6|8:10.
#   - for 2|4|6|8, that will happen when there were no plots within that 
#     PECA-domain for that MCMC sample with the sp present at t1 (stat_id's
#     2|4 for mortality) or at either t1|t2 (stat_id's 6|8 for net-growth).
#   - for 9|10, this will happen if gross growth at the population-level 
#     was 0, which can occur either b/c net-growth and mortality are equally
#     balanced or because the species did not occur within that domain..
#   Either way, we can't actually make estimates for these MCMC sample's 
#   PECA-domains. Therefore check where this happens and convert these NaN's 
#   to NA.
# - Case 2: if the denominator of the ratio-of-means estimator (i.e., mean 
#   prop of area within the domain) is very close to, but not quite, zero, 
#   it can happen that we can have a non-NaN estimate, but the SEM will be 
#   NaN. This would be rare, but did happen for 16 {spp x mcmc x peca-domain
#   x stat_id x scenario} as shown below.

# - investigate when this happens ---------------------------------------- ####
# first case (when est & SEM are both NaN, will also show case 2 for 9&10)
# - by n_mcmc & perc_mcmc where this was the case (note that for abpr, there 
#   was only 1 mcmc so 100% for RCP 8.5 CCSM4)
# - note, b/c rel-ng (9&10) calc, they're more likely to result in NaN's b/c
#   of machine precision issues etc. with gross growth in the denom
est_peca_cna %>% 
  map(function(sp.x) {
    sp.x %>%
      group_by(scenario, peca_domain, stat_id) %>%
      mutate(mcmc_tot = n()) %>%
      ungroup %>% 
      filter(is.nan(est) & stat_id %in% paste0('stat_', c((1:4)*2, 9:10))) %>% 
      group_by(scenario, peca_domain, stat_id) %>% 
      summarize(n_mcmc = n(),
                perc_mcmc = n_mcmc / unique(mcmc_tot)) %>%
      ungroup %>%
      pivot_longer(c("n_mcmc", "perc_mcmc"),
                   names_to = "stat",
                   values_to = "value") %>%
      split(.$stat) %>%
      map(~.x %>% pivot_wider(names_from = 'stat_id', values_from = 'value'))
  }) %>% 
  keep(~length(.x)>0) %>%
  transpose %>% 
  map(~.x %>% bind_rows(.id = "spp"))

# less common case (when SEM is NaN, but est is not)
# - not shown here but falls into this category will be some/all of the stat 
#   9/10 as mentioned above...
est_peca_cna %>% 
  map(function(sp.x) {
    sp.x %>%
      group_by(scenario, peca_domain, stat_id) %>%
      mutate(mcmc_tot = n()) %>%
      ungroup %>% 
      filter(stat_id %in% paste0('stat_', c((1:4)*2, 9:10))) %>%
      filter(!is.nan(est) & is.nan(SEM)) %>% 
      group_by(scenario, peca_domain, stat_id) %>% 
      summarize(n_mcmc = n(),
                perc_mcmc = n_mcmc / unique(mcmc_tot)) %>%
      ungroup %>%
      pivot_longer(c("n_mcmc", "perc_mcmc"),
                   names_to = "stat",
                   values_to = "value") %>%
      split(.$stat) %>%
      map(~.x %>% pivot_wider(names_from = 'stat_id', values_from = 'value'))
  }) %>% 
  keep(~length(.x)>0) %>%
  transpose %>% 
  map(~.x %>% bind_rows(.id = "spp"))

# - what other stats for this {spp x mcmc x peca-domain} looked like:
est_peca_cna$qudo %>% filter(mcmc_n == 32 & peca_domain == "C")

# percent of MCMC samples x PECA-domains x scenarios where this was an issue...
# ... here only reporting for stat_ 2|4|6|8
est_peca_cna %>% 
  bind_rows(.id = 'spp') %>% 
  filter(stat_id %in% paste0("stat_", (1:4)*2)) %>%
  group_by(spp, scenario, mcmc_n, peca_domain) %>%
  summarize(TF_issue = as.numeric(sum(is.nan(est) | is.na(SEM)) > 0)) %>% 
  ungroup %>%
  group_by(scenario) %>%
  mutate(tot_perc_issue = sum(TF_issue)/n()) %>%
  group_by(scenario, tot_perc_issue, spp) %>%
  summarize(spp_perc_issue = sum(TF_issue)/n()) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x, 3))) %>%
  pivot_longer(cols = c("tot_perc_issue", "spp_perc_issue")) %>%
  pivot_wider(names_from = scenario,
              values_from = value) %>%
  mutate(spp = if_else(name == "tot_perc_issue", "ALL", spp)) %>%
  distinct(spp, .keep_all = T) %>%
  arrange(desc(name))

# ... stats 9&10
est_peca_cna %>% 
  bind_rows(.id = 'spp') %>% 
  filter(stat_id %in% paste0("stat_", 9:10)) %>%
  group_by(spp, scenario, mcmc_n, peca_domain) %>%
  summarize(TF_issue = as.numeric(sum(is.nan(est) | is.na(SEM)) > 0)) %>% 
  ungroup %>%
  group_by(scenario) %>%
  mutate(tot_perc_issue = sum(TF_issue)/n()) %>%
  group_by(scenario, tot_perc_issue, spp) %>%
  summarize(spp_perc_issue = sum(TF_issue)/n()) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x, 3))) %>%
  pivot_longer(cols = c("tot_perc_issue", "spp_perc_issue")) %>%
  pivot_wider(names_from = scenario,
              values_from = value) %>%
  mutate(spp = if_else(name == "tot_perc_issue", "ALL", spp)) %>%
  distinct(spp, .keep_all = T) %>%
  arrange(desc(name))

# ... just stat 9
est_peca_cna %>% 
  bind_rows(.id = 'spp') %>% 
  filter(stat_id %in% paste0("stat_", 9)) %>%
  group_by(spp, scenario, mcmc_n, peca_domain) %>%
  summarize(TF_issue = as.numeric(sum(is.nan(est) | is.na(SEM)) > 0)) %>% 
  ungroup %>%
  group_by(scenario) %>%
  mutate(tot_perc_issue = sum(TF_issue)/n()) %>%
  group_by(scenario, tot_perc_issue, spp) %>%
  summarize(spp_perc_issue = sum(TF_issue)/n()) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x, 3))) %>%
  pivot_longer(cols = c("tot_perc_issue", "spp_perc_issue")) %>%
  pivot_wider(names_from = scenario,
              values_from = value) %>%
  mutate(spp = if_else(name == "tot_perc_issue", "ALL", spp)) %>%
  distinct(spp, .keep_all = T) %>%
  arrange(desc(name))

# ... just stat 10
est_peca_cna %>% 
  bind_rows(.id = 'spp') %>% 
  filter(stat_id %in% paste0("stat_", 10)) %>%
  group_by(spp, scenario, mcmc_n, peca_domain) %>%
  summarize(TF_issue = as.numeric(sum(is.nan(est) | is.na(SEM)) > 0)) %>% 
  ungroup %>%
  group_by(scenario) %>%
  mutate(tot_perc_issue = sum(TF_issue)/n()) %>%
  group_by(scenario, tot_perc_issue, spp) %>%
  summarize(spp_perc_issue = sum(TF_issue)/n()) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x, 3))) %>%
  pivot_longer(cols = c("tot_perc_issue", "spp_perc_issue")) %>%
  pivot_wider(names_from = scenario,
              values_from = value) %>%
  mutate(spp = if_else(name == "tot_perc_issue", "ALL", spp)) %>%
  distinct(spp, .keep_all = T) %>%
  arrange(desc(name))

#.... stats 2|4|6|8:10
est_peca_cna %>% 
  bind_rows(.id = 'spp') %>% 
  filter(stat_id %in% paste0("stat_", c((1:4)*2, 9, 10))) %>%
  group_by(spp, scenario, mcmc_n, peca_domain) %>%
  summarize(TF_issue = as.numeric(sum(is.nan(est) | is.na(SEM)) > 0)) %>% 
  ungroup %>%
  group_by(scenario) %>%
  mutate(tot_perc_issue = sum(TF_issue)/n()) %>%
  group_by(scenario, tot_perc_issue, spp) %>%
  summarize(spp_perc_issue = sum(TF_issue)/n()) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x, 3))) %>%
  pivot_longer(cols = c("tot_perc_issue", "spp_perc_issue")) %>%
  pivot_wider(names_from = scenario,
              values_from = value) %>%
  mutate(spp = if_else(name == "tot_perc_issue", "ALL", spp)) %>%
  distinct(spp, .keep_all = T) %>%
  arrange(desc(name))

# - convert these NaN's to NA -------------------------------------------- ####
# Note: there will still be NaN's for SE_perc when the SEM was 0 we want to 
# make sure we can differentiate these from this type of NA case... 
est_peca_cna <- est_peca_cna %>% map(function(sp.x) {
  sp.x %>%
    mutate(SEM_perc = if_else(is.nan(SEM), NA_real_, SEM_perc),
           est = if_else(is.nan(est), NA_real_, est),
           SEM = if_else(is.nan(SEM), NA_real_, SEM))
})

# - Q: how many mcmc were we able to make estimates for? ----------------- ####
est_peca_cna %>%
  bind_rows(.id = "spp") %>%
  filter(!is.na(est) & 
           !is.na(SEM) & 
           !stat_id %in% paste0("stat_", c((1:4)*2 -1))) %>%
  group_by(spp, scenario, peca_domain, var_y) %>%
  summarize(n_mcmc = n()) %>%
  mutate(scenario_short = case_when(
    scenario == "RCP 4.5 CCSM4" ~ "c45",
    scenario == "RCP 8.5 CCSM4" ~ "c85",
    scenario == "RCP 4.5 HadGEM2-ES" ~ "h45",
    scenario == "RCP 8.5 HadGEM2-ES" ~ "h85")) %>%
  ungroup %>%
  select(-scenario) %>%
  pivot_wider(names_from = scenario_short,
              values_from = n_mcmc,
              values_fill = 0) %>%
  mutate(nMCMC_c45_h45_c85_h85 = paste(c45, h45, c85, h85, 
                                       sep = "|")) %>%
  select(-c45, -h45, -c85, -h85) %>%
  pivot_wider(names_from = var_y,
              values_from = nMCMC_c45_h45_c85_h85)

# == Overall estimates =================================================== ####
# - overall (range-wide): table ------------------------------------------ ####
est_overall %>%
  # for each spp, calculate 90% CI for each estimate and classify 
  map(function(sp.x, z_crit2 = qnorm(0.95)) {
    sp.x %>%
      filter(!stat_id %in% paste0("stat_", c((1:4)*2 - 1))) %>%
      separate(col = var_y,
               into = c("stat", "units"),
               sep = "_") %>%
      select(stat, units, est, SEM) %>%
      mutate(stat = case_when(grepl("mort", stat) ~ "m",
                              grepl("rel", stat) ~ "rg",
                              TRUE ~ "ng")) %>%
      filter(!is.na(est) & !is.na(SEM)) %>%
      mutate(cil = est - SEM * z_crit2,
             ciu = est + SEM * z_crit2) %>%
      mutate(est = if_else(stat == "rg", round(est*100), round(est, 2)),
             cil = if_else(stat == "rg", round(cil*100), round(cil, 2)),
             ciu = if_else(stat == "rg", round(ciu*100), round(ciu, 2))) %>%
      mutate(value = paste0(est,  " (", cil,  ", ", ciu,  ")")) %>%
      select(units, stat, value) %>%
      pivot_wider(names_from = 'stat',
                  values_from = 'value')
  }) %>%
  bind_rows(.id = "spp") %>%
  select(spp, units, ng, m, rg)

# == naive-divisions estimates =========================================== ####
# - prevalence heatmap & table ------------------------------------------- ####
# Within each naive-division, what was the number of plots on which the 
# species was present at t1? the prevalence? 
# - the idea here is to see which estimates have very prevalence going into
#   a statistics estimate so that we can have a better understanding of why 
#   some estimates have very wide ci's than others...
elevlat %>% 
  map(function(sp.x) {
    sp.x %>%
      filter(!is.na(el_group)) %>%
      group_by(el_group) %>%
      mutate(p_a = if_else(is.na(p_a), 0, p_a)) %>%
      summarise(n_p_plots = sum(p_a),
                prev = (sum(p_a)/n() * 100) %>% round(1)) %>%
      pivot_longer(cols = c("n_p_plots", "prev"),
                   names_to = "t1_stat",
                   values_to = "value") %>%
      pivot_wider(names_from = el_group,
                  values_from = value) %>%
      select(t1_stat, LL, LM, LH, ML, MM, MH, HL, HM, HH)
  }) %>% 
  bind_rows(.id = "spp") %>% 
  split(.$t1_stat)

# make tiles with no p-plots red
# for easier interpretation... show as a heat map
p_heatmap <- elevlat %>% 
  imap(function(sp.x, sp_name) {
    df.x <- sp.x %>%
      filter(!is.na(el_group)) %>%
      group_by(el_group) %>%
      mutate(p_a = if_else(is.na(p_a), 0, p_a)) %>%
      summarise(n_p_plots = sum(p_a),
                prevalence = (sum(p_a)/n() * 100) %>% round(1)) %>%
      ungroup %>%
      mutate(elev_p = substr(el_group, 1, 1),
             lat_p  = substr(el_group, 2, 2),
             elev_p = factor(elev_p, levels = c("L", "M", "H")),
             lat_p = factor(lat_p, levels = c("L", "M", "H")),
             spp = sp_name) %>%
      mutate(spp_common = case_when(spp == "ABPR" ~ "noble fir",
                                    spp == "PSMEM" ~ "coastal Douglas-fir",
                                    spp == "QUDO" ~ "blue oak",
                                    spp == "QUGA4" ~ "white oak",
                                    spp == "QUKE" ~ "black oak"))
    
    ggplot(df.x, aes(x = elev_p, y = lat_p, fill = prevalence)) +
      geom_tile() +
      geom_tile(data = df.x %>% filter(n_p_plots == 0),
                fill = "#cc0000") +
      facet_grid(. ~ spp_common) +
      scale_fill_gradient2(low = "white", high = "black") +
      theme_bw() +
      theme(legend.position = "bottom",
            aspect.ratio = 1,
            legend.box.margin = margin(-10, 0, 0, 0)) +
      labs(x = "elevation",
           y = "latitude") 
  }) %>% 
  imap(function(p.x, sp_name) {
    if (sp_name == "ABPR") {
      p.x + 
        theme(plot.margin = margin(1, 1, 0, 1),
              panel.grid = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = FF_font_size),
              strip.text.x = element_text(size = FF_font_size),
              axis.text.x  = element_text(size = FF_font_size),
              axis.text.y  = element_text(size = FF_font_size),
              legend.title = element_text(size = FF_font_size),
              legend.text  = element_text(size = FF_font_size)) +
        guides(fill = guide_colourbar(barwidth = 4.5, 
                                      title.position="top",
                                      barheight = 0.75,
                                      frame.colour = "black"))
    } else if (sp_name == "QUDO") {
      p.x + 
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = margin(1, 1, 0, 0),
              panel.grid = element_blank(),
              axis.title.x = element_text(size = FF_font_size),
              strip.text.x = element_text(size = FF_font_size),
              axis.text.x  = element_text(size = FF_font_size),
              legend.text  = element_text(size = FF_font_size)) +
        guides(fill = guide_colorbar(title = "",
                                     barwidth = 4.5, 
                                     title.position="top",
                                     barheight = 0.75,
                                     frame.colour = "black"))
    } else {
      p.x + 
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = margin(1, 1, 0, 0),
              panel.grid = element_blank(),
              strip.text.x = element_text(size = FF_font_size),
              axis.text.x  = element_text(size = FF_font_size),
              legend.text  = element_text(size = FF_font_size))+
        guides(fill = guide_colorbar(title = "",
                                     barwidth = 4.5, 
                                     title.position="top",
                                     barheight = 0.75,
                                     frame.colour = "black"))
    }
  }) %>%
  reduce(., `|`) 

ggsave(paste0(path_images, "heatmap_naive-division-prev-at-t1_all-spp.jpg"),
       plot = p_heatmap,
       height = 70, 
       width = FF_fig_width_2col, 
       units = FF_width_unit, 
       dpi = FF_min_dpi)

p_heatmap +
  plot_annotation(title = "Spp prevalence at t1 within naive-divisions",
                  caption = "Red if no presence-plots within division")

# - map of just naive-divisions (same color pallet as scatter) ----------- ####
 # join spatial info and plot/map
elevlat %>%
  # subset elevlat to remove plots without an assigned plot_id (this just 
  # for plotting - using Voronoi polygons dev in ch2 for vizualization)
  map(function(sp.x) {
    sp.x %>%
      select(plot_id, el_group) %>%
      filter(!is.na(plot_id)) %>%
      mutate(el_group = factor(el_group,
                               levels = c("LL", "LM", "LH",
                                          "ML", "MM", "MH",
                                          "HL", "HM", "HH")))
  }) %>% 
  # join with (public) spatial coords
  list(., plot_polygons_public.sf) %>%
  pmap(function(df.x, sf.x) {
    sf.x %>%
      select(-NAME, -p_a) %>%
      inner_join(df.x, by = "plot_id")
  }) %>%
  # make the base plots
  imap(function(sp.x, sp_name) {
    ggplot(sp.x %>% mutate(spp = sp_name), aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, 
              fill = "grey50") +
      geom_sf(aes(fill = el_group), color = NA) +
      facet_grid(. ~ spp) +
      scale_fill_viridis_d("naïve-divisions",
                           labels = c("LL", "LM", "LH",
                                      "ML", "MM", "MH",
                                      "HL", "HM", "HH")) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw()
  }) %>%
  # adjust the plots depending on name
  imap(function(plot.x, sp_name) {
    if (sp_name == "ABPR") {
      plot.x + 
        theme(plot.margin = margin(1, 1, 0, 1))
    } else {
      plot.x + 
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = margin(1, 1, 0, 0))
    }
  }) %>%
  reduce(., `|`) +
  # plot_annotation(title = "naïve-divisions") + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") &
  guides(fill = guide_legend(nrow = 3, byrow = F))


# - plot of nd by elev and lat ------------------------------------------- ####
((elevlat %>% 
    map(function(sp.x) {
      sp.x %>%
        select(el_group, plot_id, lat_actual, elev_actual) %>%
        filter(!is.na(plot_id)) %>%
        mutate(el_group = factor(el_group,
                                 levels = c("LL", "LM", "LH",
                                            "ML", "MM", "MH",
                                            "HL", "HM", "HH")),
               el_code = case_when(el_group == "LL" ~ 1,
                                   el_group == "LM" ~ 2,
                                   el_group == "LH" ~ 3,
                                   el_group == "ML" ~ 4,
                                   el_group == "MM" ~ 5,
                                   el_group == "MH" ~ 6,
                                   el_group == "HL" ~ 7,
                                   el_group == "HM" ~ 8,
                                   el_group == "HH" ~ 9),
               elev_km = elev_actual/1000) %>%
        select(-elev_actual)
    }) %>% 
    bind_rows(.id = "spp")) %>%
   ggplot() +
   stat_summary_hex(aes(y = lat_actual, x = elev_km, z = el_code),
                    # this function grabs the most common division in a hex
                    fun = function(x) {names(table(x))[which.max(table(x))]}) +
   facet_wrap(spp ~ ., scales = "free") +
   scale_fill_viridis_d(
     name = "naïve-division",
     breaks = paste(1:9),
     labels = c("LL", "LM", "LH",
                "ML", "MM", "MH",
                "HL", "HM", "HH")) +
   theme_bw() +
   theme(aspect.ratio = 1,
         panel.grid = element_blank(),
         panel.background = element_rect(fill = "grey50")) + 
   labs(x = "Elevation (km)",
        y = "Latitude (°N)",
        title = "Most common naïve-division") +
   guides(fill = guide_legend(nrow = 3,
                              title.position="top",
                              title.hjust = 0.5))) %>%
  lemon::reposition_legend(.,
                           panel = "panel-3-2",
                           just = c(0, -0.5),
                           x = 0, y = 0)

# - where est were `NA` (ng/m) or `-Inf` (rg) ---------------------------- ####
# What percent of the stats we're interested were 'NaN' (meaning that an est
# in that sp x el_group combo was not possible)? Only return sp x el_group
# combos where perc > 0%
est_elevlat %>% map(function(sp.x) {
  sp.x %>% map(function(df.x) {
    df.x %>%
      filter(!stat_id %in% paste0("stat_", c((1:4)*2 -1))) %>%
      summarize(perc_est_NaN = sum(is.nan(est))/ n() * 100)
  }) %>%
    bind_rows(.id = "el_group") %>%
    filter(perc_est_NaN > 0)
}) %>% bind_rows(.id = "spp")
# - looks like the number of plots on which these species were present at t1 
#   was 0 for the respective naive-divisions with NaN ests:
elevlat %>% 
  map(function(sp.x) {
    sp.x %>%
      filter(!is.na(el_group)) %>%
      group_by(el_group) %>%
      mutate(p_a = if_else(is.na(p_a), 0, p_a)) %>%
      summarise(n_p_plots = sum(p_a),
                prev = (sum(p_a)/n() * 100) %>% round(1))
  }) %>% 
  bind_rows(.id = "spp") %>% 
  filter((spp == "ABPR" & el_group == "LL") |
           spp == "QUDO" & el_group == "HH" |
           spp == "QUKE" & el_group == "HH")

# Relative-growth & '-Inf' est 
# - This happened for 4 species x naive-division combinations (shown below)
#   when TPH/yr-est was not possible b/c gross-growth was 0 TPH/yr; that is,
#   all trees died and were not replaced:
est_elevlat %>% map(function(sp.x) {
  sp.x %>% 
    map(~ .x %>%
          filter(!stat_id %in% paste0("stat_", c((1:4)*2 - 1))) %>%
          filter(is.infinite(est))) %>%
    bind_rows
}) %>% bind_rows(.id = "spp")

# - estimates table ------------------------------------------------------ ####
# reminder of what those naive-divisions are for each species
est_elevlat %>% 
  map_depth(2, ~.x %>% select(el_group_chr)) %>%
  map(~.x %>% bind_rows(.id = "el_group")) %>%
  bind_rows(.id = "spp") %>%
  distinct() %>%
  mutate(el_group_chr = gsub("m - ", "-", el_group_chr),
         el_group_chr = gsub("° - ", "-", el_group_chr)) %>%
  pivot_wider(names_from = spp,
              values_from = el_group_chr)

# Big table: 60 rows!!
est_elevlat %>%
  # for each spp & division, calculate 90% CI for each estimate and classify 
  map_depth(2, function(df.x, z_crit2 = qnorm(0.95)) {
    df.x %>%
      filter(!stat_id %in% paste0("stat_", c((1:4)*2 - 1))) %>%
      separate(col = var_y,
               into = c("stat", "units"),
               sep = "_") %>%
      select(el_group, el_group_chr, units, stat, est, SEM) %>%
      mutate(stat = case_when(grepl("mort", stat) ~ "m",
                              grepl("rel", stat) ~ "rg",
                              TRUE ~ "ng")) %>%
      mutate(cil = est - SEM * z_crit2,
             ciu = est + SEM * z_crit2) %>%
      mutate(est = if_else(stat == "rg", round(est*100), round(est, 2)),
             cil = if_else(stat == "rg", round(cil*100), round(cil, 2)),
             ciu = if_else(stat == "rg", round(ciu*100), round(ciu, 2))) %>%
      mutate(value = paste0(est,  " (", cil,  ", ", ciu,  ")")) %>%
      select(el_group, el_group_chr, units, stat, value) %>%
      pivot_wider(names_from = 'stat',
                  values_from = 'value') %>%
      mutate(el_group_chr = gsub("m - ", "-", el_group_chr),
             el_group_chr = gsub("° - ", "-", el_group_chr),
             el_group = factor(el_group,
                               levels = paste0(rep(c("L", "M", "H"), 
                                                   each = 3), 
                                               c("L", "M", "H"))))
  }) %>%
  # pull all naive-divisions together for a species into one df
  map(~.x %>% bind_rows(.id = "el_group")) %>%
  # now to one df for all spp
  bind_rows(.id = "spp") %>%
  select(spp, contains("group"), units, ng, m, rg) %>%
  arrange(spp, factor(el_group,
                      levels = paste0(rep(c("L", "M", "H"), 
                                          each = 3), 
                                      c("L", "M", "H"))), desc(units)) 
# %>%
# # un-comment this to save as csv (easier porting to manuscript tables)
#   write.csv(.,
#             file = paste0(path_output, "table_est_elevlat.csv"),
#             row.names = FALSE)



# - descriptions for manuscript ------------------------------------------ ####
# we'll use this format also for plotting
est_summary_el <- est_elevlat %>% 
  map(function(sp.x, z_crit2 = qnorm(0.95)) {
    sp.x %>%
      map(function(df.x) {
        df.x %>%
          # only our stats of interest...
          filter(!stat_id %in% paste0("stat_", c((1:4)*2 - 1))) %>%
          # pull apart var_y names for faceting
          separate(col = var_y,
                   into = c("stat", "stat_units"),
                   sep = "_") %>%
          select(el_group, stat, stat_units, est, SEM) %>%
          # drop the naive-divisions where we couldn't make the estimate
          # - note we'll need to drop some naive-divisions for rg later where
          #   we could make a BAH/yr est but not a TPH/yr est (b/c it all died)
          filter(!is.na(est) & !is.na(SEM)) %>%
          # el_group as factor for plotting
          mutate(el_group = factor(el_group,
                                   levels = paste0(rep(c("L", "M", "H"), 
                                                       each = 3), 
                                                   c("L", "M", "H")))) %>%
          # more complete names for statistics and simpler names for units
          mutate(stat = case_when(grepl("mort", stat) ~ "Mortality",
                                  grepl("rel", stat) ~ "Relative-growth",
                                  TRUE ~ "Net-growth"),
                 stat_units = if_else(grepl("B", stat_units), 
                                      "b", 
                                      "t")) %>%
          # calc bounds for error-bars in the plot...
          mutate(cil = est - SEM * z_crit2,
                 ciu = est + SEM * z_crit2) %>%
          select(-SEM) %>%
          pivot_wider(names_from = "stat_units",
                      values_from = c("est", "cil", "ciu"))
      }) %>%
      bind_rows() %>%
      # Not sure why... but didn't let me do this within the last map, but ok
      # here? Anyways --- drop the four cases where the rg TPH est was -Inf 
      # (i.e., all trees died so gross growth is zero) these pop out again 
      # after we pivot wider within that last map...
      filter(!is.na(est_t))
  }) 

# Short reference for where:
# - populations are increasing (ng BAH/yr & TPH/yr sig positive)
# - populations are decreasing (ng BAH/yr sig negative; mort unsustainable)
# - mortality not offset by ingrowth stems (ng TPH/yr sig negative)
est_summary_el %>%
  imap(function(sp.x, sp_name) {
    sp.x <- sp.x %>% filter(stat == "Net-growth")
    
    list(sp.x %>% filter(cil_b > 0 & cil_t > 0),
         sp.x %>% filter(ciu_b < 0 & ciu_t < 0),
         sp.x %>% filter(ciu_t < 0)) %>%
      setNames(c("pop_inc", "pop_dec", "mort_exceed_ingrowth")) %>%
      map(~.x %>% select(el_group) %>% mutate(yes = sp_name)) %>%
      bind_rows(.id = "stat") %>%
      mutate(elev_g = substr(el_group, 1, 1),
             lat_g = substr(el_group, 2, 2)) %>%
      select(-el_group)
  }) %>%
  bind_rows(.id = "dummy") %>%
  split(.$stat) %>%
  map(function(df.x) {
    df.x %>%
      group_by(elev_g, lat_g) %>%
      summarize(spp_T = paste0(yes, collapse = "|")) %>%
      pivot_wider(names_from = elev_g, 
                  names_prefix = "elev_",
                  values_from = spp_T) %>% 
      ungroup %>%
      rename(g = lat_g) %>%
      mutate(g = paste0("lat_", g),
             g = factor(g, levels = c("lat_L", "lat_M", "lat_H"))) %>%
      select(g, elev_L, contains("elev_M"), elev_H) %>%
      arrange(desc(factor(g)))
  }) %>%
  rev

# where was ng BAH/yr & TPH/yr both sig positive
est_summary_el %>%
  map(~ .x %>% filter(cil_b > 0 & cil_t > 0 & stat == "Net-growth")) %>%
  bind_rows(.id = "spp")

# where was mortality unsustainable (neg ng BAH/yr & TPH/yr)
est_summary_el %>%
  map(~ .x %>% filter(ciu_b < 0 & ciu_t < 0 & stat == "Net-growth")) %>%
  bind_rows(.id = "spp")

# where did ingrowth not offset the number of mortality stems?
est_summary_el %>%
  map(~ .x %>% filter(ciu_t < 0 & stat == "Net-growth")) %>%
  bind_rows(.id = "spp")

# what was the rg where the above was true (that is, if rg was able to be
# calculated...)?
est_summary_el %>% map(function(sp.x) {
    sp.x %>% 
      filter(el_group %in% 
               (sp.x %>% filter(ciu_t < 0 & stat == "Net-growth"))$el_group &
               stat == "Relative-growth") %>%
      select(el_group, contains("_t"))
  }) %>% bind_rows(.id = "spp")

# - naive-divisions: scatter-plot by ng/m/rg ----------------------------- ####
# summarize as a figure by statistic (3 figs) 
# - here, coloring by naive-division which are ordered by elevation-partition, 
#   so also use symbols to indicate the latitude-partition.
# (facet_wrap)
# - (using lemon here to place legend within facets)
est_summary_el %>%
  bind_rows(.id = "spp") %>%
  split(.$stat) %>%
  imap(function(df.x, stat_name) {
    # - set plotting arguments & prep data -------------------------- ####
    alpha_set <- 0.75
    
    x_lab <- ifelse(stat_name == "Relative-growth", 
                    "TPH/yr percent (hundreds)", 
                    "TPH/yr")
    y_lab <- ifelse(stat_name == "Relative-growth", 
                    "BAH/yr percent (hundreds)", 
                    "BAH/yr")
    
    df.x <- df.x %>% 
      mutate(spp_common = case_when(spp == "ABPR" ~ "noble fir",
                                    spp == "PSMEM" ~ "coastal Douglas-fir",
                                    spp == "QUDO" ~ "blue oak",
                                    spp == "QUGA4" ~ "white oak",
                                    spp == "QUKE" ~ "black oak"),
             spp_common = factor(spp_common,
                                 levels = c("noble fir",
                                            "coastal Douglas-fir",
                                            "blue oak",
                                            "white oak",
                                            "black oak")))
    
    # - make the plot ----------------------------------------------- ####
    p.x <- ggplot(df.x, aes(x = est_t,  y = est_b, 
                            fill = el_group, 
                            shape = el_group)) +
      geom_hline(yintercept = 0) + # was linetype = 'dotted' before...
      geom_vline(xintercept = 0) +
      geom_errorbar(aes(ymax = ciu_b, ymin = cil_b, color = el_group),
                    alpha = alpha_set,
                    size = .75) +
      geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t, color = el_group),
                     alpha = alpha_set,
                     size = .75) +
      geom_point(alpha = alpha_set, size = 2.5) +
      facet_wrap(spp_common ~ ., scales = "free", nrow = 2) +
      scale_fill_viridis_d("Naïve-divisions") +
      scale_color_viridis_d("Naïve-divisions") +
      scale_shape_manual("Naïve-divisions",
                         # values = c(16, 15, 17, 16, 15, 17, 16, 15, 17),
                         values = c(21, 22, 24, 21, 22, 24, 21, 22, 24),
                         labels = levels(df.x$el_group)) +
      labs(y = y_lab, x = x_lab) + 
      theme_bw()  +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "grey50"),
            strip.text.x = element_text(size = FF_font_size),
            axis.title.y = element_text(size = FF_font_size),
            axis.title.x = element_text(size = FF_font_size),
            axis.text.x  = element_text(size = FF_font_size),
            axis.text.y  = element_text(size = FF_font_size),
            legend.title = element_text(size = FF_font_size),
            legend.text  = element_text(size = FF_font_size)) +
      guides(fill = guide_legend(nrow = 3,
                                  title.position="top",
                                  title.hjust = 0.5))
    
    # - move legend over with lemon & ggsave ------------------------ ####
    ggsave(filename = paste0(path_images, 
                             "scatter-wrap_naive-divisions_", 
                             stat_name, 
                             ".jpg"),
           plot = lemon::reposition_legend(p.x,
                                           panel = "panel-3-2",
                                           just = c(-0.1, -0.2),
                                           x = 0, 
                                           y = 0),
           height = 100, 
           width = FF_fig_width_2col, 
           units = FF_width_unit, 
           dpi = FF_min_dpi)
  }) 


# == PECA: scatter of ng/rg by scenario, MCMC, & PECA-division =========== ####
# - Common axes by spp --------------------------------------------------- #### 
est_peca_cna %>% 
  # format the data
  # - this will give a 90% ci (1.645), a 95% would be 1.960...
  imap(function(sp.x, sp_name, z_crit2 = qnorm(0.95)) {
    sp.x <- sp.x %>%
      ungroup() %>%
      # only our net-growth stats
      filter(stat_id %in% paste0("stat_", c(6, 8:10))) %>%
      select(mcmc_n, scenario, peca_domain, var_y, est, SEM) %>%
      # pull apart var_y names for faceting
      separate(col = var_y,
               into = c("stat", "stat_units"),
               sep = "_") %>%
      mutate(spp = toupper(sp_name),
             peca_domain = factor(peca_domain,
                                  levels = c("E", "P", "C", "A")),
             stat = case_when(grepl("mort", stat) ~ "m",
                              grepl("rel", stat) ~ "rng",
                              TRUE ~ "ng"),
             stat_units = paste0(gsub("ann", "", stat_units), "/yr")) %>%
      # drop MCMC x PECA where we couldn't make the estimate
      filter(!is.na(est) & !is.na(SEM)) %>%
      # note - SEM: standard error of the mean
      mutate(cil = est - SEM * z_crit2,
             ciu = est + SEM * z_crit2) %>%
      select(-SEM) %>% 
      mutate(spp_common = case_when(spp == "ABPR" ~ "noble fir",
                                    spp == "PSMEM" ~ "coastal Douglas-fir",
                                    spp == "QUDO" ~ "blue oak",
                                    spp == "QUGA4" ~ "white oak",
                                    spp == "QUKE" ~ "black oak"),
             spp_common = factor(spp_common,
                                 levels = c("noble fir",
                                            "coastal Douglas-fir",
                                            "blue oak",
                                            "white oak",
                                            "black oak")))
    
    # split into cols of our final figure and drop MCMC's we cannot
    # plot (b/c missing x or y var)
    # - panel 1: cols 1 (TPH/yr by M v Ng) & 2 (BAH/yr by M v Ng)
    # - panel 2: col 1 rel-ng and ng (each by TPH/yr v BAH/yr)
    list(p_ng = sp.x %>%
           filter(stat == "ng") %>%
           mutate(stat = "Net-growth",
                  stat_units = if_else(grepl("B", stat_units), 
                                       "b", 
                                       "t")) %>%
           pivot_wider(names_from = "stat_units",
                       values_from = c("est", "cil", "ciu")),
         p_rg = sp.x %>% 
           filter(stat == "rng") %>%
           mutate(est = est, #est * 100,
                  cil = cil, #cil * 100,
                  ciu = ciu, #ciu * 100,
                  stat_units = if_else(grepl("B", stat_units), 
                                       "b", 
                                       "t"),
                  stat = "Relative-growth") %>%
           pivot_wider(names_from = "stat_units",
                       values_from = c("est", "cil", "ciu")))
  }) %>%
  transpose %>%
  map_depth(2, ~.x[complete.cases(.x),]) %>%
  # go by stat (ng/rng) & spp
  map_depth(2, function(sp.x) {
    # - set plotting arguments -------------------------------------- ####
    std_margin <- margin(t = 0, r = 0, b = 1, l = 0)
    alpha_set <- 0.25
    x_lab <- ifelse(sp.x$stat[1] == "Net-growth", 
                    "TPH/yr", 
                    "TPH/yr percent (hundreds)")
    y_lab <- ifelse(sp.x$stat[1] == "Net-growth", 
                    "BAH/yr", 
                    "BAH/yr percent (hundreds)")
    
    # - make the plots by spp with switch --------------------------- ####
    # Note - QUGA4 breaks problematic with Net-growth, so need set those 
    # breaks explicitly
    switch(sp.x$spp[1],
           "ABPR" = {
             # top: no axis label
             sp.x %>%
               ggplot(aes(x = est_t, y = est_b,
                          color = peca_domain, group = mcmc_n)) +
               geom_point(alpha = alpha_set) +
               geom_hline(yintercept = 0) + # was linetype = 'dotted'
               geom_vline(xintercept = 0) +
               geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t),
                              alpha = alpha_set) +
               geom_errorbar(aes(ymax = ciu_b, ymin = cil_b),
                             alpha = alpha_set) +
               scale_color_manual("PECA-division",
                                  values = peca_domain_palette) +
               facet_grid(spp_common ~ scenario) +
               labs(y = y_lab, x = x_lab) + 
               theme_bw()  +
               theme(plot.margin = std_margin,
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     strip.text.x = element_text(size = FF_font_size),
                     strip.text.y = element_text(size = FF_font_size),
                     # axis.title.y = element_text(size = FF_font_size),
                     # axis.title.x = element_text(size = FF_font_size),
                     axis.text.x  = element_text(size = FF_font_size),
                     axis.text.y  = element_text(size = FF_font_size),
                     legend.title = element_text(size = FF_font_size),
                     legend.text  = element_text(size = FF_font_size)) +
               scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1)) + 
               scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1))
           },
           "PSMEM" = {
             # mid: no x-axis labels, keep x facet
             sp.x %>%
               ggplot(aes(x = est_t, y = est_b,
                          color = peca_domain, group = mcmc_n)) +
               geom_point(alpha = alpha_set) +
               geom_hline(yintercept = 0) + # was linetype = 'dotted'
               geom_vline(xintercept = 0) +
               geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t),
                              alpha = alpha_set) +
               geom_errorbar(aes(ymax = ciu_b, ymin = cil_b),
                             alpha = alpha_set) +
               scale_color_manual("PECA-division",
                                  values = peca_domain_palette) +
               facet_grid(spp_common ~ scenario) +
               labs(y = y_lab, x = x_lab) + 
               theme_bw()  +
               theme(plot.margin = std_margin,
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     strip.text.x = element_blank(),
                     # strip.text.x = element_text(size = FF_font_size),
                     strip.text.y = element_text(size = FF_font_size),
                     # axis.title.y = element_text(size = FF_font_size),
                     # axis.title.x = element_text(size = FF_font_size),
                     axis.text.x  = element_text(size = FF_font_size),
                     axis.text.y  = element_text(size = FF_font_size),
                     legend.title = element_text(size = FF_font_size),
                     legend.text  = element_text(size = FF_font_size)) + 
               scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1)) + 
               scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1))
           },
           "QUDO" = {
             # mid: no x-axis labels, no y-facet
             sp.x %>%
               ggplot(aes(x = est_t, y = est_b,
                          color = peca_domain, group = mcmc_n)) +
               geom_point(alpha = alpha_set) +
               geom_hline(yintercept = 0) + # was linetype = 'dotted'
               geom_vline(xintercept = 0) +
               geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t),
                              alpha = alpha_set) +
               geom_errorbar(aes(ymax = ciu_b, ymin = cil_b),
                             alpha = alpha_set) +
               scale_color_manual("PECA-division",
                                  values = peca_domain_palette) +
               facet_grid(spp_common ~ scenario) +
               labs(y = y_lab, x = x_lab) + 
               theme_bw()  +
               theme(plot.margin = std_margin,
                     axis.title.x = element_blank(),
                     strip.text.x = element_blank(),
                     # strip.text.x = element_text(size = FF_font_size),
                     strip.text.y = element_text(size = FF_font_size),
                     axis.title.y = element_text(size = FF_font_size),
                     # axis.title.x = element_text(size = FF_font_size),
                     axis.text.x  = element_text(size = FF_font_size),
                     axis.text.y  = element_text(size = FF_font_size),
                     legend.title = element_text(size = FF_font_size),
                     legend.text  = element_text(size = FF_font_size)) + 
               scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1)) + 
               scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1))
           },
           "QUGA4" = {
             if (sp.x$stat[1] == "Net-growth") {
               # mid: no x-axis labels, keep x facet
               sp.x %>%
                 ggplot(aes(x = est_t, y = est_b,
                            color = peca_domain, group = mcmc_n)) +
                 geom_point(alpha = alpha_set) +
                 geom_hline(yintercept = 0) + # was linetype = 'dotted'
                 geom_vline(xintercept = 0) +
                 geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t),
                                alpha = alpha_set) +
                 geom_errorbar(aes(ymax = ciu_b, ymin = cil_b),
                               alpha = alpha_set) +
                 scale_color_manual("PECA-division",
                                    values = peca_domain_palette) +
                 facet_grid(spp_common ~ scenario) +
                 labs(y = y_lab, x = x_lab) + 
                 theme_bw()  +
                 theme(plot.margin = std_margin,
                       axis.title.y = element_blank(),
                       axis.title.x = element_blank(),
                       strip.text.x = element_blank(),
                       # strip.text.x = element_text(size = FF_font_size),
                       strip.text.y = element_text(size = FF_font_size),
                       # axis.title.y = element_text(size = FF_font_size),
                       # axis.title.x = element_text(size = FF_font_size),
                       axis.text.x  = element_text(size = FF_font_size),
                       axis.text.y  = element_text(size = FF_font_size),
                       legend.title = element_text(size = FF_font_size),
                       legend.text  = element_text(size = FF_font_size)) + 
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                                    labels = scales::comma_format(accuracy = 1)) + 
                 scale_x_continuous(breaks = c(-20, 0, 20))
             } else {
               # mid: no x-axis labels, keep x facet
               sp.x %>%
                 ggplot(aes(x = est_t, y = est_b,
                            color = peca_domain, group = mcmc_n)) +
                 geom_point(alpha = alpha_set) +
                 geom_hline(yintercept = 0) + # was linetype = 'dotted'
                 geom_vline(xintercept = 0) +
                 geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t),
                                alpha = alpha_set) +
                 geom_errorbar(aes(ymax = ciu_b, ymin = cil_b),
                               alpha = alpha_set) +
                 scale_color_manual("PECA-division",
                                    values = peca_domain_palette) +
                 facet_grid(spp_common ~ scenario) +
                 labs(y = y_lab, x = x_lab) + 
                 theme_bw()  +
                 theme(plot.margin = std_margin,
                       axis.title.y = element_blank(),
                       axis.title.x = element_blank(),
                       strip.text.x = element_blank(),
                       # strip.text.x = element_text(size = FF_font_size),
                       strip.text.y = element_text(size = FF_font_size),
                       # axis.title.y = element_text(size = FF_font_size),
                       # axis.title.x = element_text(size = FF_font_size),
                       axis.text.x  = element_text(size = FF_font_size),
                       axis.text.y  = element_text(size = FF_font_size),
                       legend.title = element_text(size = FF_font_size),
                       legend.text  = element_text(size = FF_font_size)) + 
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                                    labels = scales::comma_format(accuracy = 1)) + 
                 scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                                    labels = scales::comma_format(accuracy = 1))
             }
           },
           "QUKE" = {
             # bottom: no y-axis label
             sp.x %>%
               ggplot(aes(x = est_t, y = est_b,
                          color = peca_domain, group = mcmc_n)) +
               geom_point(alpha = alpha_set) +
               geom_hline(yintercept = 0) + # was linetype = 'dotted'
               geom_vline(xintercept = 0) +
               geom_errorbarh(aes(xmax = ciu_t, xmin = cil_t),
                              alpha = alpha_set) +
               geom_errorbar(aes(ymax = ciu_b, ymin = cil_b),
                             alpha = alpha_set) +
               scale_color_manual("PECA-division",
                                  values = peca_domain_palette) +
               facet_grid(spp_common ~ scenario) +
               labs(y = y_lab, x = x_lab) + 
               theme_bw()  +
               theme(plot.margin = std_margin,
                     axis.title.y = element_blank(),
                     strip.text.x = element_blank(),
                     # strip.text.x = element_text(size = FF_font_size),
                     strip.text.y = element_text(size = FF_font_size),
                     # axis.title.y = element_text(size = FF_font_size),
                     axis.title.x = element_text(size = FF_font_size),
                     axis.text.x  = element_text(size = FF_font_size),
                     axis.text.y  = element_text(size = FF_font_size),
                     legend.title = element_text(size = FF_font_size),
                     legend.text  = element_text(size = FF_font_size)) + 
               scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1)) + 
               scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                                  labels = scales::comma_format(accuracy = 1))
           })
  }) %>%
  # make panels pretty
  map_depth(2, function(plot.x) {
    plot.x  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1)
  }) %>% 
  # use patchwork pkg to stitch plots together (not sure if | & / right below..)
  map(~.x %>% reduce(., `/`)) %>%
  setNames(c("net-growth", "relative-growth")) %>%
  # save figures to sorting folder (path_images)
  imap(function(plot_x, stat_name) {
    ggsave(filename = paste0(path_images, "scatter_", stat_name,
                             "_all-spp-scenarios.jpg"),
           plot = plot_x + 
             plot_layout(guides = "collect") & 
             theme(legend.margin = margin(0,0,0,0),
                   legend.spacing.x = unit(0, "mm"),
                   legend.spacing.y = unit(0, "mm")) &
             guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                              linetype = "blank",
                                                              shape = 15,
                                                              size = 4))) ,
           height = 200, 
           width = FF_fig_width_2col, 
           units = FF_width_unit, 
           dpi = FF_min_dpi)
  })

# == PECA: agreement on dir of change within divisions (ng/m) ============ ####
# The summary:
# - Based on 90% CI for the mean, classify mortality as {positive|zero} and
#   net-growth as {positive|zero|negative} for each MCMC sample. Then 
#   summarize the percent of MCMC samples within each future scenario, spp, 
#   and PECA-division for which mortality (or net-growth) was classified into
#   each of these categories.

# - calc summaries (ng/m) ------------------------------------------------ ####
dir_change <- est_peca_cna %>%
  # for each spp, calculate 90% CI for each estimate and classify each 
  # MCMC sample x PECA-domain...
  imap(function(sp.x, sp_name, z_crit2 = qnorm(0.95)) {
    sp.x %>% 
      ungroup() %>%
      filter(stat_id %in% paste0("stat_", (1:4)*2)) %>%
      select(scenario, mcmc_n, peca_domain, var_y, est, SEM) %>%
      # pull apart stat & units
      separate(col = var_y,
               into = c("stat", "stat_units"),
               sep = "_") %>%
      mutate(spp = toupper(sp_name),
             peca_domain = factor(peca_domain,
                                  levels = c("E", "P", "C", "A")),
             stat = if_else(grepl("mort", stat), "m", "ng"),
             stat_units = paste0(gsub("ann", "", stat_units), "/yr")) %>%
      # drop MCMC x PECA where we couldn't make the estimate
      filter(!is.na(est) & !is.na(SEM)) %>%
      # note - SEM: standard error of the mean
      mutate(cil = est - SEM * z_crit2,
             ciu = est + SEM * z_crit2) %>%
      select(spp, scenario, mcmc_n, peca_domain, stat, stat_units, cil, ciu) %>%
      group_by(spp, scenario, mcmc_n, peca_domain, stat_units, stat) %>%
      summarize(result_code = case_when(
        cil > 0 ~ "positive",
        ciu < 0 ~ "negative",
        cil <= 0 & ciu >= 0 ~ "zero",
        is.na(cil) & is.na(ciu) ~ "no est")) %>%
      ungroup %>%
      filter(result_code != "no est") %>%
      mutate(result_code = factor(result_code, levels = c("positive",
                                                          "zero",
                                                          "negative",
                                                          "no est")))
  }) %>%
  # summarize classes across MCMC samples for each stat X scenario x peca-domain
  map(function(sp.x) {
    sp.x %>%
      ungroup %>%
      group_by(spp, scenario, peca_domain, stat_units, stat, result_code) %>%
      tally(name = "n_mcmc") %>%
      ungroup() %>% 
      group_by(spp, scenario, peca_domain, stat_units, stat) %>%
      mutate(tot_mcmc = sum(n_mcmc),
             perc_of_tot = n_mcmc / tot_mcmc) %>%
      ungroup %>%
      mutate(peca_domain = factor(peca_domain,
                                  levels = rev(c("E", "P", "C", "A"))))
  }) %>%
  bind_rows %>%
  split(.$stat_units)

# - by mcmc sample, where there was unsustainable mortality -------------- ####
# the percent of MCMC samples with unsustainable mortality
est_peca_cna %>%
  # for each spp, calculate 90% CI for each estimate and classify each 
  # MCMC sample x PECA-domain...
  imap(function(sp.x, sp_name, z_crit2 = qnorm(0.95)) {
    sp.x %>% 
      ungroup() %>%
      filter(stat_id %in% paste0("stat_", (1:4)*2)) %>%
      select(scenario, mcmc_n, peca_domain, var_y, est, SEM) %>%
      # pull apart stat & units
      separate(col = var_y,
               into = c("stat", "stat_units"),
               sep = "_") %>%
      mutate(spp = toupper(sp_name),
             peca_domain = factor(peca_domain,
                                  levels = c("E", "P", "C", "A")),
             stat = if_else(grepl("mort", stat), "m", "ng"),
             stat_units = paste0(gsub("ann", "", stat_units), "/yr")) %>%
      # drop MCMC x PECA where we couldn't make the estimate
      filter(!is.na(est) & !is.na(SEM)) %>%
      # note - SEM: standard error of the mean
      mutate(cil = est - SEM * z_crit2,
             ciu = est + SEM * z_crit2) %>%
      select(spp, scenario, mcmc_n, peca_domain, stat, stat_units, cil, ciu) %>%
      group_by(spp, scenario, mcmc_n, peca_domain, stat_units, stat) %>%
      summarize(result_code = case_when(
        cil > 0 ~ "positive",
        ciu < 0 ~ "negative",
        cil <= 0 & ciu >= 0 ~ "zero",
        is.na(cil) & is.na(ciu) ~ "no est")) %>%
      ungroup %>%
      filter(result_code != "no est") %>%
      mutate(result_code = factor(result_code, levels = c("positive",
                                                          "zero",
                                                          "negative",
                                                          "no est")))
  }) %>%
  # summarize classes across MCMC samples for each stat X scenario x peca-domain
  # - different here: summarize if mort was sustainable or not
  map(function(sp.x) {
    sp.x %>%
      ungroup %>%
      mutate(stat_units = if_else(stat_units == "BAH/yr", 
                                  "units_b", 
                                  "units_t")) %>%
      filter(stat == "ng") %>%
      group_by(spp, scenario, mcmc_n, peca_domain) %>%
      pivot_wider(names_from = stat_units,
                  values_from = result_code) %>%
      filter(!is.na(units_b) & !is.na(units_t)) %>%
      summarize(unsus_mort = if_else(units_b == "negative" & 
                                       units_t == "negative",
                                             T,
                                             F)) %>%
      ungroup %>%
      group_by(spp, scenario, peca_domain, unsus_mort) %>%
      tally(name = "n_mcmc") %>%
      ungroup() %>% 
      group_by(spp, scenario, peca_domain) %>%
      mutate(tot_mcmc = sum(n_mcmc),
             perc_of_tot = (n_mcmc / tot_mcmc) * 100) %>%
      ungroup %>%
      mutate(peca_domain = factor(peca_domain,
                                  levels = rev(c("E", "P", "C", "A"))),
             spp_common = case_when(spp == "ABPR" ~ "noble fir",
                                    spp == "PSMEM" ~ "coastal Douglas-fir",
                                    spp == "QUDO" ~ "blue oak",
                                    spp == "QUGA4" ~ "white oak",
                                    spp == "QUKE" ~ "black oak"),
             spp_common = factor(spp_common,
                                 levels = c("noble fir",
                                            "coastal Douglas-fir",
                                            "blue oak",
                                            "white oak",
                                            "black oak")))
  }) %>%
  bind_rows %>% 
  filter(unsus_mort == T) %>%
  filter(perc_of_tot > 0) %>%
  mutate(perc_of_tot = round(perc_of_tot),
         value = paste0(perc_of_tot, "%; ", n_mcmc, " of ", tot_mcmc)) %>%
  select(spp_common, peca_domain, value, scenario) %>%
  pivot_wider(names_from = peca_domain,
              values_from = value)

# - make figure for each scenario ---------------------------------------- ####
# Create separate figures for each future scenario that summarize the 
# direction of change (percent of MCMC samples in each category) by spp & 
# statistic (net-growth or mortality) scenario
dir_change %>%
  map(~.x %>% split(.$scenario)) %>%
  transpose() %>%
  map(function(scen.x) {
    scen.x %>%
      bind_rows(.id = "units_chr") %>%
      mutate(perc_of_tot = perc_of_tot *100,
             stat = case_when(stat == "m" ~ "mortality",
                              stat == "ng" ~ "net-growth"),
             stat = paste0(stat, "\n", units_chr),
             stat = factor(stat, levels = c("net-growth\nBAH/yr",
                                            "mortality\nBAH/yr",
                                            "net-growth\nTPH/yr",
                                            "mortality\nTPH/yr")),
             spp_common = case_when(spp == "ABPR" ~ "noble fir",
                                    spp == "PSMEM" ~ "coastal Douglas-fir",
                                    spp == "QUDO" ~ "blue oak",
                                    spp == "QUGA4" ~ "white oak",
                                    spp == "QUKE" ~ "black oak"),
             spp_common = factor(spp_common,
                                 levels = c("noble fir",
                                            "coastal Douglas-fir",
                                            "blue oak",
                                            "white oak",
                                            "black oak")))
  }) %>%
  # make the plots
  map(function(df.x) {
    df.x %>%
      ggplot() +
      geom_bar(aes(y = peca_domain, 
                   x = perc_of_tot, 
                   fill = as.factor(result_code)),
               stat = 'identity',
               position = 'stack') +
      scale_fill_viridis_d("result code",
                           option = "D",
                           direction = -1) +
      scale_x_continuous(breaks = c(25, 75)) + 
      facet_grid(stat ~ spp_common) +
      theme_bw() +
      labs(x = "Percent of MCMC samples",
           y = "PECA-division") +
      theme(panel.grid = element_blank(),
            # legend.position = "bottom",
            # aspect.ratio = 1,
            plot.margin = margin(1, 1, 2, 1),
            strip.text.x = element_text(size = FF_font_size),
            strip.text.y = element_text(size = FF_font_size),
            axis.title.y = element_text(size = FF_font_size),
            axis.title.x = element_text(size = FF_font_size),
            axis.text.x  = element_text(size = FF_font_size),
            axis.text.y  = element_text(size = FF_font_size),
            legend.title = element_text(size = FF_font_size),
            legend.text  = element_text(size = FF_font_size))
  }) %>%
  # ggsave the plots
  imap(function(p.x, scen_name.x) {
    ggsave(filename = paste0(path_images, 
                             "barplot_dir-change-peca_",
                             scen_name.x,
                             ".jpg"),
           plot = p.x,
           height = 120, 
           width = FF_fig_width_2col, 
           units = FF_width_unit, 
           dpi = FF_min_dpi)
  })

# == PECA: agreement on amt of change between divisions (ng/m) =========== ####
# About: We would expect for estimates of mortality and net-growth to be 
# different between areas of P vs C vs E... More specifically, we would expect
# for mortality to be higher and net-growth to be lower in areas projected
# to {C vs P, P vs E, C vs E}. So our questions are:
# - Q: did those estimates differ?
# - Q: what was the direction of the difference? (prev section)
# - Q: across the MCMC samples, how often do these estimates differ? 
#     (this section)
#
# Notes:
# - Here, we are interested in comparing mort / net-growth at areas where 
#   the sp was present at either t1 (mortality) or t1|t2 (net-growth).
#   (that is stat_id's 2/4/6/8)
# - ASSUMPTION - making assumption for two cases when both SEM = 0...
#   When the SEMs for both estimates (e.g., C's est & P's est) are 0, that 
#   means that there was (little) to no variation within that sample, and
#   so we make the following assumptions for two cases (where a z-test
#   cannot be performed as it will result in a value of INF):
#   - case 1: both SEM=0 and the est's are equivalent
#     - Assume not sig diff (for example, this might happen when there was no 
#       mortality on all plots within that PECA-domain over that time period).
#       This will be equivalent to an observed z-stat of 0.
#   - case 2: both SEM = 0, but the estimates are not equivalent
#     - Assume the estimates were significantly different. (for example, this 
#       could happen when one PECA-domain is like those described in case 1,
#       while all plots in the other domain had, 1 tree die on every plot).
#       For this case (and purposes of visualization/calculation) we will
#       assign an observed z-stat of `Inf` (well beyond the alpha-level we're 
#       considering) - important to note that we don't actually know what that 
#       z-stat would be, but this will allow us to identify these cases later
#       on for plotting purposes (and to make sure the 1-tailed tests are coded
#       according to our assumption).
# - calc the test statistic ---------------------------------------------- ####
# Calculate the z-statistic for the difference between estimates
# - note these z-statistics could technically be used with either a 1-tailed
#   or a 2-tailed test, however we are interested in 1-tailed test here
# - On the direction of the difference:
#   in this object differences are calculated in PECA order, that is
#   {P - E; P - C; E - C} so that we are consistently calculating these
#   differences. We'll switch around the direction of the differences to test
#   what we want it to with the 1-tailed test code in the next section
dif_z_stat_cna <- est_peca_cna %>% 
  map(~.x %>% split(.$scenario)) %>%
  map_depth(2, function(sp.x) {
    sp.x %>%
      ungroup %>%
      select(mcmc_n, stat_id, peca_domain, est, SEM) %>%
      filter(stat_id %in% paste0("stat_", (1:4)*2),
             peca_domain != "A") %>%
      split(.$mcmc_n) %>%
      map(~.x %>% split(.$stat_id)) %>%
      # split into domain comparisons to test
      map_depth(2, function(mcmc_stat) {
        # PECA order always so we know the direction of the difference,
        # (we'll fix the direction later in the script, so that it is
        #  testing what we think for our 1-tailed z-tests)...
        list("domain_PvE" = mcmc_stat %>% filter(peca_domain %in% c("P", "E")),
             "domain_PvC" = mcmc_stat %>% filter(peca_domain %in% c("P", "C")),
             "domain_EvC" = mcmc_stat %>% filter(peca_domain %in% c("E", "C")))
      }) %>%
      # calculate the z-stat
      map_depth(3, function(test.x) {
        if(nrow(test.x) < 2) {
          # not possible to make the comparison b/c mcmc sample only projected
          # one (or less) of the PECA-domains
          NA
        } else if (any(is.na(test.x$SEM))) {
          # not possible to make the comparison b/c mcmc sample's est for one
          # or both of the peca-domains was not possible to compute
          NA
        } else {
          with((test.x %>% 
                  select(peca_domain, est, SEM) %>% 
                  arrange(desc(peca_domain))),
               # handle the three possible cases:
               if (all(SEM == c(0, 0) & (est[1] == est[2]))) {
                 # ests are equal & no variability in plot-level est
                 0
               } else if (all(SEM == c(0, 0))) {
                 # ests are not equal & no variability in plot-level est
                 Inf
               } else {
                 # compare our observed z-stat to a two-tailed z for alpha =.1
                 # - i.e., is this difference significantly different from zero,
                 #   are these estimates significantly different
                 (est[1] - est[2]) / sqrt((SEM[1])^2 + (SEM[2])^2)
               })
        }
      }) %>%
      map_depth(2, ~.x %>% bind_rows(.id = "tested_domains")) %>%
      map(~.x %>% bind_rows(.id = "stat_id")) %>%
      bind_rows(.id = "mcmc_n") %>%
      inner_join(sp.x %>% distinct(stat_id, stat_descr), by = "stat_id")
  }) 

# - perform the test ----------------------------------------------------- ####
# Test that the estimates are significantly different at the 0.1-level.
# - One-tailed upper test against the null z-value of qnorm(0.90)), where
#   testing that mortality was higher or net-growth lower in {PvE, CvP, CvE}
# - note: in code below, adjusting direction for some 1-tailed comparisons
# - In the returned object, [TRUE = no sig diff] & [FALSE = sig dif]
#   Assigned this way to make sense when we calculate the 'percent agreement' 
#   for the MCMC samples' ests.

result_descr_chr <- c(
  "Mortality not significantly higher and net-growth not significantly lower",
  "Mortality higher, but net-growth not significantly lower",
  "Net-growth lower, but mortality not significantly higher",
  "Mortality higher & net-growth lower")

test_1tail <- dif_z_stat_cna %>% 
  # perform the test
  map_depth(2, function(sp.x, z_crit1 = qnorm(0.90)) {
    sp.x %>%
      pivot_longer(cols = contains("domain"),
                   names_to = "test",
                   values_to = "z_test_stat") %>%
      mutate(test = gsub("domain_", "", test)) %>%
      mutate(one_tail = case_when(
        # was mortality not significantly higher in:
        # - {P vs E}
        (test == "PvE" & stat_id %in% c("stat_2", "stat_4")) ~ 
          ifelse(is.infinite(z_test_stat), FALSE,  z_test_stat < z_crit1),
        # - {C vs P; C vs E}
        (test %in% 
           c("PvC", "EvC") & stat_id %in% c("stat_2", "stat_4")) ~ 
          ifelse(is.infinite(z_test_stat), FALSE,  (-1 * z_test_stat) < z_crit1),
        # was net-growth not significantly lower in:
        # - {P vs E}
        (test == "PvE" & stat_id %in% c("stat_6", "stat_8")) ~ 
          ifelse(is.infinite(z_test_stat), FALSE,  (-1 * z_test_stat) < z_crit1),
        # - {C vs P; C vs E; A vs C; A vs E}
        (test %in% 
           c("PvC", "EvC") & stat_id %in% c("stat_6", "stat_8")) ~ 
          ifelse(is.infinite(z_test_stat), FALSE,  z_test_stat < z_crit1))) %>%
      # move comparison names around so they make sense with the 1-tailed tests
      mutate(test = case_when(test == "PvC" ~ "C_vs_P",
                              test == "EvC" ~ "C_vs_E",
                              test == "PvE" ~ "P_vs_E"))
  }) %>%
  # summarize test results
  map_depth(2, function(sp.x) {
    sp.x %>%
      mutate(stat_units = case_when(grepl("BAH", stat_descr) ~ "BAH/yr",
                                    grepl("TPH", stat_descr) ~ "TPH/yr")) %>%
      filter(!is.na(one_tail)) %>%
      group_by(mcmc_n, test, stat_units) %>%
      summarize(n_comp = n(),
                result_code = paste(one_tail, collapse = "_"),
                result_code = case_when(result_code == "TRUE_TRUE"   ~ 1,
                                        result_code == "FALSE_TRUE"  ~ 2,
                                        result_code == "TRUE_FALSE"  ~ 3,
                                        result_code == "FALSE_FALSE" ~ 4),
                result_descr = case_when(result_code == "1" ~ result_descr_chr[1],
                                         result_code == "2" ~ result_descr_chr[2],
                                         result_code == "3" ~ result_descr_chr[3],
                                         result_code == "4" ~ result_descr_chr[4]),
                result_descr = factor(result_descr,
                                      levels = rev(result_descr_chr))) %>%
      ungroup %>%
      filter(n_comp > 1) %>%
      group_by(test, stat_units, result_code, result_descr) %>%
      tally(name = "n_mcmc") %>% 
      ungroup() %>% 
      group_by(test, stat_units) %>%
      mutate(tot_mcmc = sum(n_mcmc),
             perc_agree = n_mcmc / tot_mcmc) %>%
      ungroup %>%
      arrange(test)
  })

# - make figure summarizing test results --------------------------------- ####
# stacked figures for each PECA-comparison (C-P, C-E, P-E)
# - legend on bottom (2x2 legend)
# - facet_grid of spp by stat_units (BAH/yr, TPH/yr) where y-axis is the
#   percent of MCMC samples and x-axis is four columns for future climate 
#   scenarios 
# - But! first create a chr vector with how we want the PECA-comparisons to 
#   be labeled (peca_com_chr) and new strings for legend
peca_com_chr <- c("(A) Contraction vs. Persistence",
                  "(B) Contraction vs. Expansion",
                  "(C) Persistence vs. Expansion")
result_descr_chr2 <- c(
  "Mortality not significantly\n higher & net-growth not\nsignificantly lower",
  "Mortality higher,\nbut net-growth not\nsignificantly lower",
  "Net-growth lower,\nbut mortality not\nsignificantly higher",
  "Mortality higher &\nnet-growth lower\n")

(test_1tail %>% 
  imap(function(sp.x, sp_name) {
    sp.x %>% 
      map(function(df.x) {
        df.x %>%
          mutate(test_domains = case_when(test == "C_vs_P" ~ peca_com_chr[1],
                                          test == "C_vs_E" ~ peca_com_chr[2],
                                          test == "P_vs_E" ~ peca_com_chr[3]),
                 test_domains = factor(test_domains, 
                                       levels = peca_com_chr),
                 perc_agree = perc_agree * 100) %>%
          select(stat_units, test_domains, perc_agree, 
                 result_code, result_descr, tot_mcmc)
      }) %>%
      bind_rows(.id = "scenario") %>%
      mutate(spp = toupper(sp_name),
             spp_common = case_when(spp == "ABPR" ~ "noble fir",
                                    spp == "PSMEM" ~ "coastal Douglas-fir",
                                    spp == "QUDO" ~ "blue oak",
                                    spp == "QUGA4" ~ "white oak",
                                    spp == "QUKE" ~ "black oak"),
             spp_common = factor(spp_common,
                                 levels = c("noble fir",
                                            "coastal Douglas-fir",
                                            "blue oak",
                                            "white oak",
                                            "black oak")))
  }) %>%
  bind_rows %>%
  split(.$test_domains) %>%
  # make the plots
  imap(function(df.x, peca_com_chr.x) {
    (df.x %>%
       mutate(scenario_code = case_when(scenario == "RCP 4.5 CCSM4" ~ "f1",
                                        scenario == "RCP 4.5 HadGEM2-ES" ~ "f2",
                                        scenario == "RCP 8.5 CCSM4" ~ "f3",
                                        scenario == "RCP 8.5 HadGEM2-ES" ~ "f4"),
              scenario_code = factor(scenario_code,
                                     levels = paste0("f", 1:4)),
              result_descr2 = case_when(
                result_descr == result_descr_chr[1] ~ result_descr_chr2[1],
                result_descr == result_descr_chr[2] ~ result_descr_chr2[2],
                result_descr == result_descr_chr[3] ~ result_descr_chr2[3],
                result_descr == result_descr_chr[4] ~ result_descr_chr2[4]),
              result_descr2 = factor(result_descr2, 
                                     levels = rev(result_descr_chr2)))) %>%
      ggplot() +
      geom_bar(aes(x = scenario_code, y = perc_agree, fill = as.factor(result_descr2)),
               stat = 'identity',
               position = 'stack') +
      scale_fill_viridis_d(option = "D", 
                           direction = -1) +
      facet_grid(stat_units ~ spp_common) +
      theme_bw() +
      labs(subtitle = peca_com_chr.x,
           y = "percent",
           x = "scenario") +
      theme(panel.grid.major = element_blank(),
            plot.title.position = "plot",
            # aspect.ratio = 1,
            panel.grid.minor = element_blank(),
            plot.margin = margin(1,0,0,0),
            legend.title = element_blank(),
            legend.spacing.y = unit(4.0, 'mm'),
            plot.subtitle = element_text(size = (FF_font_size + 2)),
            strip.text.x = element_text(size = FF_font_size),
            strip.text.y = element_text(size = FF_font_size),
            axis.title.y = element_text(size = FF_font_size),
            axis.title.x = element_text(size = FF_font_size),
            axis.text.x  = element_text(size = FF_font_size),
            axis.text.y  = element_text(size = FF_font_size),
            legend.text  = element_text(size = FF_font_size)) +
      guides(fill = guide_legend(override.aes = list(size = 8),
                                 byrow = TRUE))
  }) %>%
  reduce(., `/`) +
  plot_layout(guides = "collect")  & 
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.box.margin = margin(0,0,0,0)) &
      guides(fill = guide_legend(nrow = 1, byrow = T))) %>%
  ggsave(filename = paste0(path_images, "plot_1tail-ztest_all-spp-scenarios.jpg"),
         plot = .,
         height = 200, 
         width = FF_fig_width_2col, 
         units = FF_width_unit, 
         dpi = FF_min_dpi)
