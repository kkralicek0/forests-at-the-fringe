###############################################################################
## Project - Forests at the fringe                                         ----
## Script  - 4-results-projected-impacts.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Forests at the Fringe: comparing observed
## change to projected climate change impacts for five species in the Pacific 
## Northwest' by Karin Kralicek, Tara Barrett, Jay Ver Hoef, and Temesgen
## Hailemariam.
## 
## About - this script:
## - Summarize and visualize projected climate change impacts
## - Produce figures/tables for manuscript and appendix. 
##   (note - 'observed change' results are analyzed in the next (5) script)
##
## About - output:
## - FF-4-hex-crosswalk.rds
##   - location: path_output
##   - crosswalk between plots & their containing hexagons for plotting
## - FF-4-map-base.rds
##   - location: path_output
##   - Subset of the pred objects for each scenario to just mean-preds & SE,
##     and, for the four future scenarios, the difference between future and
##     current mean-preds and SEs.
###############################################################################
# data manipulation
library(magrittr)   # for divide_by() etc. (and pipess)
library(dplyr)      # for ... E V E R Y T H I N G ...
library(tidyr)
library(purrr)      # o! beloved map...
library(sf)         # spatial data for plotting

# NOTE -- may require pkg `hexbin` for `stat_summary_hex()`
library(ggplot2)    # plots

library(patchwork)  # plots
library(parallel) 

# also references the lemon package for `lemon::reposition_legend()` 
# - note: legend panels can be identified with `lemon::gtable_show_names()`

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
# - load|create spatial data & hex-crosswalk ----------------------------- ####
# For plotting:
# - mean-preds & SE maps
# - percent-assignment of plots to PECA-domains
plot_polygons_public.sf <- readRDS(
  paste0(path_output_ch2, "ch2-8-plot_polygons_public_sf.rds"))
study_area.sf <- readRDS(
  paste0(path_output_ch2, "ch2-8-study_area_sf.rds"))

# load|create crosswalk between plots & containing hexagons' centroid coords
if (file.exists(paste0(path_output, "FF-4-hex-crosswalk.rds"))) {
  hex_plot_walk <- readRDS(paste0(path_output, "FF-4-hex-crosswalk.rds"))
} else {
  plot_points_true.sf  <- readRDS(
    paste0(path_output_ch2, "ch2-8-plot_points_true_sf.rds"))
  hex_polygons_true.sf <- readRDS(
    paste0(path_output_ch2, "ch2-8-hex_polygons_true_sf.rds"))
  
  hex_plot_walk <- list(plot_points_true.sf, hex_polygons_true.sf) %>%
    pmap(function(ppoints_sf.x, hex_sf.x) {
      hex_sf.x %>%
        st_set_agr(., "constant") %>% # specify attributes as spatially constant
        st_centroid() %>%
        mutate(x = unlist(map(.$geometry,1)),
               y = unlist(map(.$geometry,2)),
               fiahex_id = fiahex_id %>% as.character) %>%
        st_drop_geometry() %>%
        inner_join(ppoints_sf.x %>% st_drop_geometry,
                   by = "fiahex_id") %>%
        select(plot_id, fiahex_id, p_a) %>%
        rename(t1_pa_1in = p_a)
    })
  
  saveRDS(hex_plot_walk, paste0(path_output, "FF-4-hex-crosswalk.rds"))
}

# - load|create mean-preds & bio* from {PRISM, CNA, MACA} ---------------- ####
# Subset the prediction objects for each data set to just mean-preds & SE, as
# well as the difference between future and current mean-preds and SEs.
if(file.exists(paste0(path_output, "FF-4-map-base.rds"))) {
  map_base <- readRDS(paste0(path_output, "FF-4-map-base.rds"))
} else {
  map_base <- with(
    readRDS(paste0(path_output, "FF-1-preds-cNA_100MCMC.rds")),
    list(
      current = readRDS(paste0(path_output, "FF-1-preds-PRISM_100MCMC.rds")),
      CCSM4_rcp45 = CCSM4_rcp45,
      HadGEM2_ES_rcp45 = HadGEM2_ES_rcp45,
      CCSM4_rcp85 = CCSM4_rcp85,
      HadGEM2_ES_rcp85 = HadGEM2_ES_rcp85)) %>%
    map_depth(2, function(sp.x) {
      sp.x$maps_1k %>% # ideally this would say 'maps_100', but typo persists
        map(~.x %>% select(plot_id, mu_hat)) %>%
        bind_rows() %>%
        group_by(plot_id) %>%
        summarize(mu_meanpred = mean(mu_hat),
                  se_mu = sd(mu_hat))
    }) %>%
    transpose() %>%
    # add scenario names to dfs
    map(~.x %>% imap(function(df.x, scenario_name) {
      df.x %>% 
        mutate(scenario = scenario_name,
               stat_type = "orig")
    })) %>%
    # calc diff between future and current preds
    map(function(sp.x) {
      with(sp.x[2:5] %>% map(function(scenario.x) {
        rbind.data.frame(
          scenario.x,
          scenario.x %>%
            inner_join(sp.x$current %>% select(-stat_type, -scenario), 
                       by = "plot_id") %>%
            mutate(mu_meanpred = mu_meanpred.x - mu_meanpred.y,
                   se_mu = se_mu.x - se_mu.y) %>%
            select(plot_id, mu_meanpred, se_mu, scenario) %>%
            mutate(stat_type = "diff"))
      }), 
      list(current = sp.x$current,
           CCSM4_rcp45 = CCSM4_rcp45,
           HadGEM2_ES_rcp45 = HadGEM2_ES_rcp45,
           CCSM4_rcp85 = CCSM4_rcp85,
           HadGEM2_ES_rcp85 = HadGEM2_ES_rcp85))
    })
  
  saveRDS(map_base, paste0(path_output, "FF-4-map-base.rds"))
}

# - load PECA-division groups -------------------------------------------- ####
# PECA
# - area projected to Persist, Expand, Contract, or remain unsuitable (A)
#   and classification of plots within each MCMC sample based on CNA preds
peca_cna <- readRDS(paste0(path_output, "FF-1-groups_peca_cna.rds"))

# - load quick elevation/latitude key ------------------------------------ ####
elev_lat_key <- readRDS(paste0(path_data_ch2, "ch2-1-elev-lat-plot-key.rds"))

# - load elev-lat groups ------------------------------------------------- ####
# Load these for some specific summaries about prevalence within naive-div
elevlat <- readRDS(paste0(path_data, "FF-3-elev-lat-groups.rds"))

# == Color palettes & other plotting args ================================ ####
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

# == Map: naive-divisions & mode of peca-divisions, by spp =============== ####
# make naive-divisions part of the figure
p_body_naive <- elevlat %>%
  # subset elevlat to remove plots without an assigned plot_id (this just 
  # for plotting - using Voronoi polygons dev in ch2 for vizualization)
  map(function(sp.x) {
    sp.x %>%
      select(plot_id, el_group) %>%
      filter(!is.na(plot_id)) %>%
      mutate(el_group = factor(el_group,
                               levels = c("LL", "LM", "LH",
                                          "ML", "MM", "MH",
                                          "HL", "HM", "HH")),
             facet_label = "Naïve-divisions\n")
  }) %>% 
  # join with (public) spatial coords
  list(., plot_polygons_public.sf) %>%
  pmap(function(df.x, sf.x) {
    sf.x %>%
      select(-NAME, -p_a) %>%
      inner_join(df.x, by = "plot_id")
  }) %>%
  # make the base plots
  map(function(sp.x) {
    ggplot(sp.x, aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, 
              fill = "grey50") +
      geom_sf(aes(fill = el_group), color = NA) +
      facet_grid(. ~ facet_label) +
      scale_fill_viridis_d("Naïve-divisions",
                           labels = c("LL", "LM", "LH",
                                      "ML", "MM", "MH",
                                      "HL", "HM", "HH")) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() + 
      theme(plot.margin = margin(1, 1, 0, 1),
            strip.text.x = element_text(size = FF_font_size),
            axis.text.x  = element_text(size = FF_font_size),
            axis.text.y  = element_text(size = FF_font_size),
            legend.title = element_text(size = FF_font_size),
            legend.text  = element_text(size = FF_font_size)) +
      guides(fill = guide_legend(nrow = 3, byrow = F, title.vjust = 1))
  }) 

# make peca-divisions part of the figure
p_body_peca <- peca_cna %>% 
  transpose %>% 
  setNames(names(elevlat)) %>%
  # get mode of PECA-assignments for each plot
  map_depth(2, function(df.x) {
    df <- df.x$peca_plots %>%
      unlist %>% 
      matrix(ncol = 100, byrow = FALSE) %>% 
      apply(1, function(x) {
        table(factor(x, levels = c("P", "E", "C", "A")))
      }) %>%
      t %>%
      as.data.frame
    df$plot_id <- df.x$plot_id
    df %>% 
      pivot_longer(cols = -plot_id,
                   names_to = "peca_domain",
                   values_to = "n_mcmc") %>%
      group_by(plot_id) %>%
      filter(n_mcmc == max(n_mcmc)) %>%
      mutate(peca_code = case_when(peca_domain == "E" ~ 1,
                                   peca_domain == "P" ~ 2,
                                   peca_domain == "C" ~ 3,
                                   peca_domain == "A" ~ 4))
  }) %>%
  # pretty scenario names and to one df for each spp
  # - also put 'PECA-divisions' on the prev line so we can easily differntiate
  imap(function(sp.x, sp_name) {
    sp.x %>% 
      bind_rows(.id = "scenario") %>%
      # pretty labels for plotting
      mutate(scenario = case_when(
        scenario == "CCSM4_rcp45" ~ "RCP 4.5 CCSM4",
        scenario == "HadGEM2_ES_rcp45" ~ "RCP 4.5 HadGEM2-ES",
        scenario == "CCSM4_rcp85" ~ "RCP 8.5 CCSM4",
        scenario == "HadGEM2_ES_rcp85" ~ "RCP 8.5 HadGEM2-ES")) %>%
      mutate(scenario = factor(scenario, 
                               levels = c("RCP 4.5 CCSM4",
                                          "RCP 4.5 HadGEM2-ES",
                                          "RCP 8.5 CCSM4",
                                          "RCP 8.5 HadGEM2-ES")),
             spp = toupper(sp_name)) %>%
      mutate(scenario2 = paste0("PECA-divisions:\n", scenario),
             scenario2 = factor(scenario2, 
                                levels = paste0("PECA-divisions:\n",
                                                c("RCP 4.5 CCSM4",
                                                  "RCP 4.5 HadGEM2-ES",
                                                  "RCP 8.5 CCSM4",
                                                  "RCP 8.5 HadGEM2-ES"))))
  }) %>%
  # join with (public spatial coords)
  list(., plot_polygons_public.sf) %>%
  pmap(function(df.x, sf.x) {
    sf.x %>%
      select(-NAME, -p_a) %>%
      inner_join(df.x, by = "plot_id")
  }) %>%
  # make the base plots
  map(function(sp.x) {
    ggplot(sp.x, aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, 
              fill = "grey50") +
      geom_sf(aes(fill = as.factor(peca_code)), color = NA) +
      facet_grid(. ~ scenario2) +
      scale_fill_manual(
        name = "PECA-divisions",
        values = unname(peca_domain_palette),
        breaks = c("1", "2", "3", "4"),
        labels = names(peca_domain_palette)) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() + 
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(1, 1, 0, 0),
            strip.text.x = element_text(size = FF_font_size),
            axis.text.x  = element_text(size = FF_font_size),
            legend.title = element_text(size = FF_font_size),
            legend.text  = element_text(size = FF_font_size)) +
      guides(fill = guide_legend(title.vjust=1))
  })

# make and save the plots
p_list <- list(p_body_naive, p_body_peca) %>%
  transpose %>%
  imap(function(sp.x, sp_name) {
    sp.x %>% 
      reduce(., `|`) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom",
            legend.justification = "left")
  })

# adjusted height until it looked right...
p_list %>% imap(~ ggsave(paste0(path_images, "map_divisions_", .y, ".jpg"),
                         plot = .x,
                         height = 110, 
                         width = FF_fig_width_2col, 
                         units = FF_width_unit, 
                         dpi = FF_min_dpi))



# == Predictions of climatic suitable habitat ============================ ####
# - plot: mean-preds & SE by elev x lat ---------------------------------- ####
# Description: 5 (spp, rows) x 5 (scenarios, cols) figure showing either
# the mean-preds (prob for current, diff prob for future scenarios) or SE
# (similarly, SE or diff SE for current or future scenarios, respectively)
map_base %>%
  # format data for plotting
  imap(function(sp.x, sp_name) {
    sp.x %>% 
      bind_rows() %>%
      # grab current (mp & se) & future (dif-mp & dif-se)
      mutate(dummy = paste0(scenario, "_", stat_type)) %>%
      filter(dummy == "current_orig" | stat_type == "diff") %>%
      # pretty labels for plotting
      mutate(scenario = case_when(
        scenario == "current" ~ "current",
        scenario == "CCSM4_rcp45" ~ "RCP 4.5 CCSM4",
        scenario == "HadGEM2_ES_rcp45" ~ "RCP 4.5 HadGEM2-ES",
        scenario == "CCSM4_rcp85" ~ "RCP 8.5 CCSM4",
        scenario == "HadGEM2_ES_rcp85" ~ "RCP 8.5 HadGEM2-ES")) %>%
      mutate(scenario = factor(scenario, 
                               levels = c("current",
                                          "RCP 4.5 CCSM4",
                                          "RCP 4.5 HadGEM2-ES",
                                          "RCP 8.5 CCSM4",
                                          "RCP 8.5 HadGEM2-ES"))) %>%
      pivot_longer(cols = c("mu_meanpred", "se_mu"),
                   names_to = "stat") %>%
      left_join(elev_lat_key %>% select(-LON_ACTUAL), by = "plot_id") %>%
      mutate(spp = toupper(sp_name),
             elev_km = Elev_m / 1000) %>%
      select(-Elev_m) %>%
      split(.$stat)
  }) %>%
  transpose %>%
  map(~.x %>% bind_rows) %>%
  # make the figure
  imap(function(df.x, stat) {
    # set color based on stat --------------------------------------- ####
    TF_mp <- stat == "mu_meanpred"
    
    # change color-scheme based on if diff is mean-preds or SE...
    c_cb_name   <- ifelse(TF_mp, "prob.", "SE")
    c_cb_option <- ifelse(TF_mp, "viridis", "magma")
    f_cb_name   <- ifelse(TF_mp, "diff prob.", "diff SE")
    f_cb_low    <- ifelse(TF_mp, "#39568cff", "#451077FF")
    f_cb_high   <- ifelse(TF_mp, "#FDD835", "#FDD835")
    title_name  <- ifelse(TF_mp, "mean-prediction", "SE of prediction")
    
    # current plot -------------------------------------------------- ####
    p_c <- ggplot(df.x %>% filter(scenario == "current")) +
      stat_summary_hex(aes(y = LAT_ACTUAL, x = elev_km, z = value),
                       fun = ~ mean(.x),
                       bins = 40) +
      facet_grid(spp ~ scenario, scales = "free") +
      scale_fill_viridis_c(name = c_cb_name,
                           option = c_cb_option) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "grey50"),
            strip.background.y = element_blank(),
            strip.text.y = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(1,0,1,1)) +
      labs(y = "Latitude (°N)")
    
    # future plots -------------------------------------------------- ####
    # plot points on grey background for visibility
    p_f <- ggplot(df.x %>% filter(scenario != "current")) +
      stat_summary_hex(aes(y = LAT_ACTUAL, x = elev_km, z = value),
                       fun = ~ mean(.x),
                       bins = 40) +
      facet_grid(spp ~ scenario, scales = "free") +
      scale_fill_gradientn(
        name = f_cb_name,
        colors = c(f_cb_low, "#e5e5e5", f_cb_high),
        values = scales::rescale(c(-1, -0.25, 0, 0.25, 1)),
        limits = function(x) {c(-1,1)*max(abs(x))}) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "grey50"),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(1,1,1,0),
            axis.title.x = element_text(hjust = 0.35)) + 
      labs(x = "Elevation (km)")
    
    # combine into one figure & return ------------------------------ ####
    p_title <- paste0("Pairwise difference between climate data sets' ",
                      title_name,
                      "/nmean displayed for plots within binned area")
    
    (p_c | p_f) + 
      plot_layout(widths = c(1.05, 4),
                  heights = c(1,1),
                  guides = "collect") & 
      theme(legend.position = "bottom")
    
  }) %>%
  # print out only mp's atm...
  pluck(1)

# == Range-size & PECA size changes ====================================== ####
# - table & boxplot: change & rel-change in range-size-------------------- ####
# create base object
rangesize_posterior_summaries <- peca_cna  %>%
  map(~.x %>% map(~.x$df_area) %>% bind_rows(.id = "spp")) %>%
  bind_rows(.id = "scenario") %>%
  mutate(spp = toupper(spp),
         # convert to thousand km^2
         kkm_change = (future_ha - current_ha)/100/1000,
         rel_change = (future_ha - current_ha)/current_ha *100) %>%
  pivot_longer(cols = contains("_change"),
               names_to = "stat_type") %>%
  mutate(stat_type = case_when(
    stat_type == "kkm_change" ~ "change in range-size (thousand km^2)",
    TRUE ~ "relative-change in range-size (percent)")) %>%
  group_by(spp, scenario, stat_type) %>%
  summarise(y_min = min(value),
            q_05 = quantile(value, .05) %>% as.numeric,
            y_mean = mean(value),
            q_95 = quantile(value, .95) %>% as.numeric,
            y_max = max(value)) %>%
  mutate(spp = factor(spp, levels = rev(c("ABPR",
                                          "PSMEM",
                                          "QUDO", 
                                          "QUGA4",
                                          "QUKE")))) %>%
  # pretty labels for plotting
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

# table view of posterior means & 90% credible intervals for change and 
# rel-change in range-size
rangesize_posterior_summaries %>% 
  mutate(chr = paste0(round(y_mean, 1), " (", 
                      round(q_05, 1), ", ", 
                      round(q_95, 1), ")")) %>% 
  select(spp, scenario, stat_type, chr) %>% 
  pivot_wider(names_from = stat_type, 
              values_from = chr) %>% 
  arrange(desc(spp), scenario)

# figure with boxplot summaries of posterior distributions
rangesize_posterior_summaries %>%
  split(.$stat_type) %>%
  map(function(df.x) {
    df.x %>%
      ggplot() +
      geom_hline(yintercept = 0, linetype = 4) +
      geom_boxplot(aes(x = spp,
                       ymin = y_min,
                       lower = q_05,
                       middle = y_mean,
                       upper = q_95,
                       ymax = y_max,
                       fill = scenario),
                   stat = "identity",
                   width = 0.75,
                   alpha = 0.7) +
      # scale_y_continuous(n.breaks = 6, position = "right") +
      scale_x_discrete(limits = c("ABPR", "PSMEM", "QUDO", "QUGA4", "QUKE")) +
      scale_fill_manual(values = clim_ds_palette[2:5] %>% 
                          setNames(c("RCP 4.5 CCSM4", "RCP 4.5 HadGEM2-ES",
                                     "RCP 8.5 CCSM4", "RCP 8.5 HadGEM2-ES"))) +
      facet_wrap(stat_type ~ .) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            plot.margin = margin(0,1,0,1))
  }) %>%
  # adjust depending on plot
  map_at(1, function(plot.x) {
    plot.x + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }) %>%
  # make the figure
  reduce(., `/`) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Climate impacts to range-size by scenario & PECA-domain",
                  subtitle = "Estimates and 90% CI shown for each of 100 MCMC samples")

# - map: freq of plot assignment to PECA-domains ------------------------- ####
# map the freq with which individual plots are classified as P, E, C, or A...
# - individual figures for each spp
peca_cna %>%
  # get PECA frequencies
  map_depth(2, function(df.x) {
    df <- df.x$peca_plots %>%
      unlist %>% 
      matrix(ncol = 100, byrow = FALSE) %>% 
      apply(1, function(x) {
        table(factor(x, levels = c("P", "E", "C", "A")))
      }) %>%
      t %>%
      divide_by(100) %>%
      as.data.frame
    df$plot_id <- df.x$plot_id
    df %>% 
      pivot_longer(cols = -plot_id,
                   names_to = "peca_domain",
                   values_to = "freq")
  }) %>%
  # list by spp
  transpose %>%
  # pretty scenario names and to one df for each spp
  map(function(sp.x) {
    sp.x %>%
      bind_rows(.id = "scenario") %>%
      # pretty labels for plotting
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
  }) %>%
  # drop 'A' - only show PEC (immediate interest)
  map(~.x %>% filter(peca_domain != "A")) %>%
  # join spatial info and plot/map
  list(plot_polygons_public.sf, .) %>%
  pmap(function(sf.x, peca.x) {
    sf.x %>%
      select(-NAME, -p_a) %>%
      inner_join(peca.x, by = "plot_id")
  }) %>%
  imap(function(sp.x, sp_name) {
    ggplot(sp.x, aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, 
              fill = "grey50") +
      geom_sf(aes(fill = freq), color = NA) +
      facet_grid(peca_domain ~ scenario) +
      scale_fill_viridis_b(option = "C",
                           breaks = seq(0, 1, .2)) +
      # scale_fill_gradient(low = "#e5e5e5", high = "#FDD835") +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() +
      guides(fill = guide_colorbar(title.vjust = 0.8)) + # b/c misaligned title o/w (bug)
      ggtitle(paste0(toupper(sp_name), 
                     ": frequency of group membership by future scenario"),
              subtitle = "Based on 100 MCMC samples (every 20th) and BKDE range estimates")
  })

# - plot: most common plot assign by elev x lat -------------------------- ####
tmp <- peca_cna %>%
  # get PECA frequencies
  map_depth(2, function(df.x) {
    df <- df.x$peca_plots %>%
      unlist %>% 
      matrix(ncol = 100, byrow = FALSE) %>% 
      apply(1, function(x) {
        table(factor(x, levels = c("P", "E", "C", "A")))
      }) %>%
      t %>%
      as.data.frame
    df$plot_id <- df.x$plot_id
    df %>% 
      pivot_longer(cols = -plot_id,
                   names_to = "peca_domain",
                   values_to = "n_mcmc") %>%
      filter(n_mcmc > 0) %>%
      mutate(peca_code = case_when(peca_domain == "E" ~ 1,
                                   peca_domain == "P" ~ 2,
                                   peca_domain == "C" ~ 3,
                                   peca_domain == "A" ~ 4)) %>%
      uncount(n_mcmc)
  }) %>%
  # list by spp
  transpose %>%
  # pretty scenario names and to one df for each spp
  imap(function(sp.x, sp_name) {
    sp.x %>%
      bind_rows(.id = "scenario") %>%
      # pretty labels for plotting
      mutate(scenario = case_when(
        scenario == "CCSM4_rcp45" ~ "RCP 4.5 CCSM4",
        scenario == "HadGEM2_ES_rcp45" ~ "RCP 4.5 HadGEM2-ES",
        scenario == "CCSM4_rcp85" ~ "RCP 8.5 CCSM4",
        scenario == "HadGEM2_ES_rcp85" ~ "RCP 8.5 HadGEM2-ES")) %>%
      mutate(scenario = factor(scenario, 
                               levels = c("RCP 4.5 CCSM4",
                                          "RCP 4.5 HadGEM2-ES",
                                          "RCP 8.5 CCSM4",
                                          "RCP 8.5 HadGEM2-ES"))) %>%
      left_join(elev_lat_key %>% select(-LON_ACTUAL), by = "plot_id") %>%
      mutate(spp = toupper(sp_name),
             elev_km = Elev_m / 1000) %>%
      select(-Elev_m) 
  }) %>%
  bind_rows()

ggplot(tmp) +
    stat_summary_hex(aes(y = LAT_ACTUAL, x = elev_km, z = peca_code),
                     # this function grabs the most common peca-division in
                     # a given hex (i.e., both across the plots in that hex  
                     # and their assignment under the MCMC samples...)
                     fun = function(x) {names(table(x))[which.max(table(x))]}) +
    facet_grid(spp ~ scenario, scales = "free") +
    scale_fill_manual(
      name = "PECA-division",
      values = unname(peca_domain_palette),
      breaks = c("1", "2", "3", "4"),
      labels = names(peca_domain_palette)) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey50"),
          axis.title.x = element_text(hjust = 0.35)) + 
    labs(x = "Elevation (km)",
         y = "Latitude (°N)",
         title = "Most common PECA-division across mcmc-samples")


# - scatter: Expansion v. Contraction relative to current range-size ----- ####
p_body <- (peca_cna %>%
             map_depth(2, ~.x$df_area) %>%
             transpose() %>%
             imap(function(sp.x, sp_name) {
               sp.x %>%
                 bind_rows(.id = "scenario") %>%
                 mutate(E_rel = E_ha / current_ha * 100,
                        C_rel = C_ha / current_ha * 100,
                        f_rel = future_ha / current_ha * 100,
                        spp = toupper(sp_name)) %>%
                 # select(spp, scenario, mcmc_n, E_rel, C_rel) %>%
                 mutate(line_1to1 = min(max(E_rel), max(C_rel))) %>%
                 # pretty scenario labels for plotting
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
             }) %>%
             bind_rows) %>%
  ggplot() +
  geom_segment(aes(xend = line_1to1, yend = line_1to1, x = 0, y = 0)) +
  geom_point(aes(x = C_rel, y = E_rel, color = scenario),
             alpha = .5) + 
  scale_color_manual(values = clim_ds_palette[2:5] %>% 
                       setNames(c("RCP 4.5 CCSM4", "RCP 4.5 HadGEM2-ES",
                                  "RCP 8.5 CCSM4", "RCP 8.5 HadGEM2-ES"))) +
  facet_wrap(. ~ spp, scales = "free") +
  expand_limits(x = 0, y = 0) +
  labs(y = "relative-expansion (percent)",
       x = "relative-contraction (percent)",
       title = "Area of expansion and contraction relative to current range-size estimate") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1, 
                                                   size = 4, 
                                                   shape = 15)))

lemon::reposition_legend(p_body,
                         panel = "panel-3-2",
                         x = 0, y = 0, just = c(-0.1, -0.1))

p_body + 
  guide_area() + 
  plot_layout(guides = 'collect')


# - post-dist & EvC: range size ------------------------------------------ ####
p_body <- peca_cna %>%
  map_depth(2, ~.x$df_area) %>%
  transpose() %>%
  imap(function(sp.x, sp_name.x) {
    # - create base data set ---------------------------------------- ####
    legend_labels <- clim_ds_palette %>%
      setNames(c("current",
                 "RCP 4.5 CCSM4",
                 "RCP 4.5 HadGEM2-ES",
                 "RCP 8.5 CCSM4",
                 "RCP 8.5 HadGEM2-ES"))
    
    df.x <- sp.x %>%
      # pretty labels for plotting
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
                                          "RCP 8.5 HadGEM2-ES"))) %>%
      # grab spp names
      mutate(spp = toupper(sp_name.x)) %>%
      # convert from ha to thousand km^2...
      mutate(across(where(is.numeric), ~.x/100/1000)) %>%
      arrange(mcmc_n) %>%
      rename(current_est = current_ha,
             future_est = future_ha,
             persistence = P_ha,
             contraction = C_ha,
             expansion = E_ha)
    
    # - make the base figures --------------------------------------- ####
    # get end points for 1-1 line segment
    line_1to1 <- min(max(df.x$contraction), max(df.x$expansion))
    std_margin <- margin(0, 0, 5.5, 0)
    abpr_margin <- margin(0, 0, 7.5, 0)
    
    # plot it!
    # - do a little differently for abpr (many near zero future ranges est)
    if (sp_name.x == "abpr") {
      p_range_est <- (df.x %>% 
                        select(current_est, mcmc_n) %>%
                        mutate(scenario = "current") %>%
                        rename(range_est = current_est) %>%
                        select(scenario, range_est, mcmc_n) %>%
                        rbind.data.frame(df.x %>% 
                                           rename(range_est = future_est) %>%
                                           select(scenario, range_est, mcmc_n)) %>%
                        ungroup %>%
                        mutate(spp = toupper(sp_name.x))) %>%
        ggplot() +
        geom_density(aes(x = log1p(range_est), y = log1p(..density..), fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ .) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_fill_manual(values = legend_labels)  +
        xlab("log1p(range-size)") +
        ylab("log1p(density)") +
        theme_bw() +
        theme(plot.margin = abpr_margin,
              axis.title.x = element_text(vjust = 3))
      
      p_exp_v_contract <- df.x %>%
        ggplot() +
        geom_segment(aes(xend = log1p(line_1to1), yend = log1p(line_1to1)), x = 0, y = 0) +
        geom_point(aes(x = log1p(contraction), y = log1p(expansion), color = scenario),
                   alpha = .5) + 
        expand_limits(x = 0, y = 0) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_y_continuous(position = "right") +
        ylab("log1p(expan.)") +
        theme_bw() +
        theme(plot.margin = abpr_margin,
              axis.title.x = element_text(vjust = 3)) 
      
    } else {
      p_range_est <- (df.x %>% 
                        select(current_est, mcmc_n) %>%
                        mutate(scenario = "current") %>%
                        rename(range_est = current_est) %>%
                        select(scenario, range_est, mcmc_n) %>%
                        rbind.data.frame(df.x %>% 
                                           rename(range_est = future_est) %>%
                                           select(scenario, range_est, mcmc_n)) %>%
                        ungroup %>%
                        mutate(spp = toupper(sp_name.x))) %>%
        ggplot() +
        geom_density(aes(x = range_est, fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ .) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_fill_manual(values = legend_labels)  +
        xlab("range-size") +
        theme_bw() +
        theme(plot.margin = std_margin)
      
      p_exp_v_contract <- df.x %>%
        ggplot() +
        geom_segment(aes(xend = line_1to1, yend = line_1to1), x = 0, y = 0) +
        geom_point(aes(x = contraction, y = expansion, color = scenario),
                   alpha = .5)  + 
        expand_limits(x = 0, y = 0) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_y_continuous(position = "right") +
        theme_bw() +
        theme(plot.margin = std_margin)
    }
    
    # - return figs & base legend for the last ----------------------------- ####
    # do differently for quke
    # - quke will be the bottom figure, so the only one that needs x title 
    # - we'll also harvest the legend from this species figure
    if (sp_name.x == "quke") {
      # add some space to the left of the legend
      p_range_est <- p_range_est + 
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(legend.position = "bottom",
              legend.box.margin = margin(3, 0, 0, 0),
              legend.title = element_blank())
      
      list(p_range_est = p_range_est, 
           p_exp_v_contract = p_exp_v_contract + 
             guides(color = "none"))
    } else if (sp_name.x == "abpr") {
      list(p_range_est = p_range_est  + 
             guides(color = "none", fill = "none"), 
           p_exp_v_contract = p_exp_v_contract  + 
             guides(color = "none"))
    } else {
      list(p_range_est = p_range_est + 
             theme(axis.title.x = element_blank()) + 
             guides(color = "none", fill = "none"), 
           p_exp_v_contract = p_exp_v_contract + 
             theme(axis.title.x = element_blank()) + 
             guides(color = "none"))
    }
  }) %>%
  map(~ .x$p_range_est + .x$p_exp_v_contract) %>%
  reduce(., `/`)

p_body + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Bivariate KDE range-size estimates",
                  subtitle = "Based on 100 MCMC samples; thousand km^2",
                  caption = "Note: abpr visualized on logp1 scale (ln(x + 1))") & 
  theme(legend.position = "bottom")

# - table: area & rel-area by PECA-domain and scenario ------------------- ####
# table view of posterior means & 90% credible intervals for area and 
# percent of area relative to current range-size estimate by PECA-domain
peca_cna %>%
  map_depth(2, ~.x$df_area) %>%
  map(~.x %>% bind_rows(.id = "spp")) %>%
  bind_rows(.id = "scenario") %>%
  pivot_longer(cols = c(P_ha, E_ha, C_ha),
               names_to = "peca_domain",
               values_to = "peca_ha") %>%
  select(-future_ha) %>%
  mutate(spp = toupper(spp),
         # convert to thousand km^2
         peca_kkm2 = peca_ha/100/1000,
         rel_area  = (peca_ha/current_ha) * 100) %>%
  select(spp, mcmc_n, peca_domain, peca_kkm2, rel_area, scenario) %>%
  pivot_longer(cols = c(peca_kkm2, rel_area),
               names_to = "stat_type") %>%
  mutate(stat_type = case_when(
    stat_type == "peca_kkm2" ~ "area (thousand km^2)",
    TRUE ~ "rel-area (percent of current estimate)")) %>%
  group_by(spp, scenario, peca_domain, stat_type) %>%
  summarise(y_min = min(value),
            q_05 = quantile(value, .05) %>% as.numeric,
            y_mean = mean(value),
            q_95 = quantile(value, .95) %>% as.numeric,
            y_max = max(value)) %>% 
  mutate(chr = paste0(round(y_mean, 1), " (", 
                      round(q_05, 1), ", ", 
                      round(q_95, 1), ")")) %>% 
  select(spp, peca_domain, scenario, stat_type, chr) %>% 
  pivot_wider(names_from = stat_type, 
              values_from = chr) %>% 
  mutate(peca_domain = gsub("_ha", "", peca_domain),
         peca_domain = factor(peca_domain,
                              levels = c("E", "P", "C"))) %>%
  arrange(spp, scenario, peca_domain) %>%
  ungroup %>%
  select(-contains("thous")) %>%
  pivot_wider(names_from = peca_domain,
              values_from = contains("rel")) %>%
  select(spp, scenario, P, E, C) %>%
  # pretty labels for plotting
  mutate(scenario = case_when(
    scenario == "CCSM4_rcp45" ~ "RCP 4.5 CCSM4",
    scenario == "HadGEM2_ES_rcp45" ~ "RCP 4.5 HadGEM2-ES",
    scenario == "CCSM4_rcp85" ~ "RCP 8.5 CCSM4",
    scenario == "HadGEM2_ES_rcp85" ~ "RCP 8.5 HadGEM2-ES")) %>%
  mutate(scenario = factor(scenario, 
                           levels = c("RCP 4.5 CCSM4",
                                      "RCP 4.5 HadGEM2-ES",
                                      "RCP 8.5 CCSM4",
                                      "RCP 8.5 HadGEM2-ES"))) %>%
  arrange(spp, scenario)

# was there any expansion / persistence for noble fir (for which posterior
# summaries from the code above round to zero)
peca_cna %>%
  map_depth(2, ~.x$df_area) %>%
  map(~.x %>% bind_rows(.id = "spp")) %>%
  bind_rows(.id = "scenario") %>%
  pivot_longer(cols = c(P_ha, E_ha, C_ha),
               names_to = "peca_domain",
               values_to = "peca_ha") %>%
  select(-future_ha) %>%
  mutate(spp = toupper(spp),
         # keep as ha instead of converting to thousand km^2
         rel_area  = (peca_ha/current_ha) * 100) %>%
  select(spp, mcmc_n, peca_domain, peca_ha, rel_area, scenario) %>%
  pivot_longer(cols = c(peca_ha, rel_area),
               names_to = "stat_type") %>%
  group_by(spp, scenario, peca_domain, stat_type) %>%
  summarise(q_05 = quantile(value, .05) %>% as.numeric,
            y_mean = mean(value),
            q_95 = quantile(value, .95) %>% as.numeric) %>% 
  ungroup %>%
  filter(spp == "ABPR") %>% 
  mutate(peca_domain = gsub("_ha", "", peca_domain),
         peca_domain = factor(peca_domain,
                              levels = c("E", "P", "C"))) %>%
  arrange(spp, scenario, stat_type, peca_domain) %>%
  ungroup

# did expansion exceed contraction for black oak (I think yes, but double
# check posterior difference here...)
peca_cna %>%
  map_depth(2, ~.x$df_area) %>%
  map(~.x %>% bind_rows(.id = "spp")) %>%
  bind_rows(.id = "scenario") %>%
  mutate(value = (E_ha - C_ha)/100/1000) %>%
  group_by(spp, scenario) %>%
  summarise(q_05 = quantile(value, .05) %>% as.numeric,
            y_mean = mean(value),
            q_95 = quantile(value, .95) %>% as.numeric) %>% 
  ungroup %>%
  filter(spp == "quke")

