###############################################################################
## Project - Forests at the fringe                                         ----
## Script  - 3-get-domain-est-for-obsv-change.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Forests at the Fringe: comparing observed
## change to projected climate change impacts for five species in the Pacific 
## Northwest' by Karin Kralicek, Tara Barrett, Jay Ver Hoef, and Temesgen
## Hailemariam.
## 
## About - this script:
## - Get estimates and SE for net-growth & mortality metrics, both overall and 
##   for the PECA groupings using domain estimation and post-stratified 
##   sampling methods.
##
## About - output:
## - FF-3-elev-lat-groups.rds
##  - location: path_data
##  - a list by species with plot (rom_x) identification of naive-division
##    groups. These are naive divisions of habitat based on one-third equal 
##    divisions of elevation and latitude (9 total). Saving here so we don't
##    have to recreate wheel each time...
## - FF-3-domain-estimates.Rdata
##  - location: path_output
##  - contains 2 objects atm: est_overall, est_peca_cna, est_elevlat
##    each object contains the estimate and SE of the estimate of mortality
##    and net-growth (both reported in annualized BAH and TPH metrics).
##    (Also currently gives values for t1 & t2 live TPH/BAH, though need to 
##    check if it is reasonable to estimate these quantities with strata wts
##    that are intended for change components.
###############################################################################
library(readxl)   # grab actual elev/lat from prev data
library(dplyr)    # for ... E V E R Y T H I N G ...
library(tidyr)
library(purrr)    # o! beloved map...

# == Paths =============================================================== ####
# path to Rdata
path_top <- paste0("C:/Users/karinkralicek/Documents/", 
                   "projects/active/SI - forest at the fringe/code/")

# paths to intermediate data
path_output <- paste0(path_top, "x - output/")
path_data  <- paste0(path_top, "x - data/")
path_output_ch2 <- paste0("C:/Users/karinkralicek/Documents/CUI/",
                          "ch2 output sample/") # for spatial data..

# == Load data =========================================================== ####
# - load PECA groups ----------------------------------------------------- ####
# PECA
# - area projected to Persist, Expand, Contract, or remain unsuitable (A)
#   and classification of plots within each MCMC sample based on preds for 
#   the four future climate scenarios
peca_cna <- readRDS(paste0(path_output, "FF-1-groups_peca_cna.rds"))

# - load observed data --------------------------------------------------- ####
# list by spp with plot-level info as well as some info used to calculate 
# wts for the ratio-of-means (rom) estimator
# - rom numerator (y data): tph/yr & bah/yr for t1/t2/net-growth/mortality
rom_y <- readRDS(paste0(path_data, "FF-2-rom-y-data.rds"))

# - rom denominator (x data): adj-prop of plot that was forested (prop_for)
#   forested with spp present at t1 (prop_for_t1; used for mortality estimates)
#   or forested with spp present at both t1 and t2 (prop_for_t1t2; used for 
#   net-growth estimates)
rom_x <- readRDS(paste0(path_data, "FF-2-rom-x-data.rds"))

# list by spp with pop/evalid/estn_unit level info used to calculate wts and
# for some variance eqs, etc.
rom_wh <- readRDS(paste0(path_data, "FF-2-rom-wh-data.rds"))

# list of plots within study area and the states they occur within
state_mem <- readRDS(paste0(path_data, "FF-2-state-membership.rds"))

# - create|load elev & latitude groups ----------------------------------- ####
# For each species, define 9 naive-divisions of habitat based on the cross-
# combination of 3 (one-third) elevation-partitions and 3 latitude-partitions
# - We defined groups based on equal interval widths, considering only plots on
#   which the species was present. We did this rather than making the groups
#   have the same number of plots, which would not be as valuable to the spp as
#   it would only speak to our sampling design (a design that is not entirely
#   spatially balanced across the entire study area, e.g., CA temporal not 
#   spatial intensification for FIA).
# - Note on latitude divisions:
#   Decided to group North-to-South by latitude. Although we could divide based
#   on UTMs (to preserve distance), this would likely only group differently by
#   ~1 mile, so shouldn't be to radically different when estimating change in 
#   mortality
if (file.exists(paste0(path_data, "FF-3-elev-lat-groups.rds"))) {
  elev_lat <- readRDS(paste0(path_data, "FF-3-elev-lat-groups.rds"))
} else {
  # - Load actual elev/lat data ------------------------------------- ####
  # this data was provided by Tara Barrett during my PhD program and 
  # contains actual (true) elevation (m) and latitude (nad83) data for
  # the plots in OR/WA/CA
  
  # load data - will give warnings about coercion... okay though...
  elev_lat <- read_xlsx(paste0(path_elevlat_actual, "TheData_True.xlsx"),
                        col_types = c("text",             # PLT_CN
                                      rep("skip", 8),     # PLT info
                                      rep("numeric", 3),  # Lat, Lon, Elev
                                      rep("skip", 4),     # Elev/Aspect/Slope
                                      rep("skip", 70),    # Norm81m vars
                                      rep("skip", 2),     # TotBAm, TotTPH
                                      rep("skip", 27), # bam-spp
                                      rep("skip", 27)))   # tph-spp
  
  # better names & convert PLT_CN to factor
  elev_lat <- elev_lat %>%
    mutate(PLT_CN = factor(PLT_CN)) %>%
    rename(t1_PLT_CN   = PLT_CN,
           lat_actual  = LAT_ACTUAL,
           lon_actual  = LON_ACTUAL,
           elev_actual = Elev_m) 
  
  # - Add our plot_id's --------------------------------------------- ####
  preds100_prism <- readRDS(paste0(path_output, 
                                   "FF-1-preds-PRISM_100MCMC.rds"))
  
  elev_lat <- elev_lat %>%
    left_join(preds100_prism %>%
                map_dfr(function(sp.x) {
                  sp.x$newdat.x %>% 
                    select(plot_id, PLT_CN) %>% 
                    rename(t1_PLT_CN = PLT_CN)
                }) %>%
                distinct(.keep_all = TRUE), by = "t1_PLT_CN")
  
  # - Identify plots where spp present at t1 (rom_y, rom_x) --------- ####
  # For each plot, use the TPH at time 1 to indicate if the species was present
  # there. We will use this measure of spp presence (5" dbh) to  develop our
  # naive habitat divisions
  elev_lat <- list(rom_x, rom_y) %>%
    pmap(function(sp.x, sp.y) {
      sp.x %>%
        select(t2_PLT_CN, t1_PLT_CN) %>%
        left_join(elev_lat, by = "t1_PLT_CN") %>%
        # get spp P/A & add elev/lat info to plots
        left_join(sp.y %>%
                    mutate(p_a = as.numeric(t1_live_TPH >= 1)) %>%
                    select(t2_PLT_CN, p_a),
                  by = "t2_PLT_CN")
    })
  
  # - define & assign groups ---------------------------------------- ####
  elev_lat <- elev_lat %>% map(function(sp.x) {
    # 3 even groups (length.out set to 4, then drop minimum)
    e_r <- sp.x %>% filter(p_a == 1) %>% pull(elev_actual) %>% range
    e_b <- seq(e_r[1], e_r[2], length.out = 4)[-1]
    
    l_r <- sp.x %>% filter(p_a == 1) %>% pull(lat_actual) %>% range
    l_b <- seq(l_r[1], l_r[2], length.out = 4)[-1]
    
    # assign groups (new col in this df)
    sp.x %>%
      mutate(
        elev_g3 = case_when(elev_actual <= e_b[1] ~ 1, 
                            elev_actual <= e_b[2] ~ 2, 
                            elev_actual >  e_b[2] ~ 3),
        elev_g  = case_when(elev_g3 == 1 ~ e_b[1], 
                            elev_g3 == 2 ~ e_b[2], 
                            elev_g3 == 3 ~ e_b[3]),
        elev_chr = case_when(elev_g3 == 1 ~ paste0(round(e_r[1]), "m - ",
                                                   round(e_b[1]), "m"), 
                             elev_g3 == 2 ~ paste0(round(e_b[1]), "m - ",
                                                   round(e_b[2]), "m"), 
                             elev_g3 == 3 ~ paste0(round(e_b[2]), "m - ",
                                                   round(e_b[3]), "m")),
        lat_g3 = case_when(lat_actual <= l_b[1] ~ 1, 
                           lat_actual <= l_b[2] ~ 2, 
                           lat_actual >  l_b[2] ~ 3),
        lat_g  = case_when(lat_g3 == 1 ~ l_b[1], 
                           lat_g3 == 2 ~ l_b[2], 
                           lat_g3 == 3 ~ l_b[3]),
        lat_chr = case_when(lat_g3 == 1 ~ paste0(round(l_r[1], 1), "° - ",
                                                 round(l_b[1], 1), "°"), 
                            lat_g3 == 2 ~ paste0(round(l_b[1], 1), "° - ",
                                                 round(l_b[2], 1), "°"), 
                            lat_g3 == 3 ~ paste0(round(l_b[2], 1), "° - ", 
                                                 round(l_b[3], 1), "°"))) %>%
      mutate(el_group = case_when(elev_g3 == 1 & lat_g3 == 1 ~ "LL",
                                  elev_g3 == 1 & lat_g3 == 2 ~ "LM",
                                  elev_g3 == 1 & lat_g3 == 3 ~ "LH",
                                  elev_g3 == 2 & lat_g3 == 1 ~ "ML",
                                  elev_g3 == 2 & lat_g3 == 2 ~ "MM",
                                  elev_g3 == 2 & lat_g3 == 3 ~ "MH",
                                  elev_g3 == 3 & lat_g3 == 1 ~ "HL",
                                  elev_g3 == 3 & lat_g3 == 2 ~ "HM",
                                  elev_g3 == 3 & lat_g3 == 3 ~ "HH"),
             el_group_chr = paste0(elev_chr, " & ", lat_chr))
  })
  
  # - save these groups --------------------------------------------- ####
  saveRDS(elev_lat, paste0(path_data, "FF-3-elev-lat-groups.rds"))
}

# how many plots fall into each group (total and plots where spp present at t1)?
# - note the last row shows how many plots we lacked elevation or lat info for 
#   in the original data set from my PhD (need to check but these might be the 
#   pre-annual inventory plots from OR?)..
elev_lat %>% 
  map(function(sp.x) {
  sp.x %>% 
    group_by(elev_g3, lat_g3) %>%
    summarize(n_plots_and_p_plots = paste0(n(), " & ", sum(p_a)))
}) %>% 
  bind_rows(.id = "spp") %>%
  pivot_wider(names_from = spp, 
              values_from = n_plots_and_p_plots)

# == Estimate change ===================================================== ####
# - Notes ---------------------------------------------------------------- ####
# ATM equation numbers refer to those of Scott et al. (2005; ch4 of the
# 'green book') and are design-based, stratified estimators assuming
# known strata weights. Calculating per-ha estimates, so using their 
# ratio-of-means (ROM) estimation equations. Similarly, for relative-growth
# estimates (net-growth as a percent of gross growth: net-growth / gross-
# growth) we'll also need to use a ROM estimator, but with a slightly
# different structure.
#
# By definition of our range-estimates and PECA-domains, plots not within
# the study area cannot be within the species range (that is, estimated
# range or PEC-domains). With respect to PECA, these plots will always be
# witin the A-domain. Handling for these plots is carreid out below.
# 
# Not sure if the t1|t2 TPH/BAH metrics are proper to expand out using
# the same strata weights as GRM evalid's... I'm not 100% on this and 
# may need to check in with folks...
#
# - add plot-level gross growth to rom_x for rel-ng estimates ------------ ####
rom_x <- list(rom_x = rom_x, rom_y = rom_y) %>%
  pmap(function(rom_x.i, rom_y.i) {
    rom_x.i %>%
      left_join(rom_y.i %>%
                  mutate(grossg_annBAH = netgrowth_annBAH + mortality_annBAH,
                         grossg_annTPH = netgrowth_annTPH + mortality_annTPH) %>%
                  select(t2_PLT_CN, grossg_annBAH, grossg_annTPH),
                by = "t2_PLT_CN") %>%
      # note, the domain indicator function will basically treat any NA
      # values as zero (they are NA b/c no gross growth for that species
      # on that plot b/c the species did not occur on that plot). Therefore
      # handle these cases now by converting these NA's to zero (as was 
      # similarly done for, e.g., prop_for_t1 if there was no forest land
      # on which the species was present at t1 for that particular plot)
      mutate(grossg_annBAH = if_else(is.na(grossg_annBAH), 0, grossg_annBAH),
             grossg_annTPH = if_else(is.na(grossg_annTPH), 0, grossg_annTPH))
  })

# - prep rom_y & rom_x data ---------------------------------------------- ####
# add wts info to x data (only needs to be in either x or y, as they'll be
# combined together with our custom fn)
rom_x <- list(rom_x, rom_wh) %>% pmap(function(x, wh) {
  x %>% 
    inner_join(wh %>% select(pop_area_used, pop_n, EVALID, ESTN_UNIT),
               by = c("EVALID", "ESTN_UNIT"))
})

# split rom_x and rom_y into lists by variable (within spp list-elements)
rom_y <- rom_y %>% map(function(sp.x) {
  sp.x %>% 
    pivot_longer(cols = contains(c("TPH" ,"BAH")),
                 names_to = "var_y",
                 values_to = "y_hi") %>%
    select(EVALID, ESTN_UNIT, STRATUMCD, t2_PLT_CN, plot_id,
           var_y, y_hi) %>%
    split(.$var_y)
})
rom_x <- rom_x %>% map(function(sp.x) {
  sp.x %>%
    pivot_longer(cols = contains(c("prop_", "grossg")),
                 names_to = "var_x",
                 values_to = "x_hi") %>%
    select(EVALID, ESTN_UNIT, STRATUMCD, 
           pop_area_used, pop_n, EXPNS, P2POINTCNT,
           t2_PLT_CN, t1_PLT_CN, plot_id,
           var_x, x_hi) %>%
    split(.$var_x)
})

# create a 'stats key' identifying the types of estimates we want to make
# - That is, for mortality metrics, add prop_for & prop_for_t1 x's, and
#   for net-growth metrics, add prop_for & prop_for_t1t2 x's
stat_key <- data.frame(
  stat_name = paste0("stat_", 1:10),
  var_y     = c(rep("mortality_annBAH", 2),
                rep("mortality_annTPH", 2),
                rep("netgrowth_annBAH", 2),
                rep("netgrowth_annTPH", 2),
                "netgrowth_annBAH", "netgrowth_annTPH"),
  var_x     = c("prop_for", "prop_for_t1",
                "prop_for", "prop_for_t1",
                "prop_for", "prop_for_t1t2", 
                "prop_for", "prop_for_t1t2", 
                "grossg_annBAH", "grossg_annTPH"))

# - fn: post-strat mean & var for per-ha estimates ----------------------- ####
fn.rom_estimates <- function(rom_y.i, rom_x.i, stat_key) {
  # About: 
  # - This function calculates the post-stratified estimates of the mean 
  #   and the standard error of the mean for the statistics specified in
  #   `stat_key` from the observed FIA data queried in script 4-...R
  # - this fn is set up to work with the rom_y and rom_x data after they 
  #   have been subsetted to the spatial domain of interest, meaning:
  #   - overall estimates: no subsetting
  #   - elevation or latitude groups
  #   - PECA-domains
  # - note: for some variables have the same value for all plots in a strata 
  #   (eg, n_h, w_h) or for all strata in the pop (eg, n ie sa_n)... dplyr
  #   can be a pain with this... so need to specify that we want just one
  #   of these values - so using `unique()`
  stat_key %>%
    split(.$stat_name) %>%
    # grab data needed to make that estimate
    imap(function(i, stat_name.i) {
      # handle case when no plots within the group had abs(y) > 0 data
      # (that is y-data popped up in the query)
      tryCatch(expr = {
        rom_x.i[[i$var_x[1]]] %>%
          left_join(rom_y.i[[i$var_y[1]]],
                    by = c("EVALID", "ESTN_UNIT", 
                           "STRATUMCD", "t2_PLT_CN")) %>%
          # if y's are NA on left_join with x, then that means the
          # y-value is 0, so changing those NA's to 0 and those NA
          # var names to the correct name:
          mutate(var_y = if_else(is.na(y_hi), i$var_y[1], var_y),
                 y_hi = if_else(is.na(y_hi), 0, y_hi))
      },
      error = {
        # make the y_hi = 0
        rom_x.i[[i$var_x[1]]] %>%
          # if y's are NA on left_join with x, then that means the
          # y-value is 0, so changing those NA's to 0 and those NA
          # var names to the correct name:
          mutate(var_y = i$var_y[1],
                 y_hi = 0)
      })
    }) %>%
    # estimate mean & SE of the mean
    map(function(df) {
      df %>% 
        # get stratum-level means, variances, & covariances
        group_by(var_y, var_x, EVALID, ESTN_UNIT, STRATUMCD) %>%
        summarize(n   = unique(pop_n),
                  A_t = unique(pop_area_used),
                  n_h = unique(P2POINTCNT), 
                  w_h = unique((EXPNS * n_h) / A_t), 
                  # eq 4.11
                  ybar_h   = sum(y_hi) / n_h,
                  # eq 4.12
                  v_ybar_h = sum((y_hi - ybar_h)^2) / (n_h * (n_h - 1)),
                  # eq 4.3
                  xbar_h   = sum(x_hi) / n_h,
                  # eq 4.4 (note, same math as 4.12 for rel-growth)
                  v_xbar_h = sum((x_hi - xbar_h)^2) / (n_h * (n_h - 1)),
                  # eq 4.19
                  cov_ybar_xbar_h = (sum((y_hi - ybar_h)*(x_hi - xbar_h)) / 
                                       (n_h * (n_h - 1)))) %>%
        ungroup %>%
        # get pop-level mean & SE
        group_by(var_y, var_x) %>%
        summarize(min_nh = min(n_h),
                  n = unique(n),
                  # get post-strat mean (R_hat) eq 4.16
                  y_pop_mean = sum(w_h * ybar_h),
                  x_pop_mean = sum(w_h * xbar_h),
                  R_hat = y_pop_mean / x_pop_mean,
                  # get post-strat SE
                  # - note: canceled out the A_t so we don't need to
                  #   carry them through the eqs below (i.e., eqs for
                  #   pop means not pop totals below)
                  # - eq 4.14
                  v_y_pop_mean =
                    (sum(w_h * n_h * v_ybar_h) +
                       sum((1 - w_h) * (n_h / n) * v_ybar_h)) / n,
                  # - eq 4.6
                  v_x_pop_mean =
                    (sum(w_h * n_h * v_xbar_h) +
                       sum((1 - w_h) * (n_h / n) * v_xbar_h)) / n,
                  # - eq 4.18
                  cov_yx_pop_mean =
                    (sum(w_h * n_h * cov_ybar_xbar_h) +
                       sum((1 - w_h) * (n_h / n) * cov_ybar_xbar_h)) / n,
                  # - eq 4.17
                  v_est = (v_y_pop_mean +
                             R_hat^2 * v_x_pop_mean -
                             2 * R_hat * cov_yx_pop_mean) / (x_pop_mean ^ 2),
                  SEM = sqrt(v_est)) %>%
        rename(est = R_hat) %>%
        select(var_y, var_x, min_nh, est, SEM)
    }) %>%
    bind_rows(.id = "stat_id") %>%
    mutate(var_x = gsub("_t1", "_spp_p_t1", var_x),
           var_x = gsub("prop_", "", var_x),
           var_y = if_else(grepl("grossg", var_x), 
                           paste0("rel", var_y),
                           var_y),
           stat_descr = if_else(grepl("grossg", var_x), 
                                var_y,
                                paste0(var_y, "_", var_x)),
           SEM_perc = SEM/est * 100) %>%
    select(stat_id, stat_descr, var_y, var_x, min_nh, est, SEM, SEM_perc)
}

# - est: overall (no grouping) ------------------------------------------- ####
# get post-stratified estimates of the mean and variance (of the mean) of 
# these variables within the sp-pop (i.e., states, not just range or 
# area of persistence, highest 1/3 elevations, etc.)
est_overall <- list(rom_y, rom_x) %>%
  pmap(~fn.rom_estimates(..1, ..2, stat_key))

# summarize (peak ahead)
est_overall %>% map(function(sp.x) {
  sp.x %>% 
    ungroup %>%
    filter(stat_id %in% paste0("stat_", c((1:4)*2, 9, 10))) %>%
    mutate(ciu = (est + SEM * qnorm(.95)) %>% round(2),
           cil = (est - SEM * qnorm(.95)) %>% round(2),
           est = est %>% round(2),
           value = paste0(est, " (", cil, ", ", ciu, ")"),
           stat_descr = gsub("_for_spp_p_t1", "", stat_descr),
           stat_descr = gsub("t2", "", stat_descr)) %>%
    select(stat_descr, value) %>%
    separate(stat_descr,
             into = c("stat", "units"),
             sep = "_ann") %>%
    pivot_wider(names_from = "stat",
                values_from = "value") %>%
    select(units, netgrowth, mortality, relnetgrowth)
}) %>% bind_rows(.id = "spp")

# - est: elev & lat thirds ----------------------------------------------- ####
# Estimates based on naive-divisions of habitat (even 1/3 breaks in elev or
# latitude) 
est_elevlat <- list(elev_lat, rom_y, rom_x) %>%
  pmap(function(sp.x, rom_y.x, rom_x.x) {
    sp.x %>%
      select(t2_PLT_CN, plot_id, el_group, el_group_chr) %>%
      split(.$el_group) %>%
      map(function(group.x) {
        fn.rom_estimates(
          rom_y.x %>% map(~.x %>% filter(t2_PLT_CN %in% group.x$t2_PLT_CN)),
          rom_x.x %>% map(~.x %>% filter(t2_PLT_CN %in% group.x$t2_PLT_CN)),
          stat_key) %>%
          mutate(el_group  = group.x$el_group[1],
                 el_group_chr = group.x$el_group_chr[1])
      })
  })

# - est: PECA groupings -------------------------------------------------- ####
# could move part of this to run in parallel, but doesn't really take too
# long for 100 MCMC samples... so, keeping this as is for now.
est_peca_cna <- peca_cna %>% map(function(peca_cna.x) {
  list(peca_cna.x, rom_y, rom_x) %>%
    pmap(function(sp.x, rom_y.x, rom_x.x) {
      # reformat groupings to return a list of PECA groups of plot_id names
      # - plots not in a group belong to the 'remain unsuitable / absent' 
      #   group by definition, so give them that label below..
      # - note: a spp could have been present at a plot at t1 but still occur
      #   within an area estimated as outside of the species range at both
      #   time points (i.e., group 'A')... similarly for the expansion group
      sp.x$peca_plots %>% map(function(mcmc.x) {
        sp.x$plot_id %>%
          # associate plots with their PECA group for this mcmc sample
          setNames(mcmc.x) %>%
          # split into list of plot_ids for each PECA group
          split(., names(.)) %>%
          lapply(unname) %>%
          # get post-stratified estimates by grouping (domain)
          imap(function(group.x, group_name.x) {
            fn.rom_estimates(
              rom_y.x %>% map(function(x) {
                x %>% filter(plot_id %in% group.x | 
                               (group_name.x == "A" & is.na(plot_id)))
              }),
              rom_x.x %>% map(function(x) {
                x %>% filter(plot_id %in% group.x | 
                               (group_name.x == "A" & is.na(plot_id)))
              }),
              stat_key)
          })
        # drop any mcmc x est_groups where filtering by plots in our 
        # change data meant that there were no plots in that group
        # filter(!map_lgl(plot_id, is.null))
      })
    })
})

# == Save ests & SEs as Rdata ============================================ ####
save(est_overall, est_peca_cna, est_elevlat,
     file = paste0(path_output, "FF-3-domain-estimates.Rdata"))
