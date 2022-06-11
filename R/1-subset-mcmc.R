###############################################################################
## Project - Forests at the fringe                                         ----
## Script  - 1-subset-mcmc-and-get-PECA-groups.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Forests at the Fringe: comparing observed
## change to projected climate change impacts for five species in the Pacific 
## Northwest' by Karin Kralicek, Tara Barrett, Jay Ver Hoef, and Temesgen
## Hailemariam.
##  
## About - this script:
## - create & save new reduced-size pred objects for this work:
##   - subsets Ch2 models to 100 MCMC samples (every 20th of 2k MCMCs)
##   - subset preds (PRISM, 4 CNA scenarios) to those 100 MCMC samples
##   - Scenarios are:
##     - current climate data from the PRISM climate data
##     - RCP4.5 CCSM4, RCP4.5 HadGEM2-ES, RCP8.5 CCSM4, & RCP8.5 HadGEM2-ES,
##       all for 2080s (2085m) from the ClimateNA downscaled climate data
## - Identify which PECA group plots belong to for each MCMC sample; these are the
##   the groups for which we will estimate observed change (i.e., net-growth and 
##   mortality) over the 10-year period
##
## About - output:
## - FF-1-preds-PRISM_100MCMC.rds, FF-1-preds-CNA_100MCMC.rds
##  - location: path_output
##  - 2k pred objects subset to 100 MCMC samples (still contain mean-pred/se 
##    based on the 2k samples)
## - FF-1-groups_peca_cna.rds
##    - location: path_output
##    - PECA classifications of plots for each of the 100 MCMC samples for the
##      predictions based on either MACA or CNA
##
###############################################################################
# data manipulation
library(magrittr)  # for set_rownames()
library(dplyr)     # for ... E V E R Y T H I N G ...
library(purrr)     
library(tidyr)
library(stringr)   # for string help (e.g. str_detect in functions)
library(sf)        # spatial data, used in (BKDE) range estimation method
library(ggplot2)

# == Paths =============================================================== ####
# path to Rdata
path_top <- paste0("C:/Users/karinkralicek/Documents/", 
                   "projects/active/SI - forest at the fringe/code/")

# paths to intermediate data
path_output <- paste0(path_top, "x - output/")
path_ch2 <- paste0("C:/Users/karinkralicek/Documents/CUI/",
                   "ch2 output sample/")

# == Load|create subset pred objs for FF analysis ======================== ####
# Identify which MCMC samples to keep 
# - keep every 20th of 2k from ch2 work, for a total of 100 MCMC samples 
#   (i.e., MCMC # 20, 40, ..., 2000)
# - this will save on computing time/space for this project...
samp_ids <- (1:100)*20

# PRISM preds, 100 subset
if (file.exists(paste0(path_output, "FF-1-preds-PRISM_100MCMC.rds"))) {
  preds100_prism <- readRDS(paste0(path_output, "FF-1-preds-PRISM_100MCMC.rds"))
} else {
  # create it if it doesn't already exist:
  preds100_prism <- readRDS(paste0(path_ch2, 
                                   "ch2-7-preds-PRISM.rds")) %>%
    pluck("preds_PRISM") %>%
    map(~.x %>% map_at("maps_1k", ~.x[samp_ids]))
  
  saveRDS(preds100_prism,
          paste0(path_output, "FF-1-preds-PRISM_100MCMC.rds"))
}

# CNA preds, 100 subset
if (file.exists(paste0(path_output, "FF-1-preds-CNA_100MCMC.rds"))) {
  preds100_cna <- readRDS(paste0(path_output, "FF-1-preds-CNA_100MCMC.rds"))
} else {
  # create it if it doesn't already exist:
  cna_files <- list.files(path_ch2) %>% grep("ch2-7-preds-CNA", ., value = T)
  cna_labs <- cna_files %>% 
    gsub("ch2-7-preds-CNA_", "", .) %>% 
    gsub("_m2085.rds", "", .)
  
  preds100_cna <- cna_files %>%
    map(function(file.x) {
      readRDS(paste0(path_ch2, file.x)) %>%
        pluck("preds_CNA") %>%
        map(~.x %>% map_at("maps_1k", ~.x[samp_ids]))
    }) %>%
    setNames(cna_labs)
  
  saveRDS(preds100_cna,
          paste0(path_output, "FF-1-preds-CNA_100MCMC.rds"))
}

# == Create PECA groups ================================================== ####
# - load spatial data ---------------------------------------------------- ####
# Spatial data from Ch 2 for BKDE (point/hex/sa) and plotting (sa/plot-poly)
plot_points_true.sf <- readRDS(paste0(path_ch2, "ch2-8-plot_points_true_sf.rds"))
plot_polygons_public.sf <- readRDS(paste0(path_ch2, "ch2-8-plot_polygons_public_sf.rds"))
hex_polygons_true.sf <- readRDS(paste0(path_ch2, "ch2-8-hex_polygons_true_sf.rds"))
study_area.sf <- readRDS(paste0(path_ch2, "ch2-8-study_area_sf.rds"))

# - BKDE-based groups ---------------------------------------------------- ####
# About BKDE groups:
# - This is the same BKDE method for range estimation used in Ch2
# - This process groups plots within an MCMC sample based on the type of 
#   habitat change projected for that data set (e.g., CNA rcp4.5 CCSM4) 
#   relative to the current (PRISM) range estimate.
# - Group types include:
#   - P: persist (remain climatically suitable)
#   - E: expand (become climatically suitable)
#   - C: contract (become climatically unsuitable)
# - A plot could also remain Absent (remain climatically unsuitable) at both
#   time points

# create or load data if already created...
if (file.exists(paste0(path_output, "FF-1-groups_peca_cna.rds"))) {
  peca_cna <- readRDS(paste0(path_output, "FF-1-groups_peca_cna.rds"))
} else {
  # - define projection (crs) -------------------------------------------- ####
  # Primary planar projection:
  # - [distance preserving] UTM Zone 11N (EPSG: 32611)
  crs_UTM11N <- 32611
  
  # Area planar projection:
  # - [area preserving] Albers Equal Area Conic (EPSG: 5070)
  # - only used temporarily for calculations, data projected back to primary
  #   projection post.
  crs_Albers <- 5070
  
  # - define functions --------------------------------------------------- ####
  # fn.kde_args: gets sp-specific global arguments for the kde fn
  fn.kde_args <- function(ppoints_sf.x, hex_sf.x,
                          study_area.sf = study_area.sf,
                          approx_cell_area_ha = 2500/5,
                          mbuffer = 40) {
    # This fn returns: ----
    # - H: bandwidth matrix for a species. Choice of H is determined based on
    #   observed presences (orig data set), aggregated to a hexagon-level to
    #   avoid over-weighting areas where FIA plots were spatially intensified.
    #   The same H is used with each of the 2k MCMC maps for each scenario.
    # - boundary.x: The sp-specific study area boundary. Intersection of the
    #   states boundary with the union of in-sample hexagons for a given spp
    # - bbox_buff: bounding box of boundary.x plus a buffer. A four element 
    #   vector describing the area overwhich the kde grid should be evaluated
    # - n_xy: A 2 element vector with the dimensions of the kde gridsize
    
    # == Define sp-specific study area & gridsize for kde ================ ####
    # boundary for this species
    # - combination of in-sample hexagons and state's boundary
    boundary.x <- hex_sf.x %>% st_union %>% st_intersection(., study_area.sf)
    
    # gridsize & min/max values for the grid with which to est the density
    # (default is a 151 x 151 grid... since our bounding boxes are not square
    #  we probably want points more evenly spaced N-S and E-W)
    # - want to add buffer to our bounding box so that full 95% density polygon
    approx_cell_len_m <- approx_cell_area_ha^.5 * 100
    bbox_buff <- st_bbox(boundary.x) %>% as.numeric
    bbox_buff <- c(bbox_buff[1:2] - approx_cell_len_m * mbuffer,
                   bbox_buff[3:4] + approx_cell_len_m * mbuffer)
    len_xy <- c(bbox_buff[3:4] - bbox_buff[1:2])
    n_xy   <- ceiling(len_xy / approx_cell_len_m)
    
    # clean-up space
    rm(approx_cell_len_m, len_xy)
    
    # == Choose bandwidth matrix (H) for all scenarios/samples =========== ####
    # Filter to presence data & aggregate to hexagon level 
    # - Working with observed presences here so that the H matrix is consistent
    #   across scenarios & MCMC samples (this way, the only uncertainty between
    #   range size polygons will be due to variability between MCMC maps & 
    #   scenarios)
    # - Because FIA plots are not spatially intensified consistently across all
    #   states and ownership, aggregating any p-observations to the hexagon 
    #   level and associating them with the hex centroid.
    ppsp_coords <- hex_sf.x %>%
      st_set_agr(., "constant") %>% # specify attributes as spatially constant
      st_centroid() %>%
      mutate(x = unlist(map(.$geometry,1)),
             y = unlist(map(.$geometry,2)),
             fiahex_id = fiahex_id %>% as.character) %>%
      st_drop_geometry() %>%
      inner_join(ppoints_sf.x %>% st_drop_geometry,
                 by = "fiahex_id") %>%
      group_by(fiahex_id) %>%
      mutate(hex_p = as.numeric(sum(p_a) > 0)) %>%
      filter(hex_p == 1) %>%
      select(-p_a, -plot_id) %>%
      ungroup
    
    # Choose bandwidth matrix
    # - Using a 2-stage plug-in bandwidth selector with unconstrained pilot 
    #   bandwidth based on methods in Chacon & Duong (2010). From what I've
    # - From what I've read (Gitzen et al., 2006, C&D 2010, and others), a
    #   two-stage approach is preferable at reducing bias, but does tend to
    #   result in larger bandwidth estiamtes and more smoothing, but is still
    #   preferable to a 0-stage plug-in selector (like the ks::Hns, normal 
    #   scale selector). Also, the unconstrained pilot bandwidth method has
    #   shown in C&D2010's simulation study to perform better with densities
    #   similar to my species distributions...
    
    # (originally I was specifying deriv.order = 2... but I think what I really
    #  wanted to specify was nstage, as I have done below...)
    H <- ks::Hpi(ppsp_coords %>% select(x, y), 
                 nstage = 2,
                 pilot = "unconstr")
    
    # == Return items ==================================================== ####
    list(H = H,
         boundary.x = boundary.x,
         bbox_buff = bbox_buff,
         n_xy = n_xy)
  }
  # fn.kde_current: sf-polygon, range-size estimate, and T/F plot within range
  fn.kde_current <- function(preds.x, kde_args.x, points_sf.x,
                             hex_plot_walk = hex_plot_walk) {
    # This fn is intended to work with the current scenario data and returns:
    # - kde95_sf: 95% kde contour polygon as a sf object for each MCMC sample
    # - kde95_area: estimate of the range size (ha) for each MCMC sample
    # - kde95_TF: T/F plot was within the each MCMC sample's estimated range
    # How it works:
    # - map by MCMC sample to get kde, kde_95_poly, and then kde_95_area and
    #   T/F for each plot
    
    preds.x %>% map(function(pred.i) {
      if (any(pred.i$pps_pa > 0)) {
        # == aggregate PPS to hex-level, keeping pps-presences ========= ####
        ppsp_coords <- pred.i %>% 
          left_join(hex_plot_walk, by = "plot_id") %>%
          group_by(fiahex_id) %>%
          mutate(hex_p = as.numeric(sum(pps_pa) > 0)) %>%
          filter(hex_p == 1) %>%
          select(-pps_pa, -plot_id) %>%
          ungroup
        
        # == Compute kde =============================================== ####
        # compute kernel density estimate based on pps presence and fixed H
        kde <- ks::kde(x = ppsp_coords %>% select(x, y), 
                       H = kde_args.x$H,
                       gridsize = kde_args.x$n_xy,
                       xmin = kde_args.x$bbox_buff[1:2],
                       xmax = kde_args.x$bbox_buff[3:4])
        
        # == Extract 95% kde polygon =================================== ####
        # extract a contour containing 95% of the estimated density
        # - note: ks::kde() seems to label these backwards, so pulling '5%'
        kde95_sf <- with(kde,
                         contourLines(x = eval.points[[1]],
                                      y = eval.points[[2]],
                                      z = estimate,
                                      levels = cont["5%"]))
        
        # remove potential holes from polygons & make sf object
        # - convert to sf object & intersect with sp-specific boundary area
        kde95_sf <- kde95_sf %>% 
          map(~with(.x, matrix(c(x, y), ncol = 2, byrow = FALSE))) %>%
          map(~rbind(.x, .x[1,])) %>% # close polygons (just in case)
          st_polygon %>% # this step to remove potential holes
          st_sfc %>%     # make into sf so we can specify correct crs
          st_set_crs(crs_UTM11N) %>% # ... this crs!
          st_make_valid %>% # (will freak out next step otherwise)
          st_intersection(., kde_args.x$boundary.x) 
        
        # == Calculate area of range size estimates (ha) =============== ####
        # here using Albers Equal Area Conic (planar to preserve area)
        kde95_area <- kde95_sf %>%
          st_transform(crs = crs_Albers) %>% 
          st_area(.) %>% 
          units::set_units(., ha) %>% 
          as.numeric
        
        # == Get T/F-membership for each plot ========================== ####
        # return vector of T/F for each plot (plot order same as points_sf.x)
        kde95_TF <- kde95_sf %>%
          st_intersects(., points_sf.x, sparse = FALSE) %>% # T/F matrix
          as.vector
        
        # == Return the items we need ================================== ####
        list(kde95_sf = kde95_sf,
             kde95_area = kde95_area,
             kde95_TF = kde95_TF)
        
      } else {
        # the case when all of the pps were absences... so the estimated
        # range is assumed to have area = 0 ha
        list(kde95_sf = NA,
             kde95_area = 0,
             kde95_TF = rep(FALSE, nrow(points_sf.x)))
      }
    })
  }
  # fn.kde_future_diff: returns area estimates and PECA plot-membership
  fn.kde_future_diff <- function(preds.x, kde_current.x, kde_args.x,
                                 points_sf.x,
                                 hex_plot_walk = hex_plot_walk) {
    # This fn is intended to work with a future scenario data set and current
    # scenario kde results to return:
    # 1. a df where rows are MCMC samples and 5 columns are area (ha) summaries 
    #    for {current range-size, future range-size, and area of P, C, and E}.
    #    Recall that:
    #     - persistence: suitable currently and suitable in the future
    #       {overlap between kde95_sf of current & future}
    #     - expansion: unsuitable currently, but suitable in the future
    #       {within future but not within current}
    #     - contraction: suitable currently, but unsuitable in the future
    #       {within current but not within future}
    # 2. list of vectors that classify each plot as P|E|C|A (where A is 
    #    remaining climatically unsuitable, i.e. not within current or future
    #    range estimates)
    # Special consideration:
    # - Because some spp might have all pps as absent (0) for some scenarios, 
    #   need to handle that case (hence, if-else statement in map fn)
    
    # How it works:
    # - map by MCMC sample to get kde, kde_95_poly, and then kde_95_area and
    #   T/F for each plot...
    list(preds.x, kde_current.x) %>% 
      pmap(function(pred.i, kde_current.i) {
        if (any(pred.i$pps_pa > 0)) {
          # == aggregate PPS to hex-level, keeping pps-presences ========= ####
          ppsp_coords <- pred.i %>% 
            left_join(hex_plot_walk, by = "plot_id") %>%
            group_by(fiahex_id) %>%
            mutate(hex_p = as.numeric(sum(pps_pa) > 0)) %>%
            filter(hex_p == 1) %>%
            select(-pps_pa, -plot_id) %>%
            ungroup
          
          # == Compute kde =============================================== ####
          # compute kernel density estimate based on pps presence and fixed H
          kde <- ks::kde(x = ppsp_coords %>% select(x, y), 
                         H = kde_args.x$H,
                         gridsize = kde_args.x$n_xy,
                         xmin = kde_args.x$bbox_buff[1:2],
                         xmax = kde_args.x$bbox_buff[3:4])
          
          # == Extract 95% kde polygon =================================== ####
          # extract a contour containing 95% of the estimated density
          # - note: ks::kde() seems to label these backwards, so pulling '5%'
          kde95_sf <- with(kde,
                           contourLines(x = eval.points[[1]],
                                        y = eval.points[[2]],
                                        z = estimate,
                                        levels = cont["5%"]))
          
          # remove potential holes from polygons & make sf object
          # - convert to sf object & intersect with sp-specific boundary area
          kde95_sf <- kde95_sf %>% 
            map(~with(.x, matrix(c(x, y), ncol = 2, byrow = FALSE))) %>%
            map(~rbind(.x, .x[1,])) %>% # close polygons (just in case)
            st_polygon %>% # this step to remove potential holes
            st_sfc %>%     # make into sf so we can specify correct crs
            st_set_crs(crs_UTM11N) %>% # ... this crs!
            st_make_valid %>% # (will freak out next step otherwise)
            st_intersection(., kde_args.x$boundary.x) 
          
          # == Calculate area of range size estimates (ha) =============== ####
          # here using Albers Equal Area Conic (planar to preserve area)
          # - future range size estimate (ha)
          kde95_area <- kde95_sf %>%
            st_transform(crs = crs_Albers) %>% 
            st_area(.) %>% 
            units::set_units(., ha) %>% 
            as.numeric
          
          # == Get area of P, E, & C, and clean up `odd` case ============ ####
          # handle the odd cases where:
          #  1. some pps_pa were `1`, but still no area is allocated to a 
          #     95% polygon under the future scenario. Fix by setting the 
          #     future_est to 0 & calc the change vars
          #  2. there is no predicted overlap between the current and future 
          #     scenario's range estimates. This happened for abpr and occurs
          #     if the intersection area is empty (e.g. no polygon to calculate
          #     the area from...)
          
          # (case 1)
          kde95_area <- ifelse(length(kde95_area) == 0, 0, kde95_area)
          
          # Area estimated to persist (intersection of kde95s)
          p_area <- kde95_sf %>%
            st_intersection(., kde_current.i$kde95_sf) %>% 
            st_transform(crs = crs_Albers) %>% 
            st_area(.) %>% 
            units::set_units(., ha) %>% 
            as.numeric 
          
          # (case 2)
          p_area <- ifelse(length(p_area) == 0, 0, p_area)
          
          # All areas in a data.frame
          df_area <- data.frame(current_ha = kde_current.i$kde95_area,
                                future_ha  = kde95_area,
                                P_ha = p_area,
                                E_ha = kde95_area - p_area,
                                C_ha = kde_current.i$kde95_area - p_area)
          
          # == Get T/F-membership for each plot & PECA assignment ======== ####
          # return vector of T/F for each plot (plot order same as points_sf.x)
          kde95_TF <- kde95_sf %>%
            st_intersects(., points_sf.x, sparse = FALSE) %>% # T/F matrix
            as.vector
          
          # habitat is...
          peca_plots <- case_when(
            # ... remaining suitable
            (kde_current.i$kde95_TF == T) & (kde95_TF == T) ~ "P",
            # ... becoming unsuitable
            (kde_current.i$kde95_TF == T) & (kde95_TF == F) ~ "C",
            # ... becoming suitable
            (kde_current.i$kde95_TF == F) & (kde95_TF == T) ~ "E",
            # ... remaining unsuitable
            (kde_current.i$kde95_TF == F) & (kde95_TF == F) ~ "A") 
          
          # == Return the items we need ================================== ####
          list(df_area = df_area,
               peca_plots = peca_plots,
               plot_id = points_sf.x$plot_id)
          
        } else {
          # the case when all of the future pps were absences... so the 
          # estimated range in the future is assumed to have area = 0 ha
          list(df_area = data.frame(current_ha = kde_current.i$kde95_area,
                                    future_ha  = 0,
                                    P_ha = 0,
                                    E_ha = 0 - 0,
                                    C_ha = kde_current.i$kde95_area - 0),
               peca_plots = case_when(
                 (kde_current.i$kde95_TF == T) ~ "C", # becoming unsuitable
                 (kde_current.i$kde95_TF == F) ~ "A"), # remaining unsuitable
               plot_id = points_sf.x$plot_id) 
        }
      })
  }
  
  # - 1. make crosswalk between plots & their containing hexagons -------- ####
  # create crosswalk between plots & containing hexagons' centroid coords
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
                   by = "fiahex_id")
    }) %>%
    bind_rows() %>%
    distinct()
  
  # - 2. choose H based on observed data --------------------------------- ####
  kde_args <- list(plot_points_true.sf, hex_polygons_true.sf) %>%
    pmap(~fn.kde_args(..1, ..2, study_area.sf,
                      approx_cell_area_ha = 2500/5,
                      mbuffer = 40))
  
  # - 3. get/save {poly/area/plot-mem} for current (prism; 8 min) -------- ####
  start_time <- Sys.time()
  kde_current <- 
    list(preds100_prism %>% map(function(sp.x) {
      sp.x$maps_1k %>% map(~.x %>% select(plot_id, pps_pa))
    }), 
    kde_args,
    plot_points_true.sf %>% map(~.x %>% select(plot_id))) %>%
    pmap(~fn.kde_current(..1, ..2, ..3, hex_plot_walk))
  run_time <- Sys.time() - start_time
  
  saveRDS(kde_current, 
          file = paste0(path_output, "FF-1-kde-current.rds"))
  
  # - 4. get/save {areas/plot-peca} for future & dif (cna; 4*9.5 min) ------ ####
  start_time <- Sys.time()
  peca_cna <- preds100_cna %>% map(function(scenario.x) {
    scenario.x %>% 
      map(function(sp.x) {
        sp.x$maps_1k %>% map(~.x %>% select(plot_id, pps_pa))
      }) %>%
      list(., 
           kde_current, 
           kde_args,
           plot_points_true.sf %>% map(~.x %>% select(plot_id))) %>%
      pmap(~fn.kde_future_diff(..1, ..2, ..3, ..4, hex_plot_walk)) %>%
      map(function(sp.x) {
        sp.x %>% 
          transpose %>%
          map_at("df_area", ~.x %>% bind_rows(.id = "mcmc_n")) %>%
          map_at("plot_id", ~.x[[1]])
      })
  })
  run_time <- Sys.time() - start_time
  
  saveRDS(peca_cna, 
          file = paste0(path_output, "FF-1-groups_peca_cna.rds"))
}

# == Peak at PECA groups ================================================= ####
# - any PECA groups that were assigned very few (<50) plots? ------------- ####
# If just grouping by PECA... some low plot counts issues?
n_thresh <- 50 # at least 50 plots within spp' study area in a PECA group
peca_cna %>%
  map_depth(2, ~.x$peca_plots %>% 
              map(function(mcmc.x) {
                mcmc.x %>%
                  table %>%
                  as.data.frame() %>%
                  set_names(c("peca", "n_plots"))
              }) %>%
              bind_rows(.id = "mcmc_n")) %>%
  transpose() %>%
  map(function(ds.x) {
    ds.x %>%
      bind_rows(.id = "data_set") %>%
      group_by(data_set, peca) %>%
      summarize(n_mcmc = n_distinct(mcmc_n),
                min = min(n_plots), 
                q_10 = quantile(n_plots, .1) %>% as.numeric,
                q_50 = median(n_plots), 
                q_90 = quantile(n_plots, .9) %>% as.numeric,
                max = max(n_plots)) %>%
      pivot_longer(cols = c(min, q_10, q_50, q_90, max),
                   names_to = "stat",
                   values_to = "n_plots") %>%
      mutate(stat = factor(stat,
                           levels = c("min", "q_10", "q_50", "q_90", "max")))
  }) %>%
  bind_rows(.id = "spp") %>%
  filter(n_plots < n_thresh) %>%
  # return highest number of MCMC samples for which this was an issue... e.g.:
  # - if 'max' returned, then this was an issue in all 100 MCMC samples
  # - if 'q_10' returned, then this was an issue in ~10 of the MCMC samples
  #   (note for example that here, the 'n_plots' is the number of plots in this
  #    group for the 10th MCMC sample when samples are sorted by the number 
  #    of plots in this group)
  arrange(desc(stat)) %>%
  distinct(spp, data_set, peca, .keep_all = TRUE) %>%
  arrange(spp, data_set, peca)
