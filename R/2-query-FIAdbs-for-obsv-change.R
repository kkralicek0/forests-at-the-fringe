###############################################################################
## Project - Forests at the fringe                                         ----
## Script  - 2-query-FIAdbs-for-obsv-change.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Forests at the Fringe: Comparing observed
## change to projected climate change impacts for five species in the Pacific 
## Northwest, United States' by Karin Kralicek, Tara M. Barrett, Jay M. Ver 
## Hoef, and Temesgen Hailemariam (Front. For. Glob. Change 5:966953. 
## doi: 10.3389/ffgc.2022.966953).
##
## About - this script:                                                    ----
## - End-goal:
##   Queries the FIA data (previously downloaded as SQLite dbs for each state 
##   from the Datamart) for summaries used in subsequent scripts (5) to 
##   estimate of observed change (net-growth and mortality) of tally-trees 
##   (>=5"dia) on different accessible forest lands in OR, WA, and CA:
##   1. {mortality, net-growth} per ha on forest land within the spp-states
##   2. {mortality} per ha on forest land  within the spp-states where the 
##      species was present at t1
##   3. {net-growth} per ha on forest land  within the spp-states where the 
##      species was present at t1 or t2.  
## - ROM-estimation:
##   To obtain per-ha estimates in subsequent scripts (3-...R), we used a ROM
##   (ratio-of-means) estimator. For clarity, in this script we refer to the
##   ROM-estimator's numerator and denominator for queries in the script. 
##   In later scripts, these plot-level stats will be combined for PECA areas 
##   (persist/expand/contract/remain-absent) using the stratum-wts that are 
##   queried and calculated at the end of this script.
## - t1 selection-rule:
##   Calculations based on t1 selection rule, such that the TPA expansion that
##   applied at t1 applies for all components of change except ingrowth, which
##   will have t2 TPA apply (that is, for ingrowth that survives to t2 - if 
##   microplot trees cross threshold but die before t2, then t1 TPA applies);
##   macroplot-ingrowth (aka ongrowth) is ignored. 
##
## About - output:                                                         ----
## - FF-2-rom-y-data.rds
##   - location: path_data
##   - queried data from FIA SQLite dbs with observed change in net-growth and
##     mortality (TPH/yr and BAH/yr, so annualized) by species. ATM, also calc 
##     conversion metrics and living TPH/BAH at t1 and t2, though need to check
##     if the same weights used in A2A estimation should apply. Here only trees
##     with dbh >= 5" at t1 (or the midpt of the remeasurement interval for 
##     ingrowth) are considered, and mortality refers to natural death (i.e., 
##     excludes harvested, etc. trees). Summaries are for trees on forest land 
##     at both t1 and t2, but per-ha values are per-ha in the state of X. Also
##     includes the stratum's EXPNS & P2POINTCNT (number of p2 ground plots);
##     which are used in stratum wt calcs in script 3 (i.e., area of the 
##     stratum = EXPNS * P2POINTCNT)
## - FF-2-rom-x-data.rds
##   - location: path_data
##   - queried data from FIA SQLite dbs with plot-level estimates of the 
##     proportion of the plot (based on macroplot area) that was within the
##     domain of interest, calculated for 3 different types of domains:
##     (1) forest land within the spp-states, (2) [#1] on which the spp was
##     present at t1, and (3) [#1] on which the spp was present at t1 or t2.
##     These are condition-based estimates that have been adjusted for non-
##     response within the stratum for macroplots. Also includes the 
##     stratum's EXPNS & P2POINTCNT (number of p2 ground plots); which are 
##     used in stratum wt calcs in script 3.
## - FF-2-rom-wh-data.rds
##   - location: path_data
##   - list of dfs by spp containing information for stratum weights such
##     as the total area of the population and number of ground plots in
##     the population (as well as the break down by state and estn_unit),
##     where population determined using the state membership object below
## - FF-2-state-membership.rds
##   - location: path_data
##   - list by species of the state in which each plot in our analysis occurs 
##     within; this was used in this script to calculate ROM_wh (wts info) by
##     determining which states are included in the 'population' for species-
##     specific estimation... the idea being we only want to include estimates
##     for forest land (#1 in ROM_x) in states that intersect our spp-specific
##     study areas. For example, qudo is endemic to CA (only present in CA), 
##     but in developing our model we included absence plots that extended 
##     slightly into OR; so our domain of interest when doing domain estimation
##     for qudo will be 'forest land in CA and OR' as opposed to OR/CA/WA. 
###############################################################################
library(dplyr) # for ...E V E R Y T H I N G ...
library(tidyr) 
library(purrr) 
library(glue)  # for complex queries (SO@r2evans: parameterized over glue_sql)
library(DBI)   # 'Database Interface' pkg, to connect to the SQLite databases

# == Paths =============================================================== ####
# path to Rdata
path_top <- paste0("C:/Users/karinkralicek/Documents/", 
                   "projects/active/SI - forest at the fringe/code/")

# paths to intermediate data
path_output <- paste0(path_top, "x - output/")
path_data  <- paste0(path_top, "x - data/")
path_ch2_data <- paste0("C:/Users/karinkralicek/Documents/CUI/",
                          "ch2 output sample/")

# path to the SQLite dbs (Datamart) for CA, OR, and WA
# - these downloaded on 5//2022
path_FIAdb <- paste0(path_top, "data - FIAdb/")
  
# == Create states' querying fn & set query args ========================= ####
# - define fn to query for states' SQLite dbs ---------------------------- ####
# Each state in a different SQLite db, so going state by state
fn.query_states <- function(states, query) {
  states %>%
    lapply(function(state.x, query.x = query) {
      # take a reference to the db
      con <- dbConnect(RSQLite::SQLite(), 
                       paste0(path_FIAdb, "FIADB_", state.x, ".db"))
      
      # execute query
      df.x <- dbGetQuery(conn = con, statement = query.x)
      
      # close connection
      dbDisconnect(con)
      
      # return a data.frame containing the data for state.x
      return(df.x)
    }) %>%
    # bind into one df & create a 'state' variable
    setNames(states) %>%
    bind_rows(.id = "state")
}

# - set query arguments and spp of interest ------------------------------ ####
query_args <- data.frame(state = c("CA", "OR", "WA"),
                         evalid = c(61903, 411903, 531903))

# species codes in order of {ABPR, PSME, QUDO, QUGA4, QUKE}
# - note PSMEM: after this data-querying step, we will only be working with
#   plots within the PSMEM variety's range, so we don't have to worry about
#   filtering to that variety here...
spp_cd <- c(22, 202, 807, 815, 818)

# == aside - peak at {eval/estn_units/strata} summaries ================== ####
# peak at evaluation descriptions
fn.query_states(query_args$state, 'SELECT * FROM POP_EVAL') %>%
  filter(GROWTH_ACCT == "Y" & END_INVYR == 2019) %>%
  select(EVALID, EVAL_DESCR)

# Q: what are desc for evalid, estn_units, and examples for strata?
desc_table <- glue_sql(
  'SELECT ps.STATECD, pe.EVALID, pe.EVAL_DESCR, 
          peu.ESTN_UNIT, peu.ESTN_UNIT_DESCR,
          ps.STRATUMCD, ps.STRATUM_DESCR
     FROM POP_EVAL AS pe
          INNER JOIN POP_ESTN_UNIT AS peu
          ON pe.EVALID = peu.EVALID
          INNER JOIN POP_STRATUM AS ps
          ON pe.EVALID = ps.EVALID
             AND peu.ESTN_UNIT = ps.ESTN_UNIT
    WHERE pe.EVALID IN ({query_args$evalid*})'
) %>% fn.query_states(query_args$state, .)

# - what are the distinct evaluations and how many estn_units are there?
desc_table %>% 
  group_by(STATECD, EVALID, EVAL_DESCR) %>%
  summarise(n_estn_units = n_distinct(ESTN_UNIT))

# - what are the distinct estn_units & how many strata do they have?
desc_table %>% 
  group_by(EVALID, ESTN_UNIT, ESTN_UNIT_DESCR) %>%
  summarise(n_strata = n_distinct(STRATUMCD)) %>%
  arrange(EVALID, ESTN_UNIT)

# - what are some examples of the types of strata in each estn_unit?
desc_table %>% 
  distinct(EVALID, ESTN_UNIT, ESTN_UNIT_DESCR,
           .keep_all = TRUE) %>%
  select(EVALID, ESTN_UNIT, ESTN_UNIT_DESCR, STRATUMCD, STRATUM_DESCR) %>%
  arrange(EVALID, ESTN_UNIT)

# == Assemble query: ROM_y, numerator ==================================== ####
# - Notes on SQLite & glue_sql ------------------------------------------- ####
# Motto: "Work w/in a stratum, w/in an estimation-unit, w/in an evaluation."
# Notes on SQLite:
# - SQL: FROM > WHERE > GROUP BY > HAVING > SELECT > ORDER BY
# - With SQLite dbs, you can SELECT columns that are omitted from GROUP BY.
#   However, this is not always true for other database types; so, if porting
#   these queries over, may need to add those columns to the GROUP BY clause.
# - No PIVOT function in SQLite; therefore, repeating calculations for each
#   species via GROUP BY 
# Notes on query assembly:
# - Instead of creating new or temporary tables in the dbs, this code works
#   by writing sub-queries. Using the glue pkg, these sub-queries are then
#   referenced in (i.e., pasted into) later sub-queries. The output from 
#   individual sub-queries can be obtained with the fn.query_states() function.
#   The naming convention for these subqueries is `sq_*`. 
# - A caution on glue_sql...
#   Do not comment on the last line of a subquery; if you do, it will cause 
#   the rest of the line to be commented out when ({*})-ed into other queries
# - sq_y_1: id remeas-plots w. accessible forestland at t1 & t2 ---------- ####
# About:
# - Get list of remeasurement plots with accessible forest-land conditions at
#   both time 1 (t1) and time 2 (t2). This evaluation will give use partially 
#   non-sampled plots, but excludes fully non-sampled plots. It would also 
#   give use non-forest plots, but we're only interested in plots where we can
#   measure change at t2; therefore filtering to plots with acc. forest-land
#   cond at both t1 and t2 (trees only sampled/measured on this type of cond).
# - Removing 33 plots measured in OR in 1999 (not actually part of the annual
#   inventory)
# 
# Steps:
# - sq_y_1_1: plots in the evaluation with an accessible forest-land cond at t2.
#   (return: post-strat wts, t2 plot-level info)
# - sq_y_1_2: plots from sq_y_1 with an accessible forest-land cond at t1. 
#   (return: t1 plot-level info)
# - sq_y_1: add info from sq_y_1_1 to sq_y_1_2's list of plots

# (the commented out vars below were used to check things)
sq_y_1_1 <- glue_sql( 
  'SELECT psm.EVALID, 
          psm.ESTN_UNIT, 
          -- peu.ESTN_UNIT_DESCR,
          -- peu.AREA_USED,
          -- peu.P1PNTCNT_EU,
          -- ... stratum-level attributes
          psm.STRATUMCD,
          -- psm.STRATUM_DESCR,
          -- psm.P1POINTCNT, 
          psm.EXPNS,
          psm.P2POINTCNT, 
          psm.ADJ_FACTOR_MACR, 
          psm.ADJ_FACTOR_SUBP, 
          psm.ADJ_FACTOR_MICR,
          -- ... plot-level attributes
          plt.CN AS t2_PLT_CN,
          plt.PREV_PLT_CN AS t1_PLT_CN,
          plt.MEASYEAR AS t2_MEASYEAR,
          plt.REMPER
     FROM POP_EVAL AS pe
          INNER JOIN POP_ESTN_UNIT AS peu
          ON pe.EVALID = peu.EVALID
          INNER JOIN POP_STRATUM AS psm
          ON pe.EVALID = psm.EVALID
             AND peu.ESTN_UNIT = psm.ESTN_UNIT
    	    INNER JOIN POP_PLOT_STRATUM_ASSGN AS ppsa
    	    ON ppsa.STRATUM_CN = psm.CN
    	    INNER JOIN PLOT AS plt
    	    ON plt.CN = ppsa.PLT_CN
    	    INNER JOIN COND AS cnd
    	    ON cnd.PLT_CN = ppsa.PLT_CN
    WHERE psm.EVALID IN ({query_args$evalid*})
          AND cnd.COND_STATUS_CD = 1 -- at least one forest cond at t2
 GROUP BY ppsa.EVALID, -- apply aggregate function to these groups
          ppsa.ESTN_UNIT,
          ppsa.STRATUMCD,
          plt.CN,
          plt.PREV_PLT_CN',
  .con = con)

sq_y_1_2 <- glue_sql( 
  'SELECT sq_y_1_1.t1_PLT_CN,
          plt.MEASYEAR AS t1_MEASYEAR
     FROM ({sq_y_1_1}) AS sq_y_1_1
          INNER JOIN PLOT AS plt
          ON plt.CN = sq_y_1_1.t1_PLT_CN
          INNER JOIN COND AS cnd
          ON cnd.PLT_CN = plt.CN
    WHERE cnd.COND_STATUS_CD = 1
          -- Drop 33 PJ plots from OR measured in 1999...     
          -- ... these were used instead of initially measuring in 2001 as 
          --     part of the annual inventory... so drop.
          AND plt.MEASYEAR >= 2000 
 GROUP BY plt.CN',
  .con = con)

sq_y_1 <- glue_sql(
  'SELECT sq_y_1_1.*, 
          sq_y_1_2.t1_MEASYEAR
     FROM ({sq_y_1_1}) AS sq_y_1_1
          INNER JOIN ({sq_y_1_2}) AS sq_y_1_2
          ON sq_y_1_2.t1_PLT_CN = sq_y_1_1.t1_PLT_CN',
  .con = con)

# - sq_y_2: pull tree records -------------------------------------------- ####
# Pull the tree records that we will feed into our calculations. To save time,
# calculate the following variables:
# - BA m^2 for each tree at t1, t2, the midpoint of the remeasurement period,
#   and for a tree with dia = 5" (threshold to calculate ingrowth and growth
#   of ingrowth {surviving, cut, or mortality} against).
# - TPH expansion factor for the tree's corresponding plot-size
# - Adjustment factor for partially non-sampled plots corresponding to the
#   tree's plot-size (same plot-size as used in the TPH expansion factor)
# Some notes on the variables these calculations are based on (all these 
# from the GRM-component table):
# - DIA_BEGIN = modeled DIA (in GRM-Begin table) or tree.PREVDIA if the
#   modeled DIA doesn't exist
# - DIA_MIDPT = modeled DIA from GRM-Midpt table (at least for cut and
#   mortality trees)
# - DIA_END = same thing as tree.DIA
# - SUBP_TPAGROW_UNADJ_AL_FOREST = tpa expansion that agrees with the t1 
#   selection rule such that the TPA that applied at t1 is here except for 
#   ingrowth that was live at t2 which has t2 tpa expansion (ingrowth that
#   died before t2 gets t1 expansion, i.e. microplot, because it could only
#   be observed on this smaller footprint).
# - SUBP_SUBPTYP_GRM_AL_FOREST = subplot type corresponding to the variable
#   above (); used to easily indicate which ADJ_FACTOR_* to apply
# - COMPONENT (SUBP_COMPONENT_AL_FOREST): 
#   - {SURVIVOR, CUT1, MORTALITY1} were in the sample and a tally-tree at t1
#   - {INGROWTH} were living ingrowth on the subplot at t2
#   - {CUT2, MORTALITY2} ingrowth that was cut or died before t2, these trees
#     were on the microplot at t1, became tally-trees by the midpt of the 
#     REMPER but then were cut or died before t2.

sq_y_2 <- glue_sql(
  'SELECT sq_y_1.t2_PLT_CN,
          sq_y_1.REMPER,
          tree.CN AS t2_tree_CN,
          tree.SPCD,
          grm.SUBP_COMPONENT_AL_FOREST AS COMPONENT,
          tree.AGENTCD,
          tree.RECONCILECD,
          tree.STATUSCD,
          tree.DIA,
          -- calc BA in m2 for t1, ingrowth dia threshold, midpt, and t2
          (grm.DIA_BEGIN * grm.DIA_BEGIN * 0.001252) AS BAm2_BEGIN,
          (5 * 5 * 0.001252) AS BAm2_5in,
          (grm.DIA_MIDPT * grm.DIA_MIDPT * 0.001252) AS BAm2_MIDPT,
          (grm.DIA_END * grm.DIA_END * 0.001252) AS BAm2_END,
          -- convert TPA expansion factor to TPH
          (COALESCE(grm.SUBP_TPAGROW_UNADJ_AL_FOREST, 0) * 2.47105) AS grm_TPH,
          -- associate correct adjustment factor
          (CASE
             WHEN COALESCE(grm.SUBP_SUBPTYP_GRM_AL_FOREST, 0) = 0 THEN (0)
             WHEN grm.SUBP_SUBPTYP_GRM_AL_FOREST = 1 THEN ADJ_FACTOR_SUBP
             WHEN grm.SUBP_SUBPTYP_GRM_AL_FOREST = 2 THEN ADJ_FACTOR_MICR
             WHEN grm.SUBP_SUBPTYP_GRM_AL_FOREST = 3 THEN ADJ_FACTOR_MACR
             ELSE (0) END) AS adj_factor
     FROM ({sq_y_1}) AS sq_y_1
          INNER JOIN COND AS cond
          ON cond.PLT_CN = sq_y_1.t2_PLT_CN
          INNER JOIN COND AS prev_cond
          ON prev_cond.PLT_CN = sq_y_1.t1_PLT_CN
          INNER JOIN TREE AS tree
          ON tree.PLT_CN = sq_y_1.t2_PLT_CN
             AND tree.CONDID = cond.CONDID
             AND tree.PREVCOND = prev_cond.CONDID
          LEFT JOIN TREE_GRM_COMPONENT AS grm
          ON grm.TRE_CN = tree.CN
    WHERE tree.SPCD IN ({spp_cd*})
          -- only trees on acc. forest land at t1 and t2
          AND cond.COND_STATUS_CD = 1
          AND prev_cond.COND_STATUS_CD = 1
          -- only tally-trees that were live at t1 or are t2-ingrowth
          AND (grm.DIA_END >= 5 OR grm.DIA_MIDPT >= 5)
          AND (tree.PREV_STATUS_CD = 1 OR tree.PREV_STATUS_CD IS NULL)
          -- exclude NOT USED components for GRM
          -- ... all with 2 are dead ingrowth and we will handle those later on
          AND grm.SUBP_SUBPTYP_GRM_AL_FOREST <> 0',
  .con = con) 

# - sq_y_3: summarize by components of change (green book) --------------- ####
# To make net-growth/mortality and other calculations easier, calculating the
# contribution of each tree to the various components of change:
# - Component definitions:
#   (I)  Ingrowth: {INGROWTH, CUT2, MORTALITY2} trees
#   (Gs) Growth of survivors: {SURVIVOR} trees
#   (Gi) Growth of ingrowth: {INGROWTH}
#   (Gm) Growth of mortality: {MORTALITY1, MORTALITY2}
#   (Gm) Growth of cut: {CUT1, CUT2}
#   (M)  Mortality: {MORTALITY1, MORTALITY2}
# Peaking ahead, the BAH equations from these components will be:
#   net-growth = (I + Gs + Gi + Gm + Gc) - M
#   mortality  = M

sq_y_3 <- glue_sql(
  "SELECT t2_PLT_CN,
          t2_tree_CN,
          SPCD,
          COMPONENT,
          REMPER,
          grm_TPH,
          adj_factor,
          BAm2_BEGIN,
          BAm2_END,
          -- calc contribution of each to component of change as per green book
          -- .. each tree will either have a value for that component or be zero
          (CASE WHEN COMPONENT IN ('INGROWTH', 'CUT2', 'MORTALITY2')
                THEN BAm2_5in ELSE (0) END) AS comp_I,
          (CASE WHEN COMPONENT = 'SURVIVOR'
                THEN (BAm2_END - BAm2_BEGIN) ELSE (0) END) AS comp_Gs,
          (CASE WHEN COMPONENT = 'INGROWTH'
                THEN (BAm2_END - BAm2_5in) ELSE (0) END) AS comp_Gi,
          (CASE WHEN COMPONENT IN ('MORTALITY1', 'CUT1')
                THEN (BAm2_MIDPT - BAm2_BEGIN) ELSE (0) END) AS comp_Gmc1,
          (CASE WHEN COMPONENT IN ('MORTALITY2', 'CUT2')
                THEN (BAm2_MIDPT - BAm2_5in) ELSE (0) END) AS comp_Gmc2,
          (CASE WHEN COMPONENT IN ('MORTALITY1', 'MORTALITY2')
                THEN BAm2_MIDPT ELSE (0) END) AS comp_M
    FROM ({sq_y_2})",
  .con = con)

# - sq_y_4: calculate plot-level metrics by spp -------------------------- ####
# Summarize net-growth and mortality to the plot-level by spp as annualized 
# BAH and TPH, as well as the BAH and TPH of living trees at t1 and t2. Note, 
# these plot-level estimates will need to be combined using stratum weights, 
# which is the # acres in the strata divided by the number of sampled plots.
# Recall - net_growth includes growth of removals, excludes natural mortality.

sq_y_4 <- glue_sql(
  "SELECT t2_PLT_CN,
          SPCD,
          -- annBAH and BAH metrics
          SUM(adj_factor * grm_TPH / REMPER *
              (comp_I + comp_Gs + comp_Gi + comp_Gmc1 + comp_Gmc2 - comp_M)) AS netgrowth_annBAH,
          SUM(adj_factor * grm_TPH / REMPER * comp_M) AS mortality_annBAH,
          SUM(adj_factor * grm_TPH *
              (CASE WHEN COMPONENT IN ('SURVIVOR', 'MORTALITY1', 'CUT1')
                    THEN BAm2_BEGIN ELSE (0) END)) AS t1_live_BAH,
          SUM(adj_factor * grm_TPH *
              (CASE WHEN COMPONENT IN ('SURVIVOR', 'INGROWTH')
                    THEN BAm2_END ELSE (0) END)) AS t2_live_BAH,
          -- annTPH and TPH metrics
          SUM(adj_factor * grm_TPH / REMPER *
              (CASE WHEN COMPONENT IN ('SURVIVOR', 'CUT1')
                    THEN (0) -- theoretically still around
                    WHEN COMPONENT IN ('INGROWTH', 'CUT2')
                    THEN (1) -- theoretically new and still around
                    WHEN COMPONENT IN ('MORTALITY1', 'MORTALITY2')
                    THEN (-1)
                    END)) AS netgrowth_annTPH,
          SUM(adj_factor * grm_TPH / REMPER *
              (CASE WHEN COMPONENT IN ('MORTALITY1', 'MORTALITY2') 
                    THEN (1) ELSE (0) END)) AS mortality_annTPH,
          SUM(adj_factor * grm_TPH *
              (CASE WHEN COMPONENT IN ('SURVIVOR', 'MORTALITY1', 'CUT1')
                    THEN (1) ELSE (0) END)) AS t1_live_TPH,
          SUM(adj_factor * grm_TPH *
              (CASE WHEN COMPONENT IN ('SURVIVOR', 'INGROWTH')
                    THEN (1) ELSE (0) END)) AS t2_live_TPH
     FROM ({sq_y_3})
 GROUP BY t2_PLT_CN,
          SPCD",
  .con = con)

# - sq_y: assemble full query -------------------------------------------- ####
sq_y <- glue_sql(
  "SELECT sq_y_1.EVALID,
          sq_y_1.ESTN_UNIT, 
          sq_y_1.STRATUMCD,
          sq_y_1.EXPNS,
          sq_y_1.P2POINTCNT,
          sq_y_1.t2_PLT_CN, 
          sq_y_1.t1_PLT_CN, 
          sq_y_1.t2_MEASYEAR, 
          sq_y_1.t1_MEASYEAR, 
          sq_y_1.REMPER,
          sq_y_4.SPCD,
          sq_y_4.netgrowth_annBAH,
          sq_y_4.mortality_annBAH,
          sq_y_4.netgrowth_annTPH,
          sq_y_4.mortality_annTPH,
          -- ... not sure if these vars below are okay to expand with 
          --     stratum weights based on this evalid...
          sq_y_4.t1_live_BAH,
          sq_y_4.t2_live_BAH,
          sq_y_4.t1_live_TPH,
          sq_y_4.t2_live_TPH
     FROM ({sq_y_1}) AS sq_y_1
          INNER JOIN ({sq_y_4}) AS sq_y_4
          ON sq_y_4.t2_PLT_CN = sq_y_1.t2_PLT_CN
GROUP BY sq_y_1.t2_PLT_CN,
         sq_y_4.SPCD",
  .con = con)

# == Assemble query: ROM_x, denominator ================================== ####
# About:
# - Mortality and net-growth will be easier (for me) to understand if presented
#   on per-ha basis with respect to:
#     a {mortality, net-growth} per ha on forest land within the spp-states
#     b {mortality} per ha on forest land  within the spp-states where the 
#       species was present at t1
#     c {net-growth} per ha on forest land  within the spp-states where the 
#       species was present at t1 or t2.  
#   That is, as opposed to reporting 'per ha' on all land within the state 
#   (including developed, etc.). Any of these types of estimate will necessitate
#   the use of a ratio-of-means estimator. Here we calculate the denominator for
#   those estimates for all plots within those populations, which will be the
#   the proportion of the plot that is within the domain of (a), (b), or (c). 
# - we are keeping the denominator and numerator separate because we will 
#   need to account for the correlation between plot-level values in our 
#   variance estimates.
# Note:
# - For simplicity, abbreviating estimates for {a, b, c} in the sub-query below
#   as {prop_for, prop_for_t1, prop_for_t1t2} 
# - sq_y_1 will have all plots that had at least one condition that was 
#   accessible forest land at both t1 and t2 that are within the EVALID.
# - CONDPROP_UNADJ is the unadjusted proportion of the plot that is in the
#   condition and is equal to either SUBPPROP_UNADJ or MACRPROP_UNADJ, depending
#   on the value of PROP_BASIS.
# - PROP_BASIS is indicating what type of fixed-size subplots were installed 
#   when this plot was sampled. Is used to toggle to the proper adjustment
#   factor for non-response (i.e., ADJ_FACTOR_MACR in the PNW, but could be
#   ADJ_FACTOR_SUBP in other areas) in the stratum in which the plot occurs.
# - ADJ_FACTOR_MACR is the adjustment factor for non-response within the 
#   stratum for that type of plot. B/c of our interest in area, the largest
#   type of plot used in the PNW to map conditions is the Macroplot, so this
#   is what will be used here.
#
# - sq_x_1: forest land (ha) --------------------------------------------- ####
# get 'prop_for'
sq_x_1 <- glue_sql(
  "SELECT sq_y_1.t2_PLT_CN,
          SUM(cond.CONDPROP_UNADJ * sq_y_1.ADJ_FACTOR_MACR) AS prop_for
     FROM ({sq_y_1}) AS sq_y_1
          INNER JOIN COND AS cond
          ON cond.PLT_CN = sq_y_1.t2_PLT_CN
          INNER JOIN COND AS prev_cond
          ON prev_cond.PLT_CN = sq_y_1.t1_PLT_CN
    WHERE 1 = 1
          -- only areas that were acc. forest land at t1 and t2
          AND cond.COND_STATUS_CD = 1
          AND prev_cond.COND_STATUS_CD = 1
 GROUP BY sq_y_1.t2_PLT_CN", 
  .con = con) 

# - sq_x_2: forest land with spp at t1 (ha; for mortality est) ----------- ####
# get 'prop_for_t1'
# - first, return only [plots x cond] where that spp was live at t1 (on land 
#   that was accessible forest land at both t1 and t2)
sq_x_2_1 <- glue_sql(
  "SELECT sq_y_1.t2_PLT_CN,
          sq_y_1.ADJ_FACTOR_MACR,
          tree.SPCD,
          cond.CONDID,
          cond.CONDPROP_UNADJ
     FROM ({sq_y_1}) AS sq_y_1
          INNER JOIN COND AS cond
          ON cond.PLT_CN = sq_y_1.t2_PLT_CN
          INNER JOIN COND AS prev_cond
          ON prev_cond.PLT_CN = sq_y_1.t1_PLT_CN
          INNER JOIN TREE AS tree
          ON tree.PLT_CN = sq_y_1.t2_PLT_CN
             AND tree.CONDID = cond.CONDID
             AND tree.PREVCOND = prev_cond.CONDID
          LEFT JOIN TREE_GRM_COMPONENT AS grm
          ON grm.TRE_CN = tree.CN
    WHERE tree.SPCD IN ({spp_cd*})
          -- only areas that were acc. forest land at t1 and t2
          AND cond.COND_STATUS_CD = 1
          AND prev_cond.COND_STATUS_CD = 1
          -- only areas where there were tally-trees of that species alive at
          -- t1, but excluding trees marked as NOT USED components for GRM
          AND grm.DIA_END >= 5
          AND tree.PREV_STATUS_CD = 1
          AND grm.SUBPTYP_BEGIN IN (1,3)
          AND grm.SUBP_SUBPTYP_GRM_AL_FOREST <> 0
 GROUP BY sq_y_1.t2_PLT_CN,
          tree.SPCD,
          cond.CONDID",
  .con = con) 

# - for this subset of [plot x cond]s, calc prop_for_t1 by species
sq_x_2 <- glue_sql(
  "SELECT t2_PLT_CN,
          SPCD,
          SUM(CONDPROP_UNADJ * ADJ_FACTOR_MACR) AS prop_for_t1
     FROM ({sq_x_2_1}) AS sq_x_2_1
 GROUP BY t2_PLT_CN,
          SPCD",
  .con = con) 

# - sq_x_3: forest land with spp at t1 or t2 (ha; for net-growth est) ---- ####
# get 'prop_for_t1t2'
# - first, return only [plots x cond] where that spp was live at t1 or ingrowth
#   at t2 (on land that was accessible forest land at both t1 and t2)
sq_x_3_1 <- glue_sql(
  "SELECT sq_y_1.t2_PLT_CN,
          sq_y_1.ADJ_FACTOR_MACR,
          tree.SPCD,
          cond.CONDID,
          cond.CONDPROP_UNADJ
     FROM ({sq_y_1}) AS sq_y_1
          INNER JOIN COND AS cond
          ON cond.PLT_CN = sq_y_1.t2_PLT_CN
          INNER JOIN COND AS prev_cond
          ON prev_cond.PLT_CN = sq_y_1.t1_PLT_CN
          INNER JOIN TREE AS tree
          ON tree.PLT_CN = sq_y_1.t2_PLT_CN
             AND tree.CONDID = cond.CONDID
             AND tree.PREVCOND = prev_cond.CONDID
          LEFT JOIN TREE_GRM_COMPONENT AS grm
          ON grm.TRE_CN = tree.CN
    WHERE tree.SPCD IN ({spp_cd*})
          -- only areas that were acc. forest land at t1 and t2
          AND cond.COND_STATUS_CD = 1
          AND prev_cond.COND_STATUS_CD = 1
          -- only areas where there were tally-trees of that species alive
          --   at t1 or had t2-ingrowth, but excluding trees marked as NOT USED
          --   components for GRM
          AND (grm.DIA_END >= 5 OR grm.DIA_MIDPT >= 5)
          AND (tree.PREV_STATUS_CD = 1 OR tree.PREV_STATUS_CD IS NULL)
          AND grm.SUBP_SUBPTYP_GRM_AL_FOREST <> 0
 GROUP BY sq_y_1.t2_PLT_CN,
          tree.SPCD,
          cond.CONDID",
  .con = con) 

# - for this subset of [plot x cond]s, calc prop_for_t1t2 by species
sq_x_3 <- glue_sql(
  "SELECT t2_PLT_CN,
          SPCD,
          SUM(CONDPROP_UNADJ * ADJ_FACTOR_MACR) AS prop_for_t1t2
     FROM ({sq_x_3_1}) AS sq_x_3_1
 GROUP BY t2_PLT_CN,
          SPCD",
  .con = con) 

# - sq_x: assemble full query -------------------------------------------- ####
# pull into one sub-query
# - note sq_x_2 is a subset of sq_x_3 which is a subset of sq_x_1
sq_x <- glue_sql(
  "SELECT sq_y_1.EVALID,
          sq_y_1.ESTN_UNIT,
          sq_y_1.STRATUMCD,
          sq_y_1.EXPNS,
          sq_y_1.P2POINTCNT,
          sq_y_1.t2_PLT_CN,
          sq_y_1.t1_PLT_CN,
          sq_x_3.SPCD,
          sq_x_1.prop_for,
          sq_x_3.prop_for_t1t2,
          sq_x_2.prop_for_t1
     FROM ({sq_x_1}) AS sq_x_1
          LEFT JOIN ({sq_y_1}) AS sq_y_1
          ON sq_y_1.t2_PLT_CN = sq_x_1.t2_PLT_CN
          LEFT JOIN ({sq_x_3}) AS sq_x_3
          ON sq_x_1.t2_PLT_CN = sq_x_3.t2_PLT_CN
          LEFT JOIN ({sq_x_2}) AS sq_x_2
          ON sq_x_2.t2_PLT_CN = sq_x_3.t2_PLT_CN
             AND sq_x_2.SPCD = sq_x_3.SPCD",
  .con = con)

# == Query dbs: ROM_y ==================================================== ####
# - query dbs, peak at data ---------------------------------------------- ####
# query SQLite dbs from Datamart (downloaded locally), returning one df with
# results from the queried state dbs:
ROM_y <- fn.query_states(query_args$state, sq_y)

# peak at the data
ROM_y %>% summary()  # summary stats for each variable
ROM_y %>% nrow()     # total plots query returned
ROM_y %>%            # plots by state
  group_by(state) %>%
  summarise(n_plots = n())

# - check: all high values / oddities make sense & are OK ---------------- ####
# grab results from query for tree records
df_sq_y <- fn.query_states(query_args$state, sq_y)

# few vars from ROM_y
cd_short <- ROM_y %>% select(t2_PLT_CN, SPCD, contains(c("BAH", "TPH")))

# netgrowth_annBAH extremes
# - looks ok, tones of ingrowth, maybe a tree came down and killed a few of
#   the trees on the plot (AGENTCD 60) and opened things up for the ingrowth
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(netgrowth_annBAH == max(netgrowth_annBAH)) %>%
                           pull(t2_PLT_CN)))
# - ok, entire plot killed by fire
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(netgrowth_annBAH == min(netgrowth_annBAH)) %>%
                           pull(t2_PLT_CN)))
# mortality_annBAH extremes
# - ok, entire plot killed by fire
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(mortality_annBAH == max(mortality_annBAH)) %>%
                           pull(t2_PLT_CN)))
# netgrowth_annTPH extremes
# - looks ok, tones of ingrowth, with the rest of the plot cut1 or cut2, so
#   likely harvest and opened things up for the ingrowth
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(netgrowth_annTPH == max(netgrowth_annTPH)) %>%
                           pull(t2_PLT_CN)))
# - ok, most of plot killed by fire, some survivors and ingrowth, and one diseased
#   tree
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(netgrowth_annTPH == min(netgrowth_annTPH)) %>%
                           pull(t2_PLT_CN)))
# mortality_annTPH extremes
# - returns the same plot as above, lots of fire kill
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(mortality_annTPH == max(mortality_annTPH)) %>%
                           pull(t2_PLT_CN)))
# t1_live_BAH & t2_live_BAH extremes
# - many survivors looks ok, one big tree 2.2 m^2 BA at t1/t2 (t2 DIA 43.0)
#   Same plot pops up for t1 & t2
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t1_live_BAH == max(t1_live_BAH)) %>%
                           pull(t2_PLT_CN)))
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t2_live_BAH == max(t2_live_BAH)) %>%
                           pull(t2_PLT_CN)))
# - Okay: ingrowth of 818 where there was none before
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t1_live_BAH == min(t1_live_BAH)) %>%
                           pull(t2_PLT_CN))[1])
# - Okay: all trees died from fire
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t2_live_BAH == min(t2_live_BAH)) %>%
                           pull(t2_PLT_CN))[1])
# t1_live_TPH & t2_live_TPH extremes
# - many survivors looks ok
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t1_live_TPH == max(t1_live_TPH)) %>%
                           pull(t2_PLT_CN)))
# - many survivors and a bunch of ingrowth looks ok
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t2_live_TPH == max(t2_live_TPH)) %>%
                           pull(t2_PLT_CN)))
# - Okay: ingrowth of 818 where there was none before
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t1_live_TPH == min(t1_live_TPH)) %>%
                           pull(t2_PLT_CN))[1])
# - Okay: all trees died from fire
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t2_live_TPH == min(t2_live_TPH)) %>%
                           pull(t2_PLT_CN))[1])

# t1 & t2 TPH 0 (same plot comes up for t1 & t2 BAH 0)
# - one tree and it was CUT2, so ! all okay but strange
df_sq_y %>%
  filter(t2_PLT_CN %in% (ROM_y %>%
                           filter(t1_live_TPH == 0 & t2_live_TPH == 0) %>%
                           pull(t2_PLT_CN)))

# - format data: split by spp -------------------------------------------- ####
ROM_y <- ROM_y %>%
  mutate(spp = case_when(SPCD == 22  ~ "ABPR",
                         SPCD == 202 ~ "PSMEM",
                         SPCD == 807 ~ "QUDO",
                         SPCD == 815 ~ "QUGA4",
                         SPCD == 818 ~ "QUKE")) %>%
  split(.$spp)

# - format data: add our plot_id & limit to sp-specific sample ----------- #### 
  
# step 1: get key between plot_id and FIA-internal PLT_CN (not datamart)
# - also add true hex membership...
load(paste0(path_ch2_data, "ch2-1-all-vars.Rdata"))
plot_crosswalk <- FIA %>%
  select(plot_id, PLT_CN, fiahex_id)
rm(sppcd, spp_names, ES_key, FIA, FIA_spp, plot_key) # keep only plot_key

# step 2: add public (datamart) PLT_CNs
# - inner_join with key between confidential and public PLT_CNs (this created 
#   in Ch2's script 8-...R)
plot_crosswalk <- readRDS(
  paste0(path_USB, 
         "Karin phd thesis files/r - raw master data/",
         "crosswalk plots - true and public/", 
         "crosswalk_true_public_PLT_CN.rds")) %>% 
  inner_join(plot_crosswalk, by = "PLT_CN") %>%
  rename(t1_PLT_CN = PLT_CN)

# - limit to plots within Ch2 sample for each species & add p_a used in models
#   (based on 1" dbh instead of 5" like change data, so change name). 
#   Note, here using elev/lat created groups, but cna and maca PECA-based
#   groups have same # of plots for each spp)
plot_crosswalk <- readRDS(paste0(path_output, "FF-1-preds-PRISM_100MCMC.rds")) %>% 
  map(function(sp.x) {
    with(sp.x$newdat.x, data.frame(plot_id, p_a)) %>% 
      rename(t1_1in_p_a = p_a) %>%
      inner_join(plot_crosswalk, by = "plot_id") %>% 
      select(plot_id, t1_PLT_CN, fiahex_id)
  })

# step 3: make PLT_CN into factor & add plot_id
ROM_y <- list(ROM_y, plot_crosswalk) %>%
  pmap(function(sp.x, key.x) {
    # make PLT_CN's into factors
    sp.x$t1_PLT_CN <- as.factor(sp.x$t1_PLT_CN)
    sp.x$t2_PLT_CN <- as.factor(sp.x$t2_PLT_CN)
    
    # inner_join by plots in each species' sample
    sp.x %>% inner_join(key.x , by = "t1_PLT_CN")
  })

# - format data: calc other vars of interest ----------------------------- ####
# About:
# - conversion: determine if tally-trees for that spp on forested land went
#   from P (t1) to A (t2), AtoP, or remained present. Note, b/c of how data 
#   were queried we'll never see plots where the species went from AtoA and 
#   every plot in the change data set where the spp was absent at t1 will be 
#   present at t2.
# - idea on hold: relativized values
#   this would require a ratio-of-means type estimation, but we don't have
#   time before defense to implement this. so keep idea on the back-burner
#   and (ATM) comment out in the code below...
#   - concept: what do net-growth and mortality look like when annualized 
#     TPH/BAH variables are relativized by the TPH/BAH at t1; idea here 
#     being to gauge GM relative to what was originally there. For example,  
#     a loss of 7 TPH where there was 7 TPH at t1 is more drastic than for a 
#     plot where there was 70 TPH at t1.
ROM_y <- ROM_y %>%
  map(function(sp.x) {
    sp.x %>%
      # calculate conversion metrics
      mutate(conv_5in = case_when((t1_live_TPH > 0 & t2_live_TPH > 0) ~ "PtoP",
                                  (t1_live_TPH > 0 & t2_live_TPH == 0) ~ "PtoA",
                                  (t1_live_TPH == 0 & t2_live_TPH > 0) ~ "AtoP")) %>%
      # order cols
      select(SPCD, spp, 
             plot_id, fiahex_id, t1_PLT_CN, t2_PLT_CN,
             state, EVALID, ESTN_UNIT, STRATUMCD,
             EXPNS, P2POINTCNT,
             t1_MEASYEAR, t2_MEASYEAR, REMPER,
             everything())
  })

# == Query dbs: wts info & state-membership ============================== ####
# - Define population: get states that overlap sp-specific study area ---- ####
# here using plot_crosswalk, that is, the public PLT_CN's we pulled around
# line ~590 above for ROM_y...
state_mem <- plot_crosswalk %>% 
  bind_rows %>% 
  distinct(t1_PLT_CN) %>%
  pull(t1_PLT_CN) %>%
  as.numeric

state_mem <- glue_sql(
  "SELECT STATECD, 
        PREV_PLT_CN AS t1_PLT_CN
   FROM PLOT
  WHERE PREV_PLT_CN IN ({state_mem*})") %>%
  fn.query_states(query_args$state, .)

# break it apart by species
state_mem <- plot_crosswalk %>%
  map(~state_mem %>% filter(t1_PLT_CN %in% .x$t1_PLT_CN))

# summarize
state_mem %>%
  map(~.x %>% group_by(STATECD) %>% summarize(n_plots = n())) %>%
  bind_rows(.id = "spp") %>%
  pivot_wider(names_from = spp,
              values_from = n_plots)
# states that overlap sp-specific study areas
state_mem %>%
  map(~.x %>% group_by(state) %>% summarize(n_plots = n())) %>%
  bind_rows(.id = "spp") %>%
  pivot_wider(names_from = spp,
              values_from = n_plots)

# - Get weight info for each {evalid x estn_unit} ------------------------ ####
# To weight each observation, we'll use:
# - [(ac of stratum h) / (# of plots in stratum h)] / (ac of population)
#   where the population is the states that overlap the spp-specific study
#   area we defined. 

ROM_wh <- glue_sql( 
  'SELECT peu.EVALID, 
          peu.ESTN_UNIT, 
          peu.AREA_USED AS eu_area_used,
          SUM(psm.P2POINTCNT) AS eu_n -- will aggr this up to pop later
     FROM POP_ESTN_UNIT AS peu
          INNER JOIN POP_STRATUM AS psm
          ON peu.EVALID = psm.EVALID
             AND peu.ESTN_UNIT = psm.ESTN_UNIT
    WHERE peu.EVALID IN ({query_args$evalid*})
 GROUP BY peu.EVALID,
          peu.ESTN_UNIT',
  .con = con) %>%
  fn.query_states(query_args$state, .)

# break ROM_wh apart by species s.t. WA is not in the df for qudo
ROM_wh <- state_mem %>% map(function(sp.x) {
  x1 <- ROM_wh %>% 
    filter(state %in% sp.x$state) %>%
    group_by(state) %>%
    mutate(state_area_used = sum(eu_area_used),
           state_n = sum(eu_n)) %>%
    ungroup
  
  x1 %>%
    distinct(state, .keep_all = T) %>%
    mutate(pop_area_used = sum(state_area_used),
           pop_n = sum(state_n)) %>%
    select(state, pop_area_used, pop_n) %>%
    inner_join(x1, by = "state")
})

# == Query dbs: denominator ============================================== ####
# - query dbs, peak at data ---------------------------------------------- ####
# query SQLite dbs from Datamart (downloaded locally), returning one df with
# results from the queried state dbs:
ROM_x <- fn.query_states(query_args$state, sq_x)

# peak at the data
ROM_x %>% summary()  # summary stats for each variable
ROM_x %>% nrow()     # total plots query returned
ROM_x %>%            # plots by state
  group_by(state, SPCD) %>%
  summarise(n_plots = n())

# - format data: split by spp -------------------------------------------- ####
# a bit different from ROM_y here, b/c we have prop_for that we need to keep
# for all plots within the sp-specific population, this might include plots
# on which the spp was absent at t1/t2.

# approach: replicate into a list of 5 (1 df per spp), then slim things down
# from there...
ROM_x <- lapply(1:5, function(x) ROM_x) %>%
  setNames(c("ABPR", "PSMEM", "QUDO", "QUGA4", "QUKE")) %>%
  imap(function(sp.x, sp_name.x) {
    sp.x %>%
      mutate(spp = case_when(SPCD == 22  ~ "ABPR",
                             SPCD == 202 ~ "PSMEM",
                             SPCD == 807 ~ "QUDO",
                             SPCD == 815 ~ "QUGA4",
                             SPCD == 818 ~ "QUKE",
                             is.na(SPCD) ~ "absent")) %>%
      # convert to zero values of prop_for_t1t2 for other spp (or NA's if now 
      # labeled 'absent' above b/c none of our spp were on that plot)
      mutate(prop_for_t1t2 = if_else(spp == sp_name.x, prop_for_t1t2, 0)) %>%
      # convert to zero values of prop_for_t1 for other spp (or NA's if now 
      # labeled 'absent' above b/c none of our spp were on that plot)
      # - extra step here b/c can be coded as that spp for prop_for_t1t2, but
      #   the spp might've not been present at t1 (only t2, and so will be NA
      #   for prop_for_t1)
      mutate(prop_for_t1 = if_else(spp == sp_name.x, prop_for_t1, 0),
             prop_for_t1 = if_else(is.na(prop_for_t1), 0, prop_for_t1)) %>%
      select(-SPCD, -spp)
  })

# - format data: limit PSMEM & limit sp-states & add plot_id ------------- ####
# make plt_cn's factors and add our plot_id if it the plot was in our spp
# list & limit to states that overlap the spp study areas ('population')
ROM_x <- list(ROM_x, 
              plot_crosswalk %>% setNames(names(ROM_x)),
              state_mem %>% setNames(names(ROM_x))) %>%
  pmap(function(x1, x2, x3) {
    # make PLT_CN's into factors
    x1$t1_PLT_CN <- as.factor(x1$t1_PLT_CN)
    x1$t2_PLT_CN <- as.factor(x1$t2_PLT_CN)
    
    # add plot_id if it was in our spp list
    x1 %>% 
      left_join(x2 %>% select(t1_PLT_CN, plot_id), by = "t1_PLT_CN") %>%
      filter(state %in% unique(x3$state))
  })

# because FIA doesn't identify PSME to the variety level, we'll use our 
# study area to define whether or not a plot contained PSMEM (and could 
# therefore have a value other than 0 for prop_for_t1t2 or prop_for_t1)
ROM_x <- ROM_x %>%
  imap(function(sp.x, sp_name.x) {
    # adjust PSMEM if on that spp's list element & for all spp add T/F
    # toggle of if the spp was P at t1 or P at either t1 or t2
    psmem_sa <- plot_crosswalk$psmem$t1_PLT_CN
    sp.x %>%
      mutate(prop_for_t1t2 = 
               if_else((sp_name.x == "PSMEM") & (!t1_PLT_CN %in% psmem_sa),
                       0, prop_for_t1t2),
             prop_for_t1   = 
               if_else((sp_name.x == "PSMEM") & (!t1_PLT_CN %in% psmem_sa),
                       0, prop_for_t1)) %>%
      mutate(spp_p_t1t2 = prop_for_t1t2 > 0,
             spp_p_t1   = prop_for_t1 > 0)
  })

# == Save data for future scripts ======================================== ####
saveRDS(ROM_y, paste0(path_data, "FF-2-rom-y-data.rds"))
saveRDS(ROM_x, paste0(path_data, "FF-2-rom-x-data.rds"))
saveRDS(ROM_wh, paste0(path_data, "FF-2-rom-wh-data.rds"))
saveRDS(state_mem, paste0(path_data, "FF-2-state-membership.rds"))
