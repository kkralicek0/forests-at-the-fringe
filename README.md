# forests-at-the-fringe
These R scripts are in support of 'Forests at the Fringe: comparing observed change to projected climate change impacts for five species in the Pacific Northwest' by Karin Kralicek, Tara Barrett, Jay Ver Hoef, and Temesgen Hailemariam. There are 5 scripts in this project, which begin with a number indicating their order in the workflow. Some objects in this workflow are dependent on results (i.e., R-objects) created in the second research chapter of K. Kralicek's PhD dissertation (Kralicek, 2022), the scripts for which can be found online at https://github.com/kkralicek0/car-sdm. This repository is archived and read-only.


Definition of scripts:

* 1-subset-mcmc.R: takes the predictions from Kralicek (2022) for the current and four future climate scenarios and subsets to 100 of 2000 MCMC samples (every 20th of 2000 to save computing time/space for this project), calculates range estimates based on bivariate kernel density estimates, and identifies if plots are within areas of projected expansion/contraction/persistence.

* 2-query-FIAdbs-for-obsv-change.R: queries the FIAdb SQLite databases (locally downloaded, but publicly available through the FIA Datamart at https://apps.fs.usda.gov/fia/datamart/datamart.html) for the states of Oregon, California, and Washington, for net-growth, mortality, and relative-growth in terms of BAH/yr and TPH/yr for the five species. Actual estimates are calculated within the next script, however this script pulls plot-level estimates as well as post-stratified weights.

* 3-get-domain-est-for-obsv-change.R: calculates estimates and SE for net-growth, mortality, and relative-growth in terms of BAH/yr and TPH/yr within 1) the species range overall, 2) naive-divisions of habitat based on elevation and latitude, and 3) areas of projected expansion/contraction/persistence. Naive-divisions are created in this script based on FIA plot data.

* 4-results-projected-impacts.R: summarizes and visualizes projected climate change impacts to provide context for the discussion of observe change within these areas. Produces table estimates for manuscript.

* 5-results-observed-change.R: summarizes and visualizes observed change estimates within naive-divisions of habitat and PECA-divisions (i.e., persistence/expansion/contraction). Produces figures/tables for manuscript and appendix.


Note on data: 

The raw data supporting the conclusions of this article are based on confidential, precise-location coordinates for FIA plots. Fuzzed (imprecise) versions of these raw data will be made available by the corresponding author upon request, without undue reservation. For access to precise-location information see Burrill et al. (2021).


Citations:

Burrill, E. A., DiTommaso, A. M., Turner, J. A., Pugh, S. A., Menlove, J., Christiansen, G., et al. (2021). The Forest Inventory and Analysis Database: database description and user guide version 9.0.1 for Phase 2. U.S. Department of Agriculture, Forest Service. Available online at: http://www.fia.fs.fed.us/library/database-documentation/

Kralicek., K. (2022). “Climate change induced shifts in suitable habitat for five tree species in the Pacific Northwest projected with spatial-Bayesian hierarchical models,” in Characterizing Uncertainty and Assessing the Impact of Rapid Climate Change on the Distribution of Important Tree Species in the Pacific Northwest. [dissertation]. [Corvallis (OR)]: Oregon State University. 36–74. Available online at: https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/4f16cb105
