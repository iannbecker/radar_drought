# Data and code for: Spatiotemporal impacts of drought on migratory stopover across a heterogeneous landscape
**Citation:** ####

## Scripts
### Analysis
These scripts are for the main analysis described in the manuscript.

| Script | Description |
|--------|-------------|
| `random_forest.r` | Utilizes random forest to identify important predictors for subsequent model input. This script can be utilized for both seasons. |
| `drought_GAMS.r` | Runs 9 competing GAMs for each migratory season. Models were built using variable importance rankings from random forest |

### Figures
These scripts are for figures in the main body of the manuscript as well as the supplemental information indicated as either Figure# (main body) or Figure S# (supplemental).

| Script | Description |
|--------|-------------|
| `Figure1_concept_plots.r` | Simulated figures used in a flowchart depicting the effects of drought on wildlife |
| `Figure2_study_area.r` | Makes study area map with ecoregions also showing radar locations and detection radii |
| `Figure3_stopover_maps.r` | Creates Panels A/B in Figure 3; stopover maps showing seasonal stopover density by radar during exemplar seasons|
| `Figure3_drought_history.r` | Creates Panel C in Figure 3; drought history across our study area through our study period (1995-2020) |
| `Figure4_best_model.r` | Forest plots showing F-statistics of predictors in selected model for each season |
| `Figure5_drought_impact.r` | Makes all parts of Figure 5: bar graph comparing effect size of drought between fall and spring, smooth curves for each drought window |
| `Figure6_land_cover.r` | Faceted figure showing relationship between SPEI and stopover for each land cover class |
| `supp_smooth_curves.r` | Used to plot any specified smooth curve for any predictor |

## Abstract
The increased prevalence of extreme weather events poses a threat for native ecosystems worldwide. As patterns of precipitation and temperature are altered due to climate change, drought, in particular, may become a more prevalent and severe event. Drought directly affects resources that migratory animals depend on across the transient and stationary phases of their life cycle. To explore the impacts of drought on migration, we focus on migratory songbirds moving through Texas and Louisiana, an area that supports one of the greatest volumes of avian movement across the North American migratory corridor. We use 26 years of weather surveillance radar data (1995-2020) to quantify the effects of drought on migratory stopover density during fall and spring migration. We tested four competing hypotheses to understand the interactive effects of drought along with climate and habitat characteristics. Our results indicate that drought conditions had a stronger effect on fall migrants relative to spring migrants. In both migratory seasons, stopover density decreased with more intense current drought conditions. However, dry conditions in the immediate previous season (Summer; Winter) were associated with increasing fall and spring stopover density. Additionally, migrants experiencing drought conditions appeared to respond primarily to local habitat type rather than variation in ecoregion. Permanent wetlands provided stable refugia for migrants under drought conditions during both fall and spring, likely reflecting the importance of adequate resources within stopover habitat. 
Our findings highlight the nonstationary impacts of severe weather events on migratory birds and  help us better understand how migratory species make decisions when experiencing the effects of global change. 


