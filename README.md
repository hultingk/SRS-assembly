# SRS plant community assembly

Contact: Katherine Hulting, hultingk@msu.edu

### Scripts
#### 00_functions.R
Script containing custom functions used throughout analysis.

#### 01_wrangling.R
Script for data cleaning and wrangling.
* Requires: access to "L0_original" folder containing data (not committed). 
* Products: srs_plant_all.csv (cleaned plant community data)

#### 02_pcoa.R
Script for visualizing PCoA of plant community trajectories 
* Requires: 01_wrangling.R 
* Products: figureS2.pdf 

#### 03_convergence.R
Script for calculating convergence/divergence between patch type communities across time, repeat within dispersal mode groups
* Requires: 00_functions.R, 02_pcoa.R
* Products: figure4.pdf, tableS10.html, tableS11.html, tableS12.html

#### 04_cta_segments.R
Script for calculating interannual trajectory distances, repeat within dispersal mode groups
* Requires: 00_functions.R, 01_wrangling.R
* Products: tableS1.html, tableS2.html, tableS3.html, figure2.pdf

#### 05_cta_directionality.R
Script for calculating trajectory directionality, repeat within dispersal mode groups
* Requires: 04_cta_segments.R
* Products: tableS7.html, tableS8.html, table S9.html, figure3.pdf

#### 06_supplement.R
Script for supplemental exploratory plots: species richness over time, proportion of dispersal mode, turnover
* Requires: 00_functions.R, 01_wrangling.R
* Products: figureS1.pdf, figureS3.pdf, figureS4.pdf

#### 06b-supplement-subset.R
Script for supplemental repeating analyses with only blocks that were active for the full time series
* Requires: 00_functions.R, 01_wrangling.R
* Products: figureS6.pdf

* #### 06c-supplement-year.R
Script for supplemental repeating analyses using calender year as predictor
* Requires: 00_functions.R, 01_wrangling.R
* Products: figureS5.pdf




