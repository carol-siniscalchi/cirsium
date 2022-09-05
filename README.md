# Cirsium
Scripts and data from "Diversification and biogeography of North American thistles (Cirsium: Carduoideae: Compositae): drivers of a rapid continent-wide radiation", in review at the International Journal of Plant Sciences.

## Directory structure

### biogeography
This directory contains the callibrated tree used in biogeographical analysis (`new.cirsium.tre`), both geographic codings used (by geographical region [`geo.areas.data`] and by biome [`geo.biome.data`]), the correspondent R scripts used to run three models for each geographical coding (`bgb_allmodels_areas.R` and `bgb_allmodels_biomes.R`) and tables with the AICc values obtained.

### environmental
This directory contains a script (`environmental.R`) that runs linear and PGLS models for all environmental variables. It also contains a file with the average values (`cirsium_environmental_allcombined_DR.csv`), the tree needed for the analysis (`cirsium.phyluce95p.ultrametric.tre`) and a file with the tip diversification rates (`tip_DR_cirsium.csv`).

### traits
This directory contains three scripts: `cirsium_pollination.R` runs all the pollination vs. environmental analysis, including plotting violin plots of the variable distribution per pollinator; `cirsium.trait.rec.R` runs two models (ER and ARD) of character state reconstruction for selected morphological traits, builds an AIC table, plots the reconstructions on the tree, and calculates D statistics for binary traits; `FiSSE.R` runs state dependent models for binary characters and calculates the p-value of each model. There are five files needed to run these scripts: the tree (`cirsium.phyluce95p.ultrametric.tre`), a table with all traits (`traits.csv`), a table with only binary traits (`cirsium_binary.csv`), a table with pollinator states (`poll.csv`) and a table containing average values of environmental variables and diversification rates (`cirsium_environmental_allcombined_DR_pollinator.csv`). 
