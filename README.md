# COIL+

This repo contains the code needed to reproduce the results from the main text of the manuscript and extends the model and sampler developed in Papadogeorgou (2023) to account for the effects of extreme taxonomic bias. In this project, we aim to impute missing mammal-plant interactions from a multi-study data set of Afrotropical frugivory. The data are provided with generic species labels for plants and frugivores in compliance with data sharing restrictions.

## Data set

The key raw data sources are are in the folder RawData:

-   Raw frugivory interactions with generic species labels and the variables study ID, study site, country, zone, focus: frug_generic.csv

-   Plant traits with generic species labels: Obs_W.data

-   Vertebrate traits with generic species labels: Obs_X.data

-   Generic taxonomy for plants: p_taxa_generic.csv

-   Generic taxonomy for vertebrates: v_taxa_generic.csv

-   Cu_phylo.dat: phylogenetic correlation matrix for vertebrates

-   Cv_phylo.dat: phylogenetic correlation matrix for plants

Plant phylogenies is acquired using the V.PhyloMaker R package and the phylogenetic correlation matrix for vertebrates is computed via the the package ape using a consensus tree obtained from VertLife and provided in this repo. For the purposes of this vignette, we supply the respective phylogenetic correlation matrices directly without

## Code

The folder R/ includes functions with substantial modifications and additions from the [BiasedNetwork](https://github.com/gpapadog/BiasedNetwork) repo that are used in the analysis code.

The code for the analysis is in the folder Analysis/. The numbers in the beginning of the file names represent the order with which the files should be used/run. In brief the content of each analysis file is as follows:

-   Anaylsis_0_data_prep.R: This code MUST be run before any subsequent analysis. Assembles key network and meta-data sources.

-   Analysis_1a_fitCOILplus.R: This code MUST be run before any subsequent plotting functions. Fits COIL+ to the Afrotropical frugivory data with prior incorporating domain knowledge.

-   Analysis_1b_fitCOIL.R: This code MUST be run before any subsequent plotting functions. Fits basic COIL to the Afrotropical frugivory data with prior incorporating domain knowledge.

-   Analysis_1c_fitCOIL_default.R: This code MUST be run before any subsequent plotting functions. Fits basic COIL to the Afrotropical frugivory data with default prior.

-   Analysis_2a_cv_COILplus.R: This code MUST be run before any subsequent plotting functions. Performs 10-fold cross validation of COIL+ with the Afrotropical frugivory data with prior incorporating domain knowledge.

-   Analysis_2b_cv_COIL.R: This code MUST be run before any subsequent plotting functions. Performs 10-fold cross validation of COIL with the Afrotropical frugivory data with prior incorporating domain knowledge.

-   Analysis_2c_cv_COIL_default.R: This code MUST be run before any subsequent plotting functions. Performs 10-fold cross validation of COIL with the Afrotropical frugivory data with default 0/1 prior occurrence prior.

-   Analysis_3a_plot_results.R: This code analyzes in and out of sample model performance for the any of the three models fit in 1a-2c.

-   Analysis_4a_trait_matching_save.R: This code MUST be run before the subsequent plotting code. Performs trait matching for the selected model fit.

-   Analysis_4b_trait_matching_plot.R: This code plots the trait matching.

### Note

For replicating the analysis results in this code, you will need to specify a directory where processed data and results can be saved. We recommend you create a folder at the same level as the Analysis/, Data/, and R/ folders that is named Results/ and ProcessedData/.
