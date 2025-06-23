# COIL+

This repo contains the code needed to reproduce the results from the main text of the manuscript and extends the model and sampler developed in Papadogeorgou (2023) to account for the effects of extreme taxonomic bias. In this project, we aim to impute missing mammal-plant interactions from a multi-study data set of Afrotropical frugivory. The data are provided with generic species labels for plants and frugivores in compliance with data sharing restrictions.

## Data set

The key raw data sources are:

-   Raw frugivory interactions with generic species labels and the variables study ID, study site, country, zone, focus: frug_generic.csv

-   Plant traits with generic species labels: Obs_W.data

-   Vertebrate traits with generic species labels: Obs_X.data

-   Generic taxonomy for plants: p_taxa_generic.csv

-   Generic taxonomy for vertebrates: v_taxa_generic.csv

-   Cu_phylo.dat: phylogenetic correlation matrix for vertebrates

-   Cv_phylo.dat: phylogenetic correlation matrix for plants

Plant phylogenies is acquired using the V.PhyloMaker R package and the phylogenetic correlation matrix for vertebrates is computed via the the package ape using a consensus tree obtained from VertLife and provided in this repo. For the purposes of this vignette, we supply the respective phylogenetic correlation matrices directly without

## Code

The folder R/ includes functions (with substantial modifications and additions from <https://github.com/gpapadog/BiasedNetwork>) that are used in the analysis code while the folder HelperScriptsOld/ includes more minor modifications for occurrence estimation within the original model for comparison purposes.

The code for the analysis is in the folder Analysis/. The numbers in the beginning of the file names represent the order with which the files should be used/run. In brief the content of each analysis file is as follows:

-   Anaylsis_0_data_prep.R: This code MUST be run before any subsequent analysis. Assembles key network and meta-data sources.

-   Analysis_1a_fitCOILplus.R: This code MUST be run before any subsequent plotting functions. Fits COIL+ to the Afrotropical frugivory data with prior incorporating domain knowledge.

-   ... clean up below, add these:

-   Analysis_1b_fitCOIL.R: This code MUST be run before any subsequent plotting functions. Fits basic COIL to the Afrotropical frugivory data with prior incorporating domain knowledge.

-   Analysis_1c_fitCOIL_default.R: This code MUST be run before any subsequent plotting functions. Fits basic COIL to the Afrotropical frugivory data with default prior.

-   Analysis_2b_compare_results.R: This code compares the performance of the the models fit in Analysis_1a - Analysis_1c.

-   Analysis_3a_perform_trait.R: This code MUST be run before the subsequent plotting code. Performs trait matching for the selected results.

-   Analysis_3b_plot_trait.R: This code plots the trait matching.

-   1a_subset_data.R: This code MUST be run before any subsequent analysis. This file subsets data to include only mammals and plants.

-   1b_phylo_verts.R: This code MUST be run before any subsequent analysis. This file generates the vertebrate phylogenetic correlation matrix.

-   1c_phylo_plants.R: This code MUST be run before any subsequent analysis. This file generate the plant phylogenetic correlation matrix.

-   2a_subsetting.R: This code MUST be run before any subsequent analysis. This file further subsets the data to exclude isolated species.

-   2b_analysis_new_defaultprior.R: perform the analysis using default prior occurrence probabilities under the new sampler.

-   2c_analysis_new_expertprior.R: perform the analysis using expert defined prior occurrence probabilities under the new sampler

-   3a_cv_new_defaultprior.R: perform cross validation using default prior occurrence probabilities under the new sampler.

-   3b_cv_new_expertprior.R: perform cross validation using expert prior occurrence probabilities under the new sampler.

-   4a_plot_results_new.R: plot the posterior heatmap and cross validation performance metrics using the new sampler.

### Note

For replicating the analysis results in this code, you will need to specify a directory where processed data and results can be saved. We recommend you create a folder at the same level as the Analysis/, Data/, and R/ folders that is named Results/ and ProcessedData/.
