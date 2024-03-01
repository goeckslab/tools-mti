#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------------
# A collection of functions that help convert between MCMICRO/sc-verse 
# file formats and those expected and handleable by R and CELESTA
# ---------------------------------------------------------------------------------

library(janitor)
library(dplyr)
library(anndata)

clean_marker_names <- function(prior_info) {

    return(prior_info_cleaned)
}

clean_cell_type_names <- function(prior_info) {

    return(prior_info_cleaned)
}

anndata_to_celesta <- function(input_adata, rep, x_centroid, y_centroid) { 

    # initialize output as dataframe from adata.obs
    celesta_input_dataframe <- data.frame(input_adata$obs)

    # rename X,Y column names to what CELESTA wants 
    names(celesta_input_dataframe)[names(celesta_input_dataframe) == x_centroid] <- 'X'
    names(celesta_input_dataframe)[names(celesta_input_dataframe) == y_centroid] <- 'Y'

    # subset to just "X" and "Y" columns
    celesta_input_dataframe <- celesta_input_dataframe %>%
        dplyr::select(X,Y)

    # merge X,Y coords with marker intensities from the specified anndata layer


    return(celesta_input_dataframe)
}

celesta_to_anndata <- function() {

    return(output_anndata)
}



