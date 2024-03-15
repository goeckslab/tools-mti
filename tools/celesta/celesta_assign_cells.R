# ---------------------------------------------------------------------------------
# The main algorithim for CELESTA cell type assignment 
# ---------------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(janitor)))
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(anndataR)))
suppressWarnings(suppressMessages(library(Rmixmod)))
suppressWarnings(suppressMessages(library(spdep)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(zeallot)))
suppressWarnings(suppressMessages(library(CELESTA)))

### Define command line arguments 
option_list = list(
  make_option(c("-i", "--imagingdata"), action="store", default=NA, type='character',
              help="Path to imaging data"),
  make_option(c("-p", "--prior"), action="store", default=NA, type='character',
              help="Path to prior marker info file"),
  make_option(c("-x", "--xcol"), action="store", default=NA, type='character',
            help="Name of column in adata.obs containing X coordinate"),
  make_option(c("-y", "--ycol"), action="store", default=NA, type='character',
            help="Name of column in adata.obs containing Y coordinate"),
  make_option(c("--filter"), action="store_true", type='logical',
              help="Boolean to filter cells or not (default: no filtering)"),
  make_option(c("--highfilter"), action="store", default=0.9, type='double',
              help="High marker threshold if filtering cells (default: 0.9)"),
  make_option(c("--lowfilter"), action="store", default=0.4, type='double',
              help="Low marker threshold if filtering cells (default: 0.4)"),
  make_option(c("--maxiteration"), action="store", default=10, type='integer',
              help="Maximum iterations allowed in the EM algorithm per round"),
  make_option(c("--changethresh"), action="store", default=0.01, type='double',
              help="Ending condition for the EM algorithm"),    
  make_option(c("--highexpthresh"), action="store", default='default_high_thresholds', type='character',
            help="Path to file specifying high expression thresholds for anchor and index cells"),
  make_option(c("--lowexpthresh"), action="store", default='default_low_thresholds', type='character',
            help="Path to file specifying low expression thresholds for anchor and index cells")
)

### Functions 
anndata_to_celesta <- function(input_adata, x_col, y_col) { 

  #' Function to convert anndata object to dataframe readable by CELESTA
  #' Coordinates columns in adata.obs are renamed to 'X' and 'Y', and placed at beginning of dataframe 
  #' Marker intensities are concatted columnwise to the X and Y coords. cols: X,Y,Marker_1,...Marker N

  # initialize output as dataframe from adata.obs
  celesta_input_dataframe <- data.frame(input_adata$obs)

  # subset to X and Y coordinates from obs only 
  celesta_input_dataframe <- celesta_input_dataframe %>% 
      dplyr::select({{x_col}},{{y_col}})

  # rename X,Y column names to what CELESTA wants 
  colnames(celesta_input_dataframe) <- c('X','Y')

  # merge X,Y coords with marker intensities from adata.X
  x_df <- data.frame(input_adata$X)
  colnames(x_df) <- input_adata$var_names
  celesta_input_dataframe <- cbind(celesta_input_dataframe,x_df)

  return(celesta_input_dataframe)
}

### Main 
# parse args 
opt = parse_args(OptionParser(option_list=option_list))

# read anndata, convert to celesta format 
adata <- read_h5ad(opt$imagingdata)
celesta_input_df <- anndata_to_celesta(adata, x_col = opt$xcol, y_col = opt$ycol)

# read prior marker info 
prior <- read.csv(opt$prior, check.names = FALSE)

# clean prior names, keeping a copy of originals for writing output 
prior_original_names <- colnames(prior)
prior <- janitor::clean_names(prior, case = "all_caps")

# clean input dataframe names, keeping a copy of originals for writing output 
celesta_input_df_original_names <- colnames(celesta_input_df)
celesta_input_df <- janitor::clean_names(celesta_input_df, case = "all_caps")

# instantiate celesta object 
CelestaObj <- CreateCelestaObject(
  project_title = "",
  prior_marker_info = prior,
  imaging_data_file = celesta_input_df
)

# if filtering is specified, filter out cells outside high and low thresholds
if (opt$filter) {
  print("filtering cells based on expression")
    CelestaObj <- FilterCells(
        CelestaObj, 
        high_marker_threshold = opt$highfilter, 
        low_marker_threshold = opt$lowfilter)
} else {
  print("Proceeding to cell type assignment without cell filtering")
}

# check for non-default expression threshold files 
if (opt$highexpthresh != 'default_high_thresholds') {
  # read high thresholds 
  high_expression_thresholds <- read.csv(opt$highexpthresh)
  high_expression_threshold_anchor <- high_expression_thresholds$anchor
  high_expression_threshold_index <- high_expression_thresholds$index
} else {
  print("Using default high expression thresholds -- this may need adjustment") 
  high_expression_threshold_anchor <- rep(0.7,length = 50)
  high_expression_threshold_index <- rep(0.5,length = 50)
}

if (opt$lowexpthresh != 'default_low_thresholds') {
  # read low thresholds 
  low_expression_thresholds <- read.csv(opt$highexpthresh)
  low_expression_threshold_anchor <- low_expression_thresholds$anchor
  low_expression_threshold_index <- low_expression_thresholds$index
} else {
  print("Using default low expression thresholds") 
  low_expression_threshold_anchor <- rep(0.9,length = 50)
  low_expression_threshold_index <- rep(1,length = 50)
}

# run cell type assignment 
CelestaObj <- AssignCells(CelestaObj,
  max_iteration = opt$maxiteration,
  cell_change_threshold = opt$changethresh,
  high_expression_threshold_anchor = high_expression_threshold_anchor,
  low_expression_threshold_anchor = low_expression_threshold_anchor,
  high_expression_threshold_index = high_expression_threshold_index,
  low_expression_threshold_index = low_expression_threshold_index,
  save_result = FALSE)

# save object as an RDS file for cell type plotting 
# for the time being, this is not exposed to Galaxy 
saveRDS(CelestaObj, file = "celestaobj.rds")

# rename celesta assignment columns so they're obvious in output anndata 
celesta_assignments <- CelestaObj@final_cell_type_assignment
celesta_assignments <- janitor::clean_names(celesta_assignments)
colnames(celesta_assignments) <- paste0("celesta_", colnames(celesta_assignments))

# merge celesta assignments into anndata object
adata$obs <- cbind(adata$obs,celesta_assignments)

# write output anndata file 
write_h5ad(adata, 'result.h5ad')