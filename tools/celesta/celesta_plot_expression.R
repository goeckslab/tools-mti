# --------------------------------------------------------------------------------------------
# Plot marker expression probabilities for cell assignment parameter optimization with CELESTA
# --------------------------------------------------------------------------------------------

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
  make_option(c("--filter"), action="store_true", type='logical', default=FALSE,
              help="Boolean to filter cells or not (default: no filtering)"),
  make_option(c("--highfilter"), action="store", default=0.9, type='double',
              help="High marker threshold if filtering cells (default: 0.9)"),
  make_option(c("--lowfilter"), action="store", default=0.4, type='double',
              help="Low marker threshold if filtering cells (default: 0.4)"),
  make_option(c("-s", "--size"), action="store", default=1, type='double',
              help="Point size for plotting"),
  make_option(c("--width"), action="store", default=5, type='integer',
              help="Width of plot (inches)"),
  make_option(c("--height"), action="store", default=4, type='integer',
              help="Height of plot (inches)")
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
  print("Proceeding to marker expression plotting without cell filtering")
}

# create output directory if it doesn't already exist 
dir.create('marker_exp_plots')

# plot expression probability
PlotExpProb(coords=CelestaObj@coords,
            marker_exp_prob=CelestaObj@marker_exp_prob,
            prior_marker_info = CelestaObj@prior_info,
            save_plot = TRUE,
            output_dir = "./marker_exp_plots")