library(optparse)
library(CELESTA)
library(Rmixmod)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)

# probably should move main args collection to a separate script. 
# same with filtering and object instantiation

# REPLACE OPTS HERE WITH THE ONES IN ASSIGN CELLS

# define command line args 
option_list = list(
  make_option(c("-t", "--title"), action="store", default=NA, type='character',
              help="Title for the CELESTA project"),
  make_option(c("-i", "--imagingdata"), action="store", default=NA, type='character',
              help="Path to imaging data"),
  make_option(c("-p", "--prior"), action="store", default=NA, type='character',
              help="Path to prior marker info file"),
  make_option(c("-f", "--filter"), action="store_false", type='logical',
              help="Boolean to filter cells or not (default: FALSE)"),
  make_option(c("-h", "--high"), action="store", default=0.9, type='double',
              help="High marker threshold if filtering cells (default: 0.9)"),
  make_option(c("-l", "--low"), action="store", default=0.4, type='double',
              help="Low marker threshold if filtering cells (default: 0.4)") 
)

# parse args 
opt = parse_args(OptionParser(option_list=option_list))

# create CELESTA object 
CelestaObj <- CreateCelestaObject(
    project_title = opt$title,
    prior_marker_info = opt$prior,
    imaging_data_file = opt$imagingdata
)

# if filtering is specified, filter out cells outside high and low thresholds
if (opt$filter) {
    CelestaObj <- FilterCells(
        CelestaObj, 
        high_marker_threshold = opt$high, 
        low_marker_threshold = opt$low)
}

# plot expression probability
PlotExpProb(coords=CelestaObj@coords,
            marker_exp_prob=CelestaObj@marker_exp_prob,
            prior_marker_info = prior_marker_info,
            save_plot = TRUE)