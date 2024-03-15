# ---------------------------------------------------------------------------------
# Plot assigned cell type combinations with CELESTA 
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

# define command line args 
option_list = list(
  make_option(c("-r", "--rds"), action="store", default='celestaobj.rds', type='character',
              help="Path to CelestaObj RDS"),
  make_option(c("-p", "--prior"), action="store", default=NA, type='character',
              help="Path to prior marker info file"),
  make_option(c("-c", "--celltypes"), action="store", default=NA, type='character',
              help="Comma-separated list of cell types to plot"),
  make_option(c("-s", "--size"), action="store", default=NA, type='character',
              help="point size for plotting")
)

# parse args 
opt = parse_args(OptionParser(option_list=option_list))

CelestaObj <- readRDS(opt$rds)
cell_types_to_plot <- strsplit(opt$celltypes, ",")[[1]]

# read prior marker info 
prior <- read.csv(opt$prior,row.names = 1)

# get indices of cell types to plot from the prior marker table 
cell_type_indices <- which(row.names(prior) %in% cell_types_to_plot) 

print(cell_types_to_plot)
print(cell_type_indices)

print(row.names(prior))

# create output directory if it doesn't already exist 
dir.create('cell_assign_plots')

# create the cell type plot 
g <- PlotCellsAnyCombination(cell_type_assignment_to_plot=CelestaObj@final_cell_type_assignment[,(CelestaObj@total_rounds+1)],
                        coords = CelestaObj@coords,
                        prior_info = prior_marker_info,
                        cell_number_to_use=cell_type_indices,
                        test_size=1,
                        save_plot = FALSE)

# create a unique output name for the plot based on the input cell types 
cell_types_cleaned <- paste(make_clean_names(cell_types_to_plot), collapse = "")
output_name <- paste(c("plot_cells_",cell_types_cleaned,".png"), collapse = "")

# save to subdir 
# FIXME: may want to expose plotting params to galaxy 
ggsave(
    path = 'cell_assign_plots',
    filename = output_name,
    plot = g,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300
)