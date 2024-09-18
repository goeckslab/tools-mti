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
option_list <- list(
  make_option(c("-r", "--rds"), action = "store", default = "celestaobj.rds", type = "character",
              help = "Path to CelestaObj RDS"),
  make_option(c("-c", "--celltypes"), action = "store", default = NA, type = "character",
              help = "Comma-separated list of cell types to plot"),
  make_option(c("-s", "--size"), action = "store", default = 1, type = "double",
              help = "Point size for plotting"),
  make_option(c("--width"), action = "store", default = 12, type = "integer",
              help = "Width of plot (inches)"),
  make_option(c("--height"), action = "store", default = 12, type = "integer",
              help = "Height of plot (inches)"),
  make_option(c("--dpi"), action = "store", default = 300, type = "integer",
              help = "DPI (dots per inch) of plot")
)

# parse args
opt <- parse_args(OptionParser(option_list = option_list))

CelestaObj <- readRDS(opt$rds)
cell_types_to_plot <- strsplit(opt$celltypes, ",")[[1]]

# get indices of cell types to plot from the prior marker table
cell_type_indices <- which(CelestaObj@prior_info[,1] %in% cell_types_to_plot)

print("Cell types selected for plotting:")
print(cell_types_to_plot)
print("Indices of cell types selected for plotting:")
print(cell_type_indices)

# create output directory if it doesn"t already exist
dir.create("cell_assign_plots")

# create the cell type plot
g <- PlotCellsAnyCombination(cell_type_assignment_to_plot = CelestaObj@final_cell_type_assignment[, (CelestaObj@total_rounds + 1)],
                             coords = CelestaObj@coords,
                             prior_info = CelestaObj@prior_info,
                             cell_number_to_use = cell_type_indices,
                             test_size = opt$size,
                             save_plot = FALSE)

# create a unique output name for the plot based on the input cell types
cell_types_cleaned <- paste(make_clean_names(cell_types_to_plot), collapse = "")
output_name <- paste(c("plot_cells_", cell_types_cleaned, ".png"), collapse = "")

# save to subdir
ggsave(
  path = "cell_assign_plots",
  filename = output_name,
  plot = g,
  width = opt$width,
  height = opt$height,
  units = "in",
  dpi = opt$dpi
)
