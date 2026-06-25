# scdesign2_ind_gen_data.R - Simulate new count matrix with scDesign2 'ind' (independent) mode.


library(Matrix)
library(scDesign2)


# Configuration
seed_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/seed/afc/HCC3N_600spot.spotxgene.afc.rds"
cell_anno_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/cell_anno/spot_anno.tsv"
out_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/afc_scdesign2_ind/simu"
out_fn_prefix = 'scDesign2_ind'
sim_method = 'ind'

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# Load Count Matrix (cells labelled with cell types)
data_mat = readRDS(seed_fn)
str(data_mat)

cell_anno = read.table(cell_anno_fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(cell_anno) = c("cell", "cell_type")
str(cell_anno)

data_mat = t(data_mat)
if (! all(colnames(data_mat) == cell_anno$cell))
    stop("some cells do not match!")
colnames(data_mat) = cell_anno$cell_type
str(data_mat)


# remove spike-in -----------------------------------------------------------------------
nonspikes <- which(!grepl("ercc", rownames(data_mat), ignore.case = TRUE))
print(paste("number of spike-ins:", nrow(data_mat)-length(nonspikes)))
data_mat <- data_mat[nonspikes, ,drop = FALSE]


# explore basic structure of data -------------------------------------------------------
print(dim(data_mat))
print(table(colnames(data_mat)))


# Simulate Count Matrix

# Since in model fitting, we will use the mclappy() function, 
# for reproducibility, we need to change the random number generator in the following way.
RNGkind("L'Ecuyer-CMRG")


# fit model -----------------------------------------------------------
set.seed(1)
model_fit_result <- fit_model_scDesign2(
    data_mat, 
    cell_type_sel = 'normal', 
    sim_method = sim_method,
    ncores = 3
)


# simulate data -----------------------------------------------------------
sim_count_result <- simulate_count_scDesign2(
    model_fit_result, 
    n_cell_new = ncol(data_mat) * 2, 
    sim_method = sim_method,
    cell_type_prop = 1
)


# save the model parameters and the simulated data --------------------------------------
saveRDS(model_fit_result, file = sprintf("%s/%s.model_fit_result.rds", out_dir, out_fn_prefix))
saveRDS(sim_count_result, file = sprintf("%s/%s.simu_count.rds", out_dir, out_fn_prefix))
