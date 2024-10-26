## run on Mccleary
lst_fits_0.5 = prepare_zheng_data(n_cores = 1, exprs_cutoff = 0.5)
simdata_1ct = sim_one_cell_type(lst_fits_0.5)
simdata_2ct_logfc1 = sim_two_cell_types(lst_fits_0.5, n_de = 20, ct_ratio = 1.0, logfc = 1, two_side_de = T, gene_cutoff_quantile = 0.8, logfc_const = T)
simdata_2ct_logfc0.5 = sim_two_cell_types(lst_fits_0.5, n_de = 20, ct_ratio = 1.0, logfc = 0.5, two_side_de = T, gene_cutoff_quantile = 0.8, logfc_const = T)
saveRDS(simdata_1ct, "~/simdata_1ct.rds")
saveRDS(simdata_2ct_logfc1, "~/simdata_2ct_logfc1.rds")
saveRDS(simdata_2ct_logfc0.5, "~/simdata_2ct_logfc0.5.rds")
## scp rds files to local

## run locally
simdata_1ct = readRDS("../rawdata/simdata_1ct.rds")
simdata_2ct = readRDS("../rawdata/simdata_2ct_logfc.rds")
# simdata_2ct_logfc1 = readRDS("../rawdata/simdata_2ct_logfc1.rds")
# simdata_2ct_logfc0.5 = readRDS("../rawdata/simdata_2ct_logfc0.5.rds")

usethis::use_data(simdata_1ct, overwrite = T)
usethis::use_data(simdata_2ct, overwrite = T)
# usethis::use_data(simdata_2ct_logfc1, overwrite = T)
# usethis::use_data(simdata_2ct_logfc0.5, overwrite = T)
