#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
PRM <- list()

PRM[["general"]] <- list("seed" = 3894568,
                          "ncores" = 4, 
                          "time_point" = "Time",
                          "time_numeric" = "TimeNumeric",
                          "part_id_col" = "ID", 
                          "seq_id_col" = "SeqID", 
                          "diet_col" = "Diet", 
                          "penothype_col" = "Phenotype",
                          "dir_out" = "out", 
                          "dir_Rdata" = "out/Rdata",
                          "libs" =  c("phyloseq", "tidyverse", "metagenomeSeq", 
                                      "qiime2R", "ggsignif", "broom", 
                                      "MicrobiomeStat", "vegan", "cowplot", 
                                      "ComplexHeatmap", "FSA", "usedist", 
                                      "MicrobiomeStat", "Maaslin2"))


PRM[["data"]] <- list("q_path" = "data/", 
                      "m_path" = "data/intervention_metadata.csv", 
                      "min_read_tax" = 20, 
                      "tax_lvls" = c("ASV", "Genus"), 
                      "new_cols" = list("DietPhenotype" = c("Diet", "Phenotype")))


PRM[["alpha"]] <- list("dir_out" = "out/alpha",
                       "measures" = c("Observed", "Shannon",
                                      "InvSimpson", "PhyloDiversity"),
                        "tax_lvl" = "ASV",
                        "norm" = "Rare",
                        "strata_cols" = c("Diet",
                                          "Phenotype",
                                          "DietPhenotype"),
                        "alpha_cut_kw" = 0.05)


PRM[["beta"]] <- list("dir_out" ="out/beta", 
                      "tax_lvl" = c("ASV"),
                      "norm" = "CSS",
                      "distances" = c("Unweighted UniFrac" = "unifrac",
                                      "Weighted UniFrac" = "wunifrac",
                                      "Jaccard" = "jaccard",
                                      "Bray-Curtis" = "bray"),
                      "n_perm" = 9999,
                      "strata_cols" = c("Diet", "DietPhenotype"),
                      "formula" = "ID + TimeNumeric",
                      "rda_plot_formula" = "ID + TimeNumeric",
                      "rda_plot_group" = "Time",
                      "rda_plot_strata" = c("Diet", "DietPhenotype"),
                      "rda_plot_distances" = c("Unweighted UniFrac",
                                               "Weighted UniFrac",
                                               "Jaccard",
                                               "Bray-Curtis"),
                       "PlotGridColSize" = 12,
                       "PlotGridRowSize" = 3,
                       "out_dir_path" = "out/beta")


PRM[["da"]] <- list("dir_out" ="out/da",
                    "tax_lvl" = c("ASV", "Genus"),
                    "norm" = c("Raw"),
                    "strata_cols" = c("Diet", "DietPhenotype"),
                    "maas_method" = "ZINB",
                    "maas_padj_method" = "BH",
                    "maas_random_col" = "ID",
                    "maas_fixed_col" = "TimeNumeric",
                    "maas_norm" = c("TMM"),
                    "maas_trans" = "NONE",
                    "col_prev_by" = "DietPhenotype",
                    "min_prev" = c(0.5),
                    "max_qval" = 0.1,
                    "plot_bar" = TRUE,
                    "plot_area" = FALSE,
                    "violin_ps_norm" = c("CSS", "Relative"),
                    "bar_summary_functions" = c("median"))


PRM[["corr"]] <- list("dir_out" ="out/correlations",
                      "tax_lvl" = c("ASV", "Genus"),
                      "norm" = c("CSS"),
                      "strata_cols" = c("Diet", "DietPhenotype"),
                      "corr_type" = c("Shift", "Baseline"),
                      "corr_method" = "spearman",
                      "qval_method" = "BH",
                      "base_col_lvl" = c("Time" = "Week 0"),
                      "shift_col_lvl" = list("Time" = c("Week 0", "Week 12")),
                      "min_prev" = c("DietPhenotype" = 0.5),
                      "max_qval" = 0.1,
                      "qval_to_use" = "qval_group", # Options: "p.value", "qval_group", "qval_all"
                      "min_est" = 0.2,
                      "corr_cols" = c("Matsuda", "HOMA_IR", 
                                      "misi","hfmm_tag_0", "CRP"),
                      "plot_lm_method" = "lm", 
                      "scatter_p_pref" = "Q")


# Write objects 
dir.create(PRM$general$dir_Rdata, 
           recursive = TRUE, 
           showWarnings = FALSE)

save(list = c("PRM"), 
     file = "PRM.Rdata")

# Clean environment 
rm(list = ls())
gc()
