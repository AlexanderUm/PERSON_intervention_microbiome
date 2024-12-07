################################################################################
# 
################################################################################
# Load data
load("PRM.Rdata")

load(paste0(PRM$general$dir_Rdata, "/data.Rdata"))

# Load libraries
for(i in PRM$general$libs) {
  
  require(i, character.only = TRUE)
  
}

# Set seed
set.seed(PRM$general$seed)


# Create directories 
for(i in c("tabs/sig", "tabs/all", "plots")) {
  
  dir.create(paste0(PRM$corr$dir_out, "/", i), 
             recursive = TRUE, showWarnings = FALSE)
  
}

# Custom functions 
source("R/phy_taxa_filter.R")
source("R/corr_multipale.R")
source("R/fix_taxa_names_for_plot.R")

format_to_print <- function(x){
                      ifelse(round(x, 3) == 0, 
                             ">0.001", 
                              paste0("=", sprintf("%.3f", round(x, 3))))}



ResCorrLs <- list()

ResCorrPlotsLs <- list()


################################################################################
# Correlations 
################################################################################

# Indexes shift 
IndexShift <- DataLs$meta %>% 
                select(all_of(c(PRM$corr$corr_cols, 
                                PRM$corr$strata_cols, 
                                PRM$general$part_id_col, 
                                names(PRM$corr$shift_col_lvl)))) %>% 
                arrange(.data[[PRM$general$part_id_col]], 
                        factor(.data[[names(PRM$corr$shift_col_lvl)]], 
                               levels = PRM$corr$shift_col_lvl[[1]])) %>% 
                mutate(across(PRM$corr$corr_cols, diff), 
                        .by = all_of(PRM$general$part_id_col)) %>% 
                filter(.data[[names(PRM$corr$shift_col_lvl)]] == 
                         PRM$corr$shift_col_lvl[[1]][1]) %>% 
                select(-all_of(names(PRM$corr$shift_col_lvl)))


# Correlations iterations 
PrmGrid <- expand.grid("Tax_lvl" = PRM$corr$tax_lvl, 
                       "Norm" = PRM$corr$norm, 
                       "Strata" = PRM$corr$strata_cols, 
                       "Type" = PRM$corr$corr_type, 
                       "Method" = PRM$corr$corr_method, 
                       stringsAsFactors = FALSE)


for(i in 1:nrow(PrmGrid)) { 
  
  iLvl <- PrmGrid[i, "Tax_lvl"]
  
  iNorm <- PrmGrid[i, "Norm"]

  iStrata <- PrmGrid[i, "Strata"]
  
  iType <- PrmGrid[i, "Type"]
  
  iMethod <- PrmGrid[i, "Method"]
  
  iNameAdd <- paste(PrmGrid[i, ], collapse = "--")
  
  
  # Extract OTU table 
  OtuTab <- DataLs$PS[[iLvl]][[iNorm]] %>% 
                phy_taxa_filter(prev_fraction = PRM$corr$min_prev, 
                                group_col = names(PRM$corr$min_prev)) %>% 
                otu_table() %>% 
                as.matrix() %>% 
                t() %>% 
                as.data.frame() %>% 
                setNames(fix_taxa_names_for_plot(names(.)))
  
  
  # Format data for baseline correlations
  if(iType == "Baseline") {
    
    OtuTabForm <- bind_cols(OtuTab, DataLs$meta) %>% 
                      filter(.data[[names(PRM$corr$base_col_lvl)]] == 
                               PRM$corr$base_col_lvl) %>% 
                      select(all_of(c(colnames(OtuTab), 
                                      PRM$general$part_id_col)))
  }
  
  # Format data for shift correlations 
  if(iType == "Shift") {
    
    OtuTabForm <- bind_cols(OtuTab, DataLs$meta) %>% 
                        arrange(.data[[PRM$general$part_id_col]], 
                                factor(.data[[names(PRM$corr$shift_col_lvl)]], 
                                       levels = PRM$corr$shift_col_lvl[[1]])) %>% 
                        summarise(across(colnames(OtuTab), diff), 
                                .by = all_of(PRM$general$part_id_col))
  }
  
  DataCorr <- left_join(OtuTabForm, 
                        IndexShift, 
                        by = PRM$general$part_id_col)
  
  
  ResCorrAllLvls <- NULL
  
  # Correlate per strata level 
  for(j in levels(IndexShift[[iStrata]])) {
    
    # Subset data 
    DataCorrFilt <- DataCorr %>% 
                  filter(.data[[iStrata]] == j)
    
    # Correlate, adjust p-values, and add metadata 
    ResCorr <- corr_multipale(x_columns = colnames(OtuTab), 
                              y_columns = PRM$corr$corr_cols, 
                              xy_data_frame = DataCorrFilt, 
                              method = iMethod) %>% 
                arrange(y, p.value) %>% 
                mutate(qval_all = p.adjust(p.value, 
                                           method = PRM$corr$qval_method)) %>% 
                mutate(qval_group = p.adjust(p.value, 
                                             method = PRM$corr$qval_method), 
                       .by = y) %>% 
                bind_cols(PrmGrid[i, ], Strata_Lvl = j)
    
    # Combine dataframe 
    ResCorrAllLvls <- bind_rows(ResCorrAllLvls, ResCorr)
    
    # Record results 
    write.csv(ResCorr, 
              file = paste0(PRM$corr$dir_out, "/tabs/all/Corr_", 
                            iNameAdd, "--", j, ".csv"), 
              row.names = FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # Visualize correlations 
  #-----------------------------------------------------------------------------
  SigTaxaTab <- ResCorrAllLvls %>% 
                  filter(.data[[PRM$corr$qval_to_use]] <= PRM$corr$max_qval, 
                         abs(estimate) >= PRM$corr$min_est) %>% 
                  arrange(abs(estimate))
  
  if(nrow(SigTaxaTab) > 1) {
    
    # Write out significant restults 
    write.csv(SigTaxaTab, 
            file = paste0(PRM$corr$dir_out, "/tabs/sig/CorrSig_", 
                          iNameAdd, ".csv"), 
            row.names = FALSE)
    
    #-----------------------------------------------------------------------------
    # Heat map 
    # Data preparation   
    HeatDfs <- list()
    
    for(j in c("estimate", PRM$corr$qval_to_use)) {
      
      HeatDfs[[j]] <- ResCorrAllLvls %>% 
                          filter(x %in% SigTaxaTab$x) %>% 
                          select(all_of(c("x", "y" , j , "Strata_Lvl"))) %>% 
                          pivot_wider(names_from = all_of(c("Strata_Lvl", "y")), 
                                      values_from = all_of(j), 
                                      names_sep = "--") %>% 
                          column_to_rownames(var = "x") %>% 
                          select(order(colnames(.)))
    }
    
    # Heat map plot 
    CorrHeatPlot <- Heatmap(HeatDfs$estimate, 
                            name = "Correlation \nestimate",
                            column_labels = gsub(".*--", "", colnames(HeatDfs$estimate)), 
                            column_split = gsub("--.*", "", colnames(HeatDfs$estimate)), 
                            row_names_side = "left", 
                            row_names_gp = gpar(fontface = "italic", 
                                                fontfamily = "serif"), 
                            column_names_rot = 45, 
                            height = unit(nrow(HeatDfs$estimate)*0.25, "in"), 
                            width = unit(ncol(HeatDfs$estimate)*0.4, "in"),
                            cluster_rows = TRUE, 
                            show_row_dend = FALSE,
                            cluster_columns = FALSE, 
                            rect_gp = gpar(col = "gray25", lwd = 0.5),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              if(HeatDfs[[PRM$corr$qval_to_use]][i, j] <= 
                                 PRM$corr$max_qval) {
                                grid.text(sprintf("%.3f", 
                                                  HeatDfs[[PRM$corr$qval_to_use]][i, j]), 
                                          x, y, gp = gpar(fontsize = 8, 
                                                          col = "white", 
                                                          fontface = "bold"))
                              } else {
                                grid.text(sprintf("%.3f", 
                                                  HeatDfs[[PRM$corr$qval_to_use]][i, j]), 
                                          x, y, gp = gpar(fontsize = 6, 
                                                          col = "gray40"))}})
    
    # Collect & write
    CorrHeatPlotLs <- list("plot" = CorrHeatPlot, 
                           "w" = ncol(HeatDfs$estimate)*0.4 + 5, 
                           "h" = nrow(HeatDfs$estimate)*0.25 + 2)
    
    png(filename =  paste0(PRM$corr$dir_out, "/plots/Heat_", iNameAdd, ".png"), 
        width = CorrHeatPlotLs$w, 
        height = CorrHeatPlotLs$h, 
        units = "in", 
        res = 600)
    
    draw(CorrHeatPlotLs$plot)
    
    dev.off()
    
    
    #---------------------------------------------------------------------------
    # Scatter plots 
    #---------------------------------------------------------------------------
    dir.create(paste0(PRM$corr$dir_out, "/plots/scatter/", iNameAdd), 
               showWarnings = FALSE, recursive = TRUE)
    
    # Text color 
    ColorForQvalue <- c("gray" = "gray", "black" = "black")
    
    ResCorrAllLvls <- ResCorrAllLvls %>%
                      mutate(across(all_of(c(PRM$corr$qval_to_use, "estimate")), 
                                    format_to_print, 
                                    .names = "text_{.col}")) %>% 
                      mutate(text = paste0(PRM$corr$scatter_p_pref, 
                                           .data[[paste0("text_", 
                                                         PRM$corr$qval_to_use)]], 
                                           " [est", text_estimate, "]"), 
                             !!iStrata := Strata_Lvl)
    
    
    for(j in 1:nrow(SigTaxaTab)) {
      
      jTaxa <- SigTaxaTab$x[j]
      
      jMet <- SigTaxaTab$y[j]
      
      # Data preparation 
      # Scatter plot data 
      DataScatter <- DataCorr %>% 
                      select(all_of(c(jTaxa, jMet, iStrata))) %>% 
                      mutate(feature = jTaxa, 
                             top_strip = paste0(jMet, ": ", .data[[iStrata]])) 
      
      # Min and max 
      MinYDf<- DataScatter %>% 
                    summarise(across(all_of(jMet), 
                                     function(x) min(x, na.rm = TRUE)), 
                              .by = all_of(iStrata))
      
      # Significance label
      TextDf <- ResCorrAllLvls %>% 
                      filter(x == jTaxa, 
                             y == jMet) %>% 
                      mutate(text_color = ifelse(.data[[PRM$corr$qval_to_use]] <= 
                                                   PRM$corr$max_qval, 
                                                  "black", "gray"), 
                             top_strip = paste0(jMet, ": ", 
                                                .data[[iStrata]]), 
                             !!jTaxa := max(DataScatter[[jTaxa]])) %>% 
                      left_join(MinYDf, by = iStrata)
      
    # Plot                
    ScatterPlot <- ggplot(DataScatter, 
                     aes(x = .data[[jMet]], 
                         y = .data[[jTaxa]])) +
                      geom_point() + 
                      geom_smooth(method = PRM$corr$plot_lm_method, 
                                  se = FALSE) + 
                      geom_text(data = TextDf, 
                                aes(label = text, 
                                    color = text_color), 
                                hjust = 0,
                                vjust = -0.5) + 
                      facet_grid(c("feature", "top_strip"), 
                                 scales = "free") + 
                      theme_bw() + 
                      theme(axis.title = element_blank(), 
                            strip.text.y = element_text(face = "italic", 
                                                        family = "serif"), 
                            legend.position = "none") + 
                      scale_color_manual(values = ColorForQvalue) + 
                      ylim(c(min(DataScatter[[jTaxa]])*1.15, 
                             max(DataScatter[[jTaxa]])*1.15))
    
    ggsave(filename = paste0(PRM$corr$dir_out, 
                             "/plots/scatter/", iNameAdd, "/",  
                             gsub(" |\\.", "_", jTaxa), "-", jMet, ".svg"), 
           plot = ScatterPlot, 
           width = nrow(TextDf)*2.5+0.5, 
           height = 3)
        
    }
  
  }
  
}

# Clean environment 
rm(list = ls())
gc()
