################################################################################
# Differential abundance
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

# Custom functions
source("R/fix_taxa_names_for_plot.R")
source("R/phy_taxa_filter.R")

# Create directories 
for(i in c("tabs", "plots")) {
  
  dir.create(paste0(PRM$da$dir_out, "/", i), 
             recursive = TRUE, showWarnings = FALSE)
  
}

ResDA <- list()

ResDAPlots <- list()

#=============================================================================
# Taxa summary 
#-----------------------------------------------------------------------------
# Will be used to construct supplementary tables 
SummTabsLs <- list()

for(i in PRM$da$tax_lvl) {
  
  OtuRelat <- DataLs$PS[[i]][[PRM$da$tax_summary_norm]] %>% 
                  otu_table() %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame() 
  
  for(j in PRM$da$strata_cols) {
    
    SummTabsLs[[i]][[j]] <- 
              bind_cols(OtuRelat, DataLs$meta) %>% 
              pivot_longer(cols = colnames(OtuRelat), 
                           names_to = "Taxa", 
                           values_to = "Abundance") %>% 
              summarise(Mean = mean(Abundance), 
                        Median = median(Abundance), 
                        SD = sd(Abundance), 
                        Iqr25 = quantile(Abundance, probs = 0.25),
                        Iqr75 = quantile(Abundance, probs = 0.75),
                        .by = all_of(c(j, PRM$general$time_point, "Taxa"))) %>% 
              mutate(across(where(is.numeric), function(x){sprintf("%.3f", round(x, 3))})) %>% 
              mutate(`Mean±SD(%)` = paste0(Mean, "±", SD), 
                     `Median[IQR:25%-75%](%)` = paste0(Median, "[", Iqr25, "-", Iqr75, "]"), 
                     NameCol = paste0("{", !!sym(PRM$general$time_point), "}")) %>% 
              select(c(Taxa, `Mean±SD(%)`, `Median[IQR:25%-75%](%)`, NameCol, !!sym(j))) %>% 
              pivot_wider(id_cols = all_of(c("Taxa", j)), 
                          names_from = NameCol, 
                          values_from = c(`Mean±SD(%)`, `Median[IQR:25%-75%](%)`), 
                          names_sep = " ")
  }

}


################################################################################
# MaAsLin2
#-------------------------------------------------------------------------------
PrmGrid <- expand.grid("Tax_lvl" = PRM$da$tax_lvl, 
                       "Normal" = PRM$da$norm, 
                       "Strata" = PRM$da$strata_cols, 
                       "MassNorm" = PRM$da$maas_norm, 
                       "Prevalence" = PRM$da$min_prev,
                       stringsAsFactors = FALSE)

for(i in 1:nrow(PrmGrid)) { 
  
  # Parameters instances 
  InstLvl <- PrmGrid[i, "Tax_lvl"]
  
  InstNorm <- PrmGrid[i, "Normal"]
  
  InstStrata <- PrmGrid[i, "Strata"]
  
  InstMaasNorm <- PrmGrid[i, "MassNorm"]
  
  InstMinPrev <- PrmGrid[i, "Prevalence"]
  
  InstOutID <- paste0(PRM$da$maas_method, "_", 
                      InstNorm, 
                      "(", InstLvl, "_" ,InstMaasNorm, "-", PRM$da$maas_trans, ")_", 
                      InstMinPrev, "(MinPrev)")
  
  # Phyloseq instance 
  InstOtu <- DataLs$PS[[InstLvl]][[InstNorm]] %>% 
                  phy_taxa_filter(., prev_fraction = InstMinPrev,
                                  group_col = PRM$da$col_prev_by) %>%
                  otu_table() %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame()
  
  InstResDf <- NULL
  
  InstSummCombDf <- NULL
  
  # Stratification
  for(j in levels(DataLs$meta[[InstStrata]])) { 
    
    # Maaslin output directory
    InstMaasDir <- paste0("out/DA/maaslin/", InstOutID, "/", j, "/", InstMaasNorm)
    
    dir.create(InstMaasDir, recursive = TRUE, showWarnings = FALSE)
    
    # Stratify data 
    InstMetaStrata <- DataLs$meta %>% 
                          filter(.data[[InstStrata]] == j)
    
    InstOtuStrata <- InstOtu[rownames(InstMetaStrata), ] 
    
    # Run MaasLin
    InstMaasRes <- Maaslin2(
                      input_data =  InstOtuStrata,
                      input_metadata = InstMetaStrata,
                      output = InstMaasDir,
                      fixed_effects = PRM$da$maas_fixed_col, 
                      random_effects = PRM$da$maas_random_col,
                      correction = PRM$da$maas_padj_method, 
                      standardize = FALSE,
                      cores = PRM$general$ncores,
                      min_abundance = 0,
                      min_prevalence = -Inf,
                      min_variance = 0,
                      normalization = InstMaasNorm,
                      transform = PRM$da$maas_trans, 
                      analysis_method = PRM$da$maas_method,
                      max_significance = 1,
                      plot_scatter = FALSE,
                      plot_heatmap = FALSE,
                      save_scatter = FALSE)
    
    
    InstResTab <- InstMaasRes$results %>% 
                      mutate(Starta = InstStrata, 
                             !!InstStrata := j, 
                             PrevalenceCutOff = InstMinPrev, 
                             PrevalenceCutOffCol = PRM$da$col_prev_by,
                             Method = PRM$da$maas_method,
                             PriorNormalization = InstNorm,
                             MaasNormalization = InstMaasNorm, 
                             MaasTransformation = PRM$da$maas_trans,
                             Prev = N.not.zero/N)
    
    InstMaasRes$results <- InstResTab
    
    ResDA[[InstLvl]][[InstNorm]][[InstMaasNorm]][[InstStrata]][[j]] <- InstMaasRes
    
    InstResDf <- bind_rows(InstResDf, InstResTab)
    
    #===========================================================================
    # Summary tables 
    InstSummCombDf <- InstResTab %>% 
                      select(feature, coef, stderr, pval, qval, !!sym(InstStrata)) %>% 
                      rename(Taxa = feature, `Coef` = coef, SE = stderr, 
                             `P-value` = pval, `Q-value` = qval) %>% 
                      left_join(SummTabsLs[[InstLvl]][[InstStrata]], 
                                by = c("Taxa", InstStrata)) %>% 
                      mutate(across(where(is.numeric), function(x){round(x, 3)})) %>% 
                      mutate(across(where(is.numeric), 
                                    function(x){ifelse(x == 0, 
                                                       "< 0.001", 
                                                       sprintf("%.3f", x))})) %>% 
                      select(all_of(c(colnames(SummTabsLs[[InstLvl]][[InstStrata]]), 
                                      "Coef", "SE","P-value", "Q-value"))) %>% 
                      bind_rows(InstSummCombDf, .)
  }
  
  # Significant taxa names 
  InstSigTaxa <-  InstResDf %>% 
                      filter(!is.na(qval), 
                             qval <= PRM$da$max_qval) %>% 
                      pull(feature) %>% 
                      unique()
  
  dir.create(paste0(PRM$da$dir_out, "/tabs/", InstStrata), 
             recursive = TRUE, showWarnings = FALSE)
  
  write.csv(InstResDf, 
            paste0(PRM$da$dir_out, "/tabs/", InstStrata, 
                   "/MaasLin--", InstOutID, ".csv"), 
            row.names = FALSE)
  
  
  ##############################################################################
  # Plot results 
  #-----------------------------------------------------------------------------
  dir.create(paste0(PRM$da$dir_out, "/plots/", InstStrata), 
             recursive = TRUE, showWarnings = FALSE)
  
  #=============================================================================
  # Overview plot (Maaslin LMM)
  #-----------------------------------------------------------------------------
  # Plot data 
  InstResDfFilt <- InstResDf %>% 
                        filter(feature %in% InstSigTaxa, 
                               !is.na(qval)) %>% 
                        arrange(coef) %>% 
                        mutate(qval_char = ifelse(round(qval, 3) == 0, "q>0.001", 
                                                  paste0("q=", round(qval, 3))), 
                               !!InstStrata := factor(.data[[InstStrata]], 
                                                   levels = levels(DataLs$meta[[InstStrata]])), 
                               feature = gsub("_", "-", fix_taxa_names_for_plot(feature)), 
                               Alpha = ifelse(qval <= PRM$da$max_qval, 
                                              paste0("Q-value", " \u2264 ", PRM$da$max_qval), 
                                              paste0("Q-value > ", PRM$da$max_qval))) %>% 
                        mutate(feature = factor(feature, levels = unique(feature)), 
                               Alpha = factor(Alpha, 
                                              levels = (c(paste0("Q-value", 
                                                                 " \u2264 ", 
                                                                 PRM$da$max_qval), 
                                                          paste0("Q-value > ", 
                                                                 PRM$da$max_qval)))))
  
  # Plot accessories
   AlphaColor <- c("gray", "black") %>% 
                    setNames(c(paste0("Q-value > ", PRM$da$max_qval), 
                               paste0("Q-value", " \u2264 ", PRM$da$max_qval)))
   
   LimX <- (max(abs(InstResDfFilt$coef)) + max(InstResDfFilt$stderr))
  
   # Plot itself 
   MaasPlot <- ggplot(InstResDfFilt, 
                           aes(x = coef, 
                               y = feature, 
                               color = Alpha)) + 
                geom_errorbar(aes(xmin = coef - stderr, 
                                  xmax = coef + stderr), 
                             size=0.75, width=0.15) + 
                geom_vline(xintercept = 0, 
                           alpha = 0.75, 
                           color = "gray", 
                           size = 0.75) + 
                geom_point(size = 3.5) + 
                geom_text(aes(label = qval_char), 
                          vjust = -1, 
                          fontface = "italic", 
                          size = 2.5, 
                          show.legend = FALSE) +
                facet_grid(as.formula(paste0("~", InstStrata)), 
                           scales = "free_y", 
                           space = "free") + 
                theme_bw() +
                theme(panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(), 
                      axis.text.y = element_text(face = "italic", 
                                                 size = 12, 
                                                 family = "serif", 
                                                 color = "black")) +
                ylab("") + 
                xlab("") +
                scale_color_manual(values = AlphaColor) + 
                theme(legend.text = element_text(size=12),
                      legend.title = element_text(size=15), 
                      axis.text.x = element_text(angle = 45, hjust = 1)) +
               scale_x_continuous(limits = c(-LimX, LimX))
   
   MaasPlotLs <- list("plot" = MaasPlot, 
                      "w" = length(levels(DataLs$meta[[InstStrata]]))*2.75 + 2.5, 
                      "h" = length(InstSigTaxa)*0.275 + 1.5)
   
   ggsave(paste0(PRM$da$dir_out, "/plots/", InstStrata, "/", "MaasLin--", 
                 InstOutID, ".png"), 
          plot = MaasPlotLs$plot, 
          width = MaasPlotLs$w, 
          height = MaasPlotLs$h, 
          limitsize = FALSE)
   
   ResDAPlots[[InstLvl]][[InstStrata]][[InstOutID]][["General"]] <- MaasPlotLs
   
   #----------------------------------------------------------------------------
   # Summary tables associated with the plot (Supplementary table 3)
   InstSummCombDf %>% 
     filter(Taxa %in% InstSigTaxa) %>% 
     mutate(Taxa = gsub("_", "-", fix_taxa_names_for_plot(Taxa))) %>% 
     arrange(!!sym(InstStrata), 
             factor(Taxa, levels = rev(levels(InstResDfFilt$feature)))) %>% 
     write.csv(file = paste0(PRM$general$dir_main_fig, 
                             "/Supp_DA_Tab3_", InstLvl, "_", InstStrata, ".csv"))
     
   
   #============================================================================
   # Violin plot
   #----------------------------------------------------------------------------
   for(ViolinNorm in PRM$da$violin_ps_norm) { 
     
   # Main data 
   InstOtuViolin <- DataLs$PS[[InstLvl]][[ViolinNorm]] %>% 
                      otu_table() %>% 
                      as.matrix() %>% 
                      t() %>% 
                      as.data.frame()
   
   InstViolinDf <- bind_cols(InstOtuViolin[, InstSigTaxa], 
                             DataLs$meta) %>% 
                    pivot_longer(cols = all_of(InstSigTaxa), 
                                 names_to = "feature", 
                                 values_to = "Abundance") %>% 
                    mutate(feature = factor(gsub("_", "-", fix_taxa_names_for_plot(feature)), 
                                            levels = rev(levels(InstResDfFilt$feature)))) %>% 
                    arrange(across(all_of(c("feature", 
                                            PRM$general$part_id_col, 
                                            PRM$general$time_point))))
   
   InstViolinDfSum <- InstViolinDf %>% 
                        summarise(AbundanceDiff = diff(Abundance), 
                                  .by = all_of(c("feature", 
                                                 PRM$general$part_id_col))) %>% 
                        mutate(AbundanceShift = ifelse(AbundanceDiff > 0, "Increase",
                                                       ifelse(AbundanceDiff < 0, 
                                                              "Decrease", 
                                                              "No Change")))
   
   InstViolinDfComb <- left_join(InstViolinDf, 
                                 InstViolinDfSum, 
                                 by = c("feature", PRM$general$part_id_col)) 
   
   # Significance data 
   InstSigDf <- InstViolinDf %>% 
                  summarise(AundanceMax = max(Abundance)*1.15, 
                            .by = all_of(c("feature", InstStrata))) %>% 
                  left_join(InstResDfFilt, ., by = c("feature", InstStrata)) %>% 
                  mutate(coef_char = ifelse(round(coef, 2) == 0, "\u03B2>0.01", 
                                              paste0("\u03B2=", round(coef, 2)))) %>% 
                  mutate(label = paste0(qval_char, 
                                        " [", coef_char, "]"), 
                         x_start = levels(InstViolinDfComb[[PRM$general$time_point]])[1], 
                         x_end = levels(InstViolinDfComb[[PRM$general$time_point]])[2])
            
   # Plot violin plot 
   ViolinPlot <- ggplot(InstViolinDfComb, aes(x = .data[[PRM$general$time_point]], 
                                y = Abundance)) + 
                     geom_line(aes(group = ID, color = AbundanceShift), 
                               alpha = 0.25, show.legend = FALSE) +
                     geom_violin(fill = NA) + 
                     facet_grid(c("feature", InstStrata), 
                                scales = "free_y") + 
                     theme_bw() +
                     theme(axis.title.x = element_text()) + 
                     guides(colour = guide_legend(order = 2), 
                            linetype = guide_legend(order = 1)) +
                     theme(axis.text.x = element_text(angle=45, 
                                                      vjust=1.1, 
                                                      hjust=1), 
                           strip.text = element_text(size = 6), 
                           axis.title.x = element_blank(), 
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) + 
                     scale_color_manual(values = c("Increase" = "forestgreen", 
                                                   "Decrease" = "firebrick", 
                                                   "No Change" = "gray")) + 
                     geom_signif(data = InstSigDf,
                                 aes(xmin = x_start,
                                     xmax = x_end,
                                     annotations = label,
                                     y_position = AundanceMax, 
                                     alpha = Alpha),
                                 textsize = 2, 
                                 vjust = -0.1,
                                 manual = TRUE,
                                 margin_top = 1, 
                                 show.legend = FALSE) +
                     geom_point(data = InstSigDf,
                                aes(x = x_start, 
                                    y = AundanceMax*1.05), 
                                x = NA)  + 
                     scale_alpha_manual(values = setNames(c(1, 0.5), 
                                                 levels(InstSigDf$Alpha))) +
                     ylab(label = paste0("Abundance (", ViolinNorm, ")"))
       
   
   ViolinPlotLs <- list("plot" = ViolinPlot, 
                        "w" = length(levels(DataLs$meta[[InstStrata]]))*1.2+1, 
                        "h" = length(InstSigTaxa)*1.5 + 1.5)
   
   ggsave(paste0(PRM$da$dir_out, "/plots/", InstStrata, "/", 
                 "Violin(", ViolinNorm, ")--",
                 InstOutID, ".png"), 
          plot = ViolinPlotLs$plot, 
          width = ViolinPlotLs$w, 
          height = ViolinPlotLs$h, 
          limitsize = FALSE)
   
   ResDAPlots[[InstLvl]][[InstStrata]][[InstOutID]][["Volin"]][[ViolinNorm]] <- 
     ViolinPlotLs
   
   
   #============================================================================
   # Bar plot
   #----------------------------------------------------------------------------
   for(k in PRM$da$bar_summary_functions) { 
     
     BarDf <- InstViolinDfComb %>% 
               summarise(AbundanceSummary = eval(call(k, Abundance)), 
                         .by = all_of(c("feature", InstStrata, 
                                        PRM$general$time_point)))
     
     BarTaxaColor <- AesLs$ColVect[2:(length(levels(BarDf$feature))+1)] %>% 
                      setNames(levels(BarDf$feature))
     
     if(PRM$da$plot_bar) {
         
       BarPlot <- ggplot(BarDf, aes(y = AbundanceSummary, 
                         x = .data[[PRM$general$time_point]], 
                         fill = feature)) + 
                     geom_bar(stat = "identity") + 
                     facet_grid(as.formula(paste0("~", InstStrata))) + 
                     scale_fill_manual(values = BarTaxaColor) + 
                     theme_classic() + 
                     ylab(label = paste0("Abundance (", ViolinNorm, ")")) + 
                     theme(legend.text = element_text(size=12,
                                                      family = "serif", 
                                                      face = "italic"), 
                           legend.position = "left", 
                           legend.title = element_blank(), 
                           legend.text.position = "left", 
                           axis.title.x = element_blank(), 
                           axis.text.x = element_text(angle = 45, hjust = 1))
       
       
       # Plot dimensions 
       BarPlotLs <- list("plot" = BarPlot, 
                         "w" = length(levels(BarDf[[InstStrata]])) * 
                               length(levels(BarDf[[PRM$general$time_point]]))*0.7 +
                               ceiling(length(levels(BarDf$feature))/20)*2.5 + 0.2, 
                         "h" = 5.5)
                            
       
       ggsave(filename = paste0(PRM$da$dir_out, "/plots/", InstStrata, "/", 
                                "BarSummary(", ViolinNorm, "_",k, ")--",
                                InstOutID, ".png"), 
              plot = BarPlotLs$plot, 
              width = BarPlotLs$w, 
              height =  BarPlotLs$h, 
              limitsize = FALSE,
              dpi = 600)
     
       
       ResDAPlots[[InstLvl]][[InstStrata]][[InstOutID]][["Bar"]][[ViolinNorm]] <- 
         BarPlotLs
      }
    }
  }
}

#-------------------------------------------------------------------------------
# Write out results tables 
#-------------------------------------------------------------------------------
save(list = c("ResDAPlots", "ResDA"), 
     file = paste0(PRM$da$dir_out, "/da.Rdata"))


################################################################################
# Combine plots into a panel for publication 
#-------------------------------------------------------------------------------
PlotMainDiet <- ResDAPlots$Genus$Diet[[1]]$General$plot + 
                  theme(legend.position = "none")

PlotMainDietPhenotype <- 
          ResDAPlots$Genus$DietPhenotype[[1]]$General$plot + 
          theme(legend.position = "none")

PlotBarRelatDietPhenotype <- 
          ResDAPlots$Genus$DietPhenotype[[1]]$Bar$Relative$plot + 
          guides(fill = guide_legend(ncol = 1))

PlotDaComb <- plot_grid(PlotMainDietPhenotype, 
                  plot_grid(PlotBarRelatDietPhenotype, NULL, PlotMainDiet, 
                            rel_widths = c(0.5, 0.005, 0.6), 
                            labels = c("B", "", "C"), 
                            nrow = 1, 
                            label_size = 22), 
                  ncol = 1, 
                  labels = c("A", ""), 
                  label_size = 22)

save_plot(plot = PlotDaComb, 
          filename = paste0(PRM$general$dir_main_fig, "/Fig3_DA.png"), 
          base_height = 15, 
          base_width = 15)

# Clean environment 
rm(list = ls())
gc()
