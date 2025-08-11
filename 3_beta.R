################################################################################
# Beta diversity
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
source("R/phy_dists_ls.R")
source("R/rda_extract_for_plot.R")
source("R/plot_extracted_rda_data.R")

# Create directories 
for(i in c("tabs", "plots")) {
  
  dir.create(paste0(PRM$beta$dir_out, "/", i), 
             recursive = TRUE, showWarnings = FALSE)
  
}

# Empty list for results
ResBeta <- list()


################################################################################
# Calculate distances for the full data set 
#-------------------------------------------------------------------------------
PrmGird <- expand.grid("Tax_lvl" = PRM$beta$tax_lvl, 
                       "Normal" = PRM$beta$norm, 
                       stringsAsFactors = FALSE)

DistsLs <- list()

for(i in 1:nrow(PrmGird)) {
  
  iLvl <- PrmGird[i, "Tax_lvl"]
  
  iNorm <- PrmGird[i, "Normal"] 
  
  iPs <- DataLs$PS[[iLvl]][[iNorm]]
  
  iDist <- phy_dist_ls(iPs, dists = PRM$beta$distances) %>% 
                  setNames(names(PRM$beta$distances))
  
  DistsLs[[iLvl]][[iNorm]] <- iDist
  
}


################################################################################
# Statistical testing 
#-------------------------------------------------------------------------------
PrmGirdShiftInTime <- expand.grid("Tax_lvl" = PRM$beta$tax_lvl, 
                                  "Normal" = PRM$beta$norm, 
                                  "Strata" = PRM$beta$strata_cols, 
                                  "TestFormula" = PRM$beta$formula,
                                  "Distance" = names(PRM$beta$distances), 
                                  stringsAsFactors = FALSE)
ResShiftInTime <- list()

ResShiftInTimeDf <- NULL

for(i in 1:nrow(PrmGirdShiftInTime)) {
  
  iLvl <- PrmGirdShiftInTime[i, "Tax_lvl"]
  
  iNorm <- PrmGirdShiftInTime[i, "Normal"]
  
  iStrata <- PrmGirdShiftInTime[i, "Strata"]
  
  iFormula <- PrmGirdShiftInTime[i, "TestFormula"]
  
  iDistName <- PrmGirdShiftInTime[i, "Distance"]
  
  
  # Run model per strata level 
  for(j in levels(DataLs$meta[[iStrata]])) { 
    
    # Subset metadata 
    jMeta <- DataLs$meta %>% 
                    filter(.data[[iStrata]] == j)
    
    # Extract and subset corresponding distance
    jDistMatrix <- dist_subset(DistsLs[[iLvl]][[iNorm]][[iDistName]], 
                            rownames(jMeta))
    
    # Run adonis2 
    jResAdon <- adonis2(formula = as.formula(paste0("jDistMatrix ~ ", iFormula)), 
                       data = jMeta, 
                       by = "terms", 
                       permutations = PRM$beta$n_perm, 
                       parallel = 4) %>% 
                  as.data.frame() %>% 
                  mutate(across(everything(), 
                                function(x){round(x, PRM$general$round_to)})) %>% 
                  rownames_to_column(var = "Term") %>% 
                  mutate(Distance = iDistName,
                         Strata = iStrata,
                         Strata_lvl = j,
                         Taxa_lvl = iLvl,
                         Norm_lvl = iNorm,
                         Permutations = PRM$beta$n_perm,
                         Formula = paste0("iDist ~ ", iFormula)) %>% 
                  add_row()
    
    # Collect and write out data
    ResShiftInTimeDf <- bind_rows(ResShiftInTimeDf, jResAdon)
    
    ResShiftInTime[[iLvl]][[iNorm]][[iStrata]][[j]][[iDistName]] <- jResAdon
    
    
    dir.create(paste0(PRM$beta$dir_out, 
                      "/tabs/", iStrata, "/", j),
               showWarnings = FALSE, recursive = TRUE)
    
    write.csv(jResAdon, 
              paste0(PRM$beta$dir_out, 
                     "/tabs/", iStrata, "/", j, "/",
                     iLvl, "-", iNorm, "-", iDistName, "-adonis2.csv"), 
              na = "")
    
  }
  
}

write.csv(ResShiftInTimeDf, 
          paste0(PRM$beta$dir_out, "/tabs/1Combined-adonis2.csv"), 
          na = "", row.names = FALSE)


################################################################################
# Visualize ordination 
#-------------------------------------------------------------------------------
PrmGirdRda <- expand.grid("Tax_lvl" = PRM$beta$tax_lvl, 
                          "Normal" = PRM$beta$norm, 
                          "Strata" = PRM$beta$strata_cols, 
                          stringsAsFactors = FALSE)

ResRdaPlotsLs <- list()

ResRdaPlotsOrdLs <- list()

# Plot text (stat results)
StatTextDf <- ResShiftInTimeDf %>% 
                  filter(Term == trimws(gsub(".*\\+", "",
                                             PRM$beta$rda_plot_formula))) %>% 
                  mutate(across(all_of(c("R2", "Pr(>F)")), 
                                function(x){sprintf("%.3f", round(x, 3))})) %>% 
                  mutate(R2_text = ifelse(R2 == "0.000", 
                                          "R^2<0.001", paste0("R^2==", R2)), 
                         p_text = ifelse(`Pr(>F)` == "0.000", 
                                         "P<0.001", paste0("P==", `Pr(>F)`))) %>% 
                  mutate(Text = paste0(p_text, "~~(", R2_text, ")")) %>% 
                  select(Distance, Strata, Strata_lvl, Taxa_lvl, Norm_lvl, Text)


# Make dbRDA plots 
for(i in 1:nrow(PrmGirdRda)) { 
  
  iLvl <- PrmGirdRda[i,"Tax_lvl"]
  
  iNorm <- PrmGirdRda[i, "Normal"]
  
  iStrata <- PrmGirdRda[i, "Strata"]
  
  for(j in levels(DataLs$meta[[iStrata]])) { 
    
    jMeta <- DataLs$meta %>% 
                filter(.data[[iStrata]] == j) %>% 
                droplevels()
    
    jDistMatrixLs <- DistsLs[[iLvl]][[iNorm]] %>% 
                      lapply(., function(x){dist_subset(x, rownames(jMeta))})
    
    jStatText <- StatTextDf %>% 
                        filter(Taxa_lvl == iLvl, 
                               Norm_lvl == iNorm, 
                               Strata == iStrata, 
                               Strata_lvl == j)
    
    # Custom function: Make RDA data  
    jRdaData <- rda_extract_for_plot(dists_ls = jDistMatrixLs, 
                                     metadata = jMeta, 
                                     form = PRM$beta$rda_plot_formula)
    
    # Plot 
    jRdaPlot <- plot_extracted_rda_data(extracted_data = jRdaData, 
                                        connect_observations = TRUE, 
                                        add_elepses = TRUE, 
                                        group_col = PRM$general$time_point,
                                        observ_time_group_col = PRM$general$part_id_col, 
                                        time_col = PRM$general$time_point, 
                                        stat_text = jStatText,
                                        color_vec = AesLs$col[[PRM$general$time_point]], 
                                        shape_vec = AesLs$shape[[PRM$general$time_point]])
    
    ResRdaPlotsLs[[iLvl]][[iNorm]][[iStrata]][[j]] <- jRdaPlot
    
    
    #-----------------------------------------------------------------------------
    # Reorder plot for the figure panel and write out individual plots 
    #-----------------------------------------------------------------------------
    dir.create(paste0(PRM$beta$dir_out, "/plots/", 
                      iStrata, "/", j), 
               recursive = TRUE, 
               showWarnings = FALSE)
    
    for (k in names(jRdaPlot$Ind)) { 
      
      kPlotOrd <- jRdaPlot$Ind[[k]] + 
                        ggtitle(label = j) + 
                        theme(legend.position="none")
      
      ResRdaPlotsOrdLs[[iLvl]][[iNorm]][[iStrata]][[k]][[j]] <- kPlotOrd
      
      # Save individual plots 
      ggsave(filename = paste0(PRM$beta$dir_out, "/plots/", 
                               iStrata, "/",j, "/",
                               "/dbRDA--", iLvl, "_", iNorm, "--", k, ".png"), 
             plot = jRdaPlot$Ind[[k]], width = 5.5, height = 3.5, dpi = 500)
      
    }
    
  }
  
}


#===============================================================================
# Main results - plot and stat tables 
#-------------------------------------------------------------------------------
PrmGrid <- expand.grid("Tax_lvl" = PRM$beta$tax_lvl, 
                       "Normal" = PRM$beta$norm, 
                       "Strata" = PRM$beta$strata_cols,
                       stringsAsFactors = FALSE)
GrobLs <- list()

for(i in 1:nrow(PrmGrid)) {
  
  iLvl <- PrmGrid[i,"Tax_lvl"]
  
  iNorm <- PrmGrid[i, "Normal"]
  
  iStrata <- PrmGrid[i, "Strata"]
  
  if(iStrata == "Diet") {
    
    iNameAdd <- "Supp_" 
    
  } else {
    
    iNameAdd <- ""
    
  }
  
  # Plots 
  GrobLs <- list()
  
  for(j in names(ResRdaPlotsOrdLs[[iLvl]][[iNorm]][[iStrata]])) {
    
    InstPlotsLs <- ResRdaPlotsOrdLs[[iLvl]][[iNorm]][[iStrata]][[j]]
    
    GrobLs[[j]] <- plot_grid(plotlist = InstPlotsLs, 
                             scale = 0.975, 
                             nrow = 1)
  }
  
  GrobComb <- plot_grid(plotlist = GrobLs, 
                        ncol = 1, 
                        labels = LETTERS[1:length(GrobLs)], 
                        label_size = 22)
  
  GrobCombLeg <- plot_grid(GrobComb, ResRdaPlotsLs[[1]][[1]][[1]]$Legend, 
                           rel_widths = c(0.9, 0.1))
  
  save_plot(filename = paste0(PRM$general$dir_main_fig, "/",
                              iPlotNameAdd, "Fig1_dbRDA_", 
                              iLvl, "_" , iNorm, ".png"), 
            plot = GrobCombLeg, 
            base_height = PRM$beta$PlotGridRowSize*length(GrobLs), 
            base_width = PRM$beta$PlotGridColSize*length(unique(DataLs$meta[[iStrata]])))
  
  # Statistical tables 
  ResShiftInTime[[iLvl]][[iNorm]][[iStrata]] %>% 
    list_flatten() %>% 
    bind_rows() %>% 
    select(-Strata) %>% 
    rename(`P-value` = `Pr(>F)`, 
           Strata = Strata_lvl, 
           `Taxonomic Level` = Taxa_lvl, 
           Normalization = Norm_lvl, 
           `N premutations` = Permutations) %>% 
    write_csv(paste0(PRM$general$dir_main_fig, "/",
                     iPlotNameAdd, "Fig1_dbRDA_", 
                     iLvl, "_" , iNorm, ".csv"), na = "")
  
}


# Clean environment 
rm(list = ls())
gc()

