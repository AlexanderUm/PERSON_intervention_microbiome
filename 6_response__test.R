################################################################################
# Response analysis
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
for(i in c("tabs/", "plots")) {
  
  dir.create(paste0(PRM$resp$dir_out, "/", i), 
             recursive = TRUE, showWarnings = FALSE)
  
}


# Custom functions
source("R/split_funs.R")
source("R/phy_taxa_filter.R")

################################################################################
# Data preparation
################################################################################
# Indexes shift 
IndexShift <- DataLs$meta %>% 
                  select(all_of(c(PRM$resp$resp_base_cols, 
                                  PRM$resp$strata_cols, 
                                  PRM$general$part_id_col, 
                                  names(PRM$resp$shift_col_lvl)))) %>% 
                  arrange(.data[[PRM$general$part_id_col]], 
                          factor(.data[[names(PRM$resp$shift_col_lvl)]], 
                                 levels = PRM$resp$shift_col_lvl[[1]])) %>% 
                  mutate(across(PRM$resp$resp_base_cols, diff), 
                         .by = all_of(PRM$general$part_id_col)) %>% 
                  filter(.data[[names(PRM$resp$shift_col_lvl)]] == 
                           PRM$resp$shift_col_lvl[[1]][1]) %>% 
                  select(-all_of(names(PRM$resp$shift_col_lvl)))



VarGrid <- expand.grid("fun" = PRM$resp$split_fun, 
                       "strata" = PRM$resp$strata_cols, 
                       stringsAsFactors = FALSE)

RespLs <- list()

for(i in 1:nrow(VarGrid)) {
  
  RespDf <- IndexShift %>% 
              mutate(across(PRM$resp$resp_base_cols, 
                            get(VarGrid$fun[i]), 
                            .names = "Response--{.col}"), 
                     .by = VarGrid$strata[i]) %>% 
              rename_with(~paste0("Index--", .x), 
                          all_of(PRM$resp$resp_base_cols))
  
  RespLs[[VarGrid$strata[i]]][[VarGrid$fun[i]]] <- RespDf %>% 
                                                    select(-starts_with("Index--"))
 
  RespDfLong <- RespDf %>% 
                  dplyr::select(all_of(c(paste0("Response--", 
                                                PRM$resp$resp_base_cols), 
                                         paste0("Index--", 
                                                PRM$resp$resp_base_cols), 
                                         PRM$general$part_id_col, 
                                         VarGrid$strata[i]))) %>% 
                  pivot_longer(cols = -c(PRM$general$part_id_col, 
                                         VarGrid$strata[i]), 
                                         names_to = c(".value", "Name"), 
                                         names_sep = "--" )
    
  RespPlot <- ggplot(RespDfLong,
                     aes(x = .data[[VarGrid$strata[i]]],
                         y = Index)) +
                geom_jitter(aes(color = Response), 
                            width = 0.2,
                            height = 0,
                            alpha = 0.9, 
                            size = 1) +
                geom_violin(fill = NA) +
                facet_wrap(. ~ Name, scales = "free", nrow = 1) +
                theme_bw() + 
                theme(axis.title.x = element_blank())
  
  PlotWidth <- length(PRM$resp$resp_base_cols)*
                length(levels(RespDf[[VarGrid$strata[i]]]))+2
  
  ggsave(filename = paste0(PRM$resp$dir_out, 
                           "/plots/Response--", 
                           VarGrid$fun[i], "--", 
                           VarGrid$strata[i], ".svg"), 
         plot = RespPlot, 
         width = PlotWidth, 
         height = 3.5)
  
}


################################################################################
# RF analysis 
################################################################################
# Correlations iterations 
PrmGrid <- expand.grid("Tax_lvl" = PRM$resp$tax_lvl, 
                       "Norm" = PRM$resp$norm, 
                       "Strata" = PRM$resp$strata_cols, 
                       "Type" = PRM$resp$data_type, 
                       "Method" = PRM$resp$split_fun, 
                       "Index" = PRM$resp$resp_base_cols,
                       stringsAsFactors = FALSE)

RocDataDf <- NULL

ResRfCombDf <- NULL

for(i in 1:nrow(PrmGrid)) { 
  
  iLvl <- PrmGrid[i, "Tax_lvl"]
  
  iNorm <- PrmGrid[i, "Norm"]
  
  iStrata <- PrmGrid[i, "Strata"]
  
  iType <- PrmGrid[i, "Type"]
  
  iMethod <- PrmGrid[i, "Method"]
  
  iIndex <- PrmGrid[i, "Index"]
  
  iNameAdd <- paste(PrmGrid[i, ], collapse = "--")
  
  
  # Extract OTU table 
  OtuTab <- DataLs$PS[[iLvl]][[iNorm]] %>% 
                phy_taxa_filter(prev_fraction = PRM$resp$min_prev, 
                                group_col = names(PRM$resp$min_prev)) %>% 
                otu_table() %>% 
                as.matrix() %>% 
                t() %>% 
                as.data.frame() 
  
  
  # Format data for baseline correlations
  if(iType == "Baseline") {
    
    OtuTabForm <- bind_cols(OtuTab, DataLs$meta) %>% 
                      filter(.data[[names(PRM$resp$base_col_lvl)]] == 
                               PRM$resp$base_col_lvl) %>% 
                      select(all_of(c(colnames(OtuTab), 
                                      PRM$general$part_id_col)))
  }
  
  # Format data for shift correlations 
  if(iType == "Shift") {
    
    OtuTabForm <- bind_cols(OtuTab, DataLs$meta) %>% 
                      arrange(.data[[PRM$general$part_id_col]], 
                              factor(.data[[names(PRM$resp$shift_col_lvl)]], 
                                     levels = PRM$resp$shift_col_lvl[[1]])) %>% 
                      summarise(across(colnames(OtuTab), diff), 
                                .by = all_of(PRM$general$part_id_col))
  }
  
  
  # Add response data 
  ColName <- setNames(paste0("Response--", iIndex), iIndex)
  
  DataRfFull <- left_join(OtuTabForm, 
                        RespLs[[iStrata]][[iMethod]][c(paste0("Response--", iIndex), 
                                                       PRM$general$part_id_col, 
                                                       iStrata)], 
                        by = PRM$general$part_id_col) %>% 
                        rename(all_of(ColName)) %>% 
                        select(-all_of(PRM$general$part_id_col)) %>% 
                        filter(!is.na(.data[[iIndex]]))
  
  for(j in levels(DataRfFull[[iStrata]])) {
    
    DataRf <- DataRfFull %>% 
                filter(.data[[iStrata]] == j) %>% 
                select(-all_of(iStrata))
    
    Nsize <- DataRf[[iIndex]] %>% 
              table() %>% 
              as.matrix() %>% 
              t() %>% 
              as.data.frame()
    
    # Vector with the true response
    RespTrue <- DataRf[[iIndex]]
    
    #-----------------------------------------------------------------------------
    # RF with CV
    #-----------------------------------------------------------------------------
    # Mtry
    SeedMtry <- round(sqrt((ncol(DataRf)-1)), 0) 
    
    MtryGrid <- expand.grid(.mtry = SeedMtry)
    
    
    # Control
    Control <- trainControl(method = "repeatedcv", 
                            repeats = PRM$resp$cv_repeats,
                            number = PRM$resp$n_cv, 
                            classProbs = TRUE, 
                            savePredictions = "all", 
                            summaryFunction = twoClassSummary) 
    
    # Model
    RfCvRes <- train(as.formula(paste0(iIndex, "~.")),
                     data = DataRf,
                     method = "rf",
                     metric = "ROC",
                     importance=TRUE, 
                     tuneGrid=MtryGrid, 
                     trControl = Control,
                     ntree = PRM$resp$n_trees)
    
    # Data for ROC
    RocDataDf <- roc(RfCvRes$pred$obs, 
                     RfCvRes$pred$Resp)[c("sensitivities", "specificities")] %>% 
                    as.data.frame() %>% 
                    mutate(Random = "True", 
                           Correspondance = 100) %>% 
                    arrange(-row_number()) %>% 
                    bind_cols(., PrmGrid[i, ], Strata_Lvl = j) %>% 
                    bind_rows(RocDataDf, .) %>% 
                    suppressMessages()
    
    # Results table
    ResRfCombDf <- RfCvRes$results %>% 
                        bind_cols(., 
                                  PrmGrid[i, ], 
                                  Strata_Lvl = j, 
                                  Nsize) %>% 
                        mutate(Random = "True", 
                               Correspondance = 100, 
                               Mtry = SeedMtry) %>% 
                        bind_rows(ResRfCombDf, .)
    
    #!!!--------------------------------------------------------------------------
    #!!! Randomize response
    #!!!--------------------------------------------------------------------------
    for(k in 1:PRM$resp$n_random) {
      
      # Randomize response with sample
      DataRf <- DataRf %>% 
                  mutate(!!iIndex := sample(.data[[iIndex]]))
      
      # Correspondence between randomized and true response 
      CoresPrec <- sum(DataRf[[iIndex]] == RespTrue)/length(RespTrue)*100 %>% 
                       round(2)
      
      # Run the model (control from above)
      RfCvRes <- train(as.formula(paste0(iIndex, "~.")),
                       data = DataRf,
                       method = "rf",
                       metric = "ROC",
                       importance=TRUE, 
                       tuneGrid=MtryGrid, 
                       trControl = Control,
                       ntree = PRM$resp$n_trees)
      
      # Add to ROC data frame 
      RocDataDf <- roc(RfCvRes$pred$obs, 
                       RfCvRes$pred$Resp)[c("sensitivities", "specificities")] %>% 
                    as.data.frame() %>% 
                    mutate(Random = paste0("Random_", k), 
                           Correspondance = CoresPrec) %>% 
                    arrange(-row_number()) %>% 
                    bind_cols(., PrmGrid[i, ], Strata_Lvl = j) %>% 
                    bind_rows(RocDataDf, .) %>% 
                    suppressMessages()
      
      # Add to results data frame 
      ResRfCombDf <- RfCvRes$results %>% 
                          bind_cols(., 
                                    PrmGrid[i, ], 
                                    Strata_Lvl = j, 
                                    Nsize) %>% 
                          mutate(Random = paste0("Random_", k), 
                                 Correspondance = CoresPrec, 
                                 Mtry = SeedMtry) %>% 
                          bind_rows(ResRfCombDf, .)
      
    }
    
  }
  
  print(paste0(i, " from ", nrow(PrmGrid)))
  
}


#-------------------------------------------------------------------------------
# Plot AUC for original and random
#-------------------------------------------------------------------------------
AucPlot <- list()

for(i in unique(ResRfCombDf[["Type"]])) {
  
  AucPlot[[i]] <- ResRfCombDf %>% 
                    filter(Type == i) %>% 
                    arrange(NonResp) %>% 
                    mutate(Color = gsub("_.*", "", Random), 
                           Strata_Lvl = factor(Strata_Lvl, 
                                               levels = unique(Strata_Lvl))) %>% 
                    ggplot(aes(y = ROC, x = Strata_Lvl)) + 
                    geom_jitter(aes(colour = Color), 
                                height = 0, width = 0.05, 
                                size = 2.5, alpha = 0.5) + 
                    facet_grid(c("Index", "Method")) + 
                    theme_bw() +
                    theme(strip.text = element_text(size = 8, face = "italic"),
                          axis.title = element_blank(), 
                          axis.text.x = element_text(angle = 45, hjust = 1), 
                          legend.position = "bottom")
  
  
  ggsave(filename = paste0(PRM$resp$dir_out, 
                           "/plots/AUC_comparison--",  i,".svg"), 
         plot = AucPlot[[i]], 
         width = 6, 
         height = 7)
  
}


