################################################################################
# Alpha diversity
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
source("R/phy_alpha.R")
source("R/miscellaneous.R")


# Create directories 
for(i in c("tabs", "plots")) {
  
  dir.create(paste0(PRM$alpha$dir_out, "/", i), 
             recursive = TRUE, showWarnings = FALSE)
  
}

# Empty list for results
ResAlpha <- list()


#-------------------------------------------------------------------------------
# Objects 
#-------------------------------------------------------------------------------
Ps <- DataLs$PS[[PRM$alpha$tax_lvl]][[PRM$alpha$norm]]


################################################################################
# Calculate alpha diversity
#-------------------------------------------------------------------------------
AlphaDf <- phy_alpha(Ps, measures = PRM$alpha$measures) %>% 
               bind_cols(., DataLs$meta) 

AlphaDfLong <- AlphaDf %>% 
                pivot_longer(cols = all_of(PRM$alpha$measures), 
                             names_to = "Index")


# Write out results 
ResAlpha[["Tables"]][["Diversity per Sample"]] <- AlphaDf


################################################################################
# Alpha diversity shift in time per group level (Paired Wilcoxon)
#-------------------------------------------------------------------------------
ShiftTestRes <- list()
  
# Grid for iterations 
IterGrid <- expand.grid("Variable" = PRM$alpha$strata_cols, 
                        "Index" = PRM$alpha$measures, 
                          stringsAsFactors = FALSE)
  
  
for(i in 1:nrow(IterGrid)) {
    
    iVar <- IterGrid[i, "Variable"]
    
    iInd <- IterGrid[i, "Index"] 
    
    ResWl <- NULL
    
    for(j in levels(DataLs$meta[[iVar]])) {
      
      # Subset data 
      jAlphaDataLs <- list()
      
      for(k in 1:2) {
        
        jAlphaDataLs[[k]] <- AlphaDf %>% 
                      filter(.[[iVar]] == j, 
                             .[[PRM$general$time_numeric]] == 
                               sort(unique(.[[PRM$general$time_numeric]]))[k]) %>% 
                      droplevels() %>% 
                      arrange(.[[PRM$general$part_id_col]]) 
      }
      
      # Test and collect results per variable level into a single data frame
      ResWl <- wilcox.test(jAlphaDataLs[[1]][[iInd]], jAlphaDataLs[[2]][[iInd]],
                           paired = TRUE, exact = FALSE) %>% 
                  tidy() %>% 
                  mutate(!!iVar := j) %>% 
                  bind_rows(ResWl, .) 
      
    }
    
    ShiftTestRes[[iVar]][[iInd]] <- ResWl %>% 
                                        mutate(Index = iInd, 
                                               Variable_col = iVar)
  }
  
ResAlpha[["Tables"]][["Stat_shift"]] <-  ShiftTestRes
  
  
#-----------------------------------------------------------------------------
# Plot the data
#-----------------------------------------------------------------------------
# Plot annotation dataframe
MaxMinDf <- AlphaDfLong %>% 
                reframe(y_min = min(value), 
                        y_max = max(value), 
                        .by = Index) 
  
ShiftPlots <- list()
  
for(i in PRM$alpha$strata_cols) { 
  
  if(i == "Diet") {
    
    iNameAdd <- "Supp_Fig2_" 
    
  } else {
    
    iNameAdd <- "Fig2_"
    
  }
  
  # Plot significance annotation
  iStatRes <- ShiftTestRes[[i]] %>% 
                    bind_rows() %>% 
                    mutate(p_short = round(p.value, 3), 
                           CompID = paste0(Index, "_", Variable_col)) %>% 
                    mutate(p_text = ifelse(p_short == 0, 
                                           "P<0.001", 
                                           paste0("P=", sprintf("%.3f", 
                                                                p_short)))) %>% 
                    left_join(., MaxMinDf, by = "Index")
  
  iSigAnnotDf <- iStatRes %>% 
                    mutate(Start = levels(DataLs$meta[[PRM$general$time_point]])[1], 
                           End = levels(DataLs$meta[[PRM$general$time_point]])[2], 
                           y_sig = y_max + (y_max - y_min) * 0.2, 
                           y_inv_point = y_max + (y_max - y_min) * 0.4)
  
  write.csv(iStatRes, 
            file = paste0(PRM$alpha$dir_out, "/tabs/", 
                          "stat_alpha_", i, ".csv"), 
            row.names = FALSE)
    
  # Plot dimensions 
  iWidth <- length(levels(DataLs$meta[[PRM$general$time_point]])) * 
                    length(levels(DataLs$meta[[i]])) + 
                    (length(levels(DataLs$meta[[i]])) - 1)*0.25
    
  iHeight <- length(PRM$alpha$measures)*1.75
    
    
  # Plot
  iPlot <- ggplot(AlphaDfLong, 
                     aes(y = value, 
                         x = .data[[PRM$general$time_point]])) + 
                  geom_line(
                    aes(group = .data[[PRM$general$part_id_col]]), 
                            alpha = 0.2) +
                  geom_point(size = 0.2) +
                  geom_violin(fill = NA) + 
                  facet_grid(c("Index", i), scales = "free") + 
                  geom_signif(data = iSigAnnotDf,
                                aes(xmin = Start,
                                    xmax = End,
                                    annotations = p_text,
                                    y_position = y_sig),
                                textsize = 3.5, 
                                vjust = -0.2,
                                manual = TRUE, 
                                margin_top = 1) + 
                  geom_point(data = iSigAnnotDf,
                             aes(x = End, 
                                 y = y_inv_point), 
                             x=NA)  + 
                  theme_bw() + 
                  theme(axis.line = element_line(color='black'),
                        plot.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), 
                        axis.title.y = element_blank(), 
                        axis.text.x = element_text(angle = 45, 
                                                   hjust = 1), 
                        axis.title.x = element_blank()) 
  
    # Collect plots 
    ShiftPlots[[i]] <- list(plot = iPlot, 
                            p.width = iWidth, 
                            p.height = iHeight)
    
    # Save plots 
    ggsave(paste0(PRM$general$dir_main_fig, "/", iNameAdd, "alpha.png"), 
           plot = iPlot, 
           width = iWidth, 
           height = iHeight)
    
    
    #===========================================================================
    # Write stat tables with summary 
    AlphaSummary <- AlphaDfLong %>% 
      summarise(Mean = mean(value), 
                Median = median(value), 
                SD = sd(value), 
                Iqr25 = quantile(value, probs = 0.25),
                Iqr75 = quantile(value, probs = 0.75),
                .by = all_of(c(i, PRM$general$time_point, "Index"))) %>% 
      mutate(across(where(is.numeric), function(x){sprintf("%.3f", round(x, 3))})) %>% 
      mutate(`Mean±SD(%)` = paste0(Mean, "±", SD), 
             `Median[IQR:25%-75%](%)` = paste0(Median, "[", Iqr25, "-", Iqr75, "]"), 
             NameCol = paste0("{", !!sym(PRM$general$time_point), "}")) %>% 
      select(c(Index, `Mean±SD(%)`, `Median[IQR:25%-75%](%)`, NameCol, !!sym(i))) %>% 
      pivot_wider(id_cols = all_of(c("Index", i)), 
                  names_from = NameCol, 
                  values_from = c(`Mean±SD(%)`, `Median[IQR:25%-75%](%)`), 
                  names_sep = " ")
      
    # Combine summary and statistical results 
    lapply(ShiftTestRes[[i]], function(x) {add_row(x)}) %>% 
      bind_rows() %>% 
      rename(`P-value` = `p.value`, 
             `Statsistics (W)` = statistic) %>% 
      select(-c(Variable_col, method, alternative)) %>% 
      mutate(across(where(is.numeric), function(x){sprintf("%.3f", round(x, 3))})) %>% 
      left_join(AlphaSummary, by = c(i, "Index")) %>% 
      select(all_of(c(colnames(AlphaSummary), colnames(.)))) %>% 
      write_csv(paste0(PRM$general$dir_main_fig, "/", iNameAdd, "alpha.csv"), 
                na = "")
    
}


# Clean environment 
rm(list = ls())
gc()
