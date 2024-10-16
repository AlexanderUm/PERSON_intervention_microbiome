#-------------------------------------------------------------------------------
# Correlate multiple columns in the same dataframe 
#-------------------------------------------------------------------------------
corr_multipale <- function(x_columns, 
                           y_columns, 
                           xy_data_frame, 
                           method = "spearman") {
  
  require(tidyverse)
  require(broom)
  
  ResOut <- NULL
  
  for(i in x_columns) {
    
    for(j in y_columns) {
      
      ResOut <- cor.test(xy_data_frame[[i]], 
                   xy_data_frame[[j]], 
                   method = method) %>% 
                  tidy() %>% 
                  mutate(x = i, y = j) %>% 
                  bind_rows(ResOut, .)
      
    }
    
  }

  return(ResOut)
}