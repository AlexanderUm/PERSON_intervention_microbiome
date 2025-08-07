#-------------------------------------------------------------------------------
# Functions to split data into responders and non-responders 
#-------------------------------------------------------------------------------
resp_median <- function(x){ifelse(median(x, na.rm = TRUE) > x,
                                  "NonResp", 
                                  ifelse(median(x, na.rm = TRUE) == x, 
                                         NA, "Resp"))}

resp_5tile <- function(x){case_match(ntile(x, 5), 
                                         1 ~ "Resp", 2 ~ "Resp",
                                         3 ~ NA,
                                         4 ~ "NonResp", 5 ~ "NonResp")}

resp_3tile <- function(x){case_match(ntile(x, 3), 
                                         1 ~ "Resp", 
                                         2 ~ NA,
                                         3 ~ "NonResp")}


