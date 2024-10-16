#-------------------------------------------------------------------------------
# Normalize ASV count with CSS
#-------------------------------------------------------------------------------
phy_css_norm <- function(phylo, 
                         log_transform = TRUE, 
                         taxa_are_rows = TRUE) {
  
  require("phyloseq")
  require("tidyverse")
  require("metagenomeSeq")
  
  otu_table(phylo) <- phylo %>% 
                        phyloseq_to_metagenomeSeq(.) %>% 
                        cumNorm(., p=cumNormStat(.)) %>% 
                        MRcounts(., norm=TRUE, log=log_transform) %>% 
                        as.data.frame() %>% 
                        otu_table(., taxa_are_rows = taxa_are_rows)
  
  return(phylo)
  
}
