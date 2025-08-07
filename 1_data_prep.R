#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

set.seed(PRM$general$seed)

for(i in PRM$general$libs) {
  
  require(i, character.only = TRUE)
  
}


#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_shorten_tax_names.R")
source("R/phy_css_norm.R")


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
DietCols <- PRM$general$diet_col

PhenCol <- PRM$general$penothype_col

TimePointCol <- PRM$general$time_point

TimeNumeric <- PRM$general$time_numeric

IdCol <- PRM$general$part_id_col

SeqIdCol <- PRM$general$seq_id_col

DirRdata <- PRM$general$dir_Rdata

# Data preparation
qPath <- PRM$data$q_path

MetaPath <- PRM$data$m_path

TaxLvls <- PRM$data$tax_lvls

MinReads <- PRM$data$min_read_tax

NewCols <- PRM$data$new_cols

# Create directories 
dir.create(DirRdata, showWarnings = FALSE, recursive = TRUE)

# Object to fill up 
DataLs <- list()


################################################################################
# Read data/ Create phyloseq
#-------------------------------------------------------------------------------
ps <- qza_to_phyloseq(features = paste0(qPath, "asv_table.qza"), 
                      tree = paste0(qPath, "tree/rooted-tree.qza"), 
                      taxonomy = paste0(qPath, "taxonomy_07.qza"))

#-------------------------------------------------------------------------------
# Filter out taxa
#-------------------------------------------------------------------------------
# filter taxa with with less than X reads in total   
ps <- prune_taxa(taxa_sums(ps) >= MinReads, ps)

# Remove ASVs: 
ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Mitochondria", ps)

ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Chloroplast", ps)

ps <- prune_taxa(!is.na(tax_table(ps)[, "Phylum"])[, "Phylum"], ps)

ps <- prune_taxa(tax_table(ps)[, "Kingdom"] %in% c("d__Bacteria", 
                                                   "d__Archaea"), ps)


#-------------------------------------------------------------------------------
# Add metadata
#-------------------------------------------------------------------------------
# Metadata 
meta <- read.csv(MetaPath) 

rownames(meta) <- meta[[SeqIdCol]]

# Overlapping samples
OverSamp <- intersect(sample_names(ps), rownames(meta))

ps <- prune_samples(OverSamp, ps)

ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Adjust metadata
meta <- meta[sample_names(ps), ] %>% 
        mutate(!!TimePointCol := as.factor(.[[TimePointCol]]), 
               !!TimeNumeric := (as.numeric(gsub("[^0-9]", 
                                                "", 
                                                .data[[TimePointCol]])) + 1),
               !!DietCols := as.factor(.[[DietCols]]), 
               !!PhenCol := as.factor(.[[PhenCol]]), 
               All_samples = as.factor("All_samples")) %>% 
        dplyr::rename(TAG = hfmm_tag_0, MISI = misi)

# Add new columns
for(i in names(NewCols)) {
  
  meta <- meta %>% 
            mutate(!!i := apply(.[NewCols[[i]]], 1, 
                   function(x){paste0(x, collapse = "-")})) %>% 
            mutate(!!i := as.factor(.[[i]]))
}


# Add metadata to phyloseq
sample_data(ps) <- meta

DataLs[["meta"]]  <- meta


################################################################################
# Prepare list with data 
# Glom to higher taxonomic level ->
# Normalize count ->
# Make strata
#-------------------------------------------------------------------------------
for(i in TaxLvls)  {
  
  if(i == "ASV") { 
    
    InstPs <- ps
    
    taxa_names(InstPs) <- phy_shorten_tax_names(InstPs) %>% 
                              as.data.frame() %>% 
                              setNames("feature") %>% 
                              mutate(feature2 = if(n() > 1) {
                                    paste0(feature, "__asv", row_number())} else {
                                    paste0(feature, "__asv")}, .by = "feature") %>% 
                              pull(feature2)
  
  } else { 
    
    InstPs <- tax_glom(ps, i)}
  
  # Adjust taxa names
  taxa_names(InstPs) <- phy_shorten_tax_names(InstPs) %>% 
                                  as.data.frame() %>% 
                                  setNames("feature") %>% 
                                  mutate(feature2 = if(n() > 1) {
                                    paste0(feature, "__", 
                                           tolower(str_sub(i, 1, 1)), 
                                           row_number())} else {
                                             feature}, .by = "feature") %>% 
                                  pull(feature2)
  
  InstPsRare <- rarefy_even_depth(InstPs, rngseed = PRM$general$seed)
  
  DataLs[["PS"]][[i]] <- list("Raw" = InstPs, 
                              "Rare" = InstPsRare, 
                              "Rare_log" = transform_sample_counts(InstPsRare, 
                                              function(x){log2(x+1)}), 
                              "CSS" = phy_css_norm(InstPs, 
                                                    log_transform = TRUE, 
                                                    taxa_are_rows = taxa_are_rows(InstPs)), 
                              "Relative" = transform_sample_counts(InstPs, 
                                              function(x){x/sum(x)*100}))
}


################################################################################
# Aesthetics 
#-------------------------------------------------------------------------------
AesLs <- list()

ColVec <-  c("#999999", "#7FC97F", "#BEAED4", "#FDC086", "#386CB0", # 5
             "#BF5B17", "#1B9E77", "#7570B3", "#66A61E", "#E6AB02", # 10
             "#A6761D", "#A6CEE3", "#33A02C", "#FB9A99", "#FDBF6F", # 15
             "#FF7F00", "#6A3D9A", "#FFFF99", "#B15928", "#CCEBC5", # 20
             "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", # 25
             "#B3E2CD", "#FDCDAC", "#CBD5E8", "#E6F5C9", "#FFF2AE", # 30
             "#F1E2CC", "#CCCCCC", "#E41A1C", "#984EA3", "#FF7F00", # 35
             "#FFFF33", "#A65628", "#66C2A5", "#FC8D62", "#8DA0CB", # 40
             "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#8DD3C7", # 45
             "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", # 50
             "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F") # 55

ShapeVec <- c(16, 17, 18, 15, 8, 1:7)

AesLs[["ColVect"]] <- ColVec

AesLs[["col"]][[TimePointCol]] <- ColVec[c(2, 5)] %>% 
                                   setNames(levels(meta[[TimePointCol]]))

AesLs[["col"]][[DietCols]] <- ColVec[c(6, 8)] %>% 
                                  setNames(levels(meta[[DietCols]]))

AesLs[["col"]][[PhenCol]] <- ColVec[c(10, 13)] %>% 
                                  setNames(levels(meta[[PhenCol]]))

AesLs[["col"]][["DietPhenotype"]] <- ColVec[c(17, 13, 49, 48)] %>% 
                                  setNames(levels(meta[["DietPhenotype"]]))

AesLs[["shape"]][[TimePointCol]] <- ShapeVec[c(1, 2)] %>% 
                                      setNames(levels(meta[[TimePointCol]]))



#-------------------------------------------------------------------------------
save(list = c("DataLs", "AesLs"), 
     file = paste0(DirRdata, "/data.Rdata"))

# Clean environment 
rm(list = ls())
gc()
