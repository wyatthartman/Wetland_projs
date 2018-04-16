##########################################################################################################################
# Import Salinity OTU data for testing plotting and corr functions
##########################################################################################################################

# For testing of graphics function using easy loading of Salinity OTU data, 

# Includes:
  
# 1) VST_OTU_CPM table (otu_V):  Variance stabilized OTU table from DESEq2, def Taxonomy, etc...
# 2) Metadata, including merge with OTU table, sort, etc.
# 3) Site colors
# 4) corrRanks output data from corr module (plot testing corr subsets)

##########################################################################################################################
# Load data for testing 
##########################################################################################################################

########################################
# 1) OTU table (DESeq VST CPM), process
########################################

# Import OTU Table
OTU_v <- read.table("SF_Sal_OTU_VSTcpm.txt", sep='\t', header=T, row.names=1)                          # dim(OTU_v); head(OTU_v)

# Sort OTU table                                                                      
otu_V <-OTU_v[order(OTU_v$Consensus.lineage),]                                                         # sort by lineage

# Make new top level plotting var (should be in PRE-PROCESS ? )
otu_V$Taxonomy <- ifelse(otu_V$Phylum == "Proteobacteria", paste(otu_V$Class), paste(otu_V$Phylum))    # head(otu_V)
otu_V <- data.frame(otu_V)

########################################
# 2) Get metadata mapping file, process
########################################

# Import Sample mapping
metaDB <-read.table("SF_sal_meta_FIX3.txt", sep="\t", header=TRUE)          # Import Metadata, keep all    
row.names(metaDB) <- metaDB$Sample                                          # Row names are samples for phyloseq             
metaDB = metaDB[,-1]           

# add any extra variables as needed                                                      # head(Meta_iTag)  #colnames(Meta_iTag)
CH4 <-metaDB$CH4_ug_m2_h
metaDB$CH4_logn1 <- log10(CH4 - (1.05*min(CH4)))                               # Get log n+1 data. since negative, add 5% extra 
#head(Meta_DB)

# Prune metadata to only iTag samples
# Get Sample names in OTU table               
OTU_samps <- data.frame('Sample'=colnames(OTU_v))                                        # OTU_samps

# Merge site order and Samples
Meta_iTag <- merge(metaDB, OTU_samps, by='Sample')                                       #colnames(metaDB)
rownames(Meta_iTag) <- Meta_iTag$Sample

# Reorder location factor
Meta_iTag$Location <-factor(Meta_iTag$Location, levels=c("Sandmound","WestPond","Mayberry","Browns","RushRanch","Joice","Goodyear","WhiteSlough","Tolay","ChinaCamp","Muzzi"))  #head(Meta_iTag)
Meta_iTag$Pl_Sp <-factor(Meta_iTag$Pl_Sp, levels=c("Cattail","Tule","ThreeSq","CattailNL","Phrag","PW","Cord"))

# Resort meta itag by index
indexer = 'EWsiteHyd_index'
Meta_iTag <- Meta_iTag[order(Meta_iTag[indexer]),]
colnames(Meta_iTag)
#head(Meta_iTag)


########################################
# 3) Get site colors table
########################################

# Import site colors
site_colors <- read.table("Sal_siteColors_testR.txt", sep='\t', header=T, row.names=1)               # site_colors
colnames(site_colors) <-c('color','Location')                                                        # site_colors
site_colours <- (site_colors$color)                                                                  # only color

# Get list from vector, move inside functions?
#site_col <- levels(site_colours)[as.numeric(site_colours)]      

########################################
# 4) Get corrRanks output - here for CH4
########################################

## B) Read EXISTING corr table: 
CH4corrRanksOTU <-read.table("CH4_TaxRankOTU_corr_tab_all_sites.txt", sep="\t", header=T)
#head(CH4corrRanksOTU)

## C) Read EXISTING Delta corr table: 
CH4corrRanksOTU_Delta <-read.table("CH4_TaxRankOTU_corr_tab_Delta_sites.txt", sep="\t", header=T)
# head(CH4corrRanksOTU_Delta)


## A) Recalculate corr METHODS for controller, using functions from correlations module:

# Only to genus level:
# CH4corrRanks <- corr_TaxRanks(otu_V, Meta_iTag, "CH4_logn1")                                  # head(CH4corrRanks)

# Append corrs at OTU level --  SLOW (10 min), need only for iTol? 
# OTU_ra <-append_OTU_corr_v(otu_V, CH4corr_ranks, "CH4_logn1")            

# For Delta samples subset
#CH4corrRanks_Delta_OTU <-append_OTU_corr_v(Delta_OTU, CH4corrRanks_Delta, "CH4_logn1")   
#head(CH4corrRanks_Delta_OTU)

# typically saved to file
# write.table(OTU_ra, "CH4_TaxRankOTU_corr_tab_all_sites.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
#write.table(CH4corrRanks_Delta_OTU, "CH4_TaxRankOTU_corr_tab_Delta_sites.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
