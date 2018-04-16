# 16s corr by Tax Rank module

##########################################################################################################################
## Tax rank correlation functions, test with heatmap

# DESCRIPTION:

# Functions to correlate counts at each Taxonomic rank with external variable.  These are nested as they are built.
# These functions correlate log2 transformed VST counts, will need work otherwise. Rank names are generic (e.g. could work for KO), 
#and distinguished as non-numeric in input OTU table.

# Inputs: OTU table with split taxonomy; metadata table & corr. var within it (e.g. CH4)

### FUNCTIONS in module

### 1) Aggregate OTUs by Tax. category
# agg_by_cat(otu_V, "Taxonomy")

### 2) Correlation vector ........................................................... [on aggregated]
# corr_vector(otu_agg, "CH4_logn1")

### 3) Aggrgate correlation vector ........................................... [combines Fxn.1 and Fxn.2]
# agg_corr_vect("Taxonomy", "CH4_logn1")

### 4) Apply correlation vector across tax ranks .................... [applies Fxn.3 to tax ranks]
# corr_byRank(otu_V)  

### 5) Correlation by ranks wrapper  ....................................... [re-writes fxn.3 with fixed corr_var and applies Fxn.4]
# corr_TaxRanks(otu_V, Meta_iTag, "CH4_logn1")

### 6) Append OTU correlations ............................................... [appends OTU level correlations to output of Fxn.5]
# append_OTU_corr_v(otu_V, corr_ranks, "CH4_logn1")

### 7) Test corr output by plotting corr heatmap................................[Test on Genus level]
# NOT FUNCTION, more elaborate function for this in corr_plots modules 

##########################################################################################################################
# get packages

# Import packagaes
library(Hmisc)
#suppressMessages(library(Hmisc))
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gplots)                       # for testing only
#suppressMessages(library(cowplot))
#theme_set(theme_grey())              # important, ggplot2 default overwrites cowplot defaults, which add unwanted items to plots

##########################################################################################################################
##########################################################################################################################
# Input test data
##########################################################################################################################
setwd("~/Desktop/SF_Sal_OTU")
getwd()


#### Get OTU table (DESeq VST CPM) ####

# Import OTU Table
OTU_v <- read.table("SF_Sal_OTU_VSTcpm.txt", sep='\t', header=T, row.names=1)                          # dim(OTU_v); head(OTU_v)

# Sort OTU table                                                                      
otu_V <-OTU_v[order(OTU_v$Consensus.lineage),]                                                         # sort by lineage

# Make new top level plotting var (should be in PRE-PROCESS ? )
otu_V$Taxonomy <- ifelse(otu_V$Phylum == "Proteobacteria", paste(otu_V$Class), paste(otu_V$Phylum))    # head(otu_V)
otu_V <- data.frame(otu_V)

#### Get metadata mapping file ####

# Import Sample mapping
metaDB <-read.table("SF_sal_meta_FIX3.txt", sep="\t", header=TRUE)          # Import Metadata, keep all    
row.names(metaDB) <- metaDB$Sample                                          # Row names are samples for phyloseq             
metaDB = metaDB[,-1]       

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
# colnames(Meta_iTag)

# add any extra variables as needed                                                      # head(Meta_iTag)  #colnames(Meta_iTag)
CH4 <-Meta_iTag$CH4_ug_m2_h
Meta_iTag$CH4_logn1 <- log10(CH4 - (1.05*min(CH4)))                               # Get log n+1 data. since negative, add 5% extra 
# head(Meta_iTag)

#### Get site colors ####
# Import site colors
# site_colors <- read.table("Sal_siteColors_testR.txt", sep='\t', header=T, row.names=1)               # site_colors
# colnames(site_colors) <-c('color','Location')                                                        # site_colors
# site_colours <- (site_colors$color)     


##########################################################################################################################
# FUNCTIONS defined below for correlations by category
##########################################################################################################################

####################################################
# 1) Aggregatate OTUs by category
####################################################      
# Aggregation by category 

agg_by_cat = function(otu_V, agg_var) {                                                                    
  count_cols <- sapply(otu_V, is.numeric)                              # Select numeric columns
  otu_counts <- data.matrix(otu_V[,count_cols])                        # Keep numeric cols, matrix for phyloseq  
  
  agg_fact <- unlist(as.list(otu_V[agg_var]))                          # Get levels for aggregation
  otu_agg <- aggregate(otu_counts, by=list(agg_fact), FUN=sum)         # Aggregate OTU table by agg_fact
  colnames(otu_agg)[1] <- agg_var                                      # rename new groups "Taxa"   
  row.names(otu_agg) <- otu_agg[,1]                                    # set rownames
  return(otu_agg)               
}

# Test aggregation function  
#otu_agg <- agg_by_cat(otu_V, "Taxonomy")        #otu_agg <- agg_by_cat(otu_V, "Phyla_cat")  #otu_agg <- agg_by_cat(otu_V, "Phylum")
#head(otu_agg) #head(otu_aggP)

####################################################
# 2) Correlations vector function 
####################################################      
# Make correlation vector for a given variable and otu table (including aggregated)

corr_vector = function(otu_t, corr_var, fdr_cut=0.05, corr_type = "spearman") {
  
  # Separate numeric data and Taxonomy in OTU_table
  count_cols <- sapply(otu_t, is.numeric)                            # Select numeric columns
  otu_d <- data.frame(otu_t[,count_cols])                            # Keep numeric cols                      # head otu_d
  Tax_data <- otu_t[(count_cols == FALSE)]                           # Save taxonomy data                     # head(Tax_data)
  
  # Fill zero counts with two prior to log2 transform
  otu_d2<- otu_d                                                     # define new var as data copy
  otu_d2[otu_d2==0]<- 1                                              # replace 0 with 1 -- now log2(2) = 0
  otu_d2t <- data.frame(t(log2(otu_d2)))                             # log2 and transpose, df for rownames 
  otu_d2t$Sample <- rownames(otu_d2t)                                # add "Sample" column for merge
  
  # select only continuous data from Meta data
  meta_d_cols <- sapply(Meta_iTag, is.numeric)                       # Select numeric columns
  meta_d <- data.frame(Meta_iTag[,meta_d_cols])                      # Keep numeric cols                      # head otu_d
  meta_d$Sample <- Meta_iTag$Sample
  
  # Merge with Metadata by sample
  otu_d2_meta <- merge(meta_d, otu_d2t, by="Sample")
  otu_d2_meta <- data.matrix(otu_d2_meta[,-1])                       # Drop Sample, data matrix makes float for rcorr   #head(otu_d2_meta)
  
  # Get correlation square matrix for r and P values                                          # can be slow at OTU scale 
  envOTUcorr_M <- rcorr(otu_d2_meta, type=c(corr_type))                                      # Spearman correlations, non-param
  envOTUcorr_m <- data.frame("r"=envOTUcorr_M$r[,corr_var], "P"=envOTUcorr_M$P[,corr_var])    # combine r and P
  envOTUcorr_m$P_fdr <- p.adjust(envOTUcorr_m$P, method="BH")                                 # correct for FDR, false discov. rate
  
  # Trim off metadata correlations by rows
  lmeta <- ncol(meta_d)                                              # length of metadata columns
  lmat <- dim(envOTUcorr_m)[1]                                       # length of correlation matrix
  envOTUcorr_mOTUs <- envOTUcorr_m[lmeta:lmat,]                      # drop metadata rows
  
  groups <- data.frame(rownames(envOTUcorr_mOTUs))                   # get group names (agg_ taxa) 
  colnames(groups) <- colnames(otu_t)[1]                             # name column after input column (agg_var)   
  
  # Join corr. r values and group (Taxon) names, if P < P fdr_cut, else 0
  groups$r_fdr <- ifelse(envOTUcorr_mOTUs$P_fdr < fdr_cut, envOTUcorr_mOTUs$r, 0)  # return only corr w. P < Pfdr_cut
  return(groups)
}
# test correlation function on aggregated dataset
#tax_vect <- corr_vector(otu_agg, "CH4_logn1")
#head(tax_vect)

####################################################
# 3) Combined aggregate and correlate function
####################################################      

agg_corr_vect = function(agg_var, otu_V = otu_V, corr_var="CH4_logn1"){
  otu_agg <- agg_by_cat(otu_V, agg_var)                              # agg_by_cat function on agg_var
  corr_vector(otu_agg, corr_var)                                     # corr_vector function on corr_var, otu_agg
}

# Test of agg_corr_vect function -- need to update to include otu table var
# Acorr_vect <- agg_corr_vect("Taxonomy", "CH4_logn1")
# Acorr_vect


####################################################
# 4) Aggregate and correlate by Tax Rank function
####################################################      
# Function for correlations to corr_var across taxonomic ranks

corr_byRank = function(otu_V) {
  
  # Get list of taxa names to iterage over, from Taxa matrix
  #otu_V$OTU <-row.names(otu_V)                                      # Should corr at OTU level, not working
  count_cols <- sapply(otu_V, is.numeric)                            # Only numeric columns  
  Tax_dat <- otu_V[(count_cols == FALSE)]                            # Get taxa matrix -- non-numeric
  Tax_names <- colnames(Tax_dat)                                     # Get only names of taxa levels
  Tax_names <- c(Tax_names[-1])                                      # Drop Consensus.lineage, don't want to agg on this.#Tax_names
  
  # Get corr vector function output for each var in Tax_names  
  AllTax_corr_v <- lapply(Tax_names, agg_corr_vect, otu_V)           # lapply tax names to agg_corr_vector for each rank    #AllTax_corr_v
  # also inputs otu_V into function           # returns list of correlation vectors for each tax rank
  
  # Loop for each Taxon. rank, merge corr vectors to Full Taxonomy 
  len_tax <-length(Tax_names)                                        # get length of taxonomy list 
  Tax_corr <- Tax_dat                                                # Copy of taxonomy data for adding corr vectors
  
  for (i in seq(1:len_tax)){                             
    y <- data.frame(AllTax_corr_v[i])                              # Extract data frame (corr_vector) for item in list
    Tax_rank2 <-colnames(y)[1]                                     # Get colname (Tax rank) to relabel r column 
    Tax_lab <- paste0(Tax_rank2,'_r')                              # Make paste-appended label w Tax rank_r
    colnames(y)[colnames(y)=='r_fdr'] <- Tax_lab                   # Relabel column
    Tax_corr <-merge(Tax_corr, y, all.x=TRUE)                  # Recursively MERGE current vector with taxonomy 
    Tax_corr <-Tax_corr[order(Tax_corr$Consensus.lineage),]        # Sort by Consensus.lineage
  }
  
  # reorder columns to match input
  tax_cols <- colnames(Tax_dat)                                      # Get input Taxonomy colnames, original order 
  r_cols <- colnames(Tax_corr[(ncol(Tax_dat)+1):ncol(Tax_corr)])     # get correlaton colnames
  
  #r_cols <- colnames(corr_test[(ncol(Tax_dat)+1):ncol(corr_test)])  # get correlaton colnames
  sort_cols <- c(tax_cols, r_cols)                                   # combined Tax and correlation, right order       #sort_cols                                     
  Tax_corrS <-Tax_corr[,sort_cols]                                   # get df in new column order 
  
  return(Tax_corrS)
}

# test function for concatenating output of agg_corr_vect
#corr_table <- corr_byRank(otu_V)                           # old version # corr_table <- concat_taxrank_corr(otu_V)
#dim(corr_table); head(corr_table)

####################################################
# 5) Correlate all tax ranks wrapper function
####################################################

# Rebuilds corr_vector(fxn.3) using fixed corr_var to enable iteration by tax rank (fxn.4)
corr_TaxRanks = function(otu_V, Meta_iTag, corr_var){
  
  cvx <- deparse(substitute(corr_var))                           # paste corr_var char, use below in function
  
  agg_corr_vect = function(agg_var, otu_V, Meta_iTag, cvx){      # Redef (fxn 3) w. fixed corr_var (appl by agg_var)
    otu_agg <- agg_by_cat(otu_V, agg_var)                      # agg_by_cat function on agg_var (fxn.1)
    vector <- corr_vector(otu_agg, Meta_iTag, cvx)             # corr_vector function on corr_var, otu_agg (fxn.2)
    return(vector)                                                 # return vector of corrs in corr_var, agg_var
  }
  
  Tax_corr<-corr_byRank(otu_V)                                   # Apply (fxn 4 to fxn 3 above) to tax ranks
  Tax_corr[is.na(Tax_corr)] <-0
  return(Tax_corr)                                               # return correlations by rank
}
# Test this function
#corr_ranks <- corr_TaxRanks(otu_V, Meta_iTag, "CH4_logn1")
#head(corr_ranks)

####################################################
# 6) Get OTU corrs function 
####################################################      

# Append OTU correlations function -- possibly 8-10 min

append_OTU_corr_v = function(otu_V, corr_ranks, corr_var){#, corr_var){
  # Copy OTU table, add OTU #s
  otu_Vo <- otu_V                                                            # Copy OTU table
  otu_Vo$OTU <- row.names(otu_Vo)                                            # Add OTU numbers #head(otu_Vo)
  
  # Aggregate by OTU, make corr_vector
  OTU_agg <- agg_by_cat(otu_Vo, "OTU")                                       # Aggregate on OTUs                    #head(OTU_agg)
  OTU_corr <- corr_vector(OTU_agg, corr_var)                              # Correlate on OTUs
  colnames(OTU_corr)[2] <- "OTU_r"                                           # Rename OTU column
  
  # Get taxonomy data 
  count_cols <- sapply(otu_Vo, is.numeric)                                   # Only numeric columns  
  Tax_dat <- otu_Vo[(count_cols == FALSE)]                                   # Get taxa matrix -- non-numeric
  
  # merge with OTU_corr, tax_ranks, corr_ranks                               # Because Corr ranks dropped OTU rownames
  Tax_OTU_corr <-merge(OTU_corr, Tax_dat, all.x =TRUE)                       # Merge OTU corr, Tax_dat, keep all OTUs   # head(Tax_OTU_corr) dim(Tax_OTU_corr)
  Corr_ranks = merge(Tax_OTU_corr, corr_ranks)                               # Merge Tax_OTU_corr, corr_test   # dim(Corr_ranks)  #head(Corr_ranks)
  Corr_ranksU <- unique(Corr_ranks)                                          # Get unique after merge
  
  # Rearrange columns and rows
  corr_ranks_names <- colnames(corr_ranks)                                   # Get only corr_test names
  Corr_ranksUr <- Corr_ranksU[,corr_ranks_names]                             # Get only cols in corr_test
  Corr_ranksUr$OTU_r <- Corr_ranksU$OTU_r                                    # append OTU_r column
  Corr_ranksUr$OTU <- Corr_ranksU$OTU                                        # append OTU column 
  Corr_ranksUr <-Corr_ranksUr[order(Corr_ranksUr$Consensus.lineage),]        # Sort by Consensus.lineage
  Corr_ranksUr[is.na(Corr_ranksUr)]<-0                                   # fill NAs with 0s
  return(Corr_ranksUr)
}

# Test function to append OTUs
#OTU_ra<-append_OTU_corr_v(otu_V, corr_ranks, "CH4_logn1")
#dim(OTU_ra); head(OTU_ra)


####################################################
# 7) Test corr output by plotting corr heatmap
####################################################   

# Note that this is not using OTU level data in 6, too many rows, will crash
#Tax_corr <- corr_table
#Tax_corr <- corr_ranks                                             # substitute input
# Tax_corrU <-unique(Tax_corr)                                         # Make unique tax-corr matrix (here at Consensus.lineage level)

# Get only numeric cols for heatmap...
#r_list <-sapply(Tax_corrU,is.numeric)                    # get only numeric columns (correlations)
#Tax_corrU_r <-data.matrix(Tax_corrU[,r_list])            # get data matrix  

#library(gplots)
#options(repr.plot.width=12, repr.plot.height=5)
#heatmap.2(t(Tax_corrU_r[,2:6]), Rowv=F, Colv=F, key=TRUE, labCol=FALSE, trace="none", density.info="none", margins = c(4, 8), col = brewer.pal(11, "PiYG")) #rev(brewer.pal(11,"RdYlBu")))





