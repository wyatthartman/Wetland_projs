
# Import packagaes
suppressMessages(library(Hmisc))
library(ggplot2)
library(RColorBrewer)
library(reshape2)
suppressMessages(library(cowplot))
theme_set(theme_grey())              # important, ggplot2 default overwrites cowplot defaults, which add unwanted items to plots

# Read import module for testing convenience here, still a bit of mucking around in pre-processing
source("Import_SalOTU_dat_Plot_test_v0.1.R")

# Reimport color legend from barplots
Taxons3d <-read.table("Test_taxonomy_color_pallete.txt", sep="\t", header = T)
# head(Taxons3d)

# Separate Delta Sites (oligo and FW)
levels(Meta_iTag$SALTgroup)
Meta_iTag_FW <- Meta_iTag[Meta_iTag$SALTgroup =="FW",]
Meta_iTag_Oligo <- Meta_iTag[Meta_iTag$SALTgroup =="Oligo",]

Meta_iTag_Delta <-rbind(Meta_iTag_FW, Meta_iTag_Oligo)
# Meta_iTag_Delta

Delta_sites <-data.frame(Meta_iTag_Delta[,"Sample"])                 # DF of Delta sites
colnames(Delta_sites) [1] <- "Sample"                                # Rename column "Sample"
Delta_sites["Samp_Index"] <- seq(1:nrow(Delta_sites))                # Make sample index for reordering post-merges
#Delta_sites

# get OTU sample subset function
Samp_subset_OTU = function(otu_V, Samp_subset) { 
# Separate data and taxonomy columns
    Count_cols <- sapply(otu_V, is.numeric)                               # are columns numeric?
    Tax_cols <- otu_V[!Count_cols]                                        # Get Taxonomy columns  dim(Tax_cols)  # head(Tax_cols)
    otu_d <- data.frame(otu_V[Count_cols])                                # Get count data (Sample) columns   

    # Subset data columns by Sample subset, using merge on transp.        # Passing column list didn't work, why?
    otu_dT <-data.frame(t(otu_V[Count_cols]))                             # transp. counts and DF 
    otu_dT["Sample"] <- row.names(otu_dT)                                 # Make sample column 
    otu_SubsetT <- merge(Samp_subset, otu_dT)                             # Merge with subset of samples 
    otu_SubsetT <- otu_SubsetT[order(otu_SubsetT["Samp_Index"]),]         # Sort by Sample Index

    # Retranspose and cbind Taxonomy data
    row.names(otu_SubsetT) <- otu_SubsetT[,"Sample"]                      # replace rownames w. Sample names (pre-transp.)
    otu_SubsetT <- otu_SubsetT[3:ncol(otu_SubsetT)]                       # Drop first two cols (Sample, Samp_Index from merge)
    otu_Subset_d <-data.frame(t(otu_SubsetT))                             # Transpose to OTU table               #dim(otu_Delta_d) #head(otu_Delta_d)
    otu_Subset <- cbind(otu_Subset_d, Tax_cols)                           # add back Taxonomy data
    return(otu_Subset)
}

# Subset Delta smaples from OTU table
Delta_OTU <- Samp_subset_OTU(otu_V, Delta_sites)
#head(Delta_OTU)

## A) Recalculate corr METHODS for controller, using functions from correlations module:
# source("Corr_ranks_module_v0.3.2_strip.R")

# Only to genus level:
# CH4corrRanks <- corr_TaxRanks(otu_V, Meta_iTag, "CH4_logn1")                                  # head(CH4corrRanks)

# Append corrs at OTU level --  SLOW (10 min), need only for iTol? 
# OTU_ra <-append_OTU_corr_v(otu_V, CH4corr_ranks, "CH4_logn1")            

# For Delta samples subset
# CH4corrRanks_Delta_OTU <-append_OTU_corr_v(Delta_OTU, CH4corrRanks_Delta, "CH4_logn1")   
# head(CH4corrRanks_Delta_OTU)

# typically saved to file
# write.table(OTU_ra, "CH4_TaxRankOTU_corr_tab_all_sites.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
# write.table(CH4corrRanks_Delta_OTU, "CH4_TaxRankOTU_corr_tab_Delta_sites.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)

## B) Read EXISTING corr table: 
# CH4corrRanksOTU <-read.table("CH4_TaxRankOTU_corr_tab_all_sites.txt", sep="\t", header=T)      #head(CH4corrRanksOTU)

## C) Read EXISTING Delta corr table: 
# CH4corrRanksOTU_Delta <-read.table("CH4_TaxRankOTU_corr_tab_Delta_sites.txt", sep="\t", header=T)  #head(CH4corrRanksOTU_Delta)


# Aggregation test function working

agg_by_cat = function(otu_V, agg_var) {                                                                    
    count_cols <- sapply(otu_V, is.numeric)                              # Select numeric columns
    otu_counts <- data.matrix(otu_V[,count_cols])                        # Keep numeric cols, matrix for phyloseq  #dim(otu_d); head(otu_d)  # head(otu_V); head(otu_d)
    
    agg_fact <- unlist(as.list(otu_V[agg_var]))                          # Get levels for aggregation
    otu_agg <- aggregate(otu_counts, by=list(agg_fact), FUN=sum)         # Aggregate OTU table by agg_fact
    colnames(otu_agg)[1] <- agg_var                                      # rename new groups "Taxa"   # head(otu_agg)
    row.names(otu_agg) <- otu_agg[,1]                                    # set rownames
    return(otu_agg)               
}

# Test aggregation function
#otu_agg <- agg_by_cat(otu_V, "Taxonomy")                                # Aggregate OTU table by taxonomy
#head(otu_agg)

# Plots color bar or bargraph per sample by color, sorted
site_colbar = function(Meta_iTag, sample_var, color_var, color_table, order_var, plot){
 
    vars <-c(sample_var, color_var, order_var)                                                # Get sample info vars
    samp_info <- Meta_iTag[,vars]                                                             # Get sample info from metadata 
    samp_info$counts <-1                                                                      # Make column of 1s for uniform bar height
    sort_name = colnames(samp_info[3])
    
    # Merge color and samp info
    samp_info_color <- merge(samp_info, color_table)                                          
    samp_info_colorO <-samp_info_color[order(samp_info_color[,order_var]),]                    # Sort by index
    
    # reorder factors by order_var values
    samp_info_colorO[,sample_var] <- reorder(samp_info_colorO[,sample_var], samp_info_colorO[,order_var]) # reorder factor by index..
    samp_info_colorO[,color_var] <- reorder(samp_info_colorO[,color_var], samp_info_colorO[,order_var])   # reorder factor by index..

    # Get data to plot 
    site_colU <-unique(as.character(samp_info_colorO$color))                                  # unique color vector
    Sample <- samp_info_colorO[,sample_var]                                                   # Samples
    Group <- samp_info_colorO[,color_var]                                                     # Group IDs
    
    # make height a flat legend or barplot graph?   
    ifelse(plot=='graph', order_var <- samp_info_colorO[,order_var], order_var <-samp_info_colorO$counts)
      
    # Make barplot
    p <-ggplot(samp_info_colorO, aes(x=Sample, y=order_var, fill=Group)) +                     # Plot sample by group
          geom_bar(stat='identity',  width = 1, alpha = 0.8) + scale_fill_manual(values=c(site_colU)) +   # Manual color fill
          theme(legend.position="none", axis.title.x=element_blank(),                         # Remove legend
             axis.text.x=element_blank(), axis.ticks.x=element_blank()) +                     # Remove all x-axis labs
             labs(y = paste(sort_name)) + #axis.title.y=element_text("test"))
                
          theme(panel.background = element_blank(),                                           # Remove panel borders and grid lines  #
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()) #+
    
    ifelse(plot=='graph', pp <- p, p <- p + theme(axis.title.y=element_blank(),               # Ifelse to hide y axis
                                            axis.text.y=element_blank(), axis.ticks.y=element_blank()))#  + #)  # Remove y labs  
    p
}

# Test site_colbar function
#options(repr.plot.width=4, repr.plot.height=0.6)
#s <-site_colbar(Meta_iTag, "Sample", "Location", site_colors, 'EWsiteHyd_index',"graph")
#sb <-site_colbar(Meta_iTag, "Sample", "Location", site_colors, 'EWsiteHyd_index',"")
#s
#sb

# More examples
# options(repr.plot.width=4, repr.plot.height=1.5)
# site_colbar(Meta_iTag,"Sample", "Location",  site_colors, 'CH4_logn1', "graph")
# options(repr.plot.width=4, repr.plot.height=0.4)
# site_colbar(Meta_iTag,"Sample", "Location",  site_colors, 'CH4_logn1', "")
# options(repr.plot.width=4, repr.plot.height=1.5)
# site_colbar(Meta_iTag,"Sample", "Location",  site_colors, 'Salinity.x', "graph")
# options(repr.plot.width=4, repr.plot.height=0.4)
# site_colbar(Meta_iTag, "Sample", "Location", site_colors, 'Salinity.x', "")

# Updated to modify y label, minimally v3.3
Taxon_bar_newC = function(otu_agg, color_set, min_abund=0.005, color_no=9, output_palette=F, ylab="CPM") {
    
    # Uses first column as agg. var, passes name 
    agg_var <-otu_agg[,1]                     
    agg_name <-colnames(otu_agg)[1]
    
    # Add total and % total columns for filtering
    otu_agg$TaxTot <- rowSums(data.matrix(otu_agg[,-1]))                                      # Get Tax total counts across samples 
    otu_agg$TotFrac<- (otu_agg$TaxTot/sum(otu_agg$TaxTot))                                    # Get % Tax totals

    # New Taxa names, only > 1% (filter)                                                      # min_abund = 0.0075
    newTax <- ifelse(otu_agg$TotFrac > min_abund, paste(otu_agg[,agg_name]), paste("Other"))  # newTax, if < min, "Other"
    otu_agg_d <- otu_agg[,2:(ncol(otu_agg)-3)]                                                # Get only data, ignore added 3 cols        # OTU_VST_respM_aggD
    otu_agg2 <- aggregate(otu_agg_d, by=list(newTax), FUN=sum)                                # Agg. by new Taxa  
    colnames(otu_agg2)[1] <-colnames(otu_agg)[1]                                              # rename groups   # OTU_VST_respM_agg2  # dim(OTU_VST_respM_agg2)  # otu_agg2
    # TODO: "Other" last factor--first drop, then rbind, then reorder factor.                 # otu_agg2<-otu_agg2[otu_agg2[1] !="Other"]
    
    # melt aggregated data 
    otu_agg2_melt <-melt(otu_agg2)          
    colnames(otu_agg2_melt)[3] <- "CPM"             # Label values from agg for plotting, here CPM: counts per million

    # Make color pallete (automated, semi-custom)
    ncolor<-nrow(otu_agg2)                                                                    # n colors from Tax2 cats                 # ncolor
    cols <-colorRampPalette(brewer.pal(color_no, color_set))(ncolor)                          # Make custom pallete with n colors         # cols

    # OUTPUT Taxa color dataframe and get list
    Taxa_colors <- data.frame(otu_agg2[,1], "color" = cols)
    colnames(Taxa_colors)[1] <-agg_name

    # Make barplot
    fill_name <-paste(colnames(otu_agg2)[1])
    bar <-ggplot(data=otu_agg2_melt, aes(x=variable, y=CPM, fill=otu_agg2_melt[1])) +labs(fill=agg_name)                       # plot, fill by Taxa  
    bar_plot <- bar + geom_bar(stat="identity", alpha=0.85, width=1)+scale_fill_manual(values=cols) +   # barplot, manual colors  
          guides(fill = guide_legend(ncol=1))                                              +   # 1 col legend
          theme(axis.title.x = element_blank(), axis.text.x=element_blank(),                   # hide x - axis labels  
                axis.ticks.x=element_blank())                                              +   
          theme(legend.key.size = unit(0.4, "cm"))                                         +   # Shrink legend
          theme(panel.background = element_blank(),                                            # Remove panel borders and grid lines  #
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())                                    +   # bar_plotNL <-bar_plot + theme(legend.position="none")                                    # Remove legend entirely
     
    labs(y=ylab)                                                                               # relabel y axis 
    
    ifelse(output_palette==T, return(Taxa_colors), return(bar_plot))                   
      
}

# Test functions 2 & 3

# aggregate taxonomy by "Taxonomy"
#otu_agg <- agg_by_cat(otu_V, "Taxonomy")                            # Aggregate OTU table by taxonomy

#bar_plot = Taxon_bar_newC(otu_agg, min_abund=0.005, "Set1", ylab="x")
#bar_plot = Taxon_bar_newC(otu_agg, min_abund=0.005, "Paired")                           #bar_plot = Taxon_bar(otu_agg, min_abund=0.01)
#bar_plot = Taxon_bar(otu_agg, min_abund=0.005, output_palette=T)

#options(repr.plot.width=7, repr.plot.height=4)
#bar_plot

# Get only color Palette instead
#Tax_palette0 <- Taxon_bar_newC(otu_agg, "Set1", min_abund=0.005, output_palette=T)
#Tax_palette0

# Test on subset Archaea
#arch <- otu_V[otu_V$Kingdom=='Archaea',]                # subset     # dim(arch); head(arch)
#Arch_agg <- agg_by_cat(arch, "Family")            
#arch_plot = Taxon_bar_newC(Arch_agg, min_abund=0.005, "Set1", ylab="x")
#arch_plot
#arch_pal <-Taxon_bar_newC(Arch_agg, min_abund=0.005, "Set1", ylab="x", output_palette=T)
# arch_pal

Sort_new_palette = function(Tax_palette, Taxonomy, sort_var, otu_t) {            # Fixed in v.3.5 for "Other"

    # First remove "other" category and prepare merge
    Taxons <- as.character(Tax_palette[,Taxonomy])                               # Get only Taxonomy vector, as. character for list appending...
    Taxons2 <- data.frame(Taxons[Taxons !="Other"])                              # drop "Other" category (will add back)
    colnames(Taxons2)[1]<-Taxonomy
    #Taxons2 <- data.frame(Taxonomy = Taxons[Taxons !="Other"])                  # drop "Other" category (will add back)
    Taxons2$Taxons2 <- Taxons2[,Taxonomy]                                        # copy col with new name to generate "NA"s in merge 

    # Get Tax_data for sorting
    count_cols <- sapply(otu_t, is.numeric)                                      # get only numeric columns  
    Tax_dat <- otu_t[(count_cols == FALSE)]                                      # return only taxonomic info

    # Merge with Tax_dat to sort by                                                         
    Taxons3 <- merge(Taxons2, Tax_dat, all.y=T)                                  # Merge Taxons 2 with all Tax Data
    Taxons3s <- Taxons3[order(Taxons3[,sort_var]),]                              # Sort by sort_var, e.g. Consensus.lineage
    Taxons3u <-unique(as.character(Taxons3s$Taxons2))                            # Get list of taxa, not levels (...)
    Taxons3u <- Taxons3u[!is.na(Taxons3u)]                                       # drop NA
    Taxons3d <- data.frame(Taxons3 = (append(Taxons3u, "Other")))                # Append "Other"
    Taxons3d$Index <- seq(1:nrow(Taxons3d))                                      # Make Order index
    
    # Now have ordered levels, but problem if no "Other" in original pal.
    Other_test <- "Other" %in% Taxons                                            # Is "Other" in original pal
    ifelse(Other_test == TRUE, Taxons3d <- Taxons3d,                             # Yes, pass DF 
                                Taxons3d <-Taxons3d[Taxons3d$Taxons3 !="Other",])# No, drop new "Other" from merge 
    
    Taxons3d$color <-Tax_palette$color                                           # Add colors from original Pal.  
    colnames(Taxons3d)[1] <- Taxonomy                                            # Relabel "Taxonomy" column
    return(Taxons3d)
 
}
    

# Test sort_new_palette function
#Taxons3d <- Sort_new_palette(Tax_palette0, "Taxonomy", "Consensus.lineage", otu_V)
#head(Taxons3d); #Taxons3d

# test on Arch; Arch_sorted
#Arch_pal_3d <- Sort_new_palette(Arch_sorted, "Family", "Consensus.lineage", arch)
#Arch_pal_3d <- Sort_new_palette(arch_pal, "Family", "Consensus.lineage", arch)
# head(Arch_pal_3d); dim(arch_pal); #Arch_pal_3d

# EXPORT color pallete
#write.table(Taxons3d, "Test_taxonomy_color_pallete.txt", sep="\t", col.names=T, row.names=F, quote=T)

ColorMap_other = function(Taxons3d, otu_V){  #, tax_sort = "Consensus.lineage"){ # Revised v.3.5
    
    rankvar <-colnames(Taxons3d)[1]                                             # Get taxonomic rank from Palette -v3.3
    
    # Drop "Other" rows and merge with OTU table
    Taxons3dx <- Taxons3d[1:nrow(Taxons3d)-1,]                                  # drop "Other" entirely, assuming last col   #Taxons3dx
    TaxCol_otu <- merge(Taxons3dx, otu_V, all.y=TRUE)                           # Merge with otu_table, while keeping all of OTU table.  #Tax_col_otu 

    # get last Tax order and color, corresponding to "Other"
    Indx <- Taxons3d[,'Index']                                                  # get color pal index
    colr <- Taxons3d[,'color']                                                  # get color pal colors  
    Taxa <- Taxons3d[,rankvar]
    last_order <- Indx[nrow(Taxons3d)]                                          # get last index number, for "Other"
    last_color <- colr[nrow(Taxons3d)]                                          # get last color, for "Other"      # last_order
    last_Taxa <- as.character(Taxa[nrow(Taxons3d)])
    
    # Test for other before NA replace
    Other_test <- "Other" %in% Taxons3d[,rankvar]                               # Is "Other in original pal
    ifelse(Other_test == TRUE, Other_var <-"Other", Other_var <- last_Taxa)     # Generate Other_var 
    
    # Replace NAs in Order, Tax, and Color, by ifelse, recode new vars: 
    TaxCol_otu["newOrder"] <- ifelse(is.na(TaxCol_otu[,"Index"] == TRUE), last_order, TaxCol_otu[,'Index']) 
    TaxCol_otu["newTax"] <- ifelse(is.na(TaxCol_otu[,"Index"] == TRUE), paste(Other_var), paste(TaxCol_otu[,rankvar]))    # newTax, if < min, "Other"
    TaxCol_otu["newColor"] <- ifelse(is.na(TaxCol_otu[,"color"] == TRUE), paste(last_color), paste(TaxCol_otu[,'color']))

    # Reorder factor levels by Index
    TaxCol_otu["newTax"] <- reorder(TaxCol_otu[,"newTax"], TaxCol_otu[,"newOrder"])       # reorder factors by index..
    TaxCol_otu["newColor"] <- reorder(TaxCol_otu[,"newColor"], TaxCol_otu[,"newOrder"])   # ...needed below...
    TaxCol_otu$newOrder <- reorder(TaxCol_otu$newOrder, TaxCol_otu$newOrder)              # index as factor ""

    TaxCol_otu <- TaxCol_otu[order(TaxCol_otu[,"Index"]),]                                # sort by index     
    return(TaxCol_otu)
}

# head(Taxons3d)
#head(Arch_pal_3d)

#cDA_pal_S

# Test function
#otuTaxCol <-ColorMap_other(Taxons3d, otu_V)              #head(otuTaxCol)  
# otuTaxCol <-ColorMap_other(Taxons3d, otu_agg)
#dim(otuTaxCol); 
#head(otuTaxCol) # colnames(otuTaxCol)

# Test on corr filtered subset (here by corr, not produced until below)
#r_filt_test <- agg_by_cat(r_Filt_OTU_pos, "Taxonomy")
#r_filt_test

# Filt_otuTaxCol <-ColorMap_other(Taxons3d, r_filt_test)
# Filt_otuTaxCol
# r_filt2

# Test on recolored archaea, not sure why problem...
#otuTaxCol_A <-ColorMap_other(Arch_pal_3d, arch)
#dim(otuTaxCol_A);  
#head(otuTaxCol_A) #colnames(otuTaxCol_A)

#ColorMap_other(cDA_pal_S, corr_Delta_Arch)

barplot_colorsIn = function(TaxCol_otu, Taxons3d, ylab="CPM") {                 # Barplot using existing colors

    # Get only count data from OTU table
    Count_cols <- sapply(TaxCol_otu, is.numeric)                                # are columns numeric?
    TaxCol_otu_d <-TaxCol_otu[Count_cols]                                       # get df with only numeric columns
    TaxCol_otu_d <-TaxCol_otu_d[,-1]                                            # drop 1st col, was added index for cats.              #TaxCol_otu_d #TaxCol_otu
    TaxCol_otu_d <-data.matrix(TaxCol_otu_d[,1:ncol(TaxCol_otu_d)-1])           # drop last col, was newOrder

    # Aggregate count data by newTax, melt
    Tax_n <-TaxCol_otu[,"newTax"]                                               # Get column to aggregate on
    Tax_agg <- aggregate(TaxCol_otu_d, by = list(Tax_n), FUN=sum)               # Aggregate by category 
    colnames(Tax_agg)[1] <- colnames(Taxons3d)[1]                               # Rename category    #Tax_agg
    Tax_agg_melt <-melt(Tax_agg)                                                # melt aggregated data    
    colnames(Tax_agg_melt)[3] = "CPM"                                           # Label values from agg for plotting, here CPM: counts per million

    # OUTPUT Taxa colors and legend name
    cols <-levels(TaxCol_otu[,"newColor"])                                      # get colors for plot
    fill_name <- colnames(Tax_agg_melt)[1]                                      # get legend title 
    fill_vals <- Tax_agg_melt[1]                                                # get legend category names

    # Make barplot
    bar <-ggplot(data=Tax_agg_melt, aes(x=variable, y=CPM, fill=fill_vals)) +labs(fill=fill_name)    # plot, fill by Taxa  
    bar_plot <- bar + geom_bar(stat="identity", alpha=0.85, width=1) +scale_fill_manual(values=cols) +   # barplot, manual colors 
            guides(fill = guide_legend(ncol=1))                                                  +   # 1 col legend
            theme(axis.title.x = element_blank(), axis.text.x=element_blank(),                       # hide x - axis labels  
                  axis.ticks.x=element_blank())                                                  +   
            theme(legend.key.size = unit(0.4, "cm"))                                             +   # Shrink legend
            theme(panel.background = element_blank(),                                                # Remove panel borders and grid lines  #
                 panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())                                            +   # bar_plotNL <-bar_plot + theme(legend.position="none")                                    # Remove legend entirely
            labs(y=ylab)                                                                             # change y label

    return(bar_plot)
} 

#options(repr.plot.width=7, repr.plot.height=4)
#barplot_colorsIn(otuTaxCol, Taxons3d) #, "newTax", "newColor")
#barplot_colorsIn(otuTaxCol_A, Arch_pal_3d)

# Test on correlation filtered subsets, not defined until below
# options(repr.plot.width=7, repr.plot.height=4)
# barplot_colorsIn(Filt_otuTaxCol) #, "newTax", "newColor")
# Filt_otuTaxCol

## Combined barplot with input colors function 
recolor_barplot = function(otu_V, Taxons3d, ylab="CPM"){
    otuTaxCol <-ColorMap_other(Taxons3d, otu_V)          # map colors to otu table
    barplot_colorsIn(otuTaxCol, Taxons3d, ylab)                          # make recolored barplot from color x OTUs 

}

# Test recolor barplot function
#options(repr.plot.width=7, repr.plot.height=4)
#recolor_barplot(otu_V, Taxons3d, ylab="perc")

# Use unsorted instead, original palette instead
#Tax_palette0$Index<-seq(1:nrow(Tax_palette0))            # need to add Index for function to work

#recolor_barplot(otu_V, Tax_palette0)                     # plot using unsorted palette

#test with Archaea
# recolor_barplot(arch, Arch_pal_3d, ylab="perc")
# recolor_barplot(arch_agg, Arch_pal_3d, ylab="perc")
#recolor_barplot(Arch_sorted, Arch_pal_3d, ylab="perc")
#Arch_sorted

# 5) Make new barplot taxa color scheme    #wraps fxns . 1-4

new_tax_palette_plot = function(otu_V, agg_var, sort_var, col_set, min_abund = 0.005, output_pal=F, ylab="CPM"){
    otu_agg <- agg_by_cat(otu_V, agg_var)                                                   # Fxn.1 
    Tax_palette <- Taxon_bar_newC(otu_agg, min_abund = 0.005, col_set, output_palette=T)   # Fxn.2  
    Taxa_colors <- Sort_new_palette(Tax_palette, agg_var, sort_var, otu_V)                        # Fxn.3
    y_lab=substitute(ylab)                                                                      # pass y-axis label
    newbar <- recolor_barplot(otu_V, Taxa_colors, ylab=y_lab)                               # Fxn.4(c)
    ifelse(output_pal==T, return(Taxa_colors), return(newbar))  
}

# Test new_tax_palette_plot function
#newpal <- new_tax_palette_plot(otu_V, "Taxonomy", "Consensus.lineage", "Set1") #, ylab="perc")
#newpal <- new_tax_palette_plot(otu_V, "Taxonomy", "Consensus.lineage", "YlOrRd")
#newpal <- new_tax_palette_plot(otu_V, "Taxonomy", "Consensus.lineage", "Paired")
#newpal <- new_tax_palette_plot(otu_V, "Taxonomy", "Consensus.lineage", "Set1", output_palette=T)

#options(repr.plot.width=7, repr.plot.height=4)
#newpal

# test on Archaea
#newpala <- new_tax_palette_plot(arch, "Family", "Consensus.lineage", "Set1") #, ylab="perc")
#newpala

sort_otu_agg_by_meta = function(otu_F, Meta_iTag, Sample, Samp_index) {     # revised in 2.4
# assumes column names in OTU table correspond to "Sample" column in MetaDB

    OTU_samps <- data.frame(Sample=colnames(otu_F))                         # Get Sample names from OTU column names          # I#OTU_samps
    Meta_OTU <- merge(Meta_iTag, OTU_samps, by=Sample)                      # Merge site order and Samples, drops non-matching                  
    rownames(Meta_OTU) <- Meta_OTU[,Sample]                                 # Set row names as Sample 

    # Get site sorting data matrix
    sortr <- c(Samp_index, Sample)                                          # site, index name list
    Site_sort <-Meta_OTU[,sortr]                                            # slice DF by list 
    colnames(Site_sort) <- c('Site_index', Sample)                          # rename columns   -- needs to be generalized #Site_sort
    dim(Site_sort)

    ### Reorder Samples in OTU table      
    otu_f <- otu_F                                                          # copy df to keep orig intact--used later 
    rownames(otu_f) <- otu_f[,1]                                            # make rownames Taxa 
    otu_f <-otu_f[,-1]                                                      # drop Taxa col
    otu_ft <- t(otu_f)                                                      # transpose OTU table
    otu_fT <- data.frame(Sample=row.names(otu_ft), otu_ft)                  # add Sample col

    otu_F_metaT <- merge(Site_sort, otu_fT, by=Sample, all=TRUE)            # Merge with site sorter, keep all
    otu_F_metaT <- otu_F_metaT[order(otu_F_metaT[,"Site_index"]),]          # Sort samples by meta index
    otu_F_metaT <- otu_F_metaT[,-2]                                         # Drop Site index column
    row.names(otu_F_metaT) <- otu_F_metaT[,Sample]                          # add back rownames
    
    otu_F_reorder <- data.frame(t(otu_F_metaT[,-1]))                        # transpose while drop new samp column
    otu_F_reorder2 <- cbind(rownames(otu_F_reorder), data.frame(otu_F_reorder, row.names=NULL))   # add rownames as col.
    colnames(otu_F_reorder2)[1]<-colnames(otu_F)[1]                         # name new col after input label from agg. 
   
    
    return(otu_F_reorder2)
    
}

## Test function
#otu_agg <- agg_by_cat(otu_V, "Taxonomy")                                    # Aggregate OTU table by taxonomy
#sorted <- sort_otu_agg_by_meta(otu_agg, Meta_iTag, "Sample", "CH4_logn1")
#sorted <- sort_otu_agg_by_meta(otu_agg, Meta_iTag, "Sample", "EWsiteHyd_index")
#sorted <- sort_otu_agg_by_meta(otu_agg2, Meta_iTag, "Sample", "CH4_ug_m2_h")
# head(sorted)

# Test on Archaea
# arch_agg <- agg_by_cat(arch, "Family")                                    # Aggregate OTU table by taxonomy
# Arch_sorted <- sort_otu_agg_by_meta(arch_agg, Meta_iTag, "Sample", "CH4_logn1")
# head(Arch_sorted)

# arch_ag2 <-agg_by_cat(Arch_sorted, "Family")
# arch_ag2

# Function for percent abundance data
perc_abund = function(otu_t) {
    count_cols <- sapply(otu_t, is.numeric)                               # are columns numeric?
    counts <- otu_t[count_cols]                                           # Get numeric columns 
    count_perc <-sweep(counts*100, 2, colSums(counts),"/")                # Get percent abund. (100* counts/colSums)
    taxes <- otu_t[!count_cols]                                           # get taxonomy data 
    count_perc_tax <- cbind(count_perc, taxes)                            # paste data, taxa.
    return(count_perc_tax)
}

# test function
V_perc <- perc_abund(otu_V)                         
# dim(otu_V); dim(V_perc); head(V_perc)

#V_nums <- sapply(V_perc, is.numeric)                               # are columns numeric?
#V_counts<-V_perc[V_nums]                                           # get counts
#max(V_counts); #colSums(V_counts)                                 # Max numeric # check colsums

# 3rd attempt v.3.3 -- make new pal and sort first, then plot with fixed colors
abund_loc_barplot = function(otu_tab, agg_var, Meta_iTag, order_var, plot="", colors, color_set, relab="CPM"){    # Revised version 0.2.6    
    
    ifelse(relab == "CPM", otu_tab <-otu_tab, otu_tab <- perc_abund(otu_tab))          # if not CPM, apply perc_abund
    y_lab=substitute(relab)
    
    otu_agg <- agg_by_cat(otu_tab, agg_var)                                            # Agg by cat (Fxn. 1) on agg_var
    otu_S <- sort_otu_agg_by_meta(otu_agg, Meta_iTag, "Sample", order_var)             # Sort agg by (Fxn. n) on agg
    
    
    # If color palette not specified, make new & sort
    ifelse(colors=="",                                                                 
           colors_p <-Taxon_bar_newC(otu_agg, min_abund=0.005, color_set, output_palette=T), 
           colors <- colors)
    ifelse(colors=="",
        colors <- Sort_new_palette(colors_p, agg_var, "Consensus.lineage", otu_tab),
                                  colors <- colors)
    
    # make recolored barplot, separate plot and legend
    bar_plot <- recolor_barplot(otu_S, colors, ylab=y_lab)                             # reuse existing palette 
    bar_plotNL <- bar_plot + theme(legend.position="none")                             # drop legend from ggplot
    legend <-get_legend(bar_plot)                                                      # get legend as separate obj.
    
    # Make site_barplot for top
    plot_s <-substitute(plot)
    sb <-site_colbar(Meta_iTag,"Sample", "Location", site_colors, order_var, plot_s)
    
    # make plot grid of site colors and abundance barplot
    ifelse(plot=="graph", proportion <- c(1/4,3/4), proportion <- c(1/8,7/8))          # set proportion for graph
    pg <- plot_grid(sb, bar_plotNL, align="v", nrow=2, rel_heights = proportion)       # plot grid, by proportion 
    pgl <- plot_grid(pg, legend, rel_widths = c(8,3))                                  # add legend to grid, better alignment later...align="h", ncol=2, 

pgl
}

# Test combined abundane and site sorted barplot function
options(repr.plot.width=6, repr.plot.height=4)

Site_sort <- abund_loc_barplot(otu_V, "Taxonomy", Meta_iTag, "EWsiteHyd_index", "", Taxons3d, relab="% Abundance")
Site_sort

CH4_sort_taxa <-abund_loc_barplot(otu_V, "Taxonomy", Meta_iTag, 'CH4_logn1', "graph", Taxons3d)
CH4_sort_taxa

# CH4_sort_taxa2 <-abund_loc_barplot(otu_V, "Taxonomy", Meta_iTag, 'CH4_logn1', "graph", "", color_set="YlGn")
# CH4_sort_taxa2

# Other examples / test plots
# Site_sort <- abund_loc_barplot(otu_V, "Taxonomy", "EWsiteHyd_index", "")
# Site_sort <- abund_loc_barplot(otu_V, "Taxonomy", "EWsiteHyd_index", "", Taxons3d)
# Site_sort

# Site_sort2 <- abund_loc_barplot(otu_V, "Phylum", "EWsiteHyd_index", "", "")
# Site_sort2

# CH4_sort_taxa <-abund_loc_barplot(otu_V, "Taxonomy", 'CH4_ug_m2_h', "graph")
# CH4_sort_taxa <-abund_loc_barplot(otu_V, "Taxonomy", 'CH4_logn1', "graph", Taxons3d)
# CH4_sort_taxa # "CH4_ug_m2_h""CH4_logn1"

## Test subset of OTU table                                          # head(otu_V)
#arch <- otu_V[otu_V$Kingdom=='Archaea',]                # subset     # dim(arch); head(arch)

#Arch_by_CH4 <-abund_loc_barplot(arch, "Family", Meta_iTag, "CH4_logn1", "graph", "", "Accent")
#Arch_by_CH4_perc <-abund_loc_barplot(arch, "Family", Meta_iTag, "EWsiteHyd_index", "", "", "Accent", relab="% Total")

# Arch_by_CH4 <-abund_loc_barplot(arch, "Family", Meta_iTag, "CH4_logn1", "graph", "", "Dark2")
# Arch_by_CH4_perc <-abund_loc_barplot(arch, "Family", Meta_iTag, "CH4_logn1", "graph", "", "Dark2", relab="% Total")

#options(repr.plot.width=7, repr.plot.height=4)
#Arch_by_CH4
#Arch_by_CH4_perc

# Function to get corr. filtered otu table from corrRanks object

OTU_corr_filt_OTU = function(otu_V, corrRanks, chunk_var, r_cut){

    # Filter corrRanks by abs(chunk_var) > r_cut 
    chunk_r <- paste0(chunk_var,"_r")                               # Get colname tax rank corr. var., here [chunk] 
    abs_chunk_r <-paste0("abs_", chunk_r)                           # Get colname for abs of corr(chunk) 
    corrRanks[abs_chunk_r] <-abs(corrRanks[,chunk_r])               # Get abs val of chunk
    corrRanks_F <- corrRanks[corrRanks[,abs_chunk_r] > r_cut,]      # Get corr filtered corrRanks data 

    # Make DF of corr. filtered Taxa, corr. r values, OTU_ID
    r_Filt_rank <-data.frame(corrRanks_F[,chunk_var])               # Taxa names 
    colnames(r_Filt_rank)[1] <- chunk_var                           # Label Taxa col
    r_Filt_rank[chunk_r]<-corrRanks_F[,chunk_r]                     # Get correlation r values
    r_Filt_rank["OTU"] <- corrRanks_F[,"OTU"]                       # Add OTU column for merge below...      # dim(r_Filt_rank); head(r_Filt_rank)

    # Merge with abundance table 
    otu_Vo <- otu_V                                                 # copy otu table
    otu_Vo["OTU"]<-row.names(otu_Vo)                                # add otu_IDs as col.
    r_Filt_OTU <- merge(otu_Vo, r_Filt_rank)                        # Merge with abundance 
    return(r_Filt_OTU)                                                                                  # dim(r_Filt_OTU); head(r_Filt_OTU)
}


# Test OTU_corr_filt function
CorrFilt_OTU <- OTU_corr_filt_OTU(otu_V, CH4corrRanksOTU, "Genus", 0.5)
# dim(Corr_filt_test_OTU); head(Corr_filt_test_OTU_Delta)

CorrFilt_OTU_Delta <- OTU_corr_filt_OTU(Delta_OTU, CH4corrRanksOTU_Delta, "Genus", 0.5)
# dim(Corr_filt_test_OTU_Delta); head(Corr_filt_test_OTU_Delta)

# map corr_filt to barplot
corrFilt2barplot = function(r_Filt_OTU, tax_cut_r, agg_var, metaData, sort_var, graph, colors, color_set) {

    chunk_r <- paste0(tax_cut_r,"_r")   # chunk_r is column name of tax cut level to cutoff on r

    # Filter by [+ / -]  corr
    r_Filt_OTU_pos <- r_Filt_OTU[r_Filt_OTU[,chunk_r] > 0,]                 # filter positive corr
    r_Filt_OTU_pos <- r_Filt_OTU_pos[,1:ncol(r_Filt_OTU_pos)-1]             # drop corr col, numeric interference below

    r_Filt_OTU_neg <- r_Filt_OTU[r_Filt_OTU[,chunk_r] < 0,]                 # filter negative corr
    r_Filt_OTU_neg <- r_Filt_OTU_neg[,1:ncol(r_Filt_OTU_neg)-1]      

    # Make separte barplots
    Site_sort_plus <- abund_loc_barplot(r_Filt_OTU_pos, agg_var, metaData, sort_var, graph, colors, color_set)
    Site_sort_minus <- abund_loc_barplot(r_Filt_OTU_neg, agg_var, metaData, sort_var, graph, colors, color_set)

    # combine barplots
    # options(repr.plot.width=11, repr.plot.height=4)
    Corr_plus_minus_pg <- plot_grid(Site_sort_plus, Site_sort_minus, rel_widths = c(1,1))   
    #Corr_plus_minus_pg
}

# CFO_p <- perc_abund(CorrFilt_OTU)

# head(CorrFilt_OTU)

test_double_bar <- corrFilt2barplot(CorrFilt_OTU, "Genus", "Family", Meta_iTag, "CH4_logn1", "graph","","Set1")
#test_double_bar <- corrFilt2barplot(Corr_filt_test_OTU_Delta, "Genus", "Taxonomy", Meta_iTag_Delta, "CH4_logn1", "graph", Taxons3d)
options(repr.plot.width=11, repr.plot.height=4)
test_double_bar

# Delta sites 
#by_Tax <- corrFilt2barplot(CorrFilt_OTU_Delta, "Genus", "Taxonomy", Meta_iTag_Delta, "CH4_logn1", "graph", Taxons3d)
#by_Class <- corrFilt2barplot(CorrFilt_OTU_Delta, "Genus", "Class", Meta_iTag_Delta, "CH4_logn1", "graph", "", "Accent")

#options(repr.plot.width=11, repr.plot.height=4)
#by_Tax
#by_Class

# Delta sites corr bar by rank
#options(repr.plot.width=14, repr.plot.height=5.25)
#by_Family <- corrFilt2barplot(CorrFilt_OTU_Delta, "Genus", "Family", Meta_iTag_Delta, "CH4_logn1", "graph", "", "Paired")
#by_Genus <- corrFilt2barplot(CorrFilt_OTU_Delta, "Genus", "Genus", Meta_iTag_Delta, "CH4_logn1", "graph", "", "Set2")

#by_Family
#by_Genus

#options(repr.plot.width=17, repr.plot.height=5.25)
#by_Class <- corrFilt2barplot(CorrFilt_OTU, "Genus", "Family", Meta_iTag, "CH4_logn1", "graph", "", "Accent")
#by_Class

# Make subset of corrFilt data
# Test subset --ARCHAEA (would function and list)
#corr_Delta_Arch <-CorrFilt_OTU_Delta[CorrFilt_OTU_Delta["Kingdom"]=="Archaea",]

#head(CorrFilt_OTU_Delta)
#head(corr_Delta_Arch)


# plot with corrFilt2barplot
#options(repr.plot.width=11, repr.plot.height=4)
#corrFilt2barplot(CorrFilt_OTU_Delta, "Genus", "Genus", Meta_iTag_Delta, "CH4_logn1", "graph", "", "Set2")
#Ar_by_Genus <- corrFilt2barplot(corr_Delta_Arch, "Genus", "Genus", Meta_iTag_Delta, "CH4_logn1", "graph", "", "Set1")
#Ar_by_Genus

# Single plot with abund loc barplot
#options(repr.plot.width=6, repr.plot.height=4)
#Arch_by_CH4_corr <-abund_loc_barplot(corr_Delta_Arch, "Genus", Meta_iTag_Delta, "CH4_logn1", "graph", "", "Accent")
#Arch_by_CH4_corr

#Arch_by_CH4_perc_corr <-abund_loc_barplot(corr_Delta_Arch, "Genus", Meta_iTag_Delta, "EWsiteHyd_index", "", "", "Accent", relab="% Total")
#Arch_by_CH4_perc_corr


