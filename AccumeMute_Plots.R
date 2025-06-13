library(tidyverse)
library(data.table)
library(ggExtra)
library(MutationalPatterns)
library(scales)
library(ggpubr)
library(patchwork)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MutationalPatterns")

filepath<-"your/filepath/here/"
filepath<-"~/Downloads/"

###
###Accumulated Consensus mutations
###


accume_mutes<-read.csv(paste0(filepath,"Supplementary_Table2_AccumeMutes_EXTERNAL.csv") )



###Annotation for consensus mutations
annot<-read.csv(paste0(filepath,"Gene_Positions.csv") )
annot<-annot[annot$gene_feature!="ORF9b",] #Exclude ORF9b, overlapping ORFs

annotate<-function(mutes,annot){
  #Adjust mutation format for conversion to nt POS
  mutes$gene<-gsub("([A-Z]+\\d*.*):.*","\\1",mutes$aaSubstitutions)
  mutes$aa_pos<-as.numeric(gsub("[A-Z]+\\d*.*:[A-Z](\\d+)\\D","\\1",mutes$aaSubstitutions))
  mutes<-mutes[mutes$gene!="ORF9b",]
  
  #AA Substitution -> nt POS Conversion
  mutes_annot<-left_join(mutes,annot,by=c("gene" = "gene_feature"))
  #Get the overall codon/nt position for each AA change
  #Subtract 1 from aa_pos because we want first nt position of codon
  mutes_annot$nt_pos<-mutes_annot$start_pos + (mutes_annot$aa_pos-1)*3
  
  return(mutes_annot)
}


accume_mutes<-annotate(accume_mutes,annot)

#Reorder levels to match genome structure
lev<-c("ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
accume_mutes$gene<-factor(accume_mutes$gene,levels=lev)



###
###Finding significantly enriched bins of mutations
###
bins<-seq(0,30000,by=500)

n_bins<-length(bins) 
tot_snps<-dim(accume_mutes)[1] #Total number of SNPs observed

#Empty p-value dataframe
sig_df<-data.frame(effect=c(),pval=c(),significant=c(),binmin=c(),binmax=c())

total_prob<-0
for(i in 1:(length(bins)-1) ){
  # assume uniform prob, adjusted by bin size
  n_snps<-sum(accume_mutes$nt_pos>=bins[i] & accume_mutes$nt_pos<bins[i+1])
  bin_size <- bins[i+1] - bins[i]
  expected_prob <- bin_size / max(bins) #For 60 bins, 1/60th of mutations
  total_prob <- total_prob + expected_prob
  
  pvalue<-binomial_test(expected_prob,tot_snps,n_snps)
  pvalue$binmin<-bins[i]
  pvalue$binmax<-bins[i+1]
  pvalue$n_snps<-n_snps
  
  sig_df<-rbind(sig_df,pvalue)
  #print(pvalue)
}

sig_df$p_adj<-p.adjust(sig_df$pval,method="bonferroni",n=length(sig_df$pval))
sig_df$significant<-ifelse(round(sig_df$p_adj,digits=3)<0.05,"*","")
sig_df



#Standalone histogram w/significance
accume_mutes %>%
  ggplot(aes(nt_pos)) + geom_histogram(fill=ifelse(sig_df$p_adj<0.05,"orange","grey"),breaks=bins,size = 0.7)


#####
#####Create final plots
#####

###Function for Lollipop Plot + marginal histogram w/significance
recur_plot_hist<-function(accume){

  #Lollipop Plot
  accum_plot<-accume %>%
    ggplot(aes(x=nt_pos,colour=gene)) +
    stat_count(geom="segment",colour="grey",yend=0) +
    stat_count(aes(size=after_stat(count)),geom="point",show.legend=T) +
    scale_size(range = c(1, 2),guide='none') +
    theme_minimal() + 
    scale_y_continuous("Number of Infections",limits = c(0,NA),expand=c(0,0,0.05,0)) + 
    scale_x_continuous("Consensus Mutation Position",limits=c(0,30000)) + #,expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3) ) +
    labs(color = "Genomic Region",tag="A") + theme(plot.tag=element_text(face="bold")) +
    scale_color_manual(values=c(hue_pal()(11),"#808080"),drop=F) + #default color pallet, but don't drop levels
    guides(linetype = "none")
  
  #Grab legend to add later, remove for now
  legendary<-get_legend(accum_plot)
  accum_plot<-accum_plot + theme(legend.position = "none")
  
  #Manually create marginal plot
  marg<-ggplot(accume,aes(x=nt_pos)) +
    geom_histogram(colour="black",fill=ifelse(sig_df$p_adj<0.05,"orange","grey"),breaks=bins,size = 0.7) +
    theme_void() + 
    scale_x_continuous(limits=c(0,30000)) + scale_size(range = c(1, 2),guide='none') + 
    geom_text(
      stat = "bin",label = ifelse(sig_df$p_adj<0.05,row.names(sig_df),""),
      color="white",size=2.5,fontface="bold",
      position=position_stack(vjust = 0.9), breaks = bins
    )
  
  full_plot<-marg + plot_spacer() + accum_plot + plot_layout(heights = c(1,-0.85,5),axes="collect_x")
  
  return(list(legendary,full_plot))
}

legendary<-recur_plot_hist(accume_mutes)[[1]]
accume_plot<-recur_plot_hist(accume_mutes)[[2]]
accume_plot



###
###iSNV plots next
###

accume_isnvs<-read.csv(paste0(filepath,"Supplementary_Table3_iSNVs_EXTERNAL.csv") )


###Convert iSNV AA position back to first nt pos of codon

#Grab ranges for gene_features, new annotation scheme
accume_isnvs<-accume_isnvs %>% left_join(annot,join_by(nt_position>=start_pos,nt_position<=end_pos)) %>% 
  select(-gene_feature.y)
#Remove extra gene columns
accume_isnvs$gene_feature<-accume_isnvs$gene_feature.x
accume_isnvs<-accume_isnvs %>% select(-gene_feature.x)


#Convert amino acid position BACK to nucleotide pos
accume_isnvs$nt_pos<-accume_isnvs$start_pos + (accume_isnvs$aa_pos-1)*3

#Reorder levels to match genome structure
lev<-c("ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
accume_isnvs$gene_feature<-factor(accume_isnvs$gene_feature,levels=lev) #levels=levels(accume_isnvs$gene_feature)[c(4:5,11,6,1:2,7:10,3)])



###
###Calculate Binomial Test for iSNVs
###
#Establish bins
bins<-seq(0,30000,by=500)


n_bins<-length(bins)-1
tot_snps<-dim(accume_isnvs)[1]
sig_df<-data.frame(effect=c(),pval=c(),significant=c(),binmin=c(),binmax=c())

total_prob<-0
for(i in 1:(length(bins)-1) ){
  # assume uniform prob, adjusted by bin size
  n_snps<-sum(accume_isnvs$nt_pos>=bins[i] & accume_isnvs$nt_pos<bins[i+1])
  bin_size <- bins[i+1] - bins[i]
  print(bin_size)
  expected_prob <- bin_size / max(bins)
  total_prob <- total_prob + expected_prob
  
  pvalue<-binomial_test(expected_prob,tot_snps,n_snps)
  print(bins[i]);print(n_snps)
  pvalue$binmin<-bins[i]
  pvalue$binmax<-bins[i+1]
  pvalue$n_snps<-n_snps
  
  sig_df<-rbind(sig_df,pvalue)
  #print(pvalue)
}

sig_df$p_adj<-p.adjust(sig_df$pval,method="bonferroni",n=length(sig_df$pval))
sig_df$significant<-ifelse(round(sig_df$p_adj,digits=3)<=0.05,"*","")
sig_df


#Standalone histogram w/significance
accume_isnvs %>%
  ggplot(aes(nt_pos)) + geom_histogram(fill=ifelse(round(sig_df$p_adj,digits=3)<=0.05,"orange","grey"),breaks=bins,size = 0.7)


###Function for Lollipop Plot + marginal histogram w/significance
recur_plot_hist<-function(isnvs){
  #Lollipop Plot
  isnv_plot<-isnvs %>%
    ggplot(aes(x=nt_pos,colour=gene_feature)) +
    stat_count(geom="segment",colour="grey",yend=0) +
    stat_count(aes(size=after_stat(count)),geom="point",show.legend=T) +
    scale_size(range = c(1, 2),guide='none') +
    theme_minimal() + 
    scale_y_continuous("Number of Infections",breaks= c(1:4),limits = c(0,4),expand=c(0,0,0.05,0)) + 
    scale_x_continuous("iSNV Position",limits=c(0,30000)) + #,expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3) ) +
    labs(color = "Genomic Region",tag="B") + theme(plot.tag=element_text(face = "bold")) +
    scale_color_manual(values=c(hue_pal()(11),"#808080"),drop=F) +
    guides(linetype = "none")
  
  legendary<-get_legend(isnv_plot)
  isnv_plot<-isnv_plot + theme(legend.position = "none")
 
  #Manually create marginal plot
  marg<-ggplot(isnvs,aes(x=nt_pos)) +
    geom_histogram(colour="black",fill=ifelse(round(sig_df$p_adj,digits=3)<=0.05,"orange","grey"),breaks=bins,size = 0.7) +
    theme_void() + 
    scale_x_continuous(limits=c(0,30000)) + scale_size(range = c(1, 2),guide='none') +
    geom_text(
      stat = "bin",label = ifelse(sig_df$p_adj<0.05,row.names(sig_df),""),
      color="white",size=2.5,fontface="bold",
      position=position_stack(vjust = 0.9), breaks = bins
    )
  
 
  full_plot<-marg + plot_spacer() + isnv_plot + plot_layout(heights = c(1,-0.85,5),axes="collect_x")
  
  
  return(list(legendary,full_plot))
}


legendary<-recur_plot_hist(accume_isnvs)[[1]]
isnv_plot<-recur_plot_hist(accume_isnvs)[[2]]
isnv_plot

#Stack the plots
plots <- ggarrange(accume_plot,isnv_plot, nrow = 2)
#Add legend to stacked plots
ggarrange(plots,legendary,widths = c(0.825, 0.175))





