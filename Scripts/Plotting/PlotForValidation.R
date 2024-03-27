#!/usr/bin/R

list.files()

library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
theme_set(theme_pubr())


setwd("/medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/Plotting/")


options(pillar.subtle=FALSE) # To fix the dark background when printing the error message
options(rlang_backtrace_on_error = "none") # Bowth are important for the dark background!


LETTERS[1:4]

scale_fill_own <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c('tan2', 'royalblue1', 'forestgreen'), c("kraken2", "diamond", "krakenuniq")), 
        ...
    )
}

Eukaryota_top10_bar<-ggplot(data=Eukaryota_top10, aes(x=Species, y=Counts_Normalized,fill=tool))+ geom_bar(stat="identity", position=position_dodge2(width=1, preserve="single")) + scale_fill_own() + ylab("Normalized Counts (%)") + coord_flip() + ggtitle("Top10 Eukaryota")


metadata<-read.table("Metadata_Validation_run1_All.txt", header=T)

#metadata<-read.table("test.txt", header=T)


i=16

for (i in 1:nrow(metadata)){
    print(i)
                                        # Start with Diamond
    diamondcounts<-metadata$Diamond[i]
    Samplename<-gsub(pattern="_CountsForplotting.txt", replacement="",tail(strsplit(diamondcounts, split="/")[[1]], n=1))
    print(Samplename)    
    diamondcounts<-read.table(diamondcounts, sep="\t", header=TRUE)
    diamondcounts$Counts_Normalized<-diamondcounts$Counts/sum(diamondcounts$Counts)*100
    diamondLinkage<-metadata$Diamond_Linkage[i]
    diamondLinkage<-read.table(diamondLinkage, sep="\t", header=TRUE)
    nrow(diamondLinkage)==nrow(diamondLinkage)
    diamond<-merge(diamondcounts, diamondLinkage, by='TaxID')[c('TaxID','Species','Counts','Counts_Normalized','Domain')]
    diamond$tool<-rep('diamond', times=nrow(diamond))
    diamond$tool<-as.factor(diamond$tool)
                                        # Kraken2
    kraken2counts<-metadata$Kraken2[i]
    kraken2counts<-read.table(kraken2counts, sep="\t", header=TRUE)
    kraken2counts<-kraken2counts[-which(kraken2counts$TaxID==9606),]
                                        # Remove human
    kraken2counts$Counts_Normalized<-kraken2counts$Counts/sum(kraken2counts$Counts)*100
    kraken2Linkage<-metadata$Kraken2_Linkage[i]
    kraken2Linkage<-read.table(kraken2Linkage, sep="\t", header=TRUE)
    kraken2Linkage<-kraken2Linkage[-which(kraken2Linkage$TaxID==9606),]  # Remove human
    nrow(kraken2counts)==nrow(kraken2Linkage)
    kraken2<-merge(kraken2counts, kraken2Linkage, by='TaxID')[c('TaxID','Species','Counts','Counts_Normalized','Domain')]
    kraken2$tool<-rep('kraken2', times=nrow(kraken2))
    kraken2$tool<-as.factor(kraken2$tool)
                                        # KrakenUniq
    krakenuniqcounts<-metadata$KrakenUniq[i]
    krakenuniqcounts<-read.table(krakenuniqcounts, sep="\t", header=TRUE)
    krakenuniqcounts<-krakenuniqcounts[-which(krakenuniqcounts$TaxID==9606),]  # Remove human  
    krakenuniqcounts$Counts_Normalized<-krakenuniqcounts$Counts/sum(krakenuniqcounts$Counts)*100
    krakenUniqLinkage<-metadata$KrakenUniq_Linkage[i]
    krakenUniqLinkage<-read.table(krakenUniqLinkage, sep="\t", header=TRUE)
    krakenUniqLinkage<-krakenUniqLinkage[-which(krakenUniqLinkage$TaxID==9606),]
    nrow(krakenuniqcounts)==nrow(krakenUniqLinkage)
    krakenuniq<-merge(krakenuniqcounts, krakenUniqLinkage, by='TaxID')[c('TaxID','Species','Counts','Counts_Normalized','Domain')]
    krakenuniq$tool<-rep('krakenuniq', times=nrow(krakenuniq))
    krakenuniq$tool<-as.factor(krakenuniq$tool)
                                        # Merge
    merged<-rbind(rbind(diamond,kraken2), krakenuniq)
    merged_sorted<-merged[order(merged$Counts_Normalized, decreasing=TRUE),]
                                        # Calc sum of normalized counts, append the domain
    merged_linkage<-merged_sorted[,c("Species","Domain")]
    merged_linkage<-merged_linkage %>% distinct() # remove duplicated rows
    sumofnormalizedCounts=aggregate(Counts_Normalized ~ Species, merged_sorted, sum) # Calc sum of normalized counts to get top 10 species
    sumofnormalizedCounts<-merge(sumofnormalizedCounts, merged_linkage, all.x=T, by='Species')
    sumofnormalizedCounts<-sumofnormalizedCounts[order(sumofnormalizedCounts$Counts_Normalized, decreasing=TRUE),]
                                        # Calc mean of normalized counts, append to domain
    meanofnormalizedCounts=aggregate(Counts_Normalized ~ Species, merged_sorted, mean) # Calc sum of normalized counts to get top 10 species
    meanofnormalizedCounts<-merge(meanofnormalizedCounts, merged_linkage, all.x=T, by='Species')
    meanofnormalizedCounts<-meanofnormalizedCounts[order(meanofnormalizedCounts$Counts_Normalized, decreasing=TRUE),]
                                        # Calc max of normalized counts, append to domain
    maxofnormalizedCounts=aggregate(Counts_Normalized ~ Species, merged_sorted, max) # Calc sum of normalized counts to get top 10 species
    maxofnormalizedCounts<-merge(maxofnormalizedCounts, merged_linkage, all.x=T, by='Species')
    maxofnormalizedCounts<-maxofnormalizedCounts[order(maxofnormalizedCounts$Counts_Normalized, decreasing=TRUE),]
                                        # Plotting
                                        #pdf("Heatmaps.pdf", height=15)
                                        # --- Bacteria ----------
                                        #Plot all bacteria, with atleast 0.01% 
    Bacteria<-merged_sorted[which(merged_sorted$Domain=="Bacteria"),]
    Bacteria<-Bacteria[which(Bacteria$Counts_Normalized >= 0.01),]
    if (nrow(Bacteria) == 0){
        Bacteria_heatmap_plot<-NULL
        Bacteria_top10_bar<-NULL
    } else {
                                        # Heatmap - To Include all!
    heat_colors<-brewer.pal(6, "YlOrRd")
    Bacteria_heatmap<-dcast(Bacteria, Species ~ tool, value.var="Counts_Normalized")
    rownames(Bacteria_heatmap)<-Bacteria_heatmap$Species
    Bacteria_heatmap$Species<-NULL
    Bacteria_heatmap$rowSum<-rowSums(Bacteria_heatmap, na.rm=TRUE) * ifelse(rowSums(is.na(Bacteria_heatmap)) == ncol(Bacteria_heatmap), NA, 1)
    Bacteria_heatmap<-Bacteria_heatmap[order(Bacteria_heatmap$rowSum,decreasing=TRUE),]    
    Bacteria_heatmap$rowSum<-NULL
    if (nrow(Bacteria_heatmap) > 1 | ncol(Bacteria_heatmap) > 1){
        Bacteria_heatmap_plot<-pheatmap(as.matrix(Bacteria_heatmap), na_col="snow2", cluster_rows=F,cluster_cols=F, color=heat_colors, main="Bacteria", cellheight=10, cellwidth=30, silent=TRUE)
    } else {
        Bacteria_heatmap_plot<-NULL
    }
                                        #dev.off()
                                        # Barplot to Include top10
    top10Bacteria<-head(maxofnormalizedCounts[which(maxofnormalizedCounts$Domain=="Bacteria"),], n=10)
    top10Bacteria<-top10Bacteria$Species
    Bacteria_top10<-Bacteria[which(Bacteria$Species %in% top10Bacteria),]
    Bacteria_top10$Species<-with(Bacteria_top10, reorder(Species, Counts_Normalized,max, decreasing=FALSE))# Order by max
    Bacteria_top10_bar<-ggplot(data=Bacteria_top10, aes(x=Species, y=Counts_Normalized,fill=tool))+ geom_bar(stat="identity", position=position_dodge2(width=1, preserve="single")) + scale_fill_own() + ylab("Normalized Counts (%)") + coord_flip() + ggtitle("Top10 Bacteria")
        }
                                        # --- Eukaryotes ----------
                                        #Plot all Eukaryota, with atleast 0.01% 
    Eukaryota<-merged_sorted[which(merged_sorted$Domain=="Eukaryota"),]
    Eukaryota<-Eukaryota[which(Eukaryota$Counts_Normalized >= 0.01),]
    if (nrow(Eukaryota) == 0){
        Eukaryota_heatmap_plot<-NULL
        Eukaryota_top10_bar<-NULL
    } else {
                                        # Heatmap - To Include all!
    heat_colors<-brewer.pal(6, "YlOrRd")
    Eukaryota_heatmap<-dcast(Eukaryota, Species ~ tool, value.var="Counts_Normalized")
    rownames(Eukaryota_heatmap)<-Eukaryota_heatmap$Species
    Eukaryota_heatmap$Species<-NULL
    Eukaryota_heatmap$rowSum<-rowSums(Eukaryota_heatmap, na.rm=TRUE) * ifelse(rowSums(is.na(Eukaryota_heatmap)) == ncol(Eukaryota_heatmap), NA, 1)
    Eukaryota_heatmap<-Eukaryota_heatmap[order(Eukaryota_heatmap$rowSum,decreasing=TRUE),]
    Eukaryota_heatmap$rowSum<-NULL
    if (nrow(Eukaryota_heatmap) > 1 | ncol(Eukaryota_heatmap) > 1){
        Eukaryota_heatmap_plot<-pheatmap(as.matrix(Eukaryota_heatmap), na_col="snow2", cluster_rows=F,cluster_cols=F, color=heat_colors, main="Eukaryotes", cellheight=10, cellwidth=30, silent=TRUE)
    } else {
        Eukaryota_heatmap_plot<-NULL
    }
                                        # Barplot to Include top10
    top10Eukaryota<-head(maxofnormalizedCounts[which(maxofnormalizedCounts$Domain=="Eukaryota"),], n=10)
    top10Eukaryota<-top10Eukaryota$Species
    Eukaryota_top10<-Eukaryota[which(Eukaryota$Species %in% top10Eukaryota),]
    Eukaryota_top10$Species<-with(Eukaryota_top10, reorder(Species, Counts_Normalized,max, decreasing=FALSE))# Order by max
    Eukaryota_top10_bar<-ggplot(data=Eukaryota_top10, aes(x=Species, y=Counts_Normalized,fill=tool))+ geom_bar(stat="identity", position=position_dodge2(width=1, preserve="single")) + scale_fill_own() + ylab("Normalized Counts (%)") + coord_flip() + ggtitle("Top10 Eukaryota")
        }
                                        # --- viral ----------
                                        #Plot all Virus
    Viral<-merged_sorted[which(merged_sorted$Domain=="Viruses"),]
    if (nrow(Viral) == 0){
        Viral_heatmap_plot<-NULL
        Viral_top10_bar<-NULL
    } else {
                                        # Heatmap - To Include all!
    heat_colors<-brewer.pal(6, "YlOrRd")
    Viral_heatmap<-dcast(Viral, Species ~ tool, value.var="Counts_Normalized")
    rownames(Viral_heatmap)<-Viral_heatmap$Species
    Viral_heatmap$Species<-NULL
    Viral_heatmap$rowSum<-rowSums(Viral_heatmap, na.rm=TRUE) * ifelse(rowSums(is.na(Viral_heatmap)) == ncol(Viral_heatmap), NA, 1)
    Viral_heatmap<-Viral_heatmap[order(Viral_heatmap$rowSum,decreasing=TRUE),]
    Viral_heatmap$rowSum<-NULL
                                        #Viral_heatmap_plot<-pheatmap(as.matrix(Viral_heatmap), na_col="snow2", cluster_rows=F,cluster_cols=F, color=heat_colors, main="Viruses", cellheight=10, cellwidth=30, silent=TRUE)
    if (nrow(Viral_heatmap) > 1 | ncol(Viral_heatmap) > 1){
        Viral_heatmap_plot<-pheatmap(as.matrix(Viral_heatmap), na_col="snow2", cluster_rows=F,cluster_cols=F, color=heat_colors, main="Viruses", cellheight=10, cellwidth=30, silent=TRUE)
    } else {
        Viral_heatmap_plot<-NULL
    }
                                        # Barplot to Include top10
    top10Viral<-head(maxofnormalizedCounts[which(maxofnormalizedCounts$Domain=="Viruses"),], n=10)
    top10Viral<-top10Viral$Species
    Viral_top10<-Viral[which(Viral$Species %in% top10Viral),]
    Viral_top10$Species<-with(Viral_top10, reorder(Species, Counts_Normalized,max, decreasing=FALSE))# Order by max
     Viral_top10_bar<-ggplot(data=Viral_top10, aes(x=Species, y=Counts_Normalized,fill=tool))+geom_bar(stat="identity", position=position_dodge2(width=1, preserve="single")) + scale_fill_own() + ylab("Normalized Counts (%)") + coord_flip() + ggtitle("Top10 Viral")
        }
    pdf(paste(Samplename, "_Barplot.pdf", sep=""), width=20, height=10)
    figure <- ggarrange(Bacteria_top10_bar, Eukaryota_top10_bar, Viral_top10_bar,labels = c("A", "B", "C"),common.legend=TRUE,ncol = 3, nrow = 1)
    plot(figure)
    dev.off()
                                        # Plot the heatmaps
    pdf(paste(Samplename, "_Heatmaps.pdf", sep=""), height=40, width=15)
    heatmaps=list()
    heatmaps[['p1']]=Bacteria_heatmap_plot[[4]]
    heatmaps[['p2']]=Eukaryota_heatmap_plot[[4]]
    heatmaps[['p3']]=Viral_heatmap_plot[[4]]
    layoutmatrix2=rbind(c('p1', 'p3'),c('p1','p2'))
    grid.arrange(grobs=heatmaps,ncol=2,  heigths=c(2,1,1), layout_matrix=layoutmatrix2)
    dev.off()
}




