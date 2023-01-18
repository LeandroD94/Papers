################# PREPARING THE ENVIRONMENT ########################

options(scipen = 100)

dir.create("Results")
setwd("Results")

{dir.create("plots")
dir.create("plots/clusters")
dir.create("plots/Abundances")
dir.create("DESeq2_plots")
}

{library("DESeq2")
library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2")
library("qiime2R")
library("readxl")
library("dplyr")
library("stringr")
library("dendextend")
library("ggpubr")
}

############### LOADING AND PREPARING DATA ##################

data<-qza_to_phyloseq(features="../Output_QIIME_Stool/table.qza", taxonomy="../Output_QIIME_Stool/assigned_taxa.qza")
# importing sample data
metadata<- read.table("../Metadata_Stool.csv", header=T)
is.character(metadata$Major_organ_involvement) # just a check
row.names(metadata)<-metadata$Sample_ID
metadata<-metadata[sample_names(data),]
identical(sample_names(data),metadata$Sample_ID) # TRUE
identical(length(!is.na(metadata$Sample_ID)), length(sample_names(data)),38) # TRUE
# adding sample data to phyloseq object
sample_data(data)<-as.data.frame(metadata)
head(sample_data(data))

############### PREPARING BASILAR OBJECTS (Stool) ####################

{data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
Taxa.genus<-as.data.frame(tax_table(data.genus))
data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
Taxa.fam<-as.data.frame(tax_table(data.fam))
data.class<- tax_glom(data, taxrank = "Class", NArm = F)
Taxa.class<-as.data.frame(tax_table(data.class))
data.order<- tax_glom(data, taxrank = "Order", NArm = F)
Taxa.order<-as.data.frame(tax_table(data.order))
data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
Taxa.phy<-as.data.frame(tax_table(data.phy))
}

{data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}

################# ASSIGNED IN SILVA DATABASE ###############

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Total","Assigned","%","Taxa")
}
write.csv2(assigned,file="Percent_assigned_taxa.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)

###################### COUNTS EXPORT ########################

dir.create("Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
}
dir.create("TSS_counts")
{write.csv(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="TSS_counts/counts_phylum.csv",quote=F)
write.csv(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="TSS_counts/counts_class.csv",quote=F)
write.csv(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="TSS_counts/counts_order.csv",quote=F)
write.csv(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="TSS_counts/counts_family.csv",quote=F)
write.csv(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="TSS_counts/counts_genus.csv",quote=F)
}

################# RAREFACTION ANALYSIS ###########################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    check<-"red"
    if(length(x) < 10) {
      points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col="cyan",pch=16,cex=1)
    } else {
      slope<-mean(diff(v[(length(v)-t):length(v)])/dx)
      if(slope<lim) {
        check<-"blue"
        sat = sat+1
      }
      cat(i,slope,check,"\n")
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

png(file="plots/rarefaction_plot.png",width=2000,height=1500, res=300)
r<-rarecurve(t(otu_table(data)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=0.6)
dev.off()

########### HIERARCHICAL CLUSTERING EPGA vs HEALTHY ###############

c<-hclust(dist(t(sqrt(otu_table(data.prop))))) # euclidean by default
c<-as.dendrogram(c)
{Color_table<-as.data.frame(cbind(sample_data(data)$Condition,cbind(sample_data(data)$Condition)))
colnames(Color_table)<-c("Gruppo","Colore")
Color_table$Colore<-gsub(pattern= "HC", replacement= "4", x= Color_table$Colore)
Color_table$Colore<-gsub(pattern= "EGPA", replacement= "2", x= Color_table$Colore)
colors_to_use <- as.numeric(Color_table$Colore)
colors_to_use <- colors_to_use[order.dendrogram(c)]
labels_colors(c)<-colors_to_use
}
png(file="plots/clusters/cluster_commmunity_structure.png",width=2800,height=1500, res=300)
plot(c,main="Community structure using Euclidean distance \n on sqrt relative abundance of ASVs")
dev.off()

################### TOP ABUNDANCES PLOTS ######################

# TOP 5 Phylum stacked bar plot
{ data.phy.na<- tax_glom(data, taxrank = "Phylum", NArm = T) # because the 5th is a NA
data.phy.na <- transform_sample_counts(data.phy.na, function(OTU) OTU/sum(OTU)*100)
top5 <- names(sort(taxa_sums(data.phy.na), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.na)
others<-taxa_names(data.phy.na)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.na)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
tabella$labels<-gsub("HC","Healthy",tabella$Condition)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(labels),scales = "free_x") +# facet grid raggruppa, scales toglie ripetizione sample inutilizzata +
  geom_bar(stat="identity", position="stack") + theme_bw(base_size = 14) +# Stacked 100% barplot
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.5, "cm"),legend.text = element_text ( size = 9 )) +  # Vertical x-axis tick labels, legend sceglie dimensioni
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(y="Relative abundance")
ggsave(file="plots/Abundances/TOP_5_phyla_EGPA_vs_Healthy_stacked.png",width=9,height=5, dpi=300)

rm(tabella, data.phy.na, top5)

# TOP 5 Generi stacked bar plot
{data.genus.prop <- transform_sample_counts(data.genus, function(OTU) OTU/sum(OTU)*100)
top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.genus.prop)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
tabella$labels<-gsub("HC","Healthy",tabella$Condition)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(labels),scales = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.5, "cm"),legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(y="Relative abundance")
ggsave(file="plots/Abundances/TOP_5_genera_EGPA_vs_Healthy_stacked.png",width=9,height=5, dpi=300)

################## ALFA DIVERSITY _EGPA vs Healthy #####################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1", "Shannon", "Observed"), nrow=1, x= "Condition", shape=NULL, color="Sample_ID")
{Chao<-dplyr::filter(pAlpha$data, variable=="Chao1")
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Codice,obs$Codice) # TRUE, stesso ordine
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating the plot
New_data<-rbind.data.frame(Chao,H,ev)
pAlpha$data<-New_data
}
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw(base_size = 20) + labs(x="") +
  guides(color=FALSE) +
  stat_compare_means(aes(label=sprintf("p=%5.2f", as.numeric(..p.format..)), group = Condition), method = "wilcox.test", label.x= 1.25, size=6, label.y.npc = "top", vjust=-0.5)
ggsave(file="plots/plot_Alfa_diversity_wilcox.png", width = 8,height =10, dpi=300)

rm(pAlpha, Chao, H, obs, ev, New_data)

############### BETA DIVERSITY EPGA vs HEALTHY #####################

# Bray curtis is sensible to different library depth --> TSS normalization
{OTU.prop<-as.data.frame(otu_table(data.prop))
OTU.genus.prop<-as.data.frame(otu_table(data.genus.prop))
OTU.fam.prop<-as.data.frame(otu_table(data.fam.prop))
OTU.class.prop<-as.data.frame(otu_table(data.class.prop))
OTU.order.prop<-as.data.frame(otu_table(data.order.prop))
OTU.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

# using PERMANOVA to test
sample_OTU<-as.data.frame(t(sqrt((OTU.prop)))) # sample on rows --> t
perm<- adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="bray")
perm
value_to_paste_ASV<-perm$aov.tab$`Pr(>F)`[1]

rm(sample_OTU, perm)

sample_OTU<-as.data.frame(t(sqrt(OTU.genus.prop))) # sample on rows --> t
perm<- adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="bray")
perm

rm(sample_OTU, perm)

sample_OTU<-as.data.frame(t(sqrt((OTU.phy.prop)))) # sample on rows --> t
perm<- adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="bray")
perm

rm(sample_OTU, perm)

# ASV in 2D
{sample<-sample_data(data.prop)
sample_data(data.prop)<-sample
data_sqrt<-transform_sample_counts(data.prop, sqrt)
DistBC = phyloseq::distance(data_sqrt, method = "bray")
ordBC = ordinate(data_sqrt, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data_sqrt, ordBC, color = "Condition") +
  geom_point(size=2.5) + theme_minimal() + stat_ellipse() + 
  labs(title="PCoA with Bray-Curtis distance on sqrt proportional ASV data", subtitle= paste0("PERMANOVA Pr(>F)=",round(value_to_paste_ASV,2)), x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), color="Group") + theme_minimal()
ggsave(file="plots/PCoA_Beta diversity_Bray_Curtis.png", width = 6, height = 4.5, dpi=300)

# ASV in 3D
{colore<-sample_data(data.prop)$Condition
colore<-gsub("EGPA","2",colore)
colore<-gsub("HC","4",colore)
Dist<- vegdist(t(otu_table(data_sqrt)), method = "bray")                                                                                                      
obj<-ecodist::pco(Dist)
matrix<-obj[["vectors"]]
row.names(matrix)<-sample_names(data_sqrt)
rgl::open3d(windowRect=c(25,25,1200,1200))
pca3d::pca3d(matrix, col=as.numeric(colore), group=sample_data(data_sqrt)$Condition, legend="topleft", radius=2, shape="sphere", show.plane = F, show.centroids = F)
rgl::snapshot3d(file="plots/PCoA_3D.png", width=1200, height=1000)
}

rm(eigval, DistBC, ordBC, data_sqrt)

############## DA WITH DESEQ2 _EGPA vs Healthy #################

# cutting out rows with less than 10 reads total, as suggested by DeSeq2 tutorial
rm(data_prune)
{data_prune<- prune_taxa(taxa_sums(data) > 10, data)
data.genus<- tax_glom(data_prune, taxrank = "Genus", NArm = F)
data.fam<- tax_glom(data_prune, taxrank = "Family", NArm = F)
data.class<- tax_glom(data_prune, taxrank = "Class", NArm = F)
data.order<- tax_glom(data_prune, taxrank = "Order", NArm = F)
data.phy<- tax_glom(data_prune, taxrank = "Phylum", NArm = F)
}

# genera
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Condition) #DESEQ2 requires row data !
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "EGPA", "HC"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res}

{res<-res[order(res$padj), ]
  res.X<-as.data.frame(res)
  res.X$ASV<-row.names(res.X)
  Taxa.genus$ASV<-row.names(Taxa.genus)
  res.X<-dplyr::left_join(res.X, Taxa.genus)
  res.X$Kingdom<-NULL
  res.X$Species<-NULL
  res.X$Order<-NULL
  res.X$Class<-NULL
  res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DA_genera_DeSeq2_EGPA_vs_Healthy.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res.X$ASV
  target<-prune_taxa(target,data.genus)
  target<- subset_taxa(target, Genus!="uncultured") # unknown taxa (not assigned in SILVA are removed from plot)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  labels<-gsub("HC","Healthy",tabella$Condition)
  tabella$labels<-labels
}
ggplot(tabella, aes(x= Genus, y=Abundance, fill= labels)) + theme_bw(base_size=16) +
  geom_boxplot(width=0.6) +
  theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0),axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=1)) +
  scale_x_discrete(expand=c(0.001, 0)) +
  labs(title= "Differently abundant genera", y="Sqrt Abundance")
ggsave(file="DESeq2_plots/Box_plot_significant_diff_genera_EGPA_vs_healthy.png",width=10,height=7, dpi=300)
dev.off()

rm(DE,res)

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "EGPA", "HC"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res}
{res<-res[order(res$padj), ]
  res.X<-as.data.frame(res)
  res.X$ASV<-row.names(res.X)
  Taxa.fam$ASV<-row.names(Taxa.fam)
  res.X<-dplyr::left_join(res.X, Taxa.fam)
  res.X$Kingdom<-NULL
  res.X$Species<-NULL
  res.X$Order<-NULL
  res.X$Class<-NULL
  res.X$Genus<-NULL
  res.X$Specie<-NULL
  res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DA_families_DeSeq2_EGPA_vs_Healthy.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res.X$ASV
  target<-prune_taxa(target,data.fam)
  target<- subset_taxa(target, Family!="uncultured") # unknown taxa (not assigned in SILVA are removed from plot)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  labels<-gsub("HC","Healthy",tabella$Condition)
  tabella$labels<-labels
}
ggplot(tabella, aes(x=Family, y=Abundance, fill= labels)) + theme_bw(base_size = 16) +
  geom_boxplot(width=0.6) +
  theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0),axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=1)) +
  scale_x_discrete(expand=c(0.001, 0)) + labs(x="") +
  labs(title= "Differently abundant families", y="Sqrt Abundance")
ggsave(file="DESeq2_plots/Box_plot_significant_diff_families_EGPA_vs_healthy.png",width=10,height=7, dpi=300)
dev.off()

rm(DE,res)

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "EGPA", "HC"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res }
# NO CLASSES

rm(DE,res)

# glom order
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "EGPA", "HC"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res
  res<-res[order(res$padj), ]
  res.X<-as.data.frame(res)
  res.X$ASV<-row.names(res.X)
  Taxa.order$ASV<-row.names(Taxa.order)
  res.X<-dplyr::left_join(res.X, Taxa.order)
  res.X$Kingdom<-NULL
  res.X$Genus<-NULL
  res.X$Specie<-NULL
  res.X$Family<-NULL
  res.X$Domain<-NULL
  res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DA_orders_DeSeq2_EGPA_vs_Healthy.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res.X$ASV
  target<-prune_taxa(target,data.order)
  target<- subset_taxa(target, Order!="uncultured") # unknown taxa (not assigned in SILVA are removed from plot)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  labels<-gsub("HC","Healthy",tabella$Condition)
  tabella$labels<-labels
}
ggplot(tabella, aes(x=Order, y=Abundance, fill= labels)) + theme_bw(base_size = 16) +
  geom_boxplot(width=0.6) +
  theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0),axis.text.x = element_text(angle = 0, hjust = 0.5, vjust= 1)) + 
  scale_x_discrete(expand=c(0.001, 0)) + labs(x="") +
  labs(title= "Differentialy abundant orders", y="Sqrt Abundance")
ggsave(file="DESeq2_plots/Box_plot_significant_diff_orders_EGPA_vs_healthy.png",width=10,height=7, dpi=300)
dev.off()

rm(DE,res)

# glom phylum
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "EGPA", "HC"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res
  res<-res[order(res$padj), ]
}
# NO PHYLA

rm(DE,res)

######### PART2: PREPARING DA WITH DESEQ2 ON OTHER FACTORS ############

dir.create("Other_factors")
setwd("Other_factors")
dir.create("DESeq2_plots")

data_sick<-subset_samples(data, Condition=="EGPA")
{data_prune<- prune_taxa(taxa_sums(data_sick) > 10, data_sick) # as suggested in DeSeq2 tutorial
data.genus<- tax_glom(data_prune, taxrank = "Genus", NArm = F)
data.fam<- tax_glom(data_prune, taxrank = "Family", NArm = F)
data.class<- tax_glom(data_prune, taxrank = "Class", NArm = F)
data.order<- tax_glom(data_prune, taxrank = "Order", NArm = F)
data.phy<- tax_glom(data_prune, taxrank = "Phylum", NArm = F)
}

############# DA ON Severe_eosinophilia FACTOR #####################
dir.create("DESeq2_plots/Severe_eosinophilia")

rm(DEseq_data, DE,res,res.X)

# genera
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Severe_eosinophilia)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Severe_eosinophilia", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Severe_eosinophilia/DA_genera_DeSeq2.csv", row.names = F, quote=F, na = "")

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Severe_eosinophilia)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Severe_eosinophilia", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Severe_eosinophilia/DA_families_DeSeq2.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Severe_eosinophilia)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Severe_eosinophilia", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
View(res.X)
}
write.csv2(res.X, file="DESeq2_plots/Severe_eosinophilia/DA_classes_DeSeq2.csv", row.names = F, quote=F, na = "")

# orders
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Severe_eosinophilia) 
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Severe_eosinophilia", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Severe_eosinophilia/DA_orders_DeSeq2.csv", row.names = F, quote=F, na = "")

# phyla
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~Severe_eosinophilia)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Severe_eosinophilia", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
#write.csv2(res.X, file="DESeq2_plots/Severe_eosinophilia/DA_phyla_DeSeq2.csv", row.names = F, quote=F, na = "")

################### DA ON Gender FACTOR ####################
dir.create("DESeq2_plots/Gender")

rm(DEseq_data, DE,res,res.X)

# glom genera
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Gender)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Gender", "F", "M"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Gender/DA_genera_DeSeq2_ ratio F vs M", row.names = F, quote=F, na = "")

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Gender)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Gender", "F", "M"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Gender/DA_families_DeSeq2_ ratio F vs M.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Gender)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Gender", "F", "M"))
res = res[order(res$padj,na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Gender/DA_classes_DeSeq2_ ratio F vs M.csv", row.names = F, quote=F, na = "")

# orders
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Gender)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Gender", "F", "M"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Gender/DA_orders_DeSeq2_ ratio F vs M.csv", row.names = F, quote=F, na = "")

# glom phylum
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~Gender)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Gender", "F", "M"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Gender/DA_phyla_DeSeq2_ ratio F vs M.csv", row.names = F, quote=F, na = "")

############ DA ON "Major_organ_involvement" FACTOR ##########
dir.create("DESeq2_plots/Major_organ_involvement")

rm(DEseq_data, DE,res,res.X)

# genera
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Major_organ_involvement)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Major_organ_involvement", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X) # unknown!
write.csv2(res.X, file="DESeq2_plots/Major_organ_involvement/DA_genera_DeSeq2", row.names = F, quote=F, na = "")

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Major_organ_involvement)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Major_organ_involvement", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Major_organ_involvement/DA_families_DeSeq2.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Major organ)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Major_organ_involvement", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Major_organ_involvement/DA_classes_DeSeq2.csv", row.names = F, quote=F, na = "")

# orders
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Major_organ_involvement)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Major_organ_involvement", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Major_organ_involvement/DA_orders_DeSeq2.csv", row.names = F, quote=F, na = "")

# glom phyla
{DEseq_data<-phyloseq_to_deseq2(data.phylum, ~Major_organ_involvement)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Major_organ_involvement", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Major_organ_involvement/DA_phyla_DeSeq2.csv", row.names = F, quote=F, na = "")

############# DA ON "Digestive symptoms" FACTORS ###############
dir.create("DESeq2_plots/Digestives_symptoms")

# genera
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Digestives_symptoms)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Digestives_symptoms", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Digestives_symptoms/DA_genera_DeSeq2", row.names = F, quote=F, na = "")

rm(DEseq_data, DE,res,res.X)

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Digestives_symptoms)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Digestives_symptoms", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Digestives_symptoms/DA_families_DeSeq2.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Digestives_symptoms)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Digestives_symptoms", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Digestives_symptoms/DA_classes_DeSeq2.csv", row.names = F, quote=F, na = "")

# orders
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Digestives_symptoms)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Digestives_symptoms", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Digestives_symptoms/DA_orders_DeSeq2.csv", row.names = F, quote=F, na = "")

# phyla
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~Digestives_symptoms)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Digestives_symptoms", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Digestives_symptoms/DA_phyla_DeSeq2.csv", row.names = F, quote=F, na = "")

################# DA ON "ANCA" FACTOR ##############
dir.create("DESeq2_plots/ANCA")

rm(DEseq_data, DE,res,res.X)

# genera
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~ANCA)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("ANCA", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/ANCA/DA_genera_DeSeq2", row.names = F, quote=F, na = "")

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~ANCA)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("ANCA", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/ANCA/DA_families_DeSeq2.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~ANCA)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("ANCA", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/ANCA/DA_classes_DeSeq2.csv", row.names = F, quote=F, na = "")

# orders
{DEseq_data<-phyloseq_to_deseq2(data.order, ~ANCA)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("ANCA", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/ANCA/DA_orders_DeSeq2.csv", row.names = F, quote=F, na = "")

# phyla
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~ANCA)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("ANCA", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/ANCA/DA_phyla_DeSeq2.csv", row.names = F, quote=F, na = "")

############### DA ON "Immunosuppression" FACTOR ###############
dir.create("DESeq2_plots/Immunosuppression")

rm(DEseq_data, DE,res,res.X)

# genus
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Immunosuppression)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Immunosuppression", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Immunosuppression/DA_genera_DeSeq2", row.names = F, quote=F, na = "")

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Immunosuppression)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Immunosuppression", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Immunosuppression/DA_families_DeSeq2.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Immunosuppression)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Immunosuppression", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
View(res.X)
}
write.csv2(res.X, file="DESeq2_plots/Immunosuppression/DA_classes_DeSeq2.csv", row.names = F, quote=F, na = "")

# glom order
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Immunosuppression)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Immunosuppression", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Species<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Immunosuppression/DA_orders_DeSeq2.csv", row.names = F, quote=F, na = "")

# phyla
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~Immunosuppression)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Immunosuppression", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Immunosuppression/DA_phyla_DeSeq2.csv", row.names = F, quote=F, na = "")

##################### DA ON "Active_disease" FACTOR ###################
dir.create("DESeq2_plots/Active_disease")

rm(DEseq_data, DE,res,res.X)

# glom genere
{DEseq_data<-phyloseq_to_deseq2(data.genus, ~Active_disease)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Active_disease", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.genus$ASV<-row.names(Taxa.genus)
res.X<-dplyr::left_join(res.X, Taxa.genus)
res.X$Kingdom<-NULL
res.X$Specie<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.XFamily<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Active_disease/DA_genera_DeSeq2", row.names = F, quote=F, na = "")

# families
{DEseq_data<-phyloseq_to_deseq2(data.fam, ~Active_disease)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Active_disease", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.fam$ASV<-row.names(Taxa.fam)
res.X<-dplyr::left_join(res.X, Taxa.fam)
res.X$Kingdom<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Class<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Active_disease/DA_families_DeSeq2.csv", row.names = F, quote=F, na = "")

# classes
{DEseq_data<-phyloseq_to_deseq2(data.class, ~Active_disease)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Active_disease", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Active_disease/DA_classes_DeSeq2.csv", row.names = F, quote=F, na = "")

# orders
{DEseq_data<-phyloseq_to_deseq2(data.order, ~Active_disease)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Active_disease", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Species<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Active_disease/DA_orders_DeSeq2.csv", row.names = F, quote=F, na = "")

# phyla
{DEseq_data<-phyloseq_to_deseq2(data.phy, ~Active_disease)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Active_disease", "1", "0"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]
res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.phy$ASV<-row.names(Taxa.phy)
res.X<-dplyr::left_join(res.X, Taxa.phy)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Specie<-NULL
res.X$Family<-NULL
res.X$Domain<-NULL
res.X$Class<-NULL
res.X$Order<-NULL
}
View(res.X)
write.csv2(res.X, file="DESeq2_plots/Active_disease/DA_phyla_DeSeq2.csv", row.names = F, quote=F, na = "")

setwd("..")

################ PART 3: MUCOSAL BIOPSY MICROBIOTA ################

dir.create("Biopsy")
setwd("Biopsy")

dir.create("Abundances")

# importing biopsy data
remove_old<-ls(pattern="data")
remove_old<-remove_old[!remove_old=="metadata"]
rm(list=remove_old)
data<-qza_to_phyloseq(features="../../Output_QIIME_Biopsy/table.qza", taxonomy="../../Output_QIIME_Biopsy/assigned_taxa.qza")

######### PART3: PREPARING BASILAR OBJECT OF BIOPSY ###########

{data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
OTU.genus<-as.data.frame(otu_table(data.genus))
Taxa.genus<-as.data.frame(tax_table(data.genus))
data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
Taxa.fam<-as.data.frame(tax_table(data.fam))
data.class<- tax_glom(data, taxrank = "Class", NArm = F)
Taxa.class<-as.data.frame(tax_table(data.class))
data.order<- tax_glom(data, taxrank = "Order", NArm = F)
Taxa.order<-as.data.frame(tax_table(data.order))
data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
Taxa.phy<-as.data.frame(tax_table(data.phy))
}
{data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}

############## PART 3: ASSIGNED IN SILVA DATABASE ##############

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assegnati<-rbind.data.frame(a,b,c,d,e)
colnames(assegnati)<-c("Total","Assigned","%","Taxa")
}
write.csv2(assegnati,file="Percent_assigned_taxa.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assegnati)

################ PART 3: BIOPSY COUNTS EXPORT #################

dir.create("Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
}

dir.create("TSS_counts")
{write.csv(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="TSS_counts/counts_phylum.csv",quote=F)
write.csv(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="TSS_counts/counts_class.csv",quote=F)
write.csv(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="TSS_counts/counts_order.csv",quote=F)
write.csv(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="TSS_counts/counts_family.csv",quote=F)
write.csv(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="TSS_counts/counts_genus.csv",quote=F)
}

############## PART3: BIOPSY TOP ABUNDANCES PLOTS ##########################

# TOP 5 Phyla stacked bar plot
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.prop)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") + theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.5, "cm"),legend.text = element_text ( size = 9 )) +  
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + labs(y="Relative abundance")
ggsave(file="Abundances/TOP_5_phyla_Biopsy.png",width=7,height=5, dpi=300)

# TOP 5 Genera stacked bar plot
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.genus.prop)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
tabella_top<-psmelt(prune.dat_top) 
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus))  +
  geom_bar(stat="identity", position="stack") + theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.5, "cm"),legend.text = element_text ( size = 9 )) +  
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + labs(y="Relative abundance")
ggsave(file="Abundances/TOP_5_genera_Biopsy.png",width=7,height=5, dpi=300)

################ PART 3: BIOPSY CORRELATIONS ###########################

library(Hmisc)
Leucocytes<- read_excel("../../Leucocytes_Abundances.xlsx")
head(Leucocytes)
colnames(Leucocytes)[1]<-"Sample"
dir.create("Correlations")

# VS top 5 genera
{top5 <- names(sort(taxa_sums(data.genus), decreasing=TRUE))[1:5]
data.genus.prop <- transform_sample_counts(data.genus, function(OTU) OTU/sum(OTU))
prune.dat_top5 <- prune_taxa(top5,data.genus.prop)
target<-as.data.frame(t(otu_table(prune.dat_top5)))
row.names(target)<-gsub("B2","",row.names(target))
target<-target[Leucocytes$Sample,]
taxa<-as.data.frame(tax_table(prune.dat_top5))
colnames(target)<-taxa$Genus
}

{x<-cbind(target,Leucocytes[,2:13])
x<-as.matrix(x)                        
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)                                                               
data_corr<-as.data.frame(as.table(r$r)) 
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])
data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("Genus","Leuco","Corr","pvalue")
data_corr<-subset(data_corr, Genus %in% taxa$Genus)
data_corr<-subset(data_corr, Leuco %in% colnames(Leucocytes))
data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
write.csv2(data_corr,file="Correlations/Leucoc_Genus.csv", row.names = F, quote = F)
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$Sign<0.05]<-"*"
data_corr$Sign[data_corr$Sign>0.05]<-""
ggplot(data_corr, aes(x = data_corr[,1], y = data_corr[,2], fill = Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +                   
  scale_fill_gradient2(low="blue", mid = "white", high="red") +          
  theme_bw(base_size=5) +                                              
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 10), axis.text.y = element_text(size = 12)) +                       
  guides(fill= guide_colourbar(title ="Correlation value")) +          
  geom_text(aes(label= Sign), color= "black", size =8) +
  labs(y="",x="", caption = "p-value<0.05 = *") +   
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  +  
  labs(title = "Genera Spearman correlations") +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) +  
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Correlations/Heatmap_correlations_genera.png", height=7,width=7, dpi = 300)

rm(x, r, correlation, top5, target)

# VS top 5 families
{top5 <- names(sort(taxa_sums(data.fam), decreasing=TRUE))[1:5]
data.fam.prop <- transform_sample_counts(data.fam, function(OTU) OTU/sum(OTU))
prune.dat_top5 <- prune_taxa(top5, data.fam.prop)
target<-as.data.frame(t(otu_table(prune.dat_top5)))
row.names(target)<-gsub("B2","",row.names(target))
target<-target[Leucocytes$Sample,]
taxa<-as.data.frame(tax_table(prune.dat_top5))
colnames(target)<-taxa$Family
}

{x<-cbind(target,Leucocytes[,2:13])
x<-as.matrix(x)       
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)           
data_corr<-as.data.frame(as.table(r$r))   
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])
data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("fam","Leuco","Corr","pvalue")
data_corr<-subset(data_corr, fam %in% taxa$Family)
data_corr<-subset(data_corr, Leuco %in% colnames(Leucocytes))
data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
write.csv2(data_corr,file="Correlations/Leucoc_families.csv", row.names = F, quote = F)
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$Sign<0.05]<-"*"
data_corr$Sign[data_corr$Sign>0.05]<-""
ggplot(data_corr, aes(x = data_corr[,1], y = data_corr[,2], fill = Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +                       
  scale_fill_gradient2(low="blue", mid = "white", high="red") +              
  theme_bw(base_size=5) +                                           
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 10), axis.text.y = element_text(size = 12)) +                            #(ruota labels dei nomi sull'asse x)
  guides(fill= guide_colourbar(title ="Correlation value")) +       
  geom_text(aes(label= Sign), color= "black", size =8) +
  labs(y="",x="", caption = "p-value<0.05 = *") +  
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  +   
  labs(title = "Families Spearman correlations") +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) + 
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Correlations/Heatmap_correlations_families.png", height=7,width=7, dpi = 300)

rm(x, r, correlation, top5, target)

# VS top 5 orders
{top5 <- names(sort(taxa_sums(data.order), decreasing=TRUE))[1:5]
data.order.prop <- transform_sample_counts(data.order, function(OTU) OTU/sum(OTU))
prune.dat_top5 <- prune_taxa(top5,data.order.prop)
target<-as.data.frame(t(otu_table(prune.dat_top5)))
row.names(target)<-gsub("B2","",row.names(target))
target<-target[Leucocytes$Sample,] 
taxa<-as.data.frame(tax_table(prune.dat_top5))
colnames(target)<-taxa$Order
}

{x<-cbind(target,Leucocytes[,2:13])
x<-as.matrix(x)                                                                                                                    #(coercizione necessaria)
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)                                                              
data_corr<-as.data.frame(as.table(r$r))                                        
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])                                   
data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("order","Leuco","Corr","pvalue")
data_corr<-subset(data_corr, order %in% taxa$Order)
data_corr<-subset(data_corr, Leuco %in% colnames(Leucocytes))
data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
write.csv2(data_corr,file="Correlations/Leucoc_order.csv", row.names = F, quote = F)
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$Sign<0.05]<-"*"
data_corr$Sign[data_corr$Sign>0.05]<-""
ggplot(data_corr, aes(x = data_corr[,1], y = data_corr[,2], fill = Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +              
  scale_fill_gradient2(low="blue", mid = "white", high="red") +      
  theme_bw(base_size=5) +                                            
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 10), axis.text.y = element_text(size = 12)) +                            #(ruota labels dei nomi sull'asse x)
  guides(fill= guide_colourbar(title ="Correlation value")) +         
  geom_text(aes(label= Sign), color= "black", size =8) +
  labs(y="",x="", caption = "p-value<0.05 = *") +     
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  +
  labs(title = "Orders Spearman correlations") +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) + 
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Correlations/Heatmap_correlations_orders.png", height=7,width=7, dpi = 300)

rm(x, r, correlation, top5, target)

# vs top 5 classes
{top5 <- names(sort(taxa_sums(data.class), decreasing=TRUE))[1:5]
data.class.prop <- transform_sample_counts(data.class, function(OTU) OTU/sum(OTU))
prune.dat_top5 <- prune_taxa(top5,data.class.prop)
target<-as.data.frame(t(otu_table(prune.dat_top5)))
row.names(target)<-gsub("B2","",row.names(target))
target<-target[Leucocytes$Sample,] 
taxa<-as.data.frame(tax_table(prune.dat_top5))
colnames(target)<-taxa$Class
}

{x<-cbind(target,Leucocytes[,2:13])
x<-as.matrix(x)                  
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)                    
data_corr<-as.data.frame(as.table(r$r))            
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])       
data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("class","Leuco","Corr","pvalue")
data_corr<-subset(data_corr, class %in% taxa$Class)
data_corr<-subset(data_corr, Leuco %in% colnames(Leucocytes))
data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
write.csv2(data_corr,file="Correlations/Leucoc_classes.csv", row.names = F, quote = F)
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$Sign<0.05]<-"*"
data_corr$Sign[data_corr$Sign>0.05]<-""
ggplot(data_corr, aes(x = data_corr[,1], y = data_corr[,2], fill = Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +                 
  scale_fill_gradient2(low="blue", mid = "white", high="red") +        
  theme_bw(base_size=5) +                                              
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 10), axis.text.y = element_text(size = 12)) +                            #(ruota labels dei nomi sull'asse x)
  guides(fill= guide_colourbar(title ="Correlation value")) +          
  geom_text(aes(label= Sign), color= "black", size =8) +
  labs(y="",x="", caption = "p-value<0.05 = *") +     
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  + 
  labs(title = "Classes Spearman correlations") +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) + 
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Correlations/Heatmap_correlations_classes.png", height=7,width=7, dpi = 300)

rm(x, r, correlation, top5, target)

# vs top 5 phyla
{top5 <- names(sort(taxa_sums(data.phy), decreasing=TRUE))[1:5]
data.phy.prop <- transform_sample_counts(data.phy, function(OTU) OTU/sum(OTU))
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
target<-as.data.frame(t(otu_table(prune.dat_top5)))
row.names(target)<-gsub("B2","",row.names(target))
target<-target[Leucocytes$Sample,] 
taxa<-as.data.frame(tax_table(prune.dat_top5))
colnames(target)<-taxa$Phylum
}
{x<-cbind(target,Leucocytes[,2:13])
x<-as.matrix(x)                                                                                                                    #(coercizione necessaria)
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)                            
data_corr<-as.data.frame(as.table(r$r))                    
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])               
data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("phy","Leuco","Corr","pvalue")
data_corr<-subset(data_corr, phy %in% taxa$Phylum)
data_corr<-subset(data_corr, Leuco %in% colnames(Leucocytes))
data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
write.csv2(data_corr,file="Correlations/Leucoc_phyla.csv", row.names = F, quote = F)
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$Sign<0.05]<-"*"
data_corr$Sign[data_corr$Sign>0.05]<-""
ggplot(data_corr, aes(x = data_corr[,1], y = data_corr[,2], fill = Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +             
  scale_fill_gradient2(low="blue", mid = "white", high="red") +    
  theme_bw(base_size=5) +                                          
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 10), axis.text.y = element_text(size = 12)) +                            #(ruota labels dei nomi sull'asse x)
  guides(fill= guide_colourbar(title ="Correlation value")) +      
  geom_text(aes(label= Sign), color= "black", size =8) +
  labs(y="",x="", caption = "p-value<0.05 = *") +     
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  + 
  labs(title = "Phyla Spearman correlations") +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) + 
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Correlations/Heatmap_correlations_phyla.png", height=7,width=7, dpi = 300)

rm(x, r, correlation, top5, target)

# eosinophils vs Th
x<-as.matrix(Leucocytes[,2:13])                                                                                                                    #(coercizione necessaria)
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)        
data_corr<-as.data.frame(as.table(r$r))
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])

{data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("Parameters","Th","Corr","pvalue")
data_corr<-subset(data_corr, Parameters %in% colnames(Leucocytes)[12:13])
data_corr<-subset(data_corr, Th %in% colnames(Leucocytes)[2:11])
data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
write.csv2(data_corr,file="Correlations/Leucoc_Nuovi parametri.csv", row.names = F, quote = F)
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$Sign<0.05]<-"*"
data_corr$Sign[data_corr$Sign>0.05]<-""
}
ggplot(data_corr, aes(x = data_corr[,2], y = data_corr[,1], fill = Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +     
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=5) +                                      
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 10), axis.text.y = element_text(size = 12)) +                            #(ruota labels dei nomi sull'asse x)
  guides(fill= guide_colourbar(title ="Correlation value")) +    
  geom_text(aes(label= Sign), color= "black", size =8) +
  labs(y="",x="", caption = "p-value<0.05 = *") +    
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  + 
  labs(title = "") +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) + 
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Correlations/Heatmap_correlations_eosinoph_vs_Th.png", height=7,width=7, dpi = 300)
