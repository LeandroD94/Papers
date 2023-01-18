{library("phyloseq")
library("ggplot2")
library("vegan")
library("dplyr")
library("DESeq2")
library("qiime2R")
library("readxl")
}

############## IMPORTING FROM QIIME2 AND PREPARING WORKSPACE ##################

data<-qza_to_phyloseq(features="QIIME2/table.qza", taxonomy="QIIME2/assigned_taxa.qza", metadata = "metadata.tsv")
sample_data(data)

######## updating sample names
sample_names(data)<-sample_data(data)$Sample

###### checking reads assigned in Prokaryote kingdom
{ Unass<-tax_glom(data, taxrank = "Kingdom")
Unass.prop<-transform_sample_counts(Unass, function (x) (x/sum(x)*100) )
a<-cbind( apply(otu_table(Unass), sum, MARGIN = 1) )
total<-sum(a)
b<-cbind( (a/total)*100, apply(otu_table(Unass), min, MARGIN = 1) , apply(otu_table(Unass), max, MARGIN = 1 ) )
c<-otu_table(Unass)
c_a<-NULL
for(x in 1:length(row.names(c))) {
  c_a<-c(c_a,colnames(c)[which.min(c[x,])]) }
c_b<-NULL
for(x in 1:length(row.names(c))) {
  c_b<-c(c_b,colnames(c)[which.max(c[x,])]) }
d<-as(tax_table(Unass),"matrix")
e<-cbind.data.frame(a,b, c_a, c_b, d[,"Kingdom"])
colnames(e)<-c("total_ASV_count","Percent_proportion","Minor_count_among_samples","Major_count_among_samples", "FASTQ_with_less_count","FASTQ_with_more_count","Kingdom")
e$Kingdom<-gsub("d__","",e$Kingdom)
row.names(e)<-e$Kingdom
}
options(scipen=100)
e
write.csv2(e[,colnames(e)!="Kingdom"], file="QIIME2/Unassigned_domain_checking.csv", row.names = T, quote = F)
####### filtering out every Unassigned ASV
data<- subset_taxa(data, Kingdom!="Unassigned")
head(tax_table(data))
write.csv2(tax_table(data), file="Every_filtered_ASV_and_taxonomic_assignment.csv", row.names = T)

################### STARTING THE ANALYSIS ##################

dir.create("Results")
setwd("Results")

{ OTU<-as.data.frame(otu_table(data))
Taxa<-as.data.frame(tax_table(data))

data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
OTU.genus<-as.data.frame(otu_table(data.genus))
Taxa.genus<-as.data.frame(tax_table(data.genus))

data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
OTU.fam<-as.data.frame(otu_table(data.fam))
Taxa.fam<-as.data.frame(tax_table(data.fam))

data.class<- tax_glom(data, taxrank = "Class", NArm = F)
OTU.class<-as.data.frame(otu_table(data.class))
Taxa.class<-as.data.frame(tax_table(data.class))

data.order<- tax_glom(data, taxrank = "Order", NArm = F)
OTU.order<-as.data.frame(otu_table(data.order))
Taxa.order<-as.data.frame(tax_table(data.order))

data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
OTU.phy<-as.data.frame(otu_table(data.phy))
Taxa.phy<-as.data.frame(tax_table(data.phy))
}
  
{data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}


############## % ASSIGNED IN SILVA DATABASE #####################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assegnati<-rbind.data.frame(a,b,c,d,e)
colnames(assegnati)<-c("Totale","Assegnati","%","Taxa")
}
assegnati
write.csv2(assegnati,file="Taxa_assigned_in_taxonomic_database.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assegnati)

################## COUNTS EXPORT ########################

{dir.create("Raw_counts")
write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
}

# exporting as relative abundances
options(scipen=100)
{dir.create("Proportional_counts")
write.csv2(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Proportional_counts/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Proportional_counts/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class.prop),"matrix")),file="Proportional_counts/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order.prop),"matrix")),file="Proportional_counts/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam.prop),"matrix")),file="Proportional_counts/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Proportional_counts/counts_genus.csv",quote=F)
}

################## RAREFACTION CURVE #################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    check<-"red"
    if(length(x) < 10) {
      points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col="cyan",pch=16,cex=1)
    } else {
      #the slope is estimated (average) at the last b rarefaction points
      slope<-mean(diff(v[(length(v)-t):length(v)])/dx)
      if(slope<lim) {
        check<-"blue"
        sat = sat+1
      }
      cat(i,slope,check,"\n")
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  legend("bottomright",paste(sat,"saturated samples"),bty="n")}

png(file="Rarefaction_curve.png",width=2100,height=1800, res=300)
r<-rarecurve(t(otu_table(data)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)  # if sample names is red then it's not saturated!
dev.off()


################ HIERARCHICAL CLUSTERING ###############

# euclidean distance with VST normalization
library(dendextend)
DEseq_data<-phyloseq_to_deseq2(data, ~Nazionality) # the normalization is inside DeSeq2 package
vst<-varianceStabilizingTransformation(DEseq_data)
{c<-hclust(dist(t(assay(vst))))
o<-c$order
c<-as.dendrogram(c)
tabella_colore<-as.data.frame(cbind(sample_data(data)$Nazionality,cbind(sample_data(data)$Nazionality)))
colnames(tabella_colore)<-c("Gruppo","Colore")
tabella_colore$Colore<-gsub(pattern= "Italy", replacement= "5", x= tabella_colore$Colore)
tabella_colore$Colore<-gsub(pattern= "Germany", replacement= "2", x= tabella_colore$Colore)
colors_to_use <- as.numeric(tabella_colore$Colore)
colors_to_use <- colors_to_use[order.dendrogram(c)]
labels_colors(c)<-colors_to_use
}
png(file="Hierarchical_clustering_Euclidean_distance_VST.png",width=2500,height=1800,res=300)
plot(c,main="Hierarchical clustering using Euclidean distance \n on ASV after VST normalization")
dev.off()

################# TOP ABUNDANCES PLOT #####################

dir.create("Abundances")

# TOP 5 Phylum stacked bar plot
{top5 <- names(sort(taxa_sums(data.phy), decreasing=TRUE))[1:5]
data.phy.prop <- transform_sample_counts(data.phy, function(OTU) OTU/sum(OTU))
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
prune.dat_top5 <- transform_sample_counts(prune.dat_top5, function(OTU) OTU/sum(OTU)*100)
tabella<-psmelt(prune.dat_top5)
}
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Nazionality),scales = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 8 )) + 
  theme(legend.position="bottom", strip.text.x=element_blank()) + guides(fill=guide_legend(nrow=1)) + labs(y="Relative abundance", x="")
ggsave(file="Abundances/TOP5_phyla_Nazionality_stacked.pdf",width=9,height=5)

# TOP 10 Generi stacked bar plot
{top10 <- names(sort(taxa_sums(data.genus), decreasing=TRUE))[1:10]
data.genus.prop <- transform_sample_counts(data.genus, function(OTU) OTU/sum(OTU))
prune.dat_top10 <- prune_taxa(top10,data.genus.prop)
prune.dat_top10 <- transform_sample_counts(prune.dat_top10, function(OTU) OTU/sum(OTU)*100)
tabella<-psmelt(prune.dat_top10)
}
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Nazionality),scales = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 8 )) +
  theme(legend.position="bottom", strip.text.x=element_blank()) + guides(fill=guide_legend(nrow=2)) + labs(y="Relative abundance", x="")
ggsave(file="Abundances/TOP10_Genera_Nazionality_stacked.pdf",width=9,height=5)

################ ALFA DIVERSITY GERMANY VS ITA ##################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1", "Shannon", "Observed"), nrow=1, x= "Nazionality", shape=NULL, color="ID_nuovo") # getting values from plot
{ H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness to plot
identical(H$Codice,obs$Codice) # TRUE --> same order
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev) # new plot
New_data$Nazionality<-gsub("Germany","G",New_data$Nazionality)
New_data$Nazionality<-gsub("Italy","I",New_data$Nazionality)
}
library(ggpubr) # to add p values to plot
ggplot(data=New_data, aes(x=Nazionality, y=value)) + geom_boxplot(data=New_data, alpha=0.1) + facet_wrap(~variable, scales = "free_y") + theme_bw() + labs(x="", y="") +
  stat_compare_means(data=New_data, aes(group = Nazionality), label="p.format", method = "wilcox.test", label.x= 1.4, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.3)
ggsave(file="plot_Alfa_diversity_GERvsITA_wilcox.png", width = 6,height =5, dpi= 300) 

################ BETA DIVERSITY GERMANY vs ITA ###################

dir.create("Beta_diversity")
setwd("Beta_diversity")

{ASV.prop<-as.data.frame(otu_table(data.prop))
ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")

sample_OTU<-as.data.frame(t(ASV.prop)) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Nazionality, data=metadata, permutations = 9999, method="bray")
perm_ASV
perm_ASV_Bray<-perm_ASV # needed later for plot

sample_OTU<-as.data.frame(t(ASV.genus.prop)) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Nazionality, data=metadata, permutations = 9999, method="bray")
perm_g

sample_OTU<-as.data.frame(t(ASV.fam.prop))
perm_f<- vegan::adonis(sample_OTU ~Nazionality, data=metadata, permutations = 9999, method="bray")
perm_f

sample_OTU<-as.data.frame(t(ASV.class.prop))
perm_o<- vegan::adonis(sample_OTU ~Nazionality, data=metadata, permutations = 9999, method="bray")
perm_o

sample_OTU<-as.data.frame(t(ASV.order.prop))
perm_c<- vegan::adonis(sample_OTU ~Nazionality, data=metadata, permutations = 9999, method="bray")
perm_c

sample_OTU<-as.data.frame(t(ASV.phy.prop))
perm_p<- vegan::adonis(sample_OTU ~Nazionality, data=metadata, permutations = 9999, method="bray")
perm_p

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Beta_diversity_permanova_Bray_Ger_vs_ITA.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(ASV.prop), distance="bray")
disper<-vegan::betadisper(BC.dist,metadata$Nazionality)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Beta_dispersion_permanova_Bray_sqrt prop ASV_GERvsITA.csv",quote=F,row.names = T)

setwd("..")

################## PCoA GERMANY VS ITA PLOTS #####################

setwd("Beta_diversity")

# on ASV
{sample_data(data.prop)
DistBC = phyloseq::distance(data.prop, method = "bray")
ordBC = ordinate(data.prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.prop, ordBC, color = "Nazionality") +
  geom_point(size=2.5) + stat_ellipse() + guides(color=FALSE) +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4) + theme_minimal() +
  labs(title="PCoA Bray-Curtis distance \n calculated on ASV percent abundance", color="Nazionality", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_div_on_ASV_Bray_Curtis.png", width = 8, height = 4.5, dpi=300)

# on genera
{DistBC = phyloseq::distance(data.genus.prop, method = "bray")
ordBC = ordinate(data.genus.prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.genus.prop, ordBC, color = "Nazionality") +
  geom_point(size=2.5) + stat_ellipse() + guides(color=FALSE) +
  geom_text(aes(label = sample_names(data.genus.prop)), position = (position_nudge(y=0.04)), size=4) + theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on genera percent abundance", color="Nazionality", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_div_on_genera_Bray_Curtis.png", width = 8, height = 4.5, dpi=300)

# on ASV again but with Jaccard
{sample_data(data.prop)
DistBC = phyloseq::distance(data.prop, method = "jaccard")
ordBC = ordinate(data.prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.prop, ordBC, color = "Nazionality") +
  geom_point(size=2.5) + stat_ellipse() + guides(color=FALSE) +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4) + theme_minimal() +
  labs(title="PCoA with Jaccard distance \n calculated on ASV percent abundance", color="Nazionality", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_div_on_ASV_Jaccard.png", width = 8, height = 4.5, dpi=300)

setwd("..")

################## BRAY CURTIS ON OTHER FACTORS ##################

setwd("Beta_diversity")

dir.create("Other_factors")
setwd("Other_factors")

ASV.prop<-as.data.frame(otu_table(data.prop))
metadata<-as(sample_data(data.prop),"data.frame")
sample_OTU<-as.data.frame(t(ASV.prop)) # samples has to be on rows --> t

perm_ASV_BAV<- vegan::adonis(sample_OTU ~BAV, data=metadata, permutations = 9999, method="bray")
perm_ASV_BAV

perm_ASV_CAD<- vegan::adonis(sample_OTU ~CAD, data=metadata, permutations = 9999, method="bray")
perm_ASV_CAD
# testing its beta dispersion
BC.dist<-vegan::vegdist(t(ASV.prop), distance="bray")
disper<-vegan::betadisper(BC.dist,metadata$CAD)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
system("touch 'CAD dispersion NOT significant' ")

perm_ASV_Gender<- vegan::adonis(sample_OTU ~Gender, data=metadata, permutations = 9999, method="bray")
perm_ASV_Gender

perm_ASV_Atrial<- vegan::adonis(sample_OTU ~Atrial.Fibrillation.AF., data=metadata, permutations = 9999, method="bray")
perm_ASV_Atrial

perm_ASV_Smoke<- vegan::adonis(sample_OTU ~Smoke, data=metadata, permutations = 9999, method="bray")
perm_ASV_Smoke

perm_ASV_Calcific<- vegan::adonis(sample_OTU ~Calcification.Degree, data=metadata, permutations = 9999, method="bray")
perm_ASV_Calcific

B<-rbind(perm_ASV_BAV$aov.tab[1,],perm_ASV_CAD$aov.tab[1,],perm_ASV_Gender$aov.tab[1,],
         perm_ASV_Atrial$aov.tab[1,],perm_ASV_Smoke$aov.tab[1,],perm_ASV_Calcific$aov.tab[1,])
B
write.csv2(B,file="Permanova_BRAY_on_ASV_other_factors.csv", row.names = T, quote = F)

sample_data(data.prop)
{ DistBC = phyloseq::distance(data.prop, method = "bray")
  ordBC = ordinate(data.prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}

plot_ordination(data.prop, ordBC, color = "CAD") +
  geom_point(size=2.5) + stat_ellipse() +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4, show.legend = F) + 
  theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on ASV percent abundance", color="CAD", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_ASV_Bray_factor_CAD.png", width = 8, height = 4.5, dpi=300)

plot_ordination(data.prop, ordBC, color = "BAV") +
  geom_point(size=2.5) + stat_ellipse() +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4, show.legend = F) +
  theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on ASV percent abundance", color="BAV", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_ASV_Bray_factor_BAV.png", width = 8, height = 4.5, dpi=300)

plot_ordination(data.prop, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4, show.legend = F) +
  theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on ASV percent abundance", color="Gender", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_ASV_Bray_factor_Gender.png", width = 8, height = 4.5, dpi=300)

plot_ordination(data.prop, ordBC, color = "Smoke") +
  geom_point(size=2.5) + stat_ellipse() +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4, show.legend = F) +
  theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on ASV percent abundance", color="Smoke", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_ASV_Bray_factor_Smoke.png", width = 8, height = 4.5, dpi=300)

plot_ordination(data.prop, ordBC, color = "Atrial.Fibrillation.AF.") +
  geom_point(size=2.5) + stat_ellipse() +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4, show.legend = F) +
  theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on ASV percent abundance", color="Atrial \nfibrilation", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_ASV_Bray_factor_Atrial_fibr.png", width = 8, height = 4.5, dpi=300)

plot_ordination(data.prop, ordBC, color = "Calcification.Degree") +
  geom_point(size=2.5) + stat_ellipse() +
  geom_text(aes(label = sample_names(data.prop)), position = (position_nudge(y=0.04)), size=4, show.legend = F) +
  theme_minimal() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on ASV percent abundance", color="Calcification \ndegree", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_bw()
ggsave(file="PCoA_Beta_ASV_Bray_factor_Calcification.png", width = 8, height = 4.5, dpi=300)

setwd("..")

#################### DA THROUGH DESEQ2 ######################

# better to prune values too low, see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
# or for example those articles PMID: 34262595 e PMID: 32708445
data_prune<- prune_taxa(taxa_sums(data) > 10, data)

# updating names beforehand for DeSeq2 plots
sample_data(data_prune)$Nazionality<-gsub("Italy","I", sample_data(data_prune)$Nazionality)
sample_data(data_prune)$Nazionality<-gsub("Germany","G", sample_data(data_prune)$Nazionality)
head(sample_data(data_prune))

# then preparing new objects
{data_pruned.genus<- tax_glom(data_prune, taxrank = "Genus", NArm = F)
data_pruned.fam<- tax_glom(data_prune, taxrank = "Family", NArm = F)
data_pruned.class<- tax_glom(data_prune, taxrank = "Class", NArm = F)
data_pruned.order<- tax_glom(data_prune, taxrank = "Order", NArm = F)
data_pruned.phy<- tax_glom(data_prune, taxrank = "Phylum", NArm = F)
}
  
# end of preparations, starting!
dir.create("DA_DEseq2")
dir.create("DA_DEseq2/DESeq2_plots")

# genera ################
DEseq_data_pruned<-phyloseq_to_deseq2(data_pruned.genus, ~Nazionality)
DE<-DESeq(DEseq_data_pruned)
res<-results(DE, contrast= c("Nazionality", "G", "I")) # ratio G/I
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]

{res_genus<-as.data.frame(res)
res_genus$ASV<-row.names(res_genus)
Taxa.genus$ASV<-row.names(Taxa.genus)
res_genus<-left_join(res_genus, Taxa.genus)
res_genus$Kingdom<-NULL
res_genus$Species<-NULL
}
View(res_genus)
write.csv2(res_genus, file="DA_DEseq2/Analisi diff dei generi con DeSeq2 _ ratio Germany vs Ita.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res_genus$ASV
target<-prune_taxa(target,data_pruned.genus)
target<- subset_taxa(target, Genus!="uncultured") # removes NAs from plot
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_g<-tabella
}

rm(DE, res)

# families ############
DEseq_data_pruned<-phyloseq_to_deseq2(data_pruned.fam, ~Nazionality)
DE<-DESeq(DEseq_data_pruned)
res<-results(DE, contrast= c("Nazionality", "G", "I"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]

{res_fam<-as.data.frame(res)
res_fam$ASV<-row.names(res_fam)
Taxa.fam$ASV<-row.names(Taxa.fam)
res_fam<-left_join(res_fam, Taxa.fam)
res_fam$Kingdom<-NULL
res_fam$Genus<-NULL
res_fam$Species<-NULL
}
View(res_fam)
write.csv2(res_fam, file="DA_DEseq2/Analisi diff delle famiglie con DeSeq2 _ ratio Germany vs Italia.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res_fam$ASV
target<-prune_taxa(target,data_pruned.fam)
target<- subset_taxa(target, Family!="uncultured") # removes NAs from plot
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_f<-tabella
}

rm(DE, res)

# classes #####################
DEseq_data_pruned<-phyloseq_to_deseq2(data_pruned.class, ~Nazionality)
DE<-DESeq(DEseq_data_pruned)
res<-results(DE, contrast= c("Nazionality", "G", "I"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]

{res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.class$ASV<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
res_class<-res.X
}
View(res.X)
write.csv2(res.X, file="DA_DEseq2/Analisi diff delle classi con DeSeq2 _ ratio Germany vs Italia.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res.X$ASV
target<-prune_taxa(target,data_pruned.class)
target<- subset_taxa(target, Class!="uncultured") # removing NAs from plot
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_c<-tabella
}

rm(DE, res)

# orders ###################
DEseq_data_pruned<-phyloseq_to_deseq2(data_pruned.order, ~Nazionality)
DE<-DESeq(DEseq_data_pruned)
res<-results(DE, contrast= c("Nazionality", "G", "I"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]

{res.X<-as.data.frame(res)
res.X$ASV<-row.names(res.X)
Taxa.order$ASV<-row.names(Taxa.order)
res.X<-dplyr::left_join(res.X, Taxa.order)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Species<-NULL
res.X$Family<-NULL
}
View(res.X)
res_order<-res.X
write.csv2(res.X, file="DA_DEseq2/Analisi diff degli ordini con DeSeq2 _ ratio Germany vs Italia.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res.X$ASV
target<-prune_taxa(target,data_pruned.order)
target<- subset_taxa(target, Order!="uncultured") # removing NAs from plot
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_o<-tabella
}
  
# phyla ####################
DEseq_data_pruned<-phyloseq_to_deseq2(data_pruned.phy, ~Nazionality)
DE<-DESeq(DEseq_data_pruned)
res<-results(DE, contrast= c("Nazionality", "G", "I"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
res<-res[order(res$padj), ]

{res_phy<-as.data.frame(res)
res_phy$ASV<-row.names(res_phy)
Taxa.phy$ASV<-row.names(Taxa.phy)
res_phy<-dplyr::left_join(res_phy, Taxa.phy)
res_phy$Kingdom<-NULL
res_phy$Genus<-NULL
res_phy$Species<-NULL
res_phy$Class<-NULL
res_phy$Order<-NULL
res_phy$Family<-NULL
}
View(res_phy)
write.csv2(res_phy, file="DA_DEseq2/Analisi diff dei phyla con DeSeq2 _ ratio Germany vs Italia.csv", row.names = F, quote=F, na = "")
# box plot
{target<-res_phy$ASV
target<-prune_taxa(target,data_pruned.phy)
target<- subset_taxa(target, Phylum!="uncultured") # removing NAs from plot
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_p<-tabella
}
  
############ TO CREATE AN UNIQUE PLOT OF ALL DESEQ2 RESULTS ############

tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_f$Taxa<-"Families"
tabella_f[,c("Phylum","Order","Class")]<-NULL
colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"

tabella_o$Taxa<-"Orders"
tabella_o[,c("Phylum","Class")]<-NULL
colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"

tabella_c$Taxa<-"Classes"
tabella_c[,"Phylum"]<-NULL
colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"

tabella_p$Taxa<-"Phyla"
colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"

# manteining only NOT redundant results (eg. same genus in one family, one family in one order and so on...)
{Good_fam<-c("Bifidobacteriaceae", "Acidaminococcaceae", "Lachnospiraceae", "Ruminococcaceae", "Rikenellaceae", "Erysipelotrichaceae", "Erysipelatoclostridiaceae", "Methanomicrobiaceae", "Sutterellaceae", "Desulfovibrionaceae")
Good_ord<-c("Erysipelotrichales","Oscillospirales","Bacteroidales","Sphingomonadales")
Good_class<-c("Actinobacteria","Alphaproteobacteria","Clostridia","Bacilli")
Good_phy<-c("Actinobacteriota","Firmicutes","Proteobacteria")

tabella_f_good<-tabella_f[tabella_f$Bacteria %in% Good_fam, ]
tabella_o_good<-tabella_o[tabella_o$Bacteria %in% Good_ord, ]
tabella_c_good<-tabella_c[tabella_c$Bacteria %in% Good_class, ]
tabella_p_good<-tabella_p[tabella_p$Bacteria %in% Good_phy , ]
}

tabella_tot<-rbind.data.frame(tabella_g,tabella_f_good,tabella_c_good,tabella_o_good,tabella_p_good)
tabella_tot$Taxa<-factor(tabella_tot$Taxa, levels = c("Phyla","Classes","Orders","Families","Genera"))

library(ggh4x) # allows to use facet_grid2 function
ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill= Nazionality)) + facet_grid2(rows=vars(Taxa), independent = "x", scales = "free") +
  geom_boxplot(width=0.7) + theme_bw(base_size = 13) + theme(strip.text.y=element_text(size=15,colour="black")) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,axis.text.x = element_text(angle = -55, vjust=0.5, hjust=0, size=11), axis.text.y = element_text(size=12)) +
  scale_x_discrete(expand=c(-0.3, 0)) + theme(plot.title= element_text(size=18) ,legend.key.size=unit(0.7,"cm"), legend.text=element_text(size=14)) +
  labs(title= "Differently abundant Taxa", y="Sqrt Abundance", fill="", x="") + scale_y_sqrt( )
ggsave(filename = "DA_DEseq2/Differenze totali Nazionalit.png", width = 10, height = 13, dpi=300)
dev.off()

############### CORRELATIONS BETWEEN CD4/CD 8 AND DESEQ2 RESULTS ##################

dir.create("Spearman_genera_DESeq2_vs_CD")
setwd("Spearman_genera_DESeq2_vs_CD")

CD <- as.data.frame(readxl::read_excel("../../Raw_counts_CD.xlsx"))
row.names(CD)<-CD$...1
CD$...1<-NULL
# mantaining only samples that have CD8 too in order to normalize through TSS
No_CD8<-na.omit(CD)
CD<-No_CD8
rm(No_CD8)
CD$'ratio CD4/CD8'<-CD$`CD4+`/CD$`CD8+`

Generi<-subset(res_genus, ! Genus %in% c("uncultured","",NA))
Selection_Generi<- Generi$Genus

#Genus
rm(X)
{X<-subset_taxa(data.genus.prop, Genus %in% Selection_Generi)
a<-as.data.frame(tax_table(X))
X<-otu_table(X)
identical(row.names(a),row.names(X))
row.names(X)<-a$Genus
X<-as.data.frame(t(X))
head(X, n=2)
row.names(X)<-gsub("I","LC",row.names(X))
row.names(X)<-gsub("G","LC",row.names(X))
X<-X[row.names(CD),] 
# only DA Genera, sample with the same order of CD
}
head(X) # OK!
X<-cbind.data.frame(X,CD)

library(Hmisc)
{x<-as.matrix(X) # NB: lower case   
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)                              
data_corr<-as.data.frame(as.table(r$r))                      
data_pvalue<-as.data.frame(as.table(r$P))
identical(data_corr[,1:2],data_pvalue[,1:2])                   
Correlations<-cbind(data_corr,data_pvalue[,3])
colnames(Correlations)<-c("Bacteria","CD","Corr","pvalue")
Correlations<-subset(Correlations, Bacteria %in% Selection_Generi)                
Correlations<-subset(Correlations, ! CD %in% Selection_Generi)                    
Correlations$padj<-p.adjust(Correlations$pvalue, method = "holm") 
}
write.csv2(Correlations,file="Correlation_Generi_CD.csv",quote = F,row.names = F)

{Correlations$Sign<-Correlations$pvalue
Correlations$Sign[Correlations$Sign<=0.05]<-"*"
Correlations$Sign[Correlations$Sign>0.05]<-""
}
ggplot(Correlations, aes(x = Correlations$CD, y = Correlations$Bacteria, fill = Correlations$Corr)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) + scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) + theme(axis.text.x=element_text(angle = 0, size= 12)) +
  guides(fill= guide_colourbar(title ="Correlation value")) +  
  geom_text(aes(label= Correlations$Sign), color= "white", size =6, nudge_y = -0.4) +
  labs(title = "Spearman correlations between \n genera resulted from DeSeq2 and CD count", 
       y= "", x= "", caption="p-value lower than 0.05 are signed with *") +   
  theme(plot.title = element_text(size=18)) + theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) +            #(dimensioni label gruppi)
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))   
ggsave(file="Correlation_between_Genera_and_CD.png", width = 8, height = 7, dpi=300)

setwd("..")

############# HANDLING PICRUST OUTPUT FOR LEFSE PLOT ##################

a <- read.delim("../PICRUST2_EPA/path_abun_unstrat_descrip.tsv.gz") # Metacyc results, after description step
colnames(a)<-gsub("ID1946.16S.[0-9][0-9].","",colnames(a))
colnames(a)<-gsub("ID1946.16S.[0-9].","",colnames(a))
colnames(a)

Descriptions<-a[,c("pathway","description")]

a<-a[, ! colnames(a) %in% c("pathway","description")]
metadata<-as(sample_data(data.prop),"data.frame")
metadata$Sample<-gsub("I","LC",metadata$Sample)
metadata$Sample<-gsub("G","LC",metadata$Sample)
identical(metadata$Sample,colnames(a)) # TRUE
a<-rbind.data.frame(metadata$Nazionality,a)
head(a)

# following modifies are necessary for LEFSE Galaxy software (those characters cause error in plot)
{ Descriptions$new_description<-gsub(",","_", Descriptions$description)
Descriptions$new_description<-gsub(" ","_", Descriptions$new_description)
Descriptions$new_description<-gsub("-","_", Descriptions$new_description)
Descriptions$new_description<-gsub("(","_", Descriptions$new_description,fixed = T)
Descriptions$new_description<-gsub(")","_", Descriptions$new_description, fixed = T)
Descriptions$new_description<-gsub("'","", Descriptions$new_description, fixed = T)
Descriptions$new_description<-gsub("._","_", Descriptions$new_description, fixed = T)
Descriptions$new_description<-gsub(".","_", Descriptions$new_description, fixed = T)
Descriptions$new_description<-gsub(":","_", Descriptions$new_description, fixed = T)
Descriptions$new_description<-gsub("/","or", Descriptions$new_description, fixed = T)
Descriptions$new_description<-gsub(";","_", Descriptions$new_description)
Descriptions$new_description<-gsub("&","_", Descriptions$new_description)
Descriptions$new_description<-gsub("__","_", Descriptions$new_description)
}
for(x in 1:length(Descriptions$new_description)) {
  if(stringr::str_detect( ( strtrim(Descriptions$new_description, width=1) [x]),"_")) {Descriptions$new_description[x]<-substring(Descriptions$new_description[x], first = 2) }
}

a<-cbind.data.frame(c("Nazionality",Descriptions$new_description),a)
colnames(a)[colnames(a)=='c("Nazionality", Descriptions$new_description)']<-"Sample"

dir.create("PICRUST2_LEFSE")
setwd("PICRUST2_LEFSE")
write.table(a,file="Output_ready_for_LefSe.tsv",row.names = F,quote = F, sep="\t")

# --> LEFSe (Galaxy)
# <--
Significative_functions_LEFSE<- read.delim("Result.lefse_internal_res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]

for(x in 1:length(Significative_functions_LEFSE$Pathway)) {
  if(stringr::str_detect(Significative_functions_LEFSE$Pathway[x],"^f_[0-9]")) {Significative_functions_LEFSE$Pathway[x]<-gsub("f_","",Significative_functions_LEFSE$Pathway[x]) }
}

Significative_functions_LEFSE<-dplyr::left_join(Significative_functions_LEFSE,Descriptions[,c("pathway","new_description")], by=c("Pathway"="new_description"))
colnames(Significative_functions_LEFSE)[colnames(Significative_functions_LEFSE)=="pathway"]<-"MetaCyc_ID"
Significative_functions_LEFSE$MetaCyc_ID[is.na(Significative_functions_LEFSE$MetaCyc_ID)] # there has not to be NA here
head(Significative_functions_LEFSE)

write.csv2(Significative_functions_LEFSE,file = "Significative_functions_LEFSE.csv",na="",quote = F, row.names = F)
unlink("Output_ready_for_LefSe.tsv")
unlink("Result.lefse_internal_res")

setwd("..")


####################### PLOTTING LEFSE PERSONALLY ###################

setwd("Risultati/PICRUST2_LEFSE/")
library(readxl)
Significative_functions_LEFSE<-read_excel("Significative_functions_LEFSE.xlsx")

# modifing names for the plot
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
Significative_functions_LEFSE$Pathway<-gsub("_"," ", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# plotting results only over 3
Significative_functions_LEFSE<- Significative_functions_LEFSE[Significative_functions_LEFSE$logLDA_score>3, ]
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Germany"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Germany"]

ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log10 LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Germany",1,0), size=3.2) +
  scale_fill_manual(values=c("Italy"="deepskyblue", "Germany"="coral")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST2_LEFSE_plot_over_3.png", width = 8, height = 3.5, dpi = 300)

setwd("..")