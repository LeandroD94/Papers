################### PREPARING THE ENVIRONMENT ###############

{library("phyloseq")
library("xlsx")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("dendextend")
library("vegan")
library("DESeq2")
library("Hmisc")
library("qiime2R")
}

####################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME2/table.qza", taxonomy="QIIME2/taxonomy.qza")
metadata <- read.csv("metadata.csv")
# updating metadata
sample<-as.data.frame(sample_names(data))
colnames(sample)<-"FASTQ_code"
original_names<-sample_names(data) # needed in order to check matching errors
sample$FASTQ_code<-gsub("ID2115-[0-9][0-9]-","",sample$FASTQ_code)
sample$FASTQ_code<-gsub("ID2115-[0-9]-","",sample$FASTQ_code)
sample<-dplyr::left_join(sample, metadata)
row.names(sample)<-sample$Sample_name
sample_names(data)<-sample$Sample_name
sample_data(data)<-sample

identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names))) # no sample cutted out
head(sample_data(data))

rm(original_names,sample)

# <-- save here!
# save.image("data.RData")

############ CHECKING UNASSIGNED IN PROKARYOTE KINGDOM ######################

Unass<-tax_glom(data, taxrank = "Kingdom") # or domain
{Unass.prop<-transform_sample_counts(Unass, function (x) (x/sum(x)*100) )
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
e
options(scipen=100)
write.csv2(e[,colnames(e)!="Kingdom"], file="QIIME2/Unassigned_domain_checking.csv", row.names = T, quote = F)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)

###################### CHECKING SEQUENCING QUALITY #######################

dir.create("Checking_results_QIIME")

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
Taxa.genus<-as.data.frame(tax_table(data.genus))
Taxa.phy<-as.data.frame(tax_table(data.phy))
}
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Checking_results_QIIME/Raw_abundances_phyla.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Checking_results_QIIME/Raw_abundances_genera.csv",quote=F)

# TOP 5 Phylum 
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.prop)
table_top<-psmelt(prune.dat_top5)
table_others<-psmelt(prune.data.others)
table_others$Phylum<-"Others"
table<-rbind.data.frame(table_top,table_others)
}
ggplot(data=table, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Condition),scales = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Sample", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' comprends every phylum below rank 5 ")
ggsave(file="Checking_results_QIIME/TOP5phyla.png",width=9,height=5, dpi=300)
dev.off()

# TOP 10 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:10]
prune.dat_top <- prune_taxa(top,data.genus.prop)
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
table_top<-psmelt(prune.dat_top)
table_others<-psmelt(prune.data.others)
table_others$Genus<-"Others"
table<-rbind.data.frame(table_top,table_others)
}
ggplot(data=table, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Samples", y="Relative abundance", title = "Ten most abundant genera", caption = " 'Others' comprends every genus below rank 5 ")
ggsave(file="Checking_results_QIIME/Abundance_TOP10_Genera.png",width=9,height=5,dpi=300)
dev.off()

rm(table_others, table_top, top, tax_selected, prune.dat_top5, prune.data.others,)

# PCoA Bray
data.prop.labels<-data.prop
# not used sqrt prop normalization because otherwise S22 and S21 would have been perfectly overlapped
DistBC = phyloseq::distance(data.prop.labels, method = "bray")
ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues 
eigval<- round((eigval/sum(eigval))*100, 1) 
plot_ordination(data.prop.labels, ordBC, color = "Condition") +
  geom_point(size=4) + theme_classic(base_size = 14) + stat_ellipse() + geom_text(aes(label=sample_names(data.prop.labels)), color="black", size=2.3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on proportional ASV", color="Condition", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) +
  ggforce::facet_zoom(xlim=c(0.125,0.175), zoom.size = 0.2)
ggsave(file="Checking_results_QIIME/PCoA Beta diversity Bray Curtis.png", width = 10, height = 8, dpi=300)

# alpha diversity
pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition")
pAlpha + geom_text(aes(label= Sample_name), color= "red", size =5) + theme_bw(base_size = 14)
ggsave(filename = "Checking_results_QIIME/Alpha_div.png", width = 7, height = 5, dpi=300)

###################### PREPARATION OF THE DATA #######################

to_remove<-ls(pattern ="")
to_remove<-to_remove[! to_remove %in% c("data", "metadata")]
rm(list=to_remove)
head(sample_data(data))

# removing controls and contaminated sample
data_totale<-data
data<-subset_samples(data, Condition != "Paraffin_control")
data<-subset_samples(data, Sample_name != "HPV-6") # contaminated by human DNA
identical(as.numeric(length(sample_names(data))),19)

sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels = c("CRC+HPV","CRC"))

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}
{data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}
{Taxa.genus<-as.data.frame(tax_table(data.genus))
Taxa.fam<-as.data.frame(tax_table(data.fam))
Taxa.phy<-as.data.frame(tax_table(data.phy))
Taxa.class<-as.data.frame(tax_table(data.class))
Taxa.order<-as.data.frame(tax_table(data.order))
}
# adding informations to missing names
taxa_temp<-Taxa.genus
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
  taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
  taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
  taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
  taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
Taxa.genus.update<-taxa_temp
}
rm(taxa_temp)

dir.create("Results")
setwd("Results")

#################### % ASSIGNED IN DATABASE #########################

a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assegnati<-rbind.data.frame(a,b,c,d,e)
colnames(assegnati)<-c("Total","Assigned","%","Taxa")
assegnati
write.csv2(assegnati,file="../QIIME2/Percentuals_assigned_taxa_in database_post_HPV_6_removal.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assegnati)

########################### COUNTS EXPORT ##########################################

dir.create("Abundances")

dir.create("Abundances/Raw_abundances")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Abundances/Raw_abundances/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Raw_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Abundances/Raw_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Abundances/Raw_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Abundances/Raw_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Raw_abundances/counts_genus.csv",quote=F)
}

options(scipen = 100) # disable scientific annotation
dir.create("Abundances/Relat_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Relat_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Abundances/Relat_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Abundances/Relat_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Abundances/Relat_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Relat_abundances/counts_genus.csv",quote=F)
write.csv2(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Abundances/Relat_abundances/counts_asv.csv",quote=F)
write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Relat_abundances/counts_genus.xlsx",row.names = F,col.names = T, showNA = F)
write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Relat_abundances/counts_phylum.xlsx",row.names = F,col.names = T, showNA = F)
}

###################### ABUNDANCES BAR PLOT ##########################

dir.create("Abundances")

# TOP 5 Phylum con stacked bar plot
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.prop)
table_top<-psmelt(prune.dat_top5)
table_others<-psmelt(prune.data.others)
table_others$Phylum<-"Others"
table<-rbind.data.frame(table_top,table_others)
}
ggplot(data=table, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Condition),scales = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) + 
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' comprends every phylum below rank 5 ")
ggsave(file="Abundances/TOP5phyla.png",width=8,height=5, dpi=300)
dev.off()

rm(table)

# TOP 5 Generi
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.genus.prop)
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_selected<-Taxa.genus.update[row.names(tax_selected),]
tax_table(prune.dat_top)<-as.matrix(tax_selected)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
table_top<-psmelt(prune.dat_top)
table_others<-psmelt(prune.data.others)
table_others$Genus<-"Others"
table<-rbind.data.frame(table_top,table_others)
}
ggplot(data=table, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera", caption = " 'Others' comprends every genus below rank 5 ")
ggsave(file="Abundances/TOP5Generi.png",width=8,height=5,dpi=300)
dev.off()

rm(table_others, table_top, top, tax_selected, prune.dat_top5, prune.data.others)

######################## HIERARCHICAL CLUSTERING ###################

dir.create("Hierarchical_clustering")

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
c<-as.dendrogram(c)
table_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
colnames(table_colore)<-c("Gruppo","Colore")
colors <- gsub("CRC+HPV","coral",table_colore$Colore, fixed = T)
colors <- gsub("CRC","cyan",colors)
colors
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Hierarchical_clustering/Hierarchical cluster Euclidean on sqrt prop normalized ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure using Euclidean distance on sqrt proportional ASVs",
     sub="CRC+HPV = red     CRC = blue")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
c<-as.dendrogram(c)
table_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
colnames(table_colore)<-c("Gruppo","Colore")
colors <- gsub("CRC+HPV","coral",table_colore$Colore, fixed = T)
colors <- gsub("CRC","cyan",colors)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Hierarchical_clustering/Hierarchical cluster Bray on sqrt proportional ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure using Bray-Curtis distance on sqrt proportional ASVs",
     sub="CRC+HPV = red     CRC = blue")
dev.off()

rm(colors, c, table_colore)

########################## ALFA DIVERSITY ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition")
{H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
identical(H$Sample_name, obs$Sample_name) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating
New_data<-rbind.data.frame(obs,H,ev)
pAlpha$data<-New_data
}
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Condition", title="Alpha diversity between CRC and CRC+HPV patients") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.45)
ggsave(file="Alfa_diversity.png", width = 6,height =6, dpi=300)

# just to check the plotted p value
alphadt<- as.data.frame(pAlpha$data)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

rm(ev, H, obs, New_data, Obser_value, pAlpha, factor, alphadt)

########################### BETA DIVERSITY BRAY CURTIS #######################

dir.create("Bray_C_Beta_diversity")

{ASV.prop<-as.data.frame(otu_table(data.prop))
ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
perm_ASV_Bray<-perm_ASV # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
perm_g_Bray<-perm_g # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
beta
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
write.csv2(beta, file="Bray_C_Beta_diversity/Beta_div_bray_CONDITION.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
disper<-vegan::betadisper(BC.dist,sample_data(data)$Condition)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Bray_C_Beta_diversity/Beta_disper_permanova_Bray_on_sqrt_prop_ASV_CONDITION.csv",quote=F,row.names = T)
# on genera
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="bray")
disper<-vegan::betadisper(BC.dist,sample_data(data)$Condition)
disp_genus<-vegan::permutest(disper, permutations=9999)
disp_genus
write.csv2(disp_genus$tab, file="Bray_C_Beta_diversity/Beta_disper_permanova_Bray_on_sqrt_prop_GENERA_CONDITION.csv",quote=F,row.names = T)

rm(disp_ASV,disp_ASV,BC.dist,disper, beta)

########################### PCoA BRAY CURTIS #####################

dir.create("Bray_C_Beta_diversity")

### on ASV
{data.prop.labels<-data.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_Bray$aov.tab$`Pr(>F)`[1]))
ggsave(file="Bray_C_Beta_diversity/PCoA_Beta_div_Bray_CONDITION.png", width = 8, height = 6, dpi=300)
#without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  geom_point(size=3) + theme_classic(base_size = 14) + geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_Bray$aov.tab$`Pr(>F)`[1]))
ggsave(file="Bray_C_Beta_diversity/PCoA_Beta_div_Bray_no_ellipses_CONDITION.png", width = 8, height = 6, dpi=300)
#genders
plot_ordination(data.sqrt_prop, ordBC, color = "Gender") +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Gender", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Bray_C_Beta_diversity/PCoA_Beta_div_Bray_GENDER.png", width = 8, height = 6, dpi=300)
#CRC Type
plot_ordination(data.sqrt_prop, ordBC, color = "CRC_Type") +
geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="CRC Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Bray_C_Beta_diversity/PCoA_Beta_div_Bray_CRC_TYPE.png", width = 8, height = 6, dpi=300)
#Stage
plot_ordination(data.sqrt_prop, ordBC, color = "Stage") +
  geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="CRC Stage", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Bray_C_Beta_diversity/PCoA_Beta_div_Bray_CRC_STAGE.png", width = 8, height = 6, dpi=300)
#Age
sample_data(data.sqrt_prop)$Age_groups<-sample_data(data.sqrt_prop)$Age
sample_data(data.sqrt_prop)$Age_groups[sample_data(data.sqrt_prop)$Age_groups>= 44 & sample_data(data.sqrt_prop)$Age_groups<60]<-"40-59"
sample_data(data.sqrt_prop)$Age_groups[sample_data(data.sqrt_prop)$Age_groups>= 60 & sample_data(data.sqrt_prop)$Age_groups<79]<-"60-79"
sample_data(data.sqrt_prop)$Age_groups[sample_data(data.sqrt_prop)$Age_groups>= 80]<-"over 80"
plot_ordination(data.sqrt_prop, ordBC, color = "Age_groups") +
  geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Bray_C_Beta_diversity/PCoA_Beta_div_Bray_CRC_AGE.png", width = 8, height = 6, dpi=300)

rm(data.sqrt_prop, ordBC, DistBC, eigval, data.prop.labels)

################### DA WITH DESEQ2 #####################
dir.create("DA_DESeq2")

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

Table_tot<-NULL
Res_tot<-NULL

for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  if(t=="Genus"){ # updating missing names (NA and uncultured) but only for genus level
    taxa_temp<-as.data.frame(tax_table(d))
    for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
      taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
    for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
    for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
      taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
    for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
    for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
    tax_table(d)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "CRC", "CRC+HPV"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  if(length(res$log2FoldChange)>0){ # if there are results...
    cat(paste(length(res$log2FoldChange),"results for the",t,"level\n"))
    Sys.sleep(1)
    r<-as.data.frame(res)
    r$ASV<-row.names(r)
    Taxa.d<-as.data.frame(tax_table(d))
    Taxa.d$ASV<-row.names(Taxa.d)
    r<-dplyr::left_join(r, Taxa.d, by="ASV")
    r$Kingdom<-NULL
    r$Species<-NULL
    assign(paste(t,"results",sep="_"), r)
    write.csv2(r, file=paste0("DA_DESeq2/DA_",t,"_ratio_CRC_vs_CRC+HPV.csv"), row.names = F, quote=F, na = "")
    r_level<-r
    r_level[, "Taxon"]<- rep(t)
    Res_tot<-rbind.data.frame(Res_tot,r_level)
    ### single box plots
    target<-r[[t]]
    colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
    target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
    Table_DE<-psmelt(target)
    colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
    Table_DE$ASV<-NULL
    # Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
    assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
    ### appending to unique box plot
    index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t)
    index<- index[-length(index)] # removing the last index, regarding the taxa of interest
    Table_DE[,index]<-NULL
    Table_DE$Taxa<-t
    colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
    Table_tot<-rbind.data.frame(Table_tot, Table_DE)
  } else {
    cat("Any results for the",t,"level\n")
    Sys.sleep(1)
  }
}

View(Res_tot)
write.csv(Res_tot, file="DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Condition<-factor(Table_tot$Condition, levels=c("CRC+HPV","CRC"))
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")),
             scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("CRC"="cyan3","CRC+HPV"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", 
       fill="Condition", x="") +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2)))
ggsave(filename = "DA_DESeq2/DA_Condition_every_result.png", width = 12, height = 7.5, dpi=300)
dev.off()

# removing redundants
Table_tot2<-subset(Table_tot, ! Taxa %in% c("Phylum","Class","Order")) # to remove redundant results
Table_tot2<-subset(Table_tot2, ! Bacteria %in% c("Alcaligenaceae","Mycobacteriaceae","Synergistaceae"))

ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  scale_fill_manual(values=c("CRC"="cyan3","CRC+HPV"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=11), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", 
       fill="Condition", x="")
ggsave(filename = "DA_DESeq2/DA_Condition_no_redundants.png", width = 11, height = 7, dpi=300)
dev.off()

# for comparison with splsda ...
Genera.DESEQ2<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])

############# HANDLING PICRUST OUTPUT FOR LEFSE PLOT ##################

to_remove<-ls(pattern ="")
to_remove<-to_remove[! to_remove %in% c("data", "metadata","Genera.DESEQ2","Table_tot")]
rm(list=to_remove)

a <- read.delim("../PICRUST/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz")
colnames(a)<-gsub("ID2115.[0-9][0-9].","",colnames(a))
colnames(a)<-gsub("ID2115.[0-9].","",colnames(a))
colnames(a)

Descriptions<-a[,c("pathway","description")]

a<-a[,! colnames(a) %in% c("pathway","description")]
a<-a[,metadata$FASTQ_code] # same samples order
identical(metadata$FASTQ_code,colnames(a))
a<-rbind.data.frame(metadata$Condition,a)
a$MR16<-NULL # contaminated sample
a$MR22<-NULL # control
a$MR21<-NULL # control
head(a)

# following modifies are necessary for LEFSE Galaxy software (those characters cause error in plot)
{Descriptions$new_description<-gsub(",","_", Descriptions$description)
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

a<-cbind.data.frame(c("Condition",Descriptions$new_description),a)
colnames(a)[colnames(a)=='c("Condition", Descriptions$new_description)']<-"Sample"
a[1,]<- gsub("CRC+HPV","CRC_and_HPV",a[1,], fixed = T)

dir.create("PICRUST2_LEFSE")
write.table(a,file="PICRUST2_LEFSE/Output_ready_for_LefSe.tsv",row.names = F,quote = F, sep="\t")

# --> LEFSe (On Galaxy platform)
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

write.csv2(Significative_functions_LEFSE,file = "PICRUST2_LEFSE/Significative_functions_LEFSE.csv",na="",quote = F, row.names = F)
unlink("Result.lefse_internal_res")

################ ANALYSIS ABOUT LEUCOCYTE POPULATIONS #####################

if(grepl("Results",getwd())==F){setwd("Results")}
dir.create("CD_Statistical_Analysis")

####### data preparation
Leuco<-read.csv("../Leuco_percentuals.csv")

# # ch2anging levels (CD4/CD25) in numeric values
# {Leuco$CD4.CD25<-gsub("1S","2",Leuco$CD4.CD25)
# Leuco$CD4.CD25<-gsub("1","3",Leuco$CD4.CD25)
# Leuco$CD4.CD25<-gsub("0","1",Leuco$CD4.CD25)
# Leuco$CD4.CD25<-as.numeric(Leuco$CD4.CD25)
# }

Leuco_tumor<-Leuco[Leuco$Site=="Tumoral",colnames(Leuco)!="Site"]
Leuco_peritumor<-Leuco[Leuco$Site!="Tumoral",colnames(Leuco)!="Site"]
rownames(Leuco_peritumor)<-Leuco_peritumor$FASTQ_Code
rownames(Leuco_tumor)<-Leuco_tumor$FASTQ_Code
rownames(metadata)<-metadata$FASTQ_code
sample_order<-metadata[metadata$Condition!="Paraffin_control","FASTQ_code"]
sample_order<-sample_order[sample_order!="MR16"]
Leuco_meta<-metadata[sample_order,] # without the controls and contaminated sample
Leuco_tumor<-Leuco_tumor[sample_order,colnames(Leuco_tumor)!="FASTQ_Code"]
Leuco_peritumor<-Leuco_peritumor[sample_order,colnames(Leuco_peritumor)!="FASTQ_Code"]

# just to test the correctness of those object
if( identical(as.numeric(length(rownames(Leuco_tumor))),19) &
identical(rownames(Leuco_peritumor),rownames(Leuco_tumor)) &
identical(rownames(Leuco_peritumor),rownames(Leuco_meta)) &
! "HPV-6" %in% Leuco_meta$Sample_name ) {
  print("Everything is fine!") } else {cat("\nSomething is wrong here...\n", fill=T)}

######## Tumoral site
dir.create("CD_Statistical_Analysis/Tumoral_site")

# Mann Whitney between HPV+ and HPV-
test<-NULL
for(x in colnames(Leuco_tumor)){
  suppressWarnings(rm(temp))
  temp<-wilcox.test(Leuco_tumor[,x]~Leuco_meta$Condition, paired=F)
  temp<-cbind(x,temp$statistic,temp$p.value)
  test<-rbind.data.frame(test,temp)
  test
}
{colnames(test)<-c("CD","W_statistic","P_value")
test$P_adj_BH<-p.adjust(test$P_value, "BH")
test$Signif<-test$P_adj_BH
test[test$P_adj_BH<0.05,"Signif"]<-"*"
test[test$P_adj_BH>0.05,"Signif"]<-""
}
write.xlsx(test,file="CD_Statistical_Analysis/Tumoral_site/Mann_Whitney_between_Conditions_Tumoral_Site.xlsx", showNA = F, row.names = F, col.names = T)

# correlations with top5 Genera
top5 <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
taxa_top5 <- prune_taxa(top5,data.genus.prop)
if(identical(row.names(otu_table(taxa_top5)),row.names(tax_table(taxa_top5)))==F){
  cat("Something is wrong here...")
} 
corr_top5 <- t(as.data.frame(otu_table(taxa_top5)))
colnames(corr_top5)<-as.data.frame(tax_table(taxa_top5))$Genus
corr_top5<-corr_top5[Leuco_meta$Sample_name, ] #same order

identical(length(which(!is.na(rownames(corr_top5)))),
          length(rownames(Leuco_meta)),
          length(sample_names(data.genus.prop))) # TRUE --> no sample cutted out

suppressWarnings(rm(x,r,correlation, data_corr))
{x<-as.matrix(cbind(corr_top5,Leuco_tumor))
r<-rcorr(x, type = "spearman")
correlation<-as.data.frame(r$r)
data_corr<-as.data.frame(as.table(r$r))
data_pvalue<-as.data.frame(as.table(r$P))
if(identical(data_corr[,1:2],data_pvalue[,1:2])==F){print("Something is wrong here...")}
data_corr<-cbind(data_corr,data_pvalue[,3])
colnames(data_corr)<-c("Bacteria","CD","Corr","pvalue")
data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_top5), ]
data_corr<-data_corr[data_corr$CD %in% colnames(Leuco_tumor), ]
data_corr$padj_BH<-p.adjust(data_corr$pvalue, method = "BH")
}
write.xlsx(data_corr, file="CD_Statistical_Analysis/Tumoral_site/Spearman_TOP5_genera_vs_CD.xlsx", row.names = F, col.names = T,showNA = F)

data_corr$simbol<-data_corr$padj_BH
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$CD, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) + 
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 9)) +
  guides(fill= guide_colourbar(title ="Correlation value")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.4) +
  labs(title = "Spearman correlation between \n five most abundant genera and CDs in tumoral site",
       y= "", x= "", caption= "* = p-adjusted (BH) < 0.05") +
  theme(plot.title = element_text(size=15)) +
  theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename = "CD_Statistical_Analysis/Tumoral_site/Heatmap_TOP5Genera_vs_CD_tumor.png", width = 8, height = 6, dpi=300)

# correlations with DESeq2 Genera results
taxa_deseq <- subset_taxa(data.genus.prop, Genus %in% Genera.DESEQ2)
if(identical(row.names(otu_table(taxa_deseq)),row.names(tax_table(taxa_deseq)))==F){
  cat("Something is wrong here...")
} 
corr_deseq <- t(as.data.frame(otu_table(taxa_deseq)))
colnames(corr_deseq)<-as.data.frame(tax_table(taxa_deseq))$Genus
corr_deseq<-corr_deseq[Leuco_meta$Sample_name, ] #same order

identical(length(which(!is.na(rownames(corr_deseq)))),
          length(rownames(Leuco_meta)),
          length(sample_names(data.genus.prop))) # TRUE --> no sample cutted out

suppressWarnings(rm(x,r,correlation, data_corr))
{x<-as.matrix(cbind(corr_deseq,Leuco_tumor))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  if(identical(data_corr[,1:2],data_pvalue[,1:2])==F){print("Something is wrong here...")}
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Bacteria","CD","Corr","pvalue")
  data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_deseq), ]
  data_corr<-data_corr[data_corr$CD %in% colnames(Leuco_tumor), ]
  data_corr$padj_BH<-p.adjust(data_corr$pvalue, method = "BH")
}
write.xlsx(data_corr, file="CD_Statistical_Analysis/Tumoral_site/Spearman_DESEQ2_genera_vs_CD.xlsx", row.names = F, col.names = T,showNA = F)

data_corr$simbol<-data_corr$padj_BH
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$CD, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) + 
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 9)) +
  guides(fill= guide_colourbar(title ="Correlation value")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.4) +
  labs(title = "Spearman correlation between \n genera resulted from DA and CDs in tumoral site",
       y= "", x= "", caption= "* = p-adjusted (BH) < 0.05") +
  theme(plot.title = element_text(size=15)) +
  theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename = "CD_Statistical_Analysis/Tumoral_site/Heatmap_DESEQ_Genera_vs_CD_tumor.png", width = 8, height = 6, dpi=300)

######## Peritumoral site
dir.create("CD_Statistical_Analysis/Peritumoral_site")

test<-NULL
for(x in colnames(Leuco_peritumor)){
  suppressWarnings(rm(temp))
  temp<-wilcox.test(Leuco_peritumor[,x]~Leuco_meta$Condition, paired=F)
  temp<-cbind(x,temp$statistic,temp$p.value)
  test<-rbind.data.frame(test,temp)
  test
}
{colnames(test)<-c("CD","W_statistic","P_value")
  test$P_adj_BH<-p.adjust(test$P_value, "BH")
  test$Signif<-test$P_adj_BH
  test[test$P_adj_BH<0.05,"Signif"]<-"*"
  test[test$P_adj_BH>0.05,"Signif"]<-""
}
write.xlsx(test,file="CD_Statistical_Analysis/Peritumoral_site/Mann_Whitney_between_Conditions_Peritumoral_Site.xlsx", showNA = F, row.names = F, col.names = T)

# correlations with top5 Genera
top5 <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
taxa_top5 <- prune_taxa(top5,data.genus.prop)
if(identical(row.names(otu_table(taxa_top5)),row.names(tax_table(taxa_top5)))==F){
  cat("Something is wrong here...")
} 
corr_top5 <- t(as.data.frame(otu_table(taxa_top5)))
colnames(corr_top5)<-as.data.frame(tax_table(taxa_top5))$Genus
corr_top5<-corr_top5[Leuco_meta$Sample_name, ] #same order

identical(length(which(!is.na(rownames(corr_top5)))),
          length(rownames(Leuco_meta)),
          length(sample_names(data.genus.prop))) # TRUE --> no sample cutted out

suppressWarnings(rm(x,r,correlation, data_corr))
{x<-as.matrix(cbind(corr_top5,Leuco_peritumor))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  if(identical(data_corr[,1:2],data_pvalue[,1:2])==F){print("Something is wrong here...")}
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Bacteria","CD","Corr","pvalue")
  data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_top5), ]
  data_corr<-data_corr[data_corr$CD %in% colnames(Leuco_peritumor), ]
  data_corr$padj_BH<-p.adjust(data_corr$pvalue, method = "BH")
}
write.xlsx(data_corr, file="CD_Statistical_Analysis/Peritumoral_site/Spearman_TOP5_genera_vs_CD.xlsx", row.names = F, col.names = T,showNA = F)

data_corr$simbol<-data_corr$padj_BH
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$CD, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) + 
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 9)) +
  guides(fill= guide_colourbar(title ="Correlation value")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.4) +
  labs(title = "Spearman correlation between \n five most abundant genera and CDs in peritumoral site",
       y= "", x= "", caption= "* = p-adjusted (BH) < 0.05") +
  theme(plot.title = element_text(size=15)) +
  theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename = "CD_Statistical_Analysis/Peritumoral_site/Heatmap_TOP5Genera_vs_CD_peritumor.png", width = 8, height = 6, dpi=300)

# correlations with DESeq2 Genera results
taxa_deseq <- subset_taxa(data.genus.prop, Genus %in% Genera.DESEQ2)
if(identical(row.names(otu_table(taxa_deseq)),row.names(tax_table(taxa_deseq)))==F){
  cat("Something is wrong here...")
} 
corr_deseq <- t(as.data.frame(otu_table(taxa_deseq)))
colnames(corr_deseq)<-as.data.frame(tax_table(taxa_deseq))$Genus
corr_deseq<-corr_deseq[Leuco_meta$Sample_name, ] #same order

identical(length(which(!is.na(rownames(corr_deseq)))),
          length(rownames(Leuco_meta)),
          length(sample_names(data.genus.prop))) # TRUE --> no sample cutted out

suppressWarnings(rm(x,r,correlation, data_corr))
{x<-as.matrix(cbind(corr_deseq,Leuco_peritumor))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  if(identical(data_corr[,1:2],data_pvalue[,1:2])==F){print("Something is wrong here...")}
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Bacteria","CD","Corr","pvalue")
  data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_deseq), ]
  data_corr<-data_corr[data_corr$CD %in% colnames(Leuco_peritumor), ]
  data_corr$padj_BH<-p.adjust(data_corr$pvalue, method = "BH")
}
write.xlsx(data_corr, file="CD_Statistical_Analysis/Peritumoral_site/Spearman_DESEQ2_genera_vs_CD.xlsx", row.names = F, col.names = T,showNA = F)

data_corr$simbol<-data_corr$padj_BH
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$CD, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) + 
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 9)) +
  guides(fill= guide_colourbar(title ="Correlation value")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.4) +
  labs(title = "Spearman correlation between \n genera resulted from DA and CDs in peritumoral site",
       y= "", x= "", caption= "* = p-adjusted (BH) < 0.05") +
  theme(plot.title = element_text(size=15)) +
  theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename = "CD_Statistical_Analysis/Peritumoral_site/Heatmap_DESEQ_Genera_vs_CD_peritumor.png", width = 8, height = 6, dpi=300)

##################### R AND PACKAGES VERSION #########################

package<-sessionInfo()

con <- file("R_version_and_used_packages.txt")
sink(con, append = TRUE) 

cat(package$R.version$version.string)
cat("   running on", package$running)
cat("\n", "\n", fill=TRUE)
package$otherPkgs$DESeq2[c(1,4)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$phyloseq[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$ggplot2[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$vegan[c(1,3)]
cat("\n", "\n", fill=TRUE)
print("Dendextend")
packageVersion("dendextend")
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
rm(con)