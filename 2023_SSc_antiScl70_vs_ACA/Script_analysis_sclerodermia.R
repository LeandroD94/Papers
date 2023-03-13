# setwd("~/Desktop/")

{library("phyloseq")
  library("ggplot2")
  library("ggpubr")
  library("dendextend")
  library("vegan")
  library("DESeq2")
}

################## IMPORTING DATA ###################

data = import_biom("table.biom") # otu and taxa table from MICCA
colnames(tax_table(data))<-c("Domain","Phylum","Class","Order","Family","Genus")  

# building sample data from names
Sample_ID<-sample_names(data)
sample<-as.data.frame(t(as.data.frame(strsplit(Sample_ID, split = "_"))))
row.names(sample)<-Sample_ID
colnames(sample)<-c("Sample_type","ID","Type")
{sample$Sample_type<-gsub("S","Saliva",sample$Sample_type)
  sample$Sample_type<-gsub("C","Skin",sample$Sample_type)
  sample$Sample_type<-gsub("F","Feces",sample$Sample_type) 
}
identical(sample_names(data),row.names(sample)) # TRUE
write.csv2(sample,file="metadata.csv", quote=F, row.names = T)

sample_data(data)<-sample
{sample_data(data)$Ig<-sample_data(data)$Type # renaming with Ig names
  sample_data(data)$Ig<-gsub("L","ACA", sample_data(data)$Ig)
  sample_data(data)$Ig<-gsub("D","anti_Scl70", sample_data(data)$Ig)
  sample_data(data)$Ig<-factor(sample_data(data)$Ig,levels = c("anti_Scl70","ACA"))
}

identical(as.numeric(length(sample_names(data))),78) # TRUE (Expecting 78 samples)

head(otu_table(data))
head(tax_table(data))
head(sample_data(data))


data_tot<-data # needed to subset later

# save.image("data.RData)     ### <- saved here
# load("data.RData")


# After the correction of a reviewer, the following sample has been discarded and the related group has been renominated
data<- subset_samples(data, ID != "06")
{sample_data(data)$Ig<-sample_data(data)$Type # renaming with Ig names
  sample_data(data)$Ig<-gsub("L","ACA", sample_data(data)$Ig)
  sample_data(data)$Ig<-gsub("D","anti_Scl70", sample_data(data)$Ig)
  sample_data(data)$Ig<-factor(sample_data(data)$Ig,levels = c("anti_Scl70","ACA"))
}
data_tot<-data


#################### STARTING ANALYSIS ##############

dir.create("Results")
setwd("Results")


{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
 data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}


############### GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(data)==1, data) # no "1" found
length(which(taxa_sums(data)==1))
# no singletons

################# GENERAL ALFA DIVERSITY #########################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

dir.create("General_Alpha_div")
setwd("General_Alpha_div")

{pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), nrow=1, x= "Sample_type", shape=NULL, color="ID")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev)
pAlpha$data<-New_data
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Sample_type, y=value, color=NULL), alpha=0.1) +
  theme_bw() + labs(x="", color="Sample_ID") +
  stat_compare_means(aes(group = Sample_type), method = "anova", label.x= 1.5, size=3, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alfa_diversity_with_ANOVA.tiff", width = 7,height =5, dpi=300)

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Sample_type, y=value, color=NULL), alpha=0.1) +
  theme_bw() + labs(x="", color="Sample_ID") +
  stat_compare_means(aes(group = Sample_type, label = sprintf("p = %5.3f", as.numeric(..p.format..)) ),
                     method = "kruskal.test", label.x= 1.65, size=3, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alfa_diversity_with_Kruskal_rounded.tiff", width = 7,height =5, dpi=300)

setwd("..")


################### GENERAL PCoA BRAY CURTIS #####################

dir.create("General_Bray_Curtis_PCoA")
setwd("General_Bray_Curtis_PCoA")

# on phyla
data.prop.labels<-data.phy.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues 
  eigval<- round((eigval/sum(eigval))*100, 1) 
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type",  shape="Sample_type",) +
  geom_point(size=2.8) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis dissimilarity",  color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_revised.tiff", width = 8.5, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type") +
  theme_classic(base_size = 14) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis dissimilarity", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_text.tiff", width = 8, height = 5, dpi=300)

# again but on genera
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type") +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type") +
  theme_classic(base_size = 14) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_text.tiff", width = 8, height = 5, dpi=300)

setwd("..")


######################## GENERAL HIERARCHICAL CLUSTERING ###################

dir.create("General_Hierarchical_clustering")
setwd("General_Hierarchical_clustering")


tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Sample_type),as.character(sample_data(data)$Sample_type)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("Saliva","chartreuse",tabella_colore$Colore) #green
colors <- gsub("Skin","grey45",colors) 
colors <- gsub("Feces","coral",colors) #orange

c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.phy.prop))),method = "bray"))

c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
tiff(file="Hierarchical_cluster_Bray_on_sqrt_prop_phyla.tiff",width=3500,height=1800, res=300)
par(mar=c(6,2,5,0), cex.lab=0.1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure using Bray-Curtis distance on sqrt proportional phyla")
dev.off()

c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.genus.prop))),method = "bray"))

c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
tiff(file="Hierarchical_cluster_Bray_on_sqrt_prop_genera.tiff",width=3500,height=1800, res=300)
par(mar=c(6,2,5,0), cex.lab=0.1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure using Bray-Curtis distance on sqrt proportional genera")
dev.off()

setwd("..")


################### SUBSETTING: SKIN ########################

setwd("~/Desktop/Results")

dir.create("Skin")
setwd("Skin")

remove<-ls(pattern = )
remove<-remove[! remove %in% c("data_tot","sample")]
rm(list=remove)

data<-subset_samples(data_tot, Sample_type=="Skin")
sample_names(data)

sample_data(data)$Ig<-factor(sample_data(data)$Ig,levels = c("anti_Scl70","ACA"))

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}


{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}


################# COUNTS EXPORT (Skin) ############################

dir.create("Counts")
setwd("Counts")

dir.create("Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)

dir.create("TSS_counts")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="TSS_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="TSS_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="TSS_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="TSS_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="TSS_counts/counts_genus.csv",quote=F)
}

setwd("..")


################# ABUNDANCES BAR PLOT (Skin) ##########################

dir.create("Counts")
setwd("Counts")

# TOP 5 Phyla
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
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Ig),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Abundace_TOP_5_phyla_skin_NO_TITLE.tiff",width=9,height=5, dpi=300)
dev.off()

rm(tabella, prune.data.others)

# TOP 5 Generi
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Ig),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera in skin sample", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Abundace_TOP_5_Generi.tiff",width=9,height=5,dpi=300)
dev.off()

setwd("..")

################# ALFA DIVERSITY (Skin) ##################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1","Shannon","Observed"), x="Ig")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
ch<-dplyr::filter(pAlpha$data, variable=="Chao1")
# adding evanness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
{ ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(ch,H,ev)
  pAlpha$data<-New_data
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Ig, y=value, color=NULL), alpha=0.1) +
  theme_bw() + labs(x="") + theme(axis.text.x = element_text(angle = -20, size= 9, vjust = 0.5, hjust = 0.1)) +
  stat_compare_means(aes(group = Ig), method = "wilcox.test", label.x= 1.15, size=3, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alfa_diversity_in_SKIN_wilcox.tiff", width = 6,height =5, dpi=300)

################## BRAY CURTIS (Skin) #####################

dir.create("Bray_Curtis_PCoA")
setwd("Bray_Curtis_PCoA")

# on phyla
data.prop.labels<-data.phy.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues 
  eigval<- round((eigval/sum(eigval))*100, 1) 
}
sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis distance", 
       color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  theme_classic(base_size = 14) +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_text.tiff", width = 8, height = 5, dpi=300)

# again but on genera
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
  sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
}
sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  theme_classic(base_size = 14) +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_text.tiff", width = 8, height = 5, dpi=300)

setwd("..")

################## DESEQ 2 (Skin) #################

# (there are any other results beyond the class level)

DEseq_data<-phyloseq_to_deseq2(data.class, ~Ig)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Ig", "ACA", "anti_Scl70"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
{res.X<-as.data.frame(res)
res.X$OTU<-row.names(res.X)
Taxa.class<-as.data.frame(tax_table(data.class))
Taxa.class$OTU<-row.names(Taxa.class)
res.X<-dplyr::left_join(res.X, Taxa.class)
res.X$Kingdom<-NULL
res.X$Genus<-NULL
res.X$Species<-NULL
res.X$Order<-NULL
res.X$Family<-NULL
}
View(res.X)
# box plot
target<-res.X$OTU
target<-prune_taxa(target,data.class)
length(row.names(tax_table(target))) 
target<- subset_taxa(target, Class!="uncultured") 
tabella<-psmelt(target)
 #tabella$Abundance<-sqrt(tabella$Abundance)
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
tiff(file="Box_plot_diff_significative_classi_Limitata_vs_Diffusa_revisionata.tiff", width=2200, height = 1500, res=300)
ggplot(tabella, aes(x=Class, y=Abundance, fill= Ig)) +
  scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_boxplot(width=0.5) + theme_bw(base_size = 15.5) +
  theme(legend.position="right", legend.margin=margin(0, 0, 0, 0), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0)) +
  scale_y_sqrt( breaks = c(0, 100, 500, seq(1000, 5000, 1000))) +
  scale_x_discrete(expand=c(1, 0))  +
  labs(title= "Epidermal microbiota", y="Reads raw abundance", fill= "Ig")
dev.off()

tiff(file="Box_plot_diff_significative_classi_Limitata_vs_Diffusa_revisionata_NO_LEGENDA.tiff", width=1800, height = 1500, res=300)
ggplot(tabella, aes(x=Class, y=Abundance, fill= Ig)) +
  scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_boxplot(width=0.5) + theme_bw(base_size = 14) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0)) +
  scale_x_discrete(expand=c(0.5, 0))  +
  scale_y_sqrt( breaks = c(0, 100, 500, seq(1000, 5000, 1000))) +
  labs(title= "Epidermal microbiota", y="Reads raw abundance", fill= "Ig")
dev.off()



################### SUBSETTING: SALIVA ########################

setwd("~/Desktop/Results/")

dir.create("Saliva")
setwd("Saliva")

remove<-ls(pattern = )
remove<-remove[! remove %in% c("data_tot","sample")]
rm(list=remove)

data<-subset_samples(data_tot, Sample_type=="Saliva")
sample_names(data)

sample_data(data)$Ig<-factor(sample_data(data)$Ig,levels = c("anti_Scl70","ACA"))

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

################# COUNTS EXPORT (Saliva) ############################

dir.create("Counts")
setwd("Counts")

dir.create("Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)

dir.create("TSS_counts")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="TSS_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="TSS_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="TSS_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="TSS_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="TSS_counts/counts_genus.csv",quote=F)
}

setwd("..")

################# ABUNDANCES BAR PLOT (Saliva) ##########################

dir.create("Counts")
setwd("Counts")

# TOP 5 Phyla
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
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Ig),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Abundace_TOP_5_phyla_Saliva_NO_TITLE.tiff",width=9,height=5, dpi=300)
dev.off()

rm(tabella, prune.data.others)

# TOP 5 Generi
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Ig),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", 
       title = "Five most abundant genera in saliva sample of dcSSc and lcSSc patients", 
       caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Abundace_TOP_5_Genera.tiff",width=9,height=5,dpi=300)
dev.off()

setwd("..")

################# ALFA DIVERSITY (Saliva) ##################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1","Shannon","Observed"), x="Ig")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
ch<-dplyr::filter(pAlpha$data, variable=="Chao1")
# adding evanness
{ identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(ch,H,ev)
  pAlpha$data<-New_data
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Ig, y=value, color=NULL), alpha=0.1) +
  theme_bw() + labs(x="") + theme(axis.text.x = element_text(angle = -20, size= 9, vjust = 0.5, hjust = 0.1)) +
  stat_compare_means(aes(group = Ig), method = "wilcox.test", label.x= 1.15, size=3, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alfa_diversity_in_Saliva_wilcox.tiff", width = 6,height =5, dpi=300)

################## BRAY CURTIS (Saliva) #####################

dir.create("Bray_Curtis_PCoA")
setwd("Bray_Curtis_PCoA")

# on phyla
data.prop.labels<-data.phy.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  ordBC[["vectors"]][, "Axis.1"]<- - ordBC[["vectors"]][, "Axis.1"] # just an aesthetic choise
  eigval<-ordBC$values$Eigenvalues 
  eigval<- round((eigval/sum(eigval))*100, 1) 
}
sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis distance", 
       color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  theme_classic(base_size = 14) +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_text.tiff", width = 8, height = 5, dpi=300)

# again but on genera
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
  sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  theme_classic(base_size = 14) +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_text.tiff", width = 8, height = 5, dpi=300)

setwd("..")

################## DESEQ 2 (Saliva) #################

# there are no results here...


################### SUBSETTING: Feces ########################

setwd("~/Desktop/Results/")

dir.create("Feces")
setwd("Feces")

remove<-ls(pattern = )
remove<-remove[! remove %in% c("data_tot","sample")]
rm(list=remove)

data<-subset_samples(data_tot, Sample_type=="Feces")
sample_names(data)

sample_data(data)$Ig<-factor(sample_data(data)$Ig,levels = c("anti_Scl70","ACA"))

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}


{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}


################# COUNTS EXPORT (Feces) ############################

dir.create("Counts")
setwd("Counts")

dir.create("Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100) # disable scientific annotation

dir.create("TSS_counts")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="TSS_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="TSS_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="TSS_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="TSS_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="TSS_counts/counts_genus.csv",quote=F)
}

setwd("..")

################# ABUNDANCES BAR PLOT (Feces) ##########################

dir.create("Counts")
setwd("Counts")

# TOP 5 Phyla
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
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
tabella$Phylum<-factor(tabella$Phylum, levels=c("Actinobacteria","Bacteroidetes","Firmicutes","Verrucomicrobia","Others","Proteobacteria"))
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Ig),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Abundace_TOP_5_phyla_Stool_NO_TITLE.tiff",width=9,height=5, dpi=300)
dev.off()

rm(tabella, prune.data.others)

# TOP 5 Generi
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Ig),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", 
       title = "Five most abundant genera in stool samples", 
       caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Abundace_TOP_5_Generi.tiff",width=9,height=5,dpi=300)
dev.off()

setwd("..")

################# ALFA DIVERSITY (Feces) ##################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1","Shannon","Observed"), x="Ig")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
ch<-dplyr::filter(pAlpha$data, variable=="Chao1")
# adding evanness
{ identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(ch,H,ev)
  pAlpha$data<-New_data
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Ig, y=value, color=NULL), alpha=0.1) +
  theme_bw() + labs(x="") + theme(axis.text.x = element_text(angle = -20, size= 9, vjust = 0.5, hjust = 0.1)) +
  stat_compare_means(aes(group = Ig, label = sprintf("p = %5.3f", as.numeric(..p.format..)) ),
                     method = "wilcox.test", label.x= 1.15, size=3,  label.y.npc = "top", vjust=-0.5)
ggsave(file="Alfa_diversity_in_Feces_wilcox.tiff", width = 6,height =5, dpi=300)


################## BRAY CURTIS (Feces) #####################

dir.create("Bray_Curtis_PCoA")
setwd("Bray_Curtis_PCoA")

# on phyla
data.prop.labels<-data.phy.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues 
  eigval<- round((eigval/sum(eigval))*100, 1) 
}
sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis distance", 
       color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  theme_classic(base_size = 14) +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop Phyla with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_PHYLA_text.tiff", width = 8, height = 5, dpi=300)

# again but on genera
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
  sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
}
sample_data(data.sqrt_prop)$Ig<-gsub("anti_Scl70","anti Scl70", sample_data(data.sqrt_prop)$Ig)
plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_points.tiff", width = 8, height = 5, dpi=300)
# without ellipses
a<-plot_ordination(data.sqrt_prop, ordBC, color = "Ig") +
  theme_classic(base_size = 14) +
  scale_color_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), size=3.2, show.legend = FALSE) +
  labs(title="PCoA computed on sqrt prop genera with Bray Curtis distance", color="Sample_type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
a[["layers"]]<-a[["layers"]][-1]
a
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_GENERA_text.tiff", width = 8, height = 5, dpi=300)

setwd("..")

################## DESEQ 2 (Feces) #################

# on genera
DEseq_data<-phyloseq_to_deseq2(data.genus, ~Ig)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Ig", "ACA", "anti_Scl70")) # ratio ACA/anti_Scl70
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res

{ res<-res[order(res$padj), ]
res_genus<-as.data.frame(res)
res_genus$OTU<-row.names(res_genus)
library(dplyr)
Taxa.genus$OTU<-row.names(Taxa.genus)
res_genus<-left_join(res_genus, Taxa.genus)
res_genus$Kingdom<-NULL
res_genus$Species<-NULL
View(res_genus)
}
# box_plot
{target<-res_genus$OTU
target<-prune_taxa(target,data.genus)
length(row.names(tax_table(target))) 
temp<-tax_table(target)[, "Genus"]
temp[is.na(temp)] <- "NA_Acidaminococcaceae"
tax_table(target)[, "Genus"]<-temp
tabella<-psmelt(target)
#tabella$Abundance<-sqrt(tabella$Abundance)
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
}
box_g<-ggplot(tabella, aes(x=Genus, y=Abundance, fill= Ig)) +
  geom_boxplot(width=0.4, show.legend = FALSE) +
  scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0)) + 
  scale_x_discrete(expand=c(0.2, 0)) + theme_bw(base_size = 13) +
  labs(title= "Stool microbiota", y="Reads raw abundance")

# on families
DEseq_data<-phyloseq_to_deseq2(data.fam, ~Ig)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Ig", "ACA", "anti_Scl70"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
{res<-res[order(res$padj), ]
res_fam<-as.data.frame(res)
res_fam$OTU<-row.names(res_fam)
Taxa.fam$OTU<-row.names(Taxa.fam)
res_fam<-left_join(res_fam, Taxa.fam)
res_fam$Kingdom<-NULL
res_fam$Genus<-NULL
res_fam$Species<-NULL
}
# box_plot
{target<-res_fam$OTU
target<-prune_taxa(target,data.fam)
length(row.names(tax_table(target))) 
target<- subset_taxa(target, Family!="uncultured")
tabella<-psmelt(target)
#tabella$Abundance<-sqrt(tabella$Abundance)
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
}
box_f<-ggplot(tabella, aes(x=Family, y=Abundance, fill= Ig)) + 
  geom_boxplot(width=0.6, show.legend = F) +
  scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0)) + 
  scale_x_discrete(expand=c(0.2, 0)) + theme_bw(base_size = 13) +
  labs(title= "", y="")


# # on orders   # SAME AS FAMILIES BUT WITH A CLEAR FALSE POSITIVE (DUE TO AN OUTLIER) AMONG RESULTS
# DEseq_data<-phyloseq_to_deseq2(data.order, ~Ig)
# DE<-DESeq(DEseq_data)
# res<-results(DE, contrast= c("Ig", "ACA", "anti_Scl70"))
# res = res[order(res$padj, na.last=NA), ]
# res<-res[(res$padj < 0.05), ]
# res
# {res<-res[order(res$padj), ]
#   res_order<-as.data.frame(res)
#   res_order$OTU<-row.names(res_order)
#   Taxa.order$OTU<-row.names(Taxa.order)
#   res_order<-left_join(res_order, Taxa.order)
#   res_order$Kingdom<-NULL
#   res_order$Genus<-NULL
#   res_order$Species<-NULL
# }
# # box_plot
# {target<-res_fam$OTU
#   target<-prune_taxa(target,data.order)
#   length(row.names(tax_table(target))) 
#   target<- subset_taxa(target, Order!="uncultured")
#   tabella<-psmelt(target)
#   #tabella$Abundance<-sqrt(tabella$Abundance)
#   tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
#   tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
# }
# box_o<-ggplot(tabella, aes(x=Order, y=Abundance, fill= Ig)) + 
#   geom_boxplot(width=0.6, show.legend = F) +
#   scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0)) + 
#   scale_x_discrete(expand=c(0.2, 0)) + theme_bw(base_size = 13) +
#   labs(title= "", y="")


# # on class     # SAME AS PHYLUM BUT WITH A CLEAR FALSE POSITIVE AMONG THE RESULTS
# DEseq_data<-phyloseq_to_deseq2(data.class, ~Ig)
# DE<-DESeq(DEseq_data)
# res<-results(DE, contrast= c("Ig", "ACA", "anti_Scl70"))
# res = res[order(res$padj, na.last=NA), ]
# res<-res[(res$padj < 0.05), ]
# res
# {res<-res[order(res$padj), ]
#   res_class<-as.data.frame(res)
#   res_class$OTU<-row.names(res_class)
#   Taxa.class$OTU<-row.names(Taxa.class)
#   res_class<-left_join(res_class, Taxa.class)
#   res_class$Kingdom<-NULL
#   res_class$Genus<-NULL
#   res_class$Species<-NULL
# }
# # box_plot
# {target<-res_fam$OTU
#   target<-prune_taxa(target,data.class)
#   length(row.names(tax_table(target))) 
#   target<- subset_taxa(target, Class!="uncultured")
#   tabella<-psmelt(target)
#   #tabella$Abundance<-sqrt(tabella$Abundance)
#   tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
#   tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
# }
# box_c<-ggplot(tabella, aes(x=Class, y=Abundance, fill= Ig)) + 
#   geom_boxplot(width=0.6, show.legend = F) +
#   scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0)) + 
#   scale_x_discrete(expand=c(0.2, 0)) + theme_bw(base_size = 13) +
#   labs(title= "", y="")



# on phyla
DEseq_data<-phyloseq_to_deseq2(data.phy, ~Ig)
DE<-DESeq(DEseq_data)
res<-results(DE, contrast= c("Ig", "ACA", "anti_Scl70"))
res = res[order(res$padj, na.last=NA), ]
res<-res[(res$padj < 0.05), ]
res
{res<-res[order(res$padj), ]
res_phy<-as.data.frame(res)
res_phy$OTU<-row.names(res_phy)
Taxa.phy$OTU<-row.names(Taxa.phy)
res_phy<-dplyr::left_join(res_phy, Taxa.phy)
res_phy$Kingdom<-NULL
res_phy$Genus<-NULL
res_phy$Species<-NULL
res_phy$Class<-NULL
res_phy$Order<-NULL
res_phy$Family<-NULL
View(res_phy)
}
# box_plot
{target<-res_phy$OTU
target<-prune_taxa(target,data.phy)
length(row.names(tax_table(target)))
target<- subset_taxa(target, Phylum!="uncultured")
tabella<-psmelt(target)
#tabella$Abundance<-sqrt(tabella$Abundance)
tabella$Ig<-gsub("anti_Scl70","anti Scl70", tabella$Ig)
tabella$Ig<-factor(tabella$Ig, levels=c("anti Scl70","ACA"))
}
box_p<-ggplot(tabella, aes(x=Phylum, y=Abundance, fill= Ig)) + 
  geom_boxplot(width=0.6) + 
  scale_fill_manual(values=c("ACA"="deepskyblue","anti Scl70"="coral")) +
  theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0),axis.text.x = element_text(angle = 90, hjust = 1, vjust=0)) +
  scale_x_discrete(expand=c(0.2, 0)) + theme_bw(base_size = 13) +
  labs(title= "", y="")

#install.packages("patchwork")
library(patchwork) 
tiff(file="Taxa_differenti_nelle_feci_side_by_side.tiff", width = 2600, height = 1200, res = 300)
box_g + box_f + box_p # to plot them side by side
dev.off()

# there are any other results in other levels


####################### PICRUST AND LEFSE ##########################

a <- read.delim("../PICRUST/path_abun_unstrat_descrip.tsv.gz") # output of Picrust, in the starting directory
colnames(a)

sample[sample$ID=="06",] # this sample has to be removed after reviewers rightful criticism
a<- a[ , ! colnames(a) %in% c("C_06_D","F_06_D","S_06_D") ]

Descriptions<-a[,c("pathway","description")]

a<-a[ ,! colnames(a) %in% c("pathway","description")]
metadata<-as(sample_data(data_tot),"data.frame")
rownames(metadata)<-metadata$Codice
metadata<-metadata[colnames(a),]
identical(metadata$Codice,colnames(a))
a<-rbind.data.frame(as.character(metadata$Ig),a)
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

a<-cbind.data.frame(c("Ig",Descriptions$new_description),a)
colnames(a)[colnames(a)=='c("Ig", Descriptions$new_description)']<-"Sample"
head(a)

library("stringr")
feces<-a[ ,colnames(a)=="Sample" | str_starts(colnames(a), pattern="F") ]
cute<-a[ ,colnames(a)=="Sample" | str_starts(colnames(a), pattern="C") ]

dir.create("PICRUST2_LEFSE")
setwd("PICRUST2_LEFSE")
write.table(feces,file="Output_pronto_per_LefSe_stool.tsv",row.names = F,quote = F, sep="\t")
write.table(cute,file="Output_pronto_per_LefSe_skin.tsv",row.names = F,quote = F, sep="\t")


########################## SKIN

# --> LEFSe 
# <--
Significative_functions_LEFSE<- read.delim("Risultato_LEFSE_Skin.lefse_internal_res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]

for(x in 1:length(Significative_functions_LEFSE$Pathway)) {
  if(stringr::str_detect(Significative_functions_LEFSE$Pathway[x],"^f_[0-9]")) {Significative_functions_LEFSE$Pathway[x]<-gsub("f_","",Significative_functions_LEFSE$Pathway[x]) }
}

Significative_functions_LEFSE<-dplyr::left_join(Significative_functions_LEFSE,Descriptions[,c("pathway","new_description")], by=c("Pathway"="new_description"))
colnames(Significative_functions_LEFSE)[colnames(Significative_functions_LEFSE)=="pathway"]<-"MetaCyc_ID"
Significative_functions_LEFSE$MetaCyc_ID[is.na(Significative_functions_LEFSE$MetaCyc_ID)]
head(Significative_functions_LEFSE)

write.csv2(Significative_functions_LEFSE,file = "Significative_functions_LEFSE_SKIN.csv",na="",quote = F, row.names = F)

# Significative_functions_LEFSE<- read.delim("Significative_functions_LEFSE.csv", header=T, sep=";")
# Significative_functions_LEFSE$logLDA_score<- as.numeric(gsub(",",".",Significative_functions_LEFSE$logLDA_score) )

# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{ Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
  Significative_functions_LEFSE$Class_with_highest_mean<-gsub("L","ACA",Significative_functions_LEFSE$Class_with_highest_mean)
  Significative_functions_LEFSE$Class_with_highest_mean<-gsub("D","anti Scl70",Significative_functions_LEFSE$Class_with_highest_mean)
  }
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="ACA"]<- - Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="ACA"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log10 LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean!="ACA",1,0), size=3.2) +
  scale_fill_manual(values=c("ACA"="deepskyblue", "anti_Scl70"="coral2")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST2_LEFSE_plot_METACYC_SKIN.png", width = 9, height = 3, dpi = 300)


###################### STOOL

# --> LEFSe 
# <--
Significative_functions_LEFSE<- read.delim("Risultato_LEFSE_Stool.lefse_internal_res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]

for(x in 1:length(Significative_functions_LEFSE$Pathway)) {
  if(stringr::str_detect(Significative_functions_LEFSE$Pathway[x],"^f_[0-9]")) {Significative_functions_LEFSE$Pathway[x]<-gsub("f_","",Significative_functions_LEFSE$Pathway[x]) }
}

Significative_functions_LEFSE<-dplyr::left_join(Significative_functions_LEFSE,Descriptions[,c("pathway","new_description")], by=c("Pathway"="new_description"))
colnames(Significative_functions_LEFSE)[colnames(Significative_functions_LEFSE)=="pathway"]<-"MetaCyc_ID"
Significative_functions_LEFSE$MetaCyc_ID[is.na(Significative_functions_LEFSE$MetaCyc_ID)]
head(Significative_functions_LEFSE)

write.csv2(Significative_functions_LEFSE,file = "Significative_functions_LEFSE_STOOL.csv",na="",quote = F, row.names = F)


# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
Significative_functions_LEFSE$Class_with_highest_mean<-gsub("L","ACA",Significative_functions_LEFSE$Class_with_highest_mean)
Significative_functions_LEFSE$Class_with_highest_mean<-gsub("D","anti Scl70",Significative_functions_LEFSE$Class_with_highest_mean)

# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log10 LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="D",1,0), size=3.2) +
  scale_fill_manual(values=c("ACA"="deepskyblue", "anti_Scl70"="coral2")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST2_LEFSE_plot_METACYC_STOOL.png", width = 8, height = 3, dpi = 300)


#### returning to results section
setwd("..")


####################################################################

package<-sessionInfo( )
package
