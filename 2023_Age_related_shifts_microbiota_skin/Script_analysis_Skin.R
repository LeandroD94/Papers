##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
  library("ggplot2")
  library("vegan")
  library("ggpubr")
  library("FSA")
  library("ggh4x")
  library("egg")
  library("reshape2")
  library("dendextend")
  library("DESeq2")
  library("dplyr")
  library("stringr")
  library("Hmisc")
  library("qiime2R")
}

{ dir.create("Results")
  dir.create("Results/General_Abundances")
  dir.create("Results/Three_age")
  dir.create("Results/Dry_Skin")
  dir.create("Results/SNP")
}

options(scipen = 100)


######################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME2/table.qza", taxonomy="QIIME2/assigned_taxa.qza", metadata = "metadata_new.csv")

# updating sample names
sample<-sample_data(data)
head(sample)
row.names(sample)<-sample$Sample
sample_names(data)<-row.names(sample)
sample_data(data)<-sample
head(sample_data(data))

# save.image("data.RData")


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



############# CHECKING THE ASV ASSIGNED IN SILVA ###########

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assegnati<-rbind.data.frame(a,b,c,d,e)
colnames(assegnati)<-c("Totale","Assegnati","%","Taxa")
}
assegnati
write.csv2(assegnati,file="Results/Percentuals_ASV_assigned_in_SILVA.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assegnati)


########################### COUNTS EXPORT ##########################################

dir.create("Results/General_Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/General_Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/General_Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/General_Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/General_Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/General_Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_genus.csv",quote=F)
}


###################### CHECKING THE ASV SATURATION OF THE SAMPLES ################################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as gnerated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of he origina samples (and in the same order!!!)
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
      #			text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],labels=slope,col=check,cex=0.5)
      #			points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col=check,pch=16,cex=1)
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

png(file="Results/Sample_saturation_check_RAREFACTION_ANALYSIS.png",width=2500,height=2000, res=300)
r<-rarecurve(as(t(otu_table(data)),"matrix"), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()



######################## ABUNDANCES BAR PLOTS ##########################

# TOP 5 Phyla
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
{tabella$Three_ages<-gsub("A","20-35 (A)",tabella$Three_ages)
  tabella$Three_ages<-gsub("B","36-52 (B)",tabella$Three_ages)
  tabella$Three_ages<-gsub("C","53-68 (C)",tabella$Three_ages)
  tabella$Dry_Skin<-gsub("Dry","Dry Skin",tabella$Dry_Skin)
  tabella$Dry_Skin<-gsub("Not dry","Normal",tabella$Dry_Skin)
}
#### THREE AGES
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Three_ages),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subjects", y="Relative abundance", title = "Four most abundant phyla",
       caption = " 'Others' includes every phylum below the rank 4 ")
ggsave(file="Results/Three_age/TOP_5_phyla.png", width=8.5, height=6, dpi=300)
dev.off()
#### DRY SKIN
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  facet_grid(cols= vars(Dry_Skin),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subjects", y="Relative abundance", title = "Four most abundant phyla",
       caption = " 'Others' includes every phylum below the rank 4 ")
ggsave(file="Results/Dry_Skin/TOP_5_phyla.png", width=8.5, height=6, dpi=300)
dev.off()
#### general (No groups)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") + theme_minimal(base_size=14) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subjects", y="Relative abundance", title = "Four most abundant phyla",
       caption = " 'Others' includes every phylum below the rank 4 ")
ggsave(file="Results/General_Abundances/TOP_5_phyla.png", width=8.5, height=6, dpi=300, bg = "white")
dev.off()


# means of TOP5 phyla
write.csv(file = "Results/General_Abundances/TOP5_phyla_average_abundance.csv", row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, prune.dat_top,tabella, tabella_top)


################ ALFA DIVERSITY THREE AGES ####################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), 
                      x="Three_ages", color = "NULL")
pAlpha
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Codice,obs$Codice) # TRUE, same order
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev)
pAlpha$data<-New_data
}
# modifing names for the plot
{pAlpha$data$Three_ages<-gsub("A","20-35 (A)",pAlpha$data$Three_ages)
  pAlpha$data$Three_ages<-gsub("B","36-52 (B)",pAlpha$data$Three_ages)
  pAlpha$data$Three_ages<-gsub("C","53-68 (C)",pAlpha$data$Three_ages)
}
pAlpha$data$variable<-gsub("Observed","Observed Richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable, levels = unique(pAlpha$data$variable) )
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Three_ages, y=value), alpha=0.1) +
  theme_bw() + labs(x="Ages") +
  theme( axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
  guides(fill="none", color="none") + 
  stat_compare_means(aes(group = Three_ages), label="p.format",
                     method = "kruskal.test", label.x= 1.7, size=3.2, label.y.npc = "top", vjust=-0.5)
ggsave(file="Results/Three_age/Alpha_diversity.png", width = 6,height =4.5, dpi=300)

# post HOC
Dunnet<-dunnTest(H[,"value"]~as.factor(H$Three_ages), method="bh")
Dunnet$res  # only A-B
write.csv(Dunnet$res, file="Results/Three_age/Alpha_Div_DUNNET_TEST_Post_HOC.csv", row.names = F, quote = F)


################ ALFA DIVERSITY DRY SKIN ####################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), 
                      x="Dry_Skin", color = "NULL")
pAlpha
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Codice,obs$Codice) # TRUE, same order
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
}
pAlpha$data$variable<-gsub("Observed","Observed Richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable, levels = unique(pAlpha$data$variable) )
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Dry_Skin, y=value), alpha=0.1) +
  theme_bw() + labs(x="Dry Skin") +
  guides(fill="none", color="none") + 
  stat_compare_means(aes(group = Dry_Skin), label="p.format",
                     method = "wilcox.test", label.x= 1.29, size=3.2, label.y.npc = "top", vjust=-0.5)
ggsave(file="Results/Dry_Skin/Alpha_diversity.png", width = 6,height =4.5, dpi=300)


######################## BETA DIVERSITY THREE AGES  #######################

suppressWarnings(rm(ASV.prop))

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")

sample_OTU<-as.data.frame(t(ASV.prop))
perm_ASV<- vegan::adonis(sample_OTU ~Three_ages, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.genus.prop))
perm_g<- vegan::adonis(sample_OTU ~Three_ages, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.fam.prop))
perm_f<- vegan::adonis(sample_OTU ~Three_ages, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.class.prop))
perm_o<- vegan::adonis(sample_OTU ~Three_ages, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.order.prop))
perm_c<- vegan::adonis(sample_OTU ~Three_ages, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.phy.prop))
perm_p<- vegan::adonis(sample_OTU ~Three_ages, data=metadata, permutations = 9999, method="bray")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta


### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(ASV.prop))
pair_ASV<- pairwise.adonis(sample_ASV, factors=metadata$Three_ages, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV


# exporting beta diversity
suppressWarnings(rm(con))
con<-file("Results/Three_age/Beta_divers_general_and_pairwise_between_Three_ages.txt")
sink(con, append = TRUE)
cat("General beta diversity on Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity on Bray Curtis (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(sqrt(t(ASV.prop)), distance="bray")
  disper<-vegan::betadisper(BC.dist,metadata$Three_ages)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
}
disp_ASV$tab
disp_ASV$pairwise$permuted
# NB: this test has been already computed before and the resulting p-values have been already submitted
# then in the plots I will use the old p-values (write below) from the old permutations to mantain the same values
# (the permutations are randomics!)
# 0.008 (dispersion)
# A-B 0.05
# A-C 0.01

suppressWarnings(rm(con))
con<-file("Results/Three_age/Beta_dispersion_general_and_pairwise_between_Three_ages.txt")
sink(con, append = TRUE)
cat("General beta dispers on Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersione on Bray Curtis (permuted p-values) \n", fill=TRUE)
disp_ASV$pairwise$permuted
sink()
close(con)


################## PCoA BETA DIVERSITY THREE_AGES ##############

data.prop.label<-data.prop # ASV
data.prop.label<-transform_sample_counts(data.prop.label, sqrt)
sample_data(data.prop.label)$Three_ages<-gsub("A","20-35 (A)",sample_data(data.prop.label)$Three_ages)
sample_data(data.prop.label)$Three_ages<-gsub("B","36-52 (B)",sample_data(data.prop.label)$Three_ages)
sample_data(data.prop.label)$Three_ages<-gsub("C","53-68 (C)",sample_data(data.prop.label)$Three_ages)
DistBC <- phyloseq::distance(data.prop.label, method = "bray")
ordBC <- ordinate(data.prop.label, method = "PCoA", distance = DistBC)
plot_ordination(data.prop.label, ordBC, color = "Three_ages") +
  geom_point(size=3) + theme_classic() + stat_ellipse() + guides(color="none") +
  scale_color_manual(values=c("20-35 (A)"="#E69F00","36-52 (B)"="#009E73","53-68 (C)"="#0072B2")) +
  geom_text(aes(label=Sample), size=2.2, col="white", show.legend = FALSE) +
  labs(title="Beta dispersion Pr(>F): 0.008",  # already computed before the paper revision
       subtitle= "Beta dispersion group A vs B Pr(>F) = 0.005\nBeta dispersion group A vs C Pr(>F) = 0.01", 
       color="Age",
       x=paste("PC1:", round(ordBC$values$Relative_eig[1],3)*100,"% variance"),
       y=paste("PC2:", round(ordBC$values$Relative_eig[2],3)*100,"% variance") ) +
  theme_classic()
ggsave(file="Results/Three_age/PCoA_in_2D.png", width = 6, height = 4.5, dpi=300)


# IN 3D
{Dist<- vegdist(t(otu_table(data.prop.label)), method = "bray")                                                                                                      
  obj<-ecodist::pco(Dist)
  matrix<-obj[["vectors"]]
  row.names(matrix)<-sample_names(data.prop.label)
}
color_table<-as.data.frame(cbind(as.character(sample_data(data.prop.label)$Three_age),as.character(sample_data(data.prop.label)$Three_age)))
colnames(color_table)<-c("Gruppo","Colore")
unique(color_table$Gruppo)
colors <-color_table$Colore
colors <- gsub("36-52 (B)","#009E73",color_table$Colore, fixed=T) #green
colors <- gsub("20-35 (A)","#E69F00",colors, fixed=T) # red
colors <- gsub("53-68 (C)","#0072B2",colors, fixed=T) # light blue
rgl::open3d(windowRect=c(25,25,1000,1200))
pca3d::pca3d(matrix, col=colors, radius=1.1, show.shadows = T, 
             axes.color = "darkgray", show.plane = F, show.centroids = F)
rgl::rgl.viewpoint(theta = 42.8, phi = 30.8, fov = 120, zoom = 0.34)
rgl::texts3d(matrix[,1:3], texts=row.names(matrix), col="black", adj=c(1.5, 1), cex=1.4)
rgl::rgl.snapshot("temp.png")
temp<-png::readPNG("temp.png")
# install.packages("pdftools") # needs also sudo apt install libpoppler-cpp-dev
png::writePNG(temp, "Results/Three_age/PCoA_in_3D.png", dpi = 300)
unlink("temp.png")
rgl::close3d()


######################## BETA DIVERSITY DRY SKIN (PERM only) #######################

suppressWarnings(rm(ASV.prop))

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")

sample_OTU<-as.data.frame(t(ASV.prop))
perm_ASV<- vegan::adonis(sample_OTU ~Dry_Skin, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.genus.prop))
perm_g<- vegan::adonis(sample_OTU ~Dry_Skin, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.fam.prop))
perm_f<- vegan::adonis(sample_OTU ~Dry_Skin, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.class.prop))
perm_o<- vegan::adonis(sample_OTU ~Dry_Skin, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.order.prop))
perm_c<- vegan::adonis(sample_OTU ~Dry_Skin, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(ASV.phy.prop))
perm_p<- vegan::adonis(sample_OTU ~Dry_Skin, data=metadata, permutations = 9999, method="bray")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

write.csv2(beta, file="Results/Dry_Skin/Beta_div_PERM_results.csv")


######################## HIERARCHICAL CLUSTERING ###################

suppressWarnings(rm(c))

# On Three ages
{c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
  c<-as.dendrogram(c)
  color_table<-as.data.frame(cbind(as.character(sample_data(data)$Three_ages),as.character(sample_data(data)$Three_ages)))
  colnames(color_table)<-c("Gruppo","Colore")
  colors <- gsub("B","#009E73",color_table$Colore) #green
  colors <- gsub("A","#E69F00",colors) #orange
  colors <- gsub("C","#0072B2",colors) #blue
  labels_colors(c) <- colors[order.dendrogram(c)] # but sort them based on their order in dendrogram
}
png(file="Results/Three_age/Hierarchical_cluster_Hellinger_Three_ages.png",width=2500,height=1800, res=300)
par(mar=c(8,4,4,3), cex.lab=0.4, cex.main=1.4, cex.sub=1.3)
plot(c)
dev.off()

suppressWarnings(rm(c,color_table,labels_colors))

# # On Dry skin      NOTHING THERE
# {c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
# c<-as.dendrogram(c)
# color_table<-as.data.frame(cbind(as.character(sample_data(data)$Dry_Skin),as.character(sample_data(data)$Dry_Skin)))
# colnames(color_table)<-c("Gruppo","Colore")
# colors <- gsub("Not dry","#009E73",color_table$Colore) #green
# colors <- gsub("Dry","#E69F00",colors) #orange
# labels_colors(c) <- colors[order.dendrogram(c)] # but sort them based on their order in dendrogram
# }
# png(file="Results/Dry_Skin/Hierarchical_cluster_Hellinger.png",width=2500,height=1800, res=300)
# par(mar=c(8,4,4,3), cex.lab=0.4, cex.main=1.4, cex.sub=1.3)
# plot(c)
# dev.off()

suppressWarnings(rm(c,color_table,labels_colors))


##################### DA WITH DESEQ2 (THREE AGES) #######################

# Trimming under sum of 10, see DeSeq2 tutorial
data_pruned<- prune_taxa(taxa_sums(data) > 10, data)

# preparing new taxa vocabularies for plots (some ASV may change glomming after trimming)
{data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
{taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
  tax_table(data_pruned.genus.prop)<-as.matrix(taxa_temp)
  # adding informations to uncultured families
  taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
  for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
    taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  tax_table(data_pruned.fam.prop)<-as.matrix(taxa_temp)
  rm(taxa_temp)
}

# function to automatize the comparison between groups
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Three_ages)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Three_ages)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Three_ages", nume[i], deno[i]), test="Wald")
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("A","B","A")
deno<-c("B","C","C")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
unlink("DE_results.txt")
write.csv2(DE_results, file="Results/Three_age//DE_every_results.csv", quote=F, row.names = F)
DE_res_age <- DE_results # for later

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
  target<-prune_taxa(target, data_pruned.genus.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Family","ASV"]
  target<-prune_taxa(target, data_pruned.fam.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Order","ASV"]
  target<-prune_taxa(target, data_pruned.order.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Class","ASV"]
  target<-prune_taxa(target, data_pruned.class.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
  target<-prune_taxa(target, data_pruned.phy.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
  rm(tabella)
}

# unique boxplot of DeSeq2 results
{tabella_g$Taxa<-"Genera"
  tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
  colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"
}
{tabella_f$Taxa<-"Families"
  tabella_f[,c("Phylum","Order","Class")]<-NULL
  colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"
}
{tabella_o$Taxa<-"Orders"
  tabella_o[,c("Phylum","Class")]<-NULL
  colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"
}
{tabella_c$Taxa<-"Classes"
  tabella_c[,"Phylum"]<-NULL
  colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"
}
{tabella_p$Taxa<-"Phyla"
  colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"
}

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_c,tabella_o,tabella_p)

# removing redundants results or almost absent bacteria
tabella_tot2<- tabella_tot[! tabella_tot$Bacteria %in% c("Exiguobacterium","Exiguobacteraceae",
                                                         "Spirochaetaceae",
                                                         "Spirochaetia","Campylobacterales",
                                                         "Spirochaetales") ,]
tabella_tot2<- tabella_tot2[! (tabella_tot2$Bacteria == "Absconditabacteriales_(SR1)" & tabella_tot2$Taxa!="Orders"), ]

ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Three_ages)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("A"="#E69F00","B"="#009E73","C"="#0072B2")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 9 ) + 
  theme(strip.text.x=element_text(size=9,colour="black")) + 
  theme(axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=6.9), 
        axis.text.y = element_text(size=6.2),
        legend.margin=margin(-20, 0, 0, 0), legend.position="bottom") +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(y="Proportional Abundance", fill="Three_ages", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.2, 0.5 ,1, seq(2,20,2),seq(22,max(tabella_tot2$Abundance)+4,3)))
ggsave(filename = "Results/Three_age/DA_DeSeq2_Three_ages.png", width = 6.5, height = 4.5, dpi=300)
dev.off()


################### DA WITH DESEQ2 (DRY SKIN) #################

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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Dry_Skin)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Dry_Skin", "Not dry", "Dry"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res
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
    # write.csv2(r, file=paste0("Results/DA_",t,"_ratio_Not dry_vs_Dry.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Dry_Skin/DA_DESeq2.csv", row.names = F)
DE_res_Dry <- Res_tot # for later


{ Table_tot2<-Table_tot[! is.na(Table_tot$Bacteria), ]
  Table_tot2$Dry_Skin<-gsub("Not dry","normal",Table_tot2$Dry_Skin)
  Table_tot2$Dry_Skin<-gsub("Dry","Dry Skin",Table_tot2$Dry_Skin)
  Table_tot2$Taxa<-gsub("Phylum","Phyla",Table_tot2$Taxa)
  Table_tot2$Taxa<-gsub("Class","Classes",Table_tot2$Taxa)
  Table_tot2$Taxa<-gsub("Genus","Genera",Table_tot2$Taxa)
}

ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Dry_Skin)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Genera")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 9) +
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  scale_fill_manual(values=c("Dry Skin"="#0072B2","normal"="#E69F00")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=7.5), 
        axis.text.y = element_text(size=6.8), legend.title = element_text(size=10),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.2, 0.5, 1, 2, 4, seq(6,max(Table_tot$Abundance)+2,3))) +
  labs(y="Proportional Abundance", fill="Skin", x="")
ggsave(filename = "Results/Dry_Skin/Dry_Skin.png", width = 7, height = 4.8, dpi=300)
dev.off()


################### PLOTTING LEFSE RESULTS (THREE AGES) ################

a <- read.delim("picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("function.","description")]
KO_abundances<-a[ , colnames(a)!="description" ] # for later

# the descriptions include the genes names, removing them
{descr<-unlist(strsplit(Descriptions$description, split = ";"))
descr<-descr[grepl("^ ", descr)] # only the strings that start with a space are the targets
Descriptions$description<-descr
}

{
Significative_functions_LEFSE<- read.table(file="Results/PICRUST_LEFSE/LEFSE_GALAXY_risultato_Three_ages.lefse_internal_res", sep="\t")
head(Significative_functions_LEFSE, n=4)
row.names(Descriptions)<-Descriptions$function.
Descriptions<-Descriptions[Significative_functions_LEFSE$V1,] # they must have the same order
colnames(Significative_functions_LEFSE)<-c("Function","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Function<-Descriptions$description
Significative_functions_LEFSE$Function_ID<-Descriptions$function.
head(Significative_functions_LEFSE, n=4)
}
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.csv2(Significative_functions_LEFSE, file="Results/PICRUST_LEFSE/Significantly_different_KO_Functions_THREE_AGES.csv", quote = F, row.names = F)

# modifing names for the plot
Significative_functions_LEFSE$Function<-gsub("&beta;-","", Significative_functions_LEFSE$Function, fixed = T)
{Significative_functions_LEFSE$Function<-paste("",Significative_functions_LEFSE$Function,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
}
{Significative_functions_LEFSE$Function<-gsub("dehydrogenase","dehydrogen.",Significative_functions_LEFSE$Function, fixed=T)
Significative_functions_LEFSE$Function<-gsub("[","sep",Significative_functions_LEFSE$Function, fixed=T)
Significative_functions_LEFSE$Function<-gsub("]","sep",Significative_functions_LEFSE$Function, fixed=T)
Significative_functions_LEFSE$Function<-gsub("sep.*sep","",Significative_functions_LEFSE$Function) # NB: " .* " taglia il pattern in mezzo
Significative_functions_LEFSE$Function<-factor(Significative_functions_LEFSE$Function, levels = Significative_functions_LEFSE$Function) # to prevent alphabetical re-sorting
}
# renaming factors too
Significative_functions_LEFSE$Class_with_highest_mean <- gsub("B","36-52 (B)",Significative_functions_LEFSE$Class_with_highest_mean)
Significative_functions_LEFSE$Class_with_highest_mean <- gsub("A","20-35 (A)",Significative_functions_LEFSE$Class_with_highest_mean)

# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="20-35 (A)"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="20-35 (A)"]

# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Function, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) +
  labs(x="  log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Function), 
            hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="36-52 (B)",1,0), size=3.2) +
  scale_fill_manual(values=c("20-35 (A)"="#E69F00", "36-52 (B)"="#009E73")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-2,-30,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST_LEFSE/PICRUST2_LEFSE_diff_KO_THREE_AGES.png", width = 9.7, height = 4, dpi = 300)


#### bar plots
KO_abundances[ ,colnames(KO_abundances)!="function."] <- apply(KO_abundances[ , colnames(KO_abundances)!="function."] ,2, function(x) (x/sum(x))*100) # relative abundance
row.names(KO_abundances)<-KO_abundances$function.
KO_abundances<-KO_abundances[Significative_functions_LEFSE$Function_ID, ] # selecting only the significative observations

KO_abundances<-melt (KO_abundances, id="function.")
meta<-as(sample_data(data),"data.frame")
row.names(meta)<-meta$ID
groups<-meta[KO_abundances$variable , c("Sample","Three_ages") ] # "variable" is the sample ID
row.names(Significative_functions_LEFSE)<-Significative_functions_LEFSE$Function_ID
descr<-Significative_functions_LEFSE[KO_abundances$function. , "Function"]
KO_abundances<-cbind.data.frame(KO_abundances,groups,descr)
head(KO_abundances, n=3)

# shortening the descriptions
unique(KO_abundances$descr)
{KO_abundances$descr<-gsub("peptide methionine sulfoxide reductase msrA/msrB", "peptide methionine\nsulfoxide reductase",KO_abundances$descr)
KO_abundances$descr<-gsub("succinate dehydrogen. / fumarate reductase, iron-sulfur subunit", "succ dehydr fumarate reduct\niron-sulfur subunit",KO_abundances$descr)
KO_abundances$descr<-gsub("succinate dehydrogen. / fumarate reductase, cytochrome b subunit ", "succ dehydr fumarate reduct\ncytochrome b subunit ",KO_abundances$descr)
KO_abundances$descr<-gsub("penicillin-binding protein 1A", "penicillin-binding\nprotein 1A",KO_abundances$descr)
KO_abundances$descr<-gsub("molybdopterin molybdotransferase", "molyb\nmolybdotransferase",KO_abundances$descr)
}

average_value <-tapply(KO_abundances$value, paste0(KO_abundances$Three_ages,KO_abundances$function.), mean)
KO_abundances <- cbind.data.frame(KO_abundances, Average=average_value[paste0(KO_abundances$Three_ages,KO_abundances$function.)])

### to many functions to plot, splitting the plot in two
first_half<-1:length(unique(KO_abundances$descr))/2 - 0.5 # 9 --> 4.5 - 0.5 --> the first 4
second_part<- (max(first_half)+1):length(unique(KO_abundances$descr))

subset_KO<-KO_abundances[KO_abundances$descr %in% unique(KO_abundances$descr)[first_half] , ] 
plot1<-ggplot( data=subset_KO, aes(x=Sample , y=value, fill=Three_ages)) +
  geom_bar(stat="identity") +
  scale_y_sqrt(breaks= c(0.01,seq(0.05, 0.5, 0.05)) ) +
  theme_classic() +
  theme(strip.text = element_text(size=7.6),
        axis.text.x = element_text(size=7, angle = 90, vjust=0.5, hjust = 1),
        axis.text.y = element_text(size=6.5),
        axis.title.y = element_text(size=7)) +
  guides(fill="none") +
  labs(x="", y="Predicted KO relative abundance (%)") +
  scale_fill_manual(values=c("A"="#E69F00","B"="#009E73","C"="#0072B2")) +
  geom_errorbar(aes(ymin=Average, ymax=Average), size=0.5, color="gray")+
  facet_grid2( .~descr + Three_ages, scales = "free_x", 
            space="free", strip = strip_nested(size="constant"))

subset_KO<-KO_abundances[KO_abundances$descr %in% unique(KO_abundances$descr)[second_part] , ] 
plot2<-ggplot( data=subset_KO, aes(x=Sample , y=value, fill=Three_ages)) +
  geom_bar(stat="identity") +
  scale_y_sqrt(breaks= c(0.01,seq(0.05, 0.5, 0.05)) ) +
  facet_grid2( .~descr + Three_ages, scales = "free_x", 
               space="free", strip = strip_nested(size="constant") ) +
  theme_classic() +
  theme(strip.text = element_text(size=7.5),
        axis.text.x = element_text(size=6, angle = 90, vjust=0.5, hjust = 1),
        axis.text.y = element_text(size=6.5),
        axis.title.y = element_text(size=7)) +
  guides(fill="none") +
  labs(x="", y="Predicted KO relative abundance (%)", 
       caption="The gray lines represent the mean of each group") +
  scale_fill_manual(values=c("A"="#E69F00","B"="#009E73","C"="#0072B2")) +
  geom_errorbar(aes(ymin=Average, ymax=Average), size=0.5, color="gray")

png(filename = "Results/PICRUST_LEFSE/Three_ages_barplots.png", res=300, width=2050,height=1450)
ggarrange(plot1, plot2, nrow = 2)
dev.off( )


rm(Significative_functions_LEFSE, a, KO_abundances, plot1, plot2)


################### PLOTTING LEFSE RESULTS (Dry Skin) ################

a <- read.delim("picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("function.","description")]
KO_abundances<-a[ , colnames(a)!="description" ] # for later

# the descriptions include the genes names, removing them
descr<-unlist(strsplit(Descriptions$description, split = ";"))
descr<-descr[grepl("^ ", descr)] # only the strings that start with a space are the targets
# there are identical descriptions, adding to them also the KO code
descr[grepl("putative transp",descr)] <- Descriptions$description[grepl("putative transp",Descriptions$description)]
Descriptions$description<-descr


Significative_functions_LEFSE<- read.table(file="Results/PICRUST_LEFSE/LEFSE_GALAXY_risultato_DRY_SKIN.lefse_internal_res", sep="\t")
{head(Significative_functions_LEFSE, n=4)
row.names(Descriptions)<-Descriptions$function.
Descriptions<-Descriptions[Significative_functions_LEFSE$V1,] # they must have the same order
head(Descriptions, n=4) # a further check
colnames(Significative_functions_LEFSE)<-c("Function","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Function<-Descriptions$description
Significative_functions_LEFSE$Function_ID<-Descriptions$function.
}
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.csv2(Significative_functions_LEFSE, file="Results/PICRUST_LEFSE/Significantly_different_KO_Functions_Dry_Skin.csv", quote = F, row.names = F)

# modifying names for the plot
Significative_functions_LEFSE$Function<-gsub("&beta;-","", Significative_functions_LEFSE$Function, fixed = T)
{Significative_functions_LEFSE$Function<-paste("",Significative_functions_LEFSE$Function,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
}
Significative_functions_LEFSE$Function<-gsub("[","sep",Significative_functions_LEFSE$Function, fixed=T)
Significative_functions_LEFSE$Function<-gsub("]","sep",Significative_functions_LEFSE$Function, fixed=T)
Significative_functions_LEFSE$Function<-gsub("sep.*sep","",Significative_functions_LEFSE$Function) # NB: " .* " taglia il pattern in mezzo
Significative_functions_LEFSE$Function<- factor(Significative_functions_LEFSE$Function, levels = Significative_functions_LEFSE$Function) # to prevent alphabetical re-sorting

# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="si"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="si"]

# new names to the groups
Significative_functions_LEFSE$Class_with_highest_mean<-gsub("no","Normal",Significative_functions_LEFSE$Class_with_highest_mean)
Significative_functions_LEFSE$Class_with_highest_mean<-gsub("si","Dry skin",Significative_functions_LEFSE$Class_with_highest_mean)

# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Function, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) +
  labs(x="  log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Function), 
            hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Normal",1,0), size=3.2) +
  scale_fill_manual(values=c( "Dry skin"="#0072B2", "Normal"="#E69F00" )) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-2,-30,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST_LEFSE/PICRUST2_LEFSE_diff_KO_Dry_Skin.png", width = 9, height = 4, dpi = 300)


#### bar plots
KO_abundances[ ,colnames(KO_abundances)!="function."] <- apply(KO_abundances[ , colnames(KO_abundances)!="function."] ,2, function(x) (x/sum(x))*100) # relative abundance
row.names(KO_abundances)<-KO_abundances$function.
KO_abundances<-KO_abundances[Significative_functions_LEFSE$Function_ID, ] # selecting only the significative observations

KO_abundances<-melt (KO_abundances, id="function.")
meta<-as(sample_data(data),"data.frame")
row.names(meta)<-meta$ID
groups<-meta[KO_abundances$variable , c("Sample","Dry_Skin") ] # "variable" is the sample ID
row.names(Significative_functions_LEFSE)<-Significative_functions_LEFSE$Function_ID
descr<-Significative_functions_LEFSE[KO_abundances$function. , "Function"]
KO_abundances<-cbind.data.frame(KO_abundances, groups[,c("Sample","Dry_Skin")] ,descr)
head(KO_abundances, n=3)

# shortening the descriptions
unique(KO_abundances$descr)
{KO_abundances$descr<-gsub("branched-chain amino acid:cation transporter, LIVCS family ", "cation transporter\nLIVCS family ",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("poly(glycerol-phosphate) alpha-glucosyltransferase", "poly(glycerol-phosphate)\nalpha-glucosyltransf",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("chromosome partitioning protein", "chromosome partitioning \n protein",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("malate dehydrogenase (quinone)", "malate dehydr (quinone) ",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("methyl-accepting chemotaxis protein", "methyl-accepting \nchemotaxis protein",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("penicillin-binding protein 1A", "penicillin-binding\nprotein 1A",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("RNA polymerase sigma-70 factor, ECF subfamily ", "RNA pol sigma-70\nECF subfamily ",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("K07498; putative transposase", "(K07498)\nputative transposase",KO_abundances$descr, fixed = T)
  KO_abundances$descr<-gsub("K07491; putative transposase", "(K07491)\nputative transposase",KO_abundances$descr, fixed = T)
}
KO_abundances$Dry_Skin<-gsub("Not dry","Normal",KO_abundances$Dry_Skin)
KO_abundances$Dry_Skin<-gsub("Dry","Dry Skin",KO_abundances$Dry_Skin)

# average_value <-tapply(KO_abundances$value, KO_abundances$Dry_Skin, mean)
# KO_abundances <- cbind.data.frame(KO_abundances, Average=average_value[KO_abundances$Dry_Skin]) # tapply returns a factor which is ordinable by names

average_value <-tapply(KO_abundances$value, paste0(KO_abundances$Dry_Skin,KO_abundances$function.), mean)
KO_abundances <- cbind.data.frame(KO_abundances, Average=average_value[paste0(KO_abundances$Dry_Skin,KO_abundances$function.)])


### to many functions to plot, splitting the plot in two
first_half<-1:length(unique(KO_abundances$descr))/2 # 10 --> the first 5
second_part<- (max(first_half)+1):length(unique(KO_abundances$descr))

subset_KO<-KO_abundances[KO_abundances$descr %in% unique(KO_abundances$descr)[first_half] , ] 
plot1<-ggplot( data=subset_KO, aes(x=Sample , y=value, fill=Dry_Skin)) +
  geom_bar(stat="identity") +
  scale_y_sqrt(breaks= c(0.01,seq(0.05, 0.5, 0.05)) ) +
  theme_classic() +
  theme(strip.text = element_text(size=7),
        axis.text.x = element_text(size=6, angle = 90, vjust=0.5, hjust = 1),
        axis.text.y = element_text(size=5.5),
        axis.title.y = element_text(size=7)) +
  guides(fill="none") +
  labs(x="", y="Predicted KO relative abundance (%)") +
  scale_fill_manual(values=c("Dry Skin"="#0072B2","Normal"="#E69F00")) +
  geom_errorbar(aes(ymin=Average, ymax=Average), size=0.5, color="gray")+ 
  facet_grid2( .~descr + Dry_Skin, scales = "free",
               space="free", strip = strip_nested(size="constant"))

subset_KO<-KO_abundances[KO_abundances$descr %in% unique(KO_abundances$descr)[second_part] , ] 
plot2<-ggplot( data=subset_KO, aes(x=Sample , y=value, fill=Dry_Skin)) +
  geom_bar(stat="identity") + 
  scale_y_sqrt(breaks= c(0.01,seq(0.05, 0.5, 0.05)) ) +
  facet_grid2( .~descr + Dry_Skin, scales = "free", 
               space="free", strip = strip_nested(size="constant") ) +
  theme_classic() +
  theme(strip.text = element_text(size=7.2),
        axis.text.x = element_text(size=6, angle = 90, vjust=0.5, hjust = 1),
        axis.text.y = element_text(size=5.5),
        axis.title.y = element_text(size=7) ) +
  guides(fill="none") +
  labs(x="", y="Predicted KO relative abundance (%)",
       caption = "The gray lines represent the mean of each group") +
  scale_fill_manual(values=c("Dry Skin"="#0072B2","Normal"="#E69F00")) +
  geom_errorbar(aes(ymin=Average, ymax=Average), size=0.5, color="gray")

png(filename = "Results/PICRUST_LEFSE/Dry_Skin_barplots.png", res=300, width=2050,height=1450)
ggarrange(plot1, plot2, nrow = 2)
dev.off( )



rm(Significative_functions_LEFSE, a, plot1, plot2, KO_abundances)


################## CORRELATION WITH SNPs : IMPORTING THE DATASETS ####################

# importing every SNP data sets
mydir = "Samples_SNPs"    # path to SNP dataframes 
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)   # sceglie i file
myfiles # path of each one
lista_dataframe= lapply(myfiles, read.csv, sep = ";", header = TRUE) # importing each one of them in a single list --> lapply
head(lista_dataframe[[1]], n=3) # just to check the first one
myfiles<-gsub("Samples_SNPs/","", myfiles)
names(lista_dataframe)<-myfiles     # list numbers overwriten by csv names (same order)

# importing the "vocabulary" of those SNP potentially associated with missense mutation (downloaded from Illumina Infinium Support files, updated through Biomart R package)
Vocabulary <- read.csv2("SNP_missense_Vocabulary.csv")

lista_dataframe<-lapply(lista_dataframe, function(x) left_join(x, Vocabulary)) # using dplyr
head(lista_dataframe[[1]], n=3) # to check
# removing those without informations (being absent in the vocabulary OR or without known mutations)
lista_dataframe <- lapply(lista_dataframe, function(x) na.omit(x)) # removing the NAs
length(which(lista_dataframe[[1]]$Minor.allele!="")) # 32692 minor alleles on 42691 tot
lista_dataframe<-lapply(lista_dataframe, function(x) subset(x, !(Minor.allele %in% "")))  # removing SNPs without known ancestral allele 
head(lista[[1]], n=3)
identical(length(lista[[1]]$Minor.allele),length(which(lista_dataframe[[1]]$Minor.allele!=""))) # TRUE
lista_dataframe<-lista # it was a temp object to check the length afterwards
rm(lista)

# overwriting col names to have the same names among data frames
nomi<-colnames(lista_dataframe[[1]])
nomi[2:3]<-c("AlleleA","AlleleB")
lista_dataframe<-lapply(lista_dataframe, setNames, nm = nomi) # It must be used "setNames" because "colnames" does not work here
identical(colnames(lista_dataframe[[1]]),colnames(lista_dataframe[[2]])) # TRUE

#  A  or  B  Allele == minor allele
lista_minor<-lapply(lista_dataframe, function(x) filter(x, AlleleA==Minor.allele | AlleleB==Minor.allele))
lista_omo<-lapply(lista_minor, function(x) filter(x, AlleleA==AlleleB)) # homozygosis              

rm(lista_dataframe, Vocabulary, mydir,myfiles,nomi) 

# Saving all the dataframes
dir.create("Results/SNP/SNP_with_info")
dir.create("Results/SNP/SNP_with_info/Homozygosis")
nms<- names(lista_minor)
lapply(1:length(lista_minor), function(x) write.csv2(lista_minor[[x]], file = paste0("Results/SNP/SNP_with_info/",nms[x], ".csv"), row.names = FALSE, na=""))
nms<- names(lista_omo)
lapply(1:length(lista_omo), function(x) write.csv2(lista_omo[[x]], file = paste0("Results/SNP/SNP_with_info/Homozygosis/",nms[x], ".csv"), row.names = FALSE, na=""))
rm(nms)


########### CORRELATION WITH SNPs : PREPARING AND COUNTING THE TARGETS ############

head(lista_minor[[1]], n=3)

# disgenet.org --> data set Skin Wrinkling ( CUI: C0037301 )
# COG5	TALDO1	ALDH18A1	PIK3R1	LMNA	KCNJ6	EFEMP2	SLC25A24	ELN	POLR3A	FBLN5	LTBP4

Targets<-c("COG5","TALDO1","ALDH18A1","PIK3R1","LMNA","KCNJ6","EFEMP2","SLC25A24","ELN","POLR3A","FBLN5","LTBP4")

# counting the SNPs
for(i in 1:length(Targets)) {
  etero<-sapply(lista_minor, function(x) length(which(str_starts(x$Gene.associato, regex(pattern = Targets[i], ignore_case = T)))))
  omo<-sapply(lista_omo, function(x) length(which(str_starts(x$Gene.associato, regex(pattern = Targets[i], ignore_case = T)))))
  assign(paste(Targets[i]), etero+omo)  # if Homozygote, both etero and omo are 1 --> 2
}
rm(etero,omo)

# each one of those objects came from the "assign" of the loop but they maintained also the dataframe name (-->sample) from which they derive, as a named vector
Wrinkle<-rbind(COG5,TALDO1,ALDH18A1,PIK3R1,LMNA,KCNJ6,EFEMP2,SLC25A24,ELN,POLR3A,FBLN5,LTBP4)
rm(COG5,TALDO1,ALDH18A1,PIK3R1,LMNA,KCNJ6,EFEMP2,SLC25A24,ELN,POLR3A,FBLN5,LTBP4,Targets,i)

# same name between Wrinkle dataframe and phyloseq object (manually checked the the proper order of those names)
colnames(Wrinkle)<-c("A1","A2","A3","A5","A4","B1","B2","B4","B3","B5","C1","C2","C5","C4","C3")
# old names # colnames(Wrinkle)<-c("A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3","E1","E2","E3")
write.csv2(Wrinkle,file="Results/SNP/Wrinkle_SNP_from_Disgenet.csv", row.names = T)

sample_quest<-as.data.frame(sample_data(data))
Wrinkle<-as.data.frame(Wrinkle[ ,sample_quest$Sample]) # same order of the phyloseq object

Wrinkle$sum_0 <- rowSums(Wrinkle) # to check those with no SNPs in any sample
Wrinkle<-Wrinkle[! Wrinkle$sum_0==0,] # removed those with 0 SNPs counted
Wrinkle$sum_0<-NULL


############ NOW SUBSETTING THE BACTERIA

# Phyla
suppressWarnings(rm(taxa))
taxa<- c( DE_res_age[DE_res_age$Rank=="Phylum", "Phylum"] , DE_res_Dry[DE_res_Dry$Taxon=="Phylum", "Phylum"] )
taxa <- taxa[! is.na(taxa)]
data_filter<-subset_taxa(data.phy.prop, Phylum %in% taxa )
target<-psmelt(data_filter)
target<-target[ , colnames(target) %in% c("Sample","Abundance","Phylum")]
for (i in unique(taxa) ) { 
  P<- target[target$Phylum == i & ! is.na(target$Phylum), ]
  # P<-filter(target, Phylum == (na.omit(DE_res_age$Phylum))[i]) # abundances of one target only
  Pt<-t(P[,1:2])
  colnames(Pt)<-Pt[1,] # sample names as colnames
  Pt<-Pt[2,colnames(Wrinkle)] # only abundances and with the same order of Wrinckle
  assign( i , as.numeric(Pt))  
}
rm(Pt,i)
Phylum<-rbind.data.frame(Campilobacterota,Firmicutes) # Spirochaetota is a redundant results (see DESeq2)
colnames(Phylum)<-colnames(Wrinkle) # they lose the original sample names with rbind 
row.names(Phylum)<-c("Campilobacterota","Firmicutes")
rm(Campilobacterota,Firmicutes,Spirochaetota)


# Classes
suppressWarnings(rm(taxa))
taxa<- c( DE_res_age[DE_res_age$Rank=="Class", "Class"] , DE_res_Dry[DE_res_Dry$Taxon=="Class", "Class"] )
taxa <- taxa[! is.na(taxa)]
data_filter<-subset_taxa(data.class.prop, Class %in% taxa )
target<-psmelt(data_filter)
target<-target[ , colnames(target) %in% c("Sample","Abundance","Class")]
for (i in unique(taxa) ) { 
  P<- target[target$Class == i & ! is.na(target$Class), ]
  Pt<-t(P[,1:2])
  colnames(Pt)<-Pt[1,]
  Pt<-Pt[2,colnames(Wrinkle)]
  assign( i , as.numeric(Pt))  
}
rm(Pt,i)
Class<-rbind.data.frame(Actinobacteria,Alphaproteobacteria,Bacilli,Campylobacteria,Clostridia,Negativicutes) # only the NOT redundants
colnames(Class)<-colnames(Wrinkle)
row.names(Class)<-c("Actinobacteria","Alphaproteobacteria","Bacilli","Campylobacteria","Clostridia","Negativicutes")
rm(Actinobacteria,Alphaproteobacteria,Bacilli,Campylobacteria,Clostridia,Negativicutes,Spirochaetia)


# Orders
suppressWarnings(rm(taxa))
taxa<- c( DE_res_age[DE_res_age$Rank=="Order", "Order"] , DE_res_Dry[DE_res_Dry$Taxon=="Order", "Order"] )
taxa <- taxa[! is.na(taxa)]
data_filter<-subset_taxa(data.order.prop, Order %in% taxa )
target<-psmelt(data_filter)
target<-target[ , colnames(target) %in% c("Sample","Abundance","Order")]
for (i in unique(taxa) ) { 
  P<- target[target$Order == i & ! is.na(target$Order), ]
  # P<-filter(target, Phylum == (na.omit(DE_res_age$Phylum))[i]) # abundances of one target only
  Pt<-t(P[,1:2])
  colnames(Pt)<-Pt[1,] # sample names as colnames
  Pt<-Pt[2,colnames(Wrinkle)] # only abundances and with the same order of Wrinckle
  assign( i , as.numeric(Pt))  
}
rm(Pt,i)
Order<-rbind.data.frame(`Absconditabacteriales_(SR1)`,Corynebacteriales,Cytophagales,Exiguobacterales)
colnames(Order)<-colnames(Wrinkle)
row.names(Order)<-c("Absconditabacteriales_(order)","Corynebacteriales","Cytophagales","Exiguobacterales")
rm(`Absconditabacteriales_(SR1)`,Campylobacterales,Corynebacteriales,Cytophagales,Exiguobacterales,Spirochaetales)


# Families
suppressWarnings(rm(taxa))
taxa<- c( DE_res_age[DE_res_age$Rank=="Family", "Family"] , DE_res_Dry[DE_res_Dry$Taxon=="Family", "Family"] )
taxa <- taxa[! is.na(taxa)]
data_filter<-subset_taxa(data.fam.prop, Family %in% taxa )
target<-psmelt(data_filter)
target<-target[ , colnames(target) %in% c("Sample","Abundance","Family")]
for (i in unique(taxa) ) { 
  P<- target[target$Family == i & ! is.na(target$Family), ]
  # P<-filter(target, Phylum == (na.omit(DE_res_age$Phylum))[i]) # abundances of one target only
  Pt<-t(P[,1:2])
  colnames(Pt)<-Pt[1,] # sample names as colnames
  Pt<-Pt[2,colnames(Wrinkle)] # only abundances and with the same order of Wrinckle
  assign( i , as.numeric(Pt))  
}
rm(Pt,i)
Family<-rbind.data.frame(Nocardioidaceae,Spirochaetaceae)
colnames(Family)<-colnames(Wrinkle)
row.names(Family)<-c("Nocardioidaceae","Spirochaetaceae")
rm(Exiguobacteraceae,Nocardioidaceae,Spirochaetaceae)


# Genera
suppressWarnings(rm(taxa))
taxa<- c( DE_res_age[DE_res_age$Rank=="Genus", "Genus"] , DE_res_Dry[DE_res_Dry$Taxon=="Genus", "Genus"] )
taxa <- taxa[! is.na(taxa)]
data_filter<-subset_taxa(data.genus.prop, Genus %in% taxa )
target<-psmelt(data_filter)
target<-target[ , colnames(target) %in% c("Sample","Abundance","Genus")]
for (i in unique(taxa) ) { 
  P<- target[target$Genus == i & ! is.na(target$Genus), ]
  # P<-filter(target, Phylum == (na.omit(DE_res_age$Phylum))[i]) # abundances of one target only
  Pt<-t(P[,1:2])
  colnames(Pt)<-Pt[1,] # sample names as colnames
  Pt<-Pt[2,colnames(Wrinkle)] # only abundances and with the same order of Wrinckle
  assign( i , as.numeric(Pt))  
}
rm(Pt,i)
Genus<-rbind.data.frame(Abiotrophia,Lactococcus,Flavobacterium,Negativicoccus,Treponema,Peptoniphilus)
colnames(Genus)<-colnames(Wrinkle)
row.names(Genus)<-c("Abiotrophia","Lactococcus","Flavobacterium","Negativicoccus","Treponema","Peptoniphilus")
rm(Abiotrophia,Lactococcus,Exiguobacterium,Flavobacterium,Negativicoccus,Treponema,Peptoniphilus,`Absconditabacteriales_(SR1)`)


##### Unique dataframe
Bacteria<-rbind.data.frame(Phylum,Class,Order,Family,Genus)
# write.csv2(Bacteria,file="DESeq2_bacteria_prop_abundance.csv", row.names = T, quote = F)

if(identical(colnames(Wrinkle), colnames(Bacteria))) { # TRUE
  tabella<-as.data.frame(t(rbind(Wrinkle,Bacteria)))
}


###################### CORRELATION WITH SNPs : CORRELATIONS ###################

{r<-rcorr(as.matrix(tabella), type="spearman")
corr<-as.data.frame(r$r) # rho
pvalue<-as.data.frame(r$P) # idem, P corrisponde a tabella p-value
corr<-as.data.frame(as.table(as.matrix(corr)))  
pvalue<-as.data.frame(as.table(as.matrix(pvalue)))
identical(corr[,1:2],pvalue[,1:2])  # just to check
data_corr<-cbind(corr,pvalue[,3])
colnames(data_corr)<-c("Bacter","SNP","Correlation","pvalue")
nomi_batteri<-row.names(Bacteria)
nomi_SNP<-row.names(Wrinkle)
data_corr<-subset(data_corr, Bacter %in% nomi_batteri)
data_corr<-subset(data_corr, SNP %in% nomi_SNP)
}
write.csv2(data_corr, file="Results/SNP/Spearman_DESe2_Bacteria_vs_SNP.csv", row.names = F, quote = F)


######## HEATMAP
data_corr$Sign<-data_corr$pvalue
data_corr$Sign[data_corr$pvalue<0.05]<-"*"
data_corr$Sign[data_corr$pvalue<0.01]<-"***"
data_corr$Sign[data_corr$Sign>0.05]<-""

ggplot(data_corr, aes(x = data_corr$Bacter, y = data_corr$SNP, 
                      fill = data_corr$Correlation)) + 
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +
  scale_fill_gradient2(high="blue", mid = "white", low="red") + 
  theme_bw(base_size=6) +
  theme(axis.text.x=element_text(angle = -30, hjust = 0, size = 9.5), 
        axis.text.y = element_text(size = 9.5)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= data_corr$Sign), color= "white", size =5.2) +
  labs(y="",x="", caption = "p-value<0.05 = *      p-value<0.01 =***") + 
  theme(plot.caption = element_text(color = "black", face = "italic", size = 12))  +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=9), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(filename="Results/SNP/Heatmap_Spearm_Bacteria_DA_versus_SNPs.png", height=5.7,width=8.5, dpi = 300)
