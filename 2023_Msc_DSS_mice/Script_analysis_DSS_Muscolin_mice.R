##################### PREPARING THE ENVIRONMENT ##############

{
library("phyloseq")
library("ggplot2")
library("ggh4x")  
library("egg")
library("ggpubr")
library("vegan")
library("DESeq2")
library("qiime2R")
library("dplyr")
library("dendextend")
library("dplyr")
}

####################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample
sample<-substring(sample, first = 15) 
sample<-gsub("-","",sample)
sample<-paste0("M",sample)
sample_names(data)<-sample
Metadata <- as.data.frame(readxl::read_excel("Metadata.xlsx"))
row.names(Metadata)<-paste0("M",as.character(Metadata$FASTQ_ID))
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID)),as.numeric(original_length))
sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))
sample_data(data)$Time<-factor(sample_data(data)$Time, levels = c("T0","T1"))

write.csv2(cbind(Metadata,original_names),file="Metadata_with_FASTQ_names.csv",row.names = F, quote =F, na="")
rm(original_length,original_names,sample)

############ CHECKING UNASSIGNED IN PROKARYOTE KINGDOM ######################

{Unass<-tax_glom(data, taxrank = "Kingdom") # or domain
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
e
options(scipen=100)
write.csv2(e[,colnames(e)!="Kingdom"], file="QIIME/Unassigned_domain_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)

######################### STARTING THE ANALYSIS #######################

dir.create("Results")
setwd("Results")

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

{taxa_temp<-Taxa.fam
for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
for( x in 1: length(which(is.na(taxa_temp$Family))) ) {
  taxa_temp$Family[ which(is.na(taxa_temp$Family))[1] ]<-paste("NA_ o",taxa_temp[which(is.na(taxa_temp$Family))[1],"Order"])}
for( x in 1: length(which(taxa_temp=="NA_ o NA")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="NA_ o NA")[1] ]<-paste("NA_ c",taxa_temp[which(taxa_temp$Family=="NA_ o NA")[1],"Class"])}
for( x in 1: length(which(duplicated(taxa_temp$Family[taxa_temp$Family=="NA_ c NA"]))) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="NA_ c NA")[1] ]<-paste("NA_ c NA",x+1) }
Taxa.fam.update<-taxa_temp
}

rm(taxa_temp)

#################### % ASSIGNED IN SILVA #########################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Total","Assigned","%","Taxa")
}
assigned
write.csv2(assigned,file="../QIIME/Percentuali_taxa_assegn_in_database_tass.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)

######################## GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(data)==1, data)
length(which(taxa_sums(data)==1)) # if zero there are no singletons
{n<-as.data.frame(otu_table(filter))
  N<-as.data.frame(otu_table(data))
  G<-1-(colSums(n)/colSums(N))
}
con<-file("../QIIME/Percentuale_di_singletons_Good's_coverage.txt")
sink(con)
cat("GOOD'S COVERAGE ESTIMATOR \n", fill=TRUE)
cat("1-(n/N) for each sample, where n is number of singletons and N is Total ASV \n \n", fill=TRUE)
G
sink()
close(con)
rm(con, filter)

########################### COUNTS EXPORT ##########################################

dir.create("Abundances")
setwd("Abundances")

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

###################### ABUNDANCES BAR PLOT ##########################

dir.create("Abundances")
setwd("Abundances")

# TOP 5 Phylum con stacked bar plot
rm(top5, others, tabella)
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
tabella$Time<-gsub("T1","Post treatment",tabella$Time)
tabella$Time<-gsub("T0","Pre treatment",tabella$Time)
tabella$Time<-factor(tabella$Time, levels = c("Pre treatment","Post treatment"))
tabella$Condition<-gsub("KO","Msc-/-",tabella$Condition)
ggplot(data=tabella, aes(y=Abundance, x=Sample_ID, fill=Phylum)) +
  theme_classic(base_size =14) + 
  facet_grid2(.~Condition+Time, scales = "free", space="free", 
              strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Sample", y="Relative abundance", 
       title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP5_phyla_abundances.png",width=8,height=6, dpi=300) 
dev.off()

# TOP 5 Generi
rm(top, others, tabella)
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.genus.prop)
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_selected<-Taxa.genus.update[row.names(tax_selected),]
tax_table(prune.dat_top)<-as.matrix(tax_selected)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
tabella$Time<-gsub("T0","Pre treatment",tabella$Time)
tabella$Time<-gsub("T1","Post treatment",tabella$Time)
tabella$Time<-factor(tabella$Time, levels = c("Pre treatment","Post treatment"))
tabella$Condition<-gsub("KO","Msc-/-",tabella$Condition)
ggplot(data=tabella, aes(y=Abundance, x=Sample_ID, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(.~Condition+Time, scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Sample", x="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="TOP5_genera_abundances.png",width=8,height=6, dpi=300) 
dev.off()

setwd("..")

######################## SETTING GROUP COLORS #########################

# same function down there will search the colors from here
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors_cond <- gsub("WT","chartreuse",tabella_colore$Colore) #green
colors_cond <- gsub("KO","coral",colors_cond) #orange

###################### HIERARCHICAL CLUSTERING #####################

dir.create("Hierarchical_clustering")
setwd("Hierarchical_clustering")

data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))

# euclidean (labels are time)
c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
if(identical(c$labels,sample_names(data.prop))){
  c$labels<-sample_data(data.prop)$Time
}
c<-as.dendrogram(c)
labels_colors(c) <- colors_cond[order.dendrogram(c)]
png(file="Hierarchical_cluster_TIME_Euclidean_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure using Euclidean distance on sqrt proportional ASVs",
     sub="WT = green     Msc-/- = orange")
dev.off()

# euclidean (labels are names)
c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
if(identical(c$labels,sample_names(data.prop))){
  c$labels<-sample_data(data.prop)$Sample_Name
}
c<-as.dendrogram(c)
labels_colors(c) <- colors_cond[order.dendrogram(c)]
png(file="Hierarchical_cluster_NAME_Euclidean_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure using Euclidean distance on sqrt proportional ASVs",
     sub="WT = green     Msc-/- = orange")
dev.off()

rm(c)

# Bray Curtis (labels are time)
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
if(identical(c$labels,sample_names(data.prop))){
  c$labels<-sample_data(data.prop)$Time
}
c<-as.dendrogram(c)
labels_colors(c) <- colors_cond[order.dendrogram(c)]
png(file="Hierarchical_cluster_TIME_Bray_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure using Bray-Curtis distance on sqrt proportional ASVs",
     sub="WT = green     T1 = orange")
dev.off()

# Bray Curtis (labels are names)
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
if(identical(c$labels,sample_names(data.prop))){
  c$labels<-sample_data(data.prop)$Sample_Names
}
c<-as.dendrogram(c)
labels_colors(c) <- colors_cond[order.dendrogram(c)]
png(file="Hierarchical_cluster_Names_Bray_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure using Bray-Curtis distance on sqrt proportional ASVs",
     sub="WT = green     Msc-/- = orange")
dev.off()

rm(c)

setwd("..")

########################## ALFA DIVERSITY _ WT ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data, Condition=="WT")
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Time, New_data$Sample_ID),]
pAlpha$data<-New_data
}
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
pAlpha + theme_bw() + 
  geom_line(aes(group = pAlpha$data$Sample_ID),col="grey",size=0.15) +
  labs(x="Time", title="Alpha diversity between pre and post treatment in wild tipe mices") +
  guides(fill=FALSE, color=FALSE) + 
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11),
        title = element_text(size=9)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", 
                     paired = T, label.x= 1.5, size=3.5, label.y.npc = "top", 
                     vjust=-0.5, hjust=0.35)
ggsave(file="Alfa_diversity_WT_pre_post_paired_wilcoxon.png", width = 6,height =5.5, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
Obser_value<-Obser_value[order(Obser_value$Time, Obser_value$Sample_ID),]
factor<-Obser_value$Time
factor
Obser_value$Sample_ID # to re-check if they have the same order
wilcox.test(Obser_value$value~factor, paired=T)

########################## ALFA DIVERSITY _ KO ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

rm(data.temp)
data.temp<-subset_samples(data, Condition!="WT")
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Time, New_data$Sample_ID),]
  pAlpha$data<-New_data
}
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
pAlpha + theme_bw() + 
  geom_line(aes(group = pAlpha$data$Sample_ID),col="grey",size=0.15) +
  labs(x="Time", title="Alpha diversity between pre and post treatment in Msc-/- mices") +
  guides(fill=FALSE, color=FALSE) + 
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11),
        title = element_text(size=10)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", 
                     paired = T, label.x= 1.5, size=3.5, label.y.npc = "top",
                     vjust=-0.5, hjust=0.35)
ggsave(file="Alfa_diversity_KO_pre_post_paired_wilcoxon.png", width = 6,height =5.5, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
Obser_value<-Obser_value[order(Obser_value$Time, Obser_value$Sample_ID),]
factor<-Obser_value$Time
factor
Obser_value$Sample_ID # to re-check if they have the same order
wilcox.test(Obser_value$value~factor, paired=T)

########################## ALFA DIVERSITY _ WT vs KO at T0 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

rm(data.temp)
data.temp<-subset_samples(data, Time!="T1")
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
}
pAlpha$data$Condition<-gsub("KO","Msc-/-",pAlpha$data$Condition)
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Condition", title="Alpha diversity between wild type and Msc-/- (before treatment)") +
  guides(fill=FALSE, color=FALSE) +
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11),
        title = element_text(size=10)) +
  stat_compare_means(aes(group = Condition), label="p.format",
                     method = "wilcox.test", label.x= 0.9, size=3.2, label.y.npc = "top",
                     vjust=-0.5, hjust=-0.35)
ggsave(file="Alfa_diversity_WT_vs_KO_T0_Mann_Withn_Wilcox.png", width = 6,height =5.5, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

########################## ALFA DIVERSITY _ WT vs KO at T1 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

rm(data.temp)
data.temp<-subset_samples(data, Time=="T1")
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
}
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
pAlpha$data$Condition<-gsub("KO","Msc-/-",pAlpha$data$Condition)
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Condition", title="Alpha diversity between wild type and Msc-/- (after treatment)") +
  guides(fill=FALSE, color=FALSE) + theme(
    axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11),
    title = element_text(size=10)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test",
                     label.x= 0.9, size=3.2, label.y.npc = "top", vjust=-0.5, hjust=-0.35)
ggsave(file="Alfa_diversity_WT_vs_KO_T1_Mann_Withn_Wilcox.png", width = 6,height =5.5, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

##################### BETA DIVERSITY BRAY CURTIS #######################

dir.create("Bray_Curtis_Beta_diversity")
setwd("Bray_Curtis_Beta_diversity")

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Cond_Time, data=metadata, permutations = 9999, method="bray")
  perm_ASV
  perm_ASV_Bray<-perm_ASV # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Cond_Time, data=metadata, permutations = 9999, method="bray")
  perm_g
  perm_g_Bray<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Cond_Time, data=metadata, permutations = 9999, method="bray")
  perm_f
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Cond_Time, data=metadata, permutations = 9999, method="bray")
  perm_o
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Cond_Time, data=metadata, permutations = 9999, method="bray")
  perm_c
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Cond_Time, data=metadata, permutations = 9999, method="bray")
  perm_p
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

### pairwise beta diversity
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=metadata$Cond_Time, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV

# exporting beta diversity
rm(con)
con<-file("Beta_diversity_general_and_pairwise_between_Cond_Time_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)
rm(beta, pair_ASV)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
  disper<-vegan::betadisper(BC.dist,metadata$Cond_Time)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
rm(con)
con<-file("Beta_dispersion_General_and_Pairwise_between_Cond_Time_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)
rm(disp_ASV,a)

######################## PCoA BRAY CURTIS

# on ASV
data.prop.labels<-data.prop
sample_data(data.prop.labels)$Time<-gsub("T0","Before",sample_data(data.prop.labels)$Time)
sample_data(data.prop.labels)$Time<-gsub("T1","After",sample_data(data.prop.labels)$Time)
sample_data(data.prop.labels)$Condition<-gsub("KO","Msc-/-",sample_data(data.prop.labels)$Condition)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Time") +
  geom_line(aes(group=Sample_ID), col="grey", size=0.15)+
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", shape="Treament", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_on_ASV.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  geom_point(size=3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", shape="Treament", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_ASV_no_ellipse.png", width = 8, height = 6, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", shape="Treament", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_on_ASV_points.png", width = 8, height = 6, dpi=300)

# again but with jaccard
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "jaccard")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Jaccard distance \n computed on sqrt proportional ASV", color="Condition", shape="Treament", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Jaccard_on_ASV_points.png", width = 8, height = 6, dpi=300)

# again but in 3D
{Dist<- vegdist(t(otu_table(data.sqrt_prop)), method = "bray")                                                                                                      
  obj<-ecodist::pco(Dist)
  matrix<-obj[["vectors"]]
  row.names(matrix)<-sample_names(data.sqrt_prop)
}
# install.packages("pca3d")
rgl::open3d(windowRect=c(25,25,1200,1200),useNULL = F)
pca3d::pca3d(matrix, col=colors_cond, radius=1.5, show.shadows = T, show.plane = F, show.centroids = F)
rgl::rgl.viewpoint(theta = -38.8, phi = 30.8, fov = 150, zoom = 0.215)
rgl::legend3d("topleft", c("Msc-/-_Before treatment", "WT_Before treatment","Msc-/-_After treatment", "WT_After treatment" ),
              col=c(2,4,2,4), pch = c(2,2,1,1))
rgl::rgl.snapshot("temp.png")
temp<-png::readPNG("temp.png")
png::writePNG(temp, "PCoA_Bray_Curtis_3D_on_ASV.png", dpi = 600)
rgl::close3d()

unlink("temp.png")

# again but on genera
{data.prop.labels<-data.genus.prop
  sample_data(data.prop.labels)$Time<-gsub("T0","Before",sample_data(data.prop.labels)$Time)
  sample_data(data.prop.labels)$Time<-gsub("T1","After",sample_data(data.prop.labels)$Time)
  sample_data(data.prop.labels)$Condition<-gsub("KO","Msc-/-",sample_data(data.prop.labels)$Condition)
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Time") +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional genera", color="Condition", shape="Treatment", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_SUI_GENERI.png", width = 8, height = 6, dpi=300)

setwd("..")

################## DAs WITH DESEQ2 _ T0 vs T1 in WT ##############################

rm(data.temp)
data.temp<-subset_samples(data, Condition=="WT")

# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
rm(data_pruned, data.genus_pruned)
data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
{ data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
  Taxa.genus_pruned<-taxa_temp
  rm(taxa_temp)
}
# adding informations to uncultured families
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
{for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  Taxa.fam_pruned<-taxa_temp
  rm(taxa_temp)
}
# other taxonomic vocabularies
{ Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
  Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
  Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

dir.create("Diff_abundance_DESeq2_WT_Treatment")
system(" echo 'In case of pair analysis, the log2foldchange displayes the change from T0 to T1, see https://support.bioconductor.org/p/105981/' > Diff_abundance_DESeq2_WT_Treatment/Nota_bene.txt ")  

# genus ##############################
rm(res, DE, target, res_genus)
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  rm(res)
  write.csv2(res_genus, file="Diff_abundance_DESeq2_WT_Treatment/Analisi diff dei generi_ratio T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)  
  tax_table(data_pruned.genus.prop)<-as.matrix(Taxa.genus_pruned)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  ##tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family ##############################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_fam<-as.data.frame(res)
  res_fam$ASV<-row.names(res_fam)
  Taxa.fam_pruned$ASV<-row.names(Taxa.fam_pruned)
  res_fam<-dplyr::left_join(res_fam, Taxa.fam_pruned, by="ASV")
  res_fam$Kingdom<-NULL
  res_fam$Genus<-NULL
  res_fam$Species<-NULL
  View(res_fam)
  rm(res)
  write.csv2(res_fam, file="Diff_abundance_DESeq2_WT_Treatment/Analisi diff delle famiglie con DeSeq2 _ ratio_T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_fam$Family
  target<-na.omit(target)
  tax_table(data_pruned.fam.prop)<-as.matrix(Taxa.fam_pruned)
  target<-subset_taxa(data_pruned.fam.prop, Family %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  ##tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
}

# class ###################################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_C<-as.data.frame(res)
  res_C$ASV<-row.names(res_C)
  Taxa.class_pruned$ASV<-row.names(Taxa.class_pruned)
  res_C<-dplyr::left_join(res_C, Taxa.class_pruned, by="ASV")
  res_C$Kingdom<-NULL
  res_C$Genus<-NULL
  res_C$Species<-NULL
  res_C$Order<-NULL
  res_C$Family<-NULL
  View(res_C)
  rm(res)
  write.csv2(res_C, file="Diff_abundance_DESeq2_WT_Treatment/Analisi diff delle classi_ ratio _T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_C$Class
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.class.prop, Class %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
}

# order ############################################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_O<-as.data.frame(res)
  res_O$ASV<-row.names(res_O)
  Taxa.order_pruned$ASV<-row.names(Taxa.order_pruned)
  res_O<-dplyr::left_join(res_O, Taxa.order_pruned, by="ASV")
  res_O$Kingdom<-NULL
  res_O$Genus<-NULL
  res_O$Species<-NULL
  res_O$Family<-NULL
  View(res_O)
  rm(res)
  write.csv2(res_O, file="Diff_abundance_DESeq2_WT_Treatment/Analisi diff degli ordini_ ratio _T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_O$Order
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.order.prop, Order %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
}

# phylum  ########################################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_phy<-as.data.frame(res)
  res_phy$ASV<-row.names(res_phy)
  Taxa.phy_pruned$ASV<-row.names(Taxa.phy_pruned)
  res_phy<-dplyr::left_join(res_phy, Taxa.phy_pruned, by="ASV")
  res_phy$Kingdom<-NULL
  res_phy$Genus<-NULL
  res_phy$Species<-NULL
  res_phy$Class<-NULL
  res_phy$Order<-NULL
  res_phy$Family<-NULL
  View(res_phy)
  rm(res)
  write.csv2(res_phy, file="Diff_abundance_DESeq2_WT_Treatment/Analisi diff dei phyla_ratio _T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_phy$Phylum
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.phy.prop, Phylum %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
}

# plots  ########################################

tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_f$Taxa<-"Families"
tabella_f[,c("Phylum","Order","Class")]<-NULL
colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"

#tabella_o$Taxa<-"Orders"
#tabella_o[,c("Phylum","Class")]<-NULL
#colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"

tabella_c$Taxa<-"Classes"
tabella_c[,"Phylum"]<-NULL
colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"

tabella_p$Taxa<-"Phyla"
colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"

tabella_2<-rbind.data.frame(tabella_f,tabella_c,tabella_p)

# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Time)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Time)


tabella_g$Bacteria<-gsub("Escherichia-Shigella","Escherichia \n -Shighella",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("[Eubacterium]","Eubacterium",tabella_g$Bacteria, fixed = T)
tabella_g$Bacteria<-gsub("_","\n",tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 10) +
  scale_color_manual(values = c("T0"="deepskyblue","T1"="red4")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, labeller = labeller(group = label_wrap_gen(width = 34)),scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time), size=3) +
  geom_line(aes(group=Sample_ID), size=0.5) +
  theme(strip.text.x=element_text(size=9.4,colour="black"), strip.switch.pad.wrap = unit(10,"line")  ) + 
  theme(axis.text.y = element_text(size=7))+
  scale_y_sqrt(breaks=c(0, 0.1, 0.5,seq(1,max(tabella_g$Abundance),1))) +
  scale_x_discrete(labels=rep(unique(levels(tabella_g$Time))),expand=c(0,0.5)) +
  theme(plot.title= element_text(size=12)) +
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
plot_1 +
  labs(title= "Differently abundant genera in WT mices before (T0) and after (T1) the treatment", y="Proportional Abundance", fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_WT_Treatment/Plot_genera_DeSeq2_Time.png", width = 10, height = 5, dpi=300)
dev.off()

plot_2<-ggplot(tabella_2, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 10) +
  scale_color_manual(values = c("T0"="deepskyblue","T1"="red4")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Time), size=2) +
  geom_line(aes(group=Sample_ID)) +
  theme(strip.text.x=element_text(size=10,colour="black"), axis.text.y = element_text(size=8)) +
  scale_x_discrete(labels=rep(unique(levels(tabella_2$Time))), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=12)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  scale_y_sqrt(breaks=c(0, 0.1, 0.5,seq(1,max(tabella_g$Abundance),1)))
#scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))
plot_2 +
  labs(title= "Differently abundant families, classes and phyla in WT mices", y="Proportional Abundance", fill="Time", x="")

ggsave(filename = "Diff_abundance_DESeq2_WT_Treatment/Plot_2_DESeq2_Time.png", width = 10, height = 5, dpi=300)
dev.off()

# now a unique plot
head_plot<-plot_1 + # refining first plot with the title
  labs(title= "Differently abundant taxa in WT mices before (T0) and after (T1) the treatment", y="Proportional Abundance", x="") + 
  theme(plot.title = element_text(size=12))
png(filename = "Diff_abundance_DESeq2_WT_Treatment/Plot_TOTAL_DESeq2_Time.png", width = 2950, height = 2200, res=300)
grid.arrange(head_plot,
             plot_2 + labs(title= "", y="Proportional Abundance", x=""))
dev.off()

rm(tabella_2, plot_2, plot_1, head_plot, tabella_g,tabella_p, tabella_c, tabella_f)

################## DAs WITH DESEQ2 _ T0 vs T1 in KO ##############################

rm(data.temp)
data.temp<-subset_samples(data, Condition=="KO")

# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
rm(data_pruned, data.genus_pruned)
data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
{ data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
  Taxa.genus_pruned<-taxa_temp
  rm(taxa_temp)
}
# adding informations to uncultured families
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
{for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  Taxa.fam_pruned<-taxa_temp
  rm(taxa_temp)
}
# other taxonomic vocabularies
{ Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
  Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
  Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

dir.create("Diff_abundance_DESeq2_KO_Treatment")
system(" echo 'In case of pair analysis, the log2foldchange displayes the change from T0 to T1, see https://support.bioconductor.org/p/105981/' > Diff_abundance_DESeq2_KO_Treatment/Nota_bene.txt ")  

# genus ##############################
rm(res, DE, target, res_genus)
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  View(res_genus)
  rm(res)
  write.csv2(res_genus, file="Diff_abundance_DESeq2_KO_Treatment/Analisi diff dei generi_ratio T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)  
  tax_table(data_pruned.genus.prop)<-as.matrix(Taxa.genus_pruned)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family ##############################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_fam<-as.data.frame(res)
  res_fam$ASV<-row.names(res_fam)
  Taxa.fam_pruned$ASV<-row.names(Taxa.fam_pruned)
  res_fam<-dplyr::left_join(res_fam, Taxa.fam_pruned, by="ASV")
  res_fam$Kingdom<-NULL
  res_fam$Genus<-NULL
  res_fam$Species<-NULL
  View(res_fam)
  rm(res)
  write.csv2(res_fam, file="Diff_abundance_DESeq2_KO_Treatment/Analisi diff delle famiglie con DeSeq2 _ ratio_T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_fam$Family
  target<-na.omit(target)
  tax_table(data_pruned.fam.prop)<-as.matrix(Taxa.fam_pruned)
  target<-subset_taxa(data_pruned.fam.prop, Family %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
}

# class ###################################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_C<-as.data.frame(res)
  res_C$ASV<-row.names(res_C)
  Taxa.class_pruned$ASV<-row.names(Taxa.class_pruned)
  res_C<-dplyr::left_join(res_C, Taxa.class_pruned, by="ASV")
  res_C$Kingdom<-NULL
  res_C$Genus<-NULL
  res_C$Species<-NULL
  res_C$Order<-NULL
  res_C$Family<-NULL
  View(res_C)
  rm(res)
  write.csv2(res_C, file="Diff_abundance_DESeq2_KO_Treatment/Analisi diff delle classi_ ratio _T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_C$Class
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.class.prop, Class %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
}

# order ############################################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_O<-as.data.frame(res)
  res_O$ASV<-row.names(res_O)
  Taxa.order_pruned$ASV<-row.names(Taxa.order_pruned)
  res_O<-dplyr::left_join(res_O, Taxa.order_pruned, by="ASV")
  res_O$Kingdom<-NULL
  res_O$Genus<-NULL
  res_O$Species<-NULL
  res_O$Family<-NULL
  View(res_O)
  rm(res)
  write.csv2(res_O, file="Diff_abundance_DESeq2_KO_Treatment/Analisi diff degli ordini_ ratio _T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_O$Order
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.order.prop, Order %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
}

# phylum  ########################################
rm(res, DE, target)
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~Sample_ID + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_phy<-as.data.frame(res)
  res_phy$ASV<-row.names(res_phy)
  Taxa.phy_pruned$ASV<-row.names(Taxa.phy_pruned)
  res_phy<-dplyr::left_join(res_phy, Taxa.phy_pruned, by="ASV")
  res_phy$Kingdom<-NULL
  res_phy$Genus<-NULL
  res_phy$Species<-NULL
  res_phy$Class<-NULL
  res_phy$Order<-NULL
  res_phy$Family<-NULL
  View(res_phy)
  rm(res)
  write.csv2(res_phy, file="Diff_abundance_DESeq2_KO_Treatment/Analisi diff dei phyla_ratio _T0 vs T1.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_phy$Phylum
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.phy.prop, Phylum %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
}

# plots  ########################################

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

tabella_2<-rbind.data.frame(tabella_f)
tabella_3<-rbind.data.frame(tabella_o)
tabella_4<-rbind.data.frame(tabella_c,tabella_p)

# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Time)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Time)
tabella_3$Xaxis<-paste0(tabella_3$Bacteria, tabella_3$Time)
tabella_4$Xaxis<-paste0(tabella_4$Bacteria, tabella_4$Time)

unique(tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("Escherichia-Shigella","Escherichia \n -Shighella",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("[Eubacterium]","Eubacterium",tabella_g$Bacteria, fixed = T)
tabella_g$Bacteria<-gsub("_","\n",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("Clostridium\nsensu\nstricto\n1","Clostridium\nsensu\nstricto_1",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("f ","of family\n",tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 12) +
  scale_color_manual(values = c("T0"="coral","T1"="red4")) +
  facet_wrap2(nrow=3, factor(Taxa,levels = "Genera") ~Bacteria, labeller = labeller(group = label_wrap_gen(width = 34)),scales = "free_x", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time), size=3) +
  geom_line(aes(group=Sample_ID), size=0.7) +
  theme(strip.text.x=element_text(size=9.45,colour="black"), strip.switch.pad.wrap = unit(5,"line")  ) + 
  theme(axis.text.y = element_text(size=9.1))+
  theme(panel.spacing.x = unit(1, "pt"))+
  scale_y_sqrt(breaks=c(0.2, 1,2,4,6,8,seq(12,max(tabella_g$Abundance),4))) +
  scale_x_discrete(labels=rep(unique(levels(tabella_g$Time))),expand=c(0,0.5)) +
  theme(plot.title= element_text(size=15)) +
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
plot_1 +
  labs(title= "Differently abundant genera in Msc-/- mices before (T0) and after (T1) the treatment", y="Proportional Abundance", fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_KO_Treatment/Plot_genera_DeSeq2_Time.png", width = 13.7, height = 11, dpi=300)
dev.off()

plot_2<-ggplot(tabella_2, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 10) +
  scale_color_manual(values = c("T0"="coral","T1"="red4")) +
  facet_wrap2(nrow=2,factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Time), size=2) +
  geom_line(aes(group=Sample_ID)) +
  theme(strip.text.x=element_text(size=10,colour="black"), 
        axis.text.y = element_text(size=8.1)) +
  scale_x_discrete(labels=rep(unique(levels(tabella_2$Time))), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=12)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,2,4,6,8,seq(12,max(tabella_g$Abundance),4)))

plot_2 +
  labs(title= "Differently abundant families in Msc-/- mices before (T0) and after (T1) the treatment", y="Proportional Abundance", fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_KO_Treatment/Plot_families_DESeq2_Time.png", width = 10.5, height = 7.5, dpi=300)
dev.off()

tabella_3$Bacteria<-gsub("Peptostreptococcales-Tissierellales","Peptostreptococcales\n-Tissierellales",tabella_3$Bacteria)
plot_3<-ggplot(tabella_3, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 10) +
  scale_color_manual(values = c("T0"="coral","T1"="red4")) +
  facet_wrap2(nrow=2,factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Time), size=2) +
  geom_line(aes(group=Sample_ID)) +
  theme(strip.text.x=element_text(size=10,colour="black"), axis.text.y = element_text(size=8)) +
  scale_x_discrete(labels=rep(unique(levels(tabella_2$Time))), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=12)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,2,4,8,seq(12,max(tabella_g$Abundance),4)))
  
plot_3 +
  labs(title= "Differently abundant orders in Msc-/- mices before (T0) and after (T1) the treatment", y="Proportional Abundance", fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_KO_Treatment/Plot_orders_DESeq2_Time.png", width = 10, height = 7.5, dpi=300)
dev.off()

tabella_4$Bacteria<-gsub("Peptostreptococcales-Tissierellales","Peptostreptococcales\n-Tissierellales",tabella_4$Bacteria)
tabella_4<-tabella_4[! tabella_4$Bacteria == "Cyanobacteria",] # almost absent and, clearly, a sequencing contaminat
plot_4<-ggplot(tabella_4, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 10) +
  scale_color_manual(values = c("T0"="coral","T1"="red4")) +
  facet_wrap2(nrow=2,factor(Taxa,levels = c("Classes","Phyla"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Time), size=2) +
  geom_line(aes(group=Sample_ID)) +
  theme(strip.text.x=element_text(size=10,colour="black"), axis.text.y = element_text(size=8)) +
  scale_x_discrete(labels=rep(unique(levels(tabella_2$Time))), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=12)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,2,4,6,8,seq(12,max(tabella_g$Abundance),4)))
  
plot_4 +
  labs(title= "Differently abundant classes and phyla in Msc-/- mices before (T0) and after (T1) the treatment", y="Proportional Abundance", fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_KO_Treatment/Plot_class_phylum_DESeq2_Time.png", width = 10, height = 7, dpi=300)
dev.off()

rm(tabella_2, plot_2, tabella_3, plot_1, plot_3, plot_4, tabella_f, tabella_o, tabella_g, tabella_c, tabella_p)


##################### DA WITH DESEQ2 _ WT vs KO at T0 #######################

rm(data.temp, data_pruned, data.genus_pruned, tabella_f, tabella_c, tabella_p, tabella_o)
data.temp<-subset_samples(data, Time=="T0")

# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
{data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
  data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
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
  Taxa.genus_pruned<-taxa_temp
}
# adding informations to uncultured families
{taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
  for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
    taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  Taxa.fam_pruned<-taxa_temp
  rm(taxa_temp)
}
# other taxonomic vocabularies
{Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
  Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
  Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

dir.create("Diff_abundance_DESeq2_WT_vs_KO_T0")

# genus
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  View(res_genus)
  rm(res)
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_fam<-as.data.frame(res)
  res_fam$ASV<-row.names(res_fam)
  Taxa.fam_pruned$ASV<-row.names(Taxa.fam_pruned)
  res_fam<-dplyr::left_join(res_fam, Taxa.fam_pruned, by="ASV")
  res_fam$Kingdom<-NULL
  res_fam$Genus<-NULL
  res_fam$Species<-NULL
  View(res_fam)
  rm(res)
  #write.csv2(res_fam, file="Diff_abundance_DESeq2_WT_vs_KO_T0/Analysis_of_diff_families_DeSeq2_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_fam$Family
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.fam.prop, Family %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
}

# class
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_C<-as.data.frame(res)
  res_C$ASV<-row.names(res_C)
  Taxa.class_pruned$ASV<-row.names(Taxa.class_pruned)
  res_C<-dplyr::left_join(res_C, Taxa.class_pruned, by="ASV")
  res_C$Kingdom<-NULL
  res_C$Genus<-NULL
  res_C$Species<-NULL
  res_C$Order<-NULL
  res_C$Family<-NULL
  View(res_C)
  rm(res)
  #write.csv2(res_C, file="Diff_abundance_DESeq2_WT_vs_KO_T0/Analysis_of_diff_classes_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_C$Class
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.class.prop, Class %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
}

# order
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_O<-as.data.frame(res)
  res_O$ASV<-row.names(res_O)
  Taxa.order_pruned$ASV<-row.names(Taxa.order_pruned)
  res_O<-dplyr::left_join(res_O, Taxa.order_pruned, by="ASV")
  res_O$Kingdom<-NULL
  res_O$Genus<-NULL
  res_O$Species<-NULL
  res_O$Family<-NULL
  View(res_O)
  rm(res)
  #write.csv2(res_O, file="Diff_abundance_DESeq2_WT_vs_KO_T0/Analysis_of_diff_orders_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_O$Order
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.order.prop, Order %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
}

# phylum
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_phy<-as.data.frame(res)
  res_phy$ASV<-row.names(res_phy)
  Taxa.phy_pruned$ASV<-row.names(Taxa.phy_pruned)
  res_phy<-dplyr::left_join(res_phy, Taxa.phy_pruned, by="ASV")
  res_phy$Kingdom<-NULL
  res_phy$Genus<-NULL
  res_phy$Species<-NULL
  res_phy$Class<-NULL
  res_phy$Order<-NULL
  res_phy$Family<-NULL
  View(res_phy)
  rm(res)
  #write.csv2(res_phy, file="Diff_abundance_DESeq2_WT_vs_KO_T0/Analysis_of_diff_phyla_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_phy$Phylum
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.phy.prop, Phylum %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
}

# plots  ########################################

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

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_c,tabella_o,tabella_p)

tabella_tot$Condition<-factor(tabella_tot$Condition, levels = c("WT","KO"))
ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("WT"="chartreuse","KO"="coral"), labels=c("WT","Msc-/-")) +
  geom_boxplot(width=0.8) + theme_bw( ) + theme(strip.text.x=element_text(size=11,colour="black")) + 
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,axis.text.x = element_text(angle = -50, vjust=0.5, hjust=0, size=10), axis.text.y = element_text(size=9)) + 
  scale_x_discrete(expand=c(-0.2, 1)) + theme(plot.title= element_text(size=14) ,
                                              legend.key.size=unit(0.7,"cm"), 
                                              legend.text=element_text(size=10)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  scale_y_sqrt(breaks=c(0.1, 1,2,4,6,8,seq(12,max(tabella_g$Abundance),4))) +
  labs(title= "Differently abundant taxa between WT and Msc-/- mices before the treatment", y="Proportional Abundance", fill="Condition", x="")
#scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))

# cutting off the redundant results
tabella_tot<-rbind.data.frame(tabella_g,tabella_f)
tabella_tot[tabella_tot$Bacteria=="Muribaculaceae" & tabella_tot$Taxa=="Families","Bacteria"]<-"Muribaculaceae_f"
tabella_tot<-tabella_tot[! tabella_tot$Bacteria %in% c("Muribaculaceae_f","Monoglobaceae","Lactobacillaceae"),]
tabella_tot$Condition<-factor(tabella_tot$Condition, levels = c("WT","KO"))
ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("WT"="chartreuse","KO"="coral"), labels=c("WT","Msc-/-")) +
  geom_boxplot(width=0.8) + 
  theme_bw( ) +
  theme(strip.text.x=element_text(size=11,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), # to rise legend position
        legend.position="bottom" ,
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=10)) + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  theme(plot.title= element_text(size=14) ,
        legend.key.size=unit(0.8, "cm"), 
        legend.title = element_text(size=13),
        legend.text=element_text(size=13)) +
  theme(panel.grid.minor.y= element_blank()) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(plot.margin = unit(c(0.1,0.1,0,0.8),"cm") )+
  scale_y_sqrt(breaks=c(0.2, 0.5, 1,2,4,6,8,seq(12,max(tabella_g$Abundance),4))) +
  labs(title= "Differently abundant taxa between WT and Msc-/- mices before the treatment", y="Proportional Abundance", fill="Condition", x="")

ggsave(filename = "Diff_abundance_DESeq2_WT_vs_KO_T0/Total_differences_DeSeq2_Condition_NO_REDUNDANTS.png", width = 10, height = 7, dpi=300)
dev.off()

##################### DA WITH DESEQ2 _ WT vs KO at T1 #######################

rm(data.temp, data_pruned, data.genus_pruned, tabella_g, tabella_f, tabella_c, tabella_p, tabella_o)
data.temp<-subset_samples(data, Time=="T1")

# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
{data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
  data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
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
  Taxa.genus_pruned<-taxa_temp
}
# adding informations to uncultured families
{taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
  for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
    taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  Taxa.fam_pruned<-taxa_temp
  rm(taxa_temp)
}
# other taxonomic vocabularies
{Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
  Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
  Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

dir.create("Diff_abundance_DESeq2_WT_vs_KO_T1")

# genus
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  View(res_genus)
  rm(res)
  write.csv2(res_genus, file="Diff_abundance_DESeq2_WT_vs_KO_T1/Analysis_of_diff_genera_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_fam<-as.data.frame(res)
  res_fam$ASV<-row.names(res_fam)
  Taxa.fam_pruned$ASV<-row.names(Taxa.fam_pruned)
  res_fam<-dplyr::left_join(res_fam, Taxa.fam_pruned, by="ASV")
  res_fam$Kingdom<-NULL
  res_fam$Genus<-NULL
  res_fam$Species<-NULL
  View(res_fam)
  rm(res)
  write.csv2(res_fam, file="Diff_abundance_DESeq2_WT_vs_KO_T1/Analysis_of_diff_families_DeSeq2_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_fam$Family
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.fam.prop, Family %in% target)
  tabella<-psmelt(target)
  tabella$ASV<-NULL
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
}

# class
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_C<-as.data.frame(res)
  res_C$ASV<-row.names(res_C)
  Taxa.class_pruned$ASV<-row.names(Taxa.class_pruned)
  res_C<-dplyr::left_join(res_C, Taxa.class_pruned, by="ASV")
  res_C$Kingdom<-NULL
  res_C$Genus<-NULL
  res_C$Species<-NULL
  res_C$Order<-NULL
  res_C$Family<-NULL
  View(res_C)
  rm(res)
  write.csv2(res_C, file="Diff_abundance_DESeq2_WT_vs_KO_T1/Analysis_of_diff_classes_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_C$Class
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.class.prop, Class %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
}

# order
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_O<-as.data.frame(res)
  res_O$ASV<-row.names(res_O)
  Taxa.order_pruned$ASV<-row.names(Taxa.order_pruned)
  res_O<-dplyr::left_join(res_O, Taxa.order_pruned, by="ASV")
  res_O$Kingdom<-NULL
  res_O$Genus<-NULL
  res_O$Species<-NULL
  res_O$Family<-NULL
  View(res_O)
  rm(res)
  write.csv2(res_O, file="Diff_abundance_DESeq2_WT_vs_KO_T1/Analysis_of_diff_orders_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_O$Order
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.order.prop, Order %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
}

# phylum
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "KO", "WT"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_phy<-as.data.frame(res)
  res_phy$ASV<-row.names(res_phy)
  Taxa.phy_pruned$ASV<-row.names(Taxa.phy_pruned)
  res_phy<-dplyr::left_join(res_phy, Taxa.phy_pruned, by="ASV")
  res_phy$Kingdom<-NULL
  res_phy$Genus<-NULL
  res_phy$Species<-NULL
  res_phy$Class<-NULL
  res_phy$Order<-NULL
  res_phy$Family<-NULL
  View(res_phy)
  rm(res)
  write.csv2(res_phy, file="Diff_abundance_DESeq2_WT_vs_KO_T1/Analysis_of_diff_phyla_ratio_KO_vs_WT.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_phy$Phylum
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.phy.prop, Phylum %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
}

# plots  ########################################

tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_g$Condition<-factor(tabella_g$Condition, levels = c("WT","KO"))
tabella_g$Bacteria<-gsub("_group","",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_","\n",tabella_g$Bacteria)
ggplot(tabella_g, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("WT"="chartreuse","KO"="coral"), labels=c("WT","Msc-/-")) +
  geom_boxplot(width=0.8) + 
  theme_bw( ) +
  theme(strip.text.x=element_text(size=11,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), # to rise legend position
        legend.position="bottom" ,
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=11), 
        axis.text.y = element_text(size=11)) + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  theme(plot.title= element_text(size=14) ,
        legend.key.size=unit(0.8, "cm"), 
        legend.title = element_text(size=13),
        legend.text=element_text(size=13)) +
  theme(panel.grid.minor.y= element_blank()) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(plot.margin = unit(c(0.1,0.1,0,1),"cm") )+
  scale_y_sqrt(breaks=c(0.2, 0.5, 1,2,4,6,8)) +
  labs(title= "Differently abundant taxa between WT and Msc-/- mices after the treatment", y="Proportional Abundance", fill="Condition", x="")

ggsave(filename = "Diff_abundance_DESeq2_WT_vs_KO_T1/Total_differences_DeSeq2_Condition_NO_REDUNDANTS.png", width = 10, height = 7, dpi=300)
dev.off()

system("touch Diff_abundance_DESeq2_WT_vs_KO_T1/'Nothing else here'")


################# PLOTTING PICRUST2_LEFSE RESULTS ######################

# this step requires file from PICRUST2_LEFSE analysis, see the processing Script

dir.create("PICRUST_LEFSE_results")
# setwd("Results") # if not already there

########### WT vs K0 at T0

rm(a)
a <- read.delim("../PICRUST2_LEFSE/picrust2_T0/pathways_out/path_abun_unstrat_descrip.tsv.gz") 

Descriptions<-a[,c("pathway","description")]
head(Descriptions, n=4)

Significative_functions_LEFSE<- read.delim("../PICRUST2_LEFSE/LEFSE/Risultato_WT_vs_KO_T0.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.csv2(Significative_functions_LEFSE, file="PICRUST_LEFSE_results/Significant_different_pathways_METACYC_WT_vs_KO_T0.csv", quote = F, row.names = F)

########### plotting results

# cutting names too much long
Significative_functions_LEFSE$Pathway<-gsub("superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation","N-acetylglucosamine, N-acetylmannosamine, N-acetylneuraminate degradation", Significative_functions_LEFSE$Pathway)
Significative_functions_LEFSE$Pathway<-gsub("UDP-N-acetylmuramoyl-pentapeptide biosynthesis I (meso-diaminopimelate containing)","UDP-N-acetylmuramoyl-pentapeptide biosynthesis I (meso-diaminopimelate)", Significative_functions_LEFSE$Pathway)
# setting the labels/bar directions
{Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="WT"]<-Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="WT"]*-1
  Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ]
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical sorting
}
# plot 1
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$logLDA_score>0,1,0))+
  scale_fill_manual(values = c("WT"="chartreuse","KO"="coral"),labels=c("WT","Msc-/-")) +
  theme_classic(base_size = 22) +
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.4)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 22)) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1,2,3))+
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_T0_ALL_RESULTS.png", width = 13.1, height = 18, dpi = 300)
# plot 2 (only effect size >3)
Significative_functions_LEFSE<-Significative_functions_LEFSE[abs(Significative_functions_LEFSE$logLDA_score)>3,]
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$logLDA_score>0,1,0), size=5)+
  scale_fill_manual(values = c("WT"="chartreuse","KO"="coral"),labels=c("WT","Msc-/-")) +
  theme_classic(base_size = 25) +
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 25)) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1,2,3))+
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_T0_ONLY_OVER_3.png", width = 15, height = 9, dpi = 300)

rm(a, Descriptions, Significative_functions_LEFSE)

########### WT vs K0 at T1

a <- read.delim("../PICRUST2_LEFSE/picrust2_T1/pathways_out/path_abun_unstrat_descrip.tsv.gz") 

Descriptions<-a[,c("pathway","description")]
head(Descriptions, n=4)

Significative_functions_LEFSE<- read.delim("../PICRUST2_LEFSE/LEFSE/Risultato_WT_vs_KO_T1.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.csv2(Significative_functions_LEFSE, file="PICRUST_LEFSE_results/Significant_different_pathways_METACYC_WT_vs_KO_T1.csv", quote = F, row.names = F)

########### plotting results

# cutting names too much long
Significative_functions_LEFSE$Pathway<-gsub("superpathway of &beta;-D-glucuronide and D-glucuronate degradation","superpathway of &beta;-D-glucuronide and D-glucuronate degradation", Significative_functions_LEFSE$Pathway)
Significative_functions_LEFSE$Pathway<-gsub("superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis (E. coli)","pyrimidine deoxyribonucleotides de novo biosynthesis (E. coli)", Significative_functions_LEFSE$Pathway, fixed=T)
Significative_functions_LEFSE$Pathway<-gsub("superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis","GDP-mannose-derived O-antigen building blocks biosynthesis", Significative_functions_LEFSE$Pathway)
# setting the labels/bar directions
{Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="WT"]<-Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="WT"]*-1
  Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ]
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical sorting
}
# plot 1
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), size=5, hjust = ifelse(Significative_functions_LEFSE$logLDA_score>0,1,0))+
  scale_fill_manual(values = c("WT"="chartreuse","KO"="coral"), labels=c("WT","Msc-/-")) +
  theme_classic(base_size = 22) +
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.4)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 22)) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1,2,3))+
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_T1_ALL_RESULTS.png", width = 13.4, height = 18, dpi = 300)
# plot 2 (only effect size >3)
Significative_functions_LEFSE<-Significative_functions_LEFSE[abs(Significative_functions_LEFSE$logLDA_score)>3,]
system("touch PICRUST_LEFSE_results/'Nessuna_oltre_3_per_i_T1' ")

rm(a, Descriptions, Significative_functions_LEFSE)

##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("R_version_and_packages.txt")
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
cat("\n", "\n", fill=TRUE)
print("egg")
packageVersion("egg")
cat("\n", "\n", fill=TRUE)
print("ggh4x")
packageVersion("ggh4x")
cat("\n", "\n", fill=TRUE)
package$otherPkgs$pca3d[c(1,4)]

sink() # restore STR OUTPUT to R console
close(con)
rm(con)