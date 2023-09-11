##################### PREPARING THE ENVIRONMENT ################

{ library("phyloseq")
  library("ggplot2")
  library("vegan")
  library("ggpubr")
  library("ggh4x")
  library("dendextend")
  library("mixOmics")
  library("DESeq2")
  library("Hmisc")
  library("egg")
  library("dplyr")
  library("xlsx")
  library("qiime2R")
}

dir.create("Data_check")
dir.create("Results")

options(scipen = 100) # disable scientific annotation


####################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree="QIIME/rooted-tree.qza")
# changing names
original_names<-sample_names(data)
original_names
Metadata_complete <- as.data.frame(read.csv(file = "metadata_stool_saliva.csv") )
# removing eventual fastq from the other sequencing (already published)
from_other_sequencing<-Metadata_complete[Metadata_complete$Sample_Type=="Stool", "FASTQ_name"]
original_names<-original_names[!original_names %in% from_other_sequencing]
Metadata<-Metadata_complete[Metadata_complete$Sample_Type=="Saliva", ]
# updating
row.names(Metadata)<-Metadata$FASTQ_name
original_length<-length(Metadata$FASTQ_name[!is.na(Metadata$FASTQ_name)])
head(Metadata)
Metadata<-Metadata[original_names,]
if(length(which(is.na(Metadata$Sample_ID)))>1){
  stop("\nError: Same obs is NA, check the rowname matchings\n")
} else { cat("Everything is OK")}
identical(length(original_names),length(Metadata$FASTQ_name)) # TRUE
identical(original_length,length(Metadata$FASTQ_name))
sample_data(data)<-Metadata # update
sample_names(data)<-sample_data(data)$Sample_ID

sample_data(data)$Treatment_time<-factor(sample_data(data)$Treatment_time, levels = c("T0","T24","No_treatm"))
sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels = c("Healthy","HIV","Treatment"))

data_tot<-data # ready to be subsetted later

rm(original_length, original_names,data)

# <-- save here!
save.image("data_10settembre_only_saliva.RData")


############ CHECKING THE DOMAIN PROPORTION ######################

{Unass<-tax_glom(data_tot, taxrank = "Kingdom", NArm =F ) # or domain
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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Domains_proportion_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


################# CHECKING % ASSIGNED IN SILVA #########################

{data_tot.phy = tax_glom(data_tot, taxrank = "Phylum", NArm = F)
data_tot.class = tax_glom(data_tot, taxrank = "Class", NArm = F)
data_tot.order = tax_glom(data_tot, taxrank = "Order", NArm = F)
data_tot.fam = tax_glom(data_tot, taxrank = "Family", NArm = F)
data_tot.genus = tax_glom(data_tot, taxrank = "Genus", NArm = F)
}

{Taxa.genus<-as.data.frame(tax_table(data_tot.genus))
  Taxa.fam<-as.data.frame(tax_table(data_tot.fam))
  Taxa.phy<-as.data.frame(tax_table(data_tot.phy))
  Taxa.class<-as.data.frame(tax_table(data_tot.class))
  Taxa.order<-as.data.frame(tax_table(data_tot.order))
}

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assegnati<-rbind.data.frame(a,b,c,d,e)
colnames(assegnati)<-c("Totale","Assegnati","%","Taxa")
}
assegnati
write.csv2(assegnati,file="Data_check/Percents_taxa_unassigned_in_SILVA_database.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assegnati)


###################### GENERAL PCoA #######################

suppressWarnings(rm(data.prop,DistBC,data.sqrt_prop))
# on ASV
data.prop <- transform_sample_counts(data_tot, function(ASV) ASV/sum(ASV)*100)
data.prop.labels<-data.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  geom_point(size=2.9, alpha=0.3) + theme_classic(base_size = 9) + stat_ellipse(size=0.13) +
  labs(title="General PCoA with Hellinger distance ", color="Condition", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="Results/PCoA_Beta_diversity_ON_EVERYTHING_euclidean.png", width = 4.5, height = 4, dpi=300)



##################### RAREFACTION/SATURATION ANALYSIS ######################

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
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

data_saliva<-subset_samples(data_tot, Sample_Type=="Saliva")
png(file="Data_check/Rarefaction_curve_Saliva.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data_saliva),"matrix")), step=100,label=F)
evalslopes(r,sample_names(data_saliva),lim=0.001,cex=1)
dev.off()
rm(r)


##################### ////// STARTING SALIVA ANALYSIS ////// ##############################

{ if( grepl("_analysis",getwd()) ){setwd("..")}
  if( ! grepl("Results",getwd()) ){setwd("Results")} # it sets wd to "Results"
}

{dir.create("Saliva_bacteriota_analysis")
  setwd("Saliva_bacteriota_analysis")
  dir.create("Abundances")
  dir.create("Responders_CD4_Analysis")
  dir.create("T0_vs_T24")
  dir.create("Healthy_vs_HIV")
}

{to_remove<-ls()
  to_remove<-to_remove[ ! to_remove %in% c("data_tot","Metadata") ]
  to_remove<- to_remove[! grepl(pattern = "DA_",to_remove, fixed = T) ]  # to mantain DESEQ results
  to_remove
  rm(list = to_remove) # reset the workspace
}

data<-subset_samples(data_tot,Sample_Type=="Saliva")
head(sample_data(data))

{data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
  Taxa.genus<-as.data.frame(tax_table(data.genus))
  data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
  data.class<- tax_glom(data, taxrank = "Class", NArm = F)
  data.order<- tax_glom(data, taxrank = "Order", NArm = F)
  data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
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


taxa_temp<-as.data.frame(tax_table(data.fam))
{for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  for( x in 1: length(which(is.na(taxa_temp$Family))) ) {
    taxa_temp$Family[ which(is.na(taxa_temp$Family))[1] ]<-paste("NA_ o",taxa_temp[which(is.na(taxa_temp$Family))[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="NA_ o NA")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="NA_ o NA")[1] ]<-paste("NA_ c",taxa_temp[which(taxa_temp$Family=="NA_ c NA")[1],"Class"])}
  for( x in 1: length(which(duplicated(taxa_temp$Family[taxa_temp$Family=="NA_ c NA"]))) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="NA_ c NA")[1] ]<-paste("NA_ c NA",x+1) }
  Taxa.fam.update<-taxa_temp
}


{data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
  data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
  data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}

system("echo 'The comparison 'Healthy vs HIV' as been performed on healthy samples and HIV samples -->at time T0<-- ' > Healthy_vs_HIV/Samples_of_this_comparison.txt")

##### PREPARING OBJECTS FOR RESPONDERS ANALYSIS (CD4/CD8 responders after treatment)
data.responders<-subset_samples(data, Responders %in% c(0,1)) 
sample_data(data.responders)$Responders<-factor(sample_data(data.responders)$Responders, levels=c("0","1"))
Meta_responders<-Metadata[Metadata$Responders %in% c(0,1),]
{data.responders.genus<- tax_glom(data.responders, taxrank = "Genus", NArm = F)
  Taxa.genus<-as.data.frame(tax_table(data.responders.genus))
  data.responders.fam<- tax_glom(data.responders, taxrank = "Family", NArm = F)
  data.responders.class<- tax_glom(data.responders, taxrank = "Class", NArm = F)
  data.responders.order<- tax_glom(data.responders, taxrank = "Order", NArm = F)
  data.responders.phy<- tax_glom(data.responders, taxrank = "Phylum", NArm = F)
  data.responders.prop <- transform_sample_counts(data.responders, function(otu) otu/sum(otu)*100)
  data.responders.phy.prop <- transform_sample_counts(data.responders.phy, function(otu) otu/sum(otu)*100)
  data.responders.class.prop <- transform_sample_counts(data.responders.class, function(otu) otu/sum(otu)*100)
  data.responders.order.prop <- transform_sample_counts(data.responders.order, function(otu) otu/sum(otu)*100)
  data.responders.fam.prop <- transform_sample_counts(data.responders.fam, function(otu) otu/sum(otu)*100)
  data.responders.genus.prop <- transform_sample_counts(data.responders.genus, function(otu) otu/sum(otu)*100)
}


########################### COUNTS EXPORT #####################################

dir.create("Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")), file="Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


#################### ABUNDANCES BAR PLOT (SALIVA) ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one

# TOP 5 Phylum con stacked bar plot
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
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Treatment_time<-gsub("No_treatm","Healthy",tabella$Treatment_time)

tabella_healthy<-tabella[tabella$Condition=="Healthy",]
plot_healthy<-ggplot(data=tabella_healthy, aes(x=Abundance, y=sample_Sample, fill=Phylum)) + theme_classic(base_size =14) + 
  facet_grid2(Treatment_time~.,
              scales = "free", space="free", strip = strip_nested(size="constant", bleed = T))+
  theme(panel.spacing.y = unit(0,"pt"))+
  scale_y_discrete (expand = c(-0.4,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9.5),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(2,2,-4,-4), "mm"),
        title = element_text(size=15)) +
  guides(fill="none") +
  labs(y="", x="", title = "Five most abundant phyla")

tabella_HIV<-tabella[tabella$Condition!="Healthy",]
plot_HIV<-ggplot(data=tabella_HIV, aes(x=Abundance, y=sample_Sample, fill=Phylum)) + theme_classic(base_size =14) + 
  facet_grid2(sample_Sample+Treatment_time~.,
              scales = "free", space="free", strip = strip_nested(size="constant", bleed = T))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10), 
        strip.text.y = element_text(angle=0, size=9),
        legend.margin = margin(-2,-5,5,-5),
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom",
        plot.margin = unit(c(-2,2, 0,-2), "mm")) +
  guides(fill=guide_legend(nrow=2)) +
  labs(y="", x="Percent abundance", caption = " 'Others' includes every phylum below rank 5")

ggarrange(plot_healthy,plot_HIV,nrow=2) %>%
  ggexport(filename ="Abundances/TOP5_phyla_abundances.png",width=1800,height=3200, res=300) 
dev.off()

suppressWarnings(rm(tabella,prune.data.others))


# TOP 5 Genera
suppressWarnings(rm(top, others, tabella, tabella_healthy, tabella_HIV, plot_healthy, plot_HIV))
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
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Treatment_time<-gsub("No_treatm","Healthy",tabella$Treatment_time)

tabella_healthy<-tabella[tabella$Condition=="Healthy",]
plot_healthy<-ggplot(data=tabella_healthy, aes(x=Abundance, y=sample_Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(Treatment_time~.,
              scales = "free", space="free", strip = strip_nested(size="constant", bleed = T))+
  theme(panel.spacing.y = unit(0,"pt"))+
  scale_y_discrete (expand = c(-0.4,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9.5),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(2,2,-4,-4), "mm"),
        title = element_text(size=15)) +
  guides(fill="none") +
  labs(y="", x="", title = "Five most abundant genera")

tabella_HIV<-tabella[tabella$Condition!="Healthy",]
plot_HIV<-ggplot(data=tabella_HIV, aes(x=Abundance, y=sample_Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(sample_Sample+Treatment_time~.,
              scales = "free", space="free", strip = strip_nested(size="constant", bleed = T))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10), 
        strip.text.y = element_text(angle=0, size=9),
        legend.margin = margin(-5,-5,0,-5),
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom",
        plot.margin = unit(c(-2,2, 0,-2), "mm")) +
  guides(fill=guide_legend(nrow=2)) +
  labs(y="", x="Percentual abundance", caption = " 'Others' includes every genus below rank 5")

ggarrange(plot_healthy,plot_HIV,nrow=2) %>%
  ggexport(filename ="Abundances/TOP5_Genera_abundances.png",width=1800,height=3200, res=300) 
dev.off()

suppressWarnings(rm(tabella,prune.data.others, fill_color_5, prune.dat_top, top, others, tabella_top, tabella_healthy, tabella_HIV, plot_HIV, plot_healthy))


# TOP 8 Genera
suppressWarnings(rm(top, others, tabella, tabella_healthy, tabella_HIV, plot_healthy, plot_HIV))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
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
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Treatment_time<-gsub("No_treatm","Healthy",tabella$Treatment_time)

tabella_healthy<-tabella[tabella$Condition=="Healthy",]
plot_healthy<-ggplot(data=tabella_healthy, aes(x=Abundance, y=sample_Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(Treatment_time~.,
              scales = "free", space="free", 
              strip = strip_nested(size="constant", bleed = T))+
  theme(panel.spacing.y = unit(0,"pt"))+
  scale_y_discrete (expand = c(-0.4,0) ) +
  scale_fill_manual(values=fill_color_8) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9.5),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(2,2,-4,-4), "mm"),
        title = element_text(size=15)) +
  guides(fill="none") +
  labs(y="", x="", title = "Eight most abundant genera")

tabella_HIV<-tabella[tabella$Condition!="Healthy",]
plot_HIV<-ggplot(data=tabella_HIV, aes(x=Abundance, y=sample_Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(sample_Sample+Treatment_time~.,
              scales = "free", space="free", strip = strip_nested(size="constant", bleed = T))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_8) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10), 
        strip.text.y = element_text(angle=0, size=9),
        legend.margin = margin(-10,-5,0,-5),
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text ( size = 11 )) + 
  theme(legend.position="bottom",
        plot.margin = unit(c(-2,2, 0,-2), "mm")) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="", x="Percent abundance", caption = " 'Others' includes every genus below rank 8")

ggarrange(plot_healthy,plot_HIV,nrow=2) %>%
  ggexport(filename ="Abundances/TOP8_Genera_abundances.png",width=1800,height=3220, res=300) 
dev.off()

suppressWarnings(rm(tabella,prune.data.others, fill_color_8, prune.dat_top, top, others, tabella_top, tabella_healthy, tabella_HIV, plot_HIV, plot_healthy))


##################### HIERARCHICAL CLUSTERING (SALIVA) ########################

#### Treatment
suppressWarnings(rm(data_sub,c,tabella_colore))
data_sub<-subset_samples(data.prop, Treatment_time %in% c("T0","T24"))
{c<-hclust(vegan::vegdist(t(sqrt(otu_table(data_sub))),method = "euclidean"))
  c$labels<-as.character(sample_data(data_sub)$Sample)
  c<-as.dendrogram(c)
  tabella_colore<-as.data.frame(cbind(as.character(sample_data(data_sub)$Treatment_time),as.character(sample_data(data_sub)$Treatment_time)))
  colnames(tabella_colore)<-c("Gruppo","Colore")
  colors <- gsub("T24","chartreuse",tabella_colore$Colore)
  colors <- gsub("T0","coral",colors)
  labels_colors(c) <- colors[order.dendrogram(c)]
}
png(file="T0_vs_T24/Hierarchical cluster euclidean on sqrt prop normalized ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,4,4), cex.lab=0.9, cex.main=1.35, cex.sub=1.2)
plot(c,main="Community structure using Hellinger distance on sqrt proportional ASVs",
     sub="T0 = orange     T24 = green")
dev.off()

#### Healthy vs HIV
suppressWarnings(rm(data_sub,c,tabella_colore))
data_sub<-subset_samples(data.prop, Condition %in% c("Healthy","HIV"))
{c<-hclust(vegan::vegdist(t(sqrt(otu_table(data_sub))),method = "euclidean"))
  c$labels<-as.character(sample_data(data_sub)$Sample)
  c<-as.dendrogram(c)
  tabella_colore<-as.data.frame(cbind(as.character(sample_data(data_sub)$Condition),as.character(sample_data(data_sub)$Condition)))
  colnames(tabella_colore)<-c("Gruppo","Colore")
  colors <- gsub("Healthy","deepskyblue",tabella_colore$Colore)
  colors <- gsub("HIV","coral",colors)
  labels_colors(c) <- colors[order.dendrogram(c)]
}
png(file="Healthy_vs_HIV/Hierarchical cluster euclidean on sqrt prop normalized ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,4,4), cex.lab=0.8, cex.main=1.35, cex.sub=1.2)
plot(c,main="Community structure using Hellinger distance on sqrt proportional ASVs",
     sub="HIV = orange     Healthy = blue")
dev.off()

#### Responding CD4
suppressWarnings(rm(data_sub,c,tabella_colore))
data_sub<-data.responders.prop
{c<-hclust(vegan::vegdist(t(sqrt(otu_table(data_sub))),method = "euclidean"))
  c$labels<-as.character(sample_data(data_sub)$Sample)
  c<-as.dendrogram(c)
  tabella_colore<-as.data.frame(cbind(as.character(sample_data(data_sub)$Responders),as.character(sample_data(data_sub)$Responders)))
  colnames(tabella_colore)<-c("Gruppo","Colore")
  colors <- gsub("1","chartreuse",tabella_colore$Colore)
  colors <- gsub("0","chartreuse4",colors)
  labels_colors(c) <- colors[order.dendrogram(c)]
}
png(file="Responders_CD4_Analysis/Hierarchical cluster euclidean on sqrt prop normalized ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,4,4), cex.lab=0.9, cex.main=1.35, cex.sub=1.2)
plot(c,main="Community structure using Hellinger distance on sqrt proportional ASVs",
     sub="Responding = light green     Not responding = dark green")
dev.off()


################ ALFA DIVERSITY ANALYSIS (SALIVA) #######################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

######## Treatment
suppressWarnings(rm(pAlpha, data_sub))
data_sub<-subset_samples(data, Treatment_time %in% c("T0","T24"))
pAlpha<-plot_richness(data_sub, measures=c("Shannon", "Observed"), x="Treatment_time", color="Treatment_time")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Treatment_time, New_data$Sample),]
pAlpha$data<-New_data
pAlpha + 
  geom_point(data=pAlpha$data, 
             aes(x=Treatment_time, y=value, color=NULL), alpha=0.1) + 
  scale_color_manual(values = c("T24"="coral","T0"="chartreuse")) +
  #geom_point(aes(color=Treatment_time), size=1, alpha=0.8) +
  theme_classic(base_size = 9.5) +
  geom_line(aes(group = pAlpha$data$Sample_code),size=0.05) +
  labs(x="Treatment_time", 
       title="A) \n\nAlpha diversity between T0 and T24 patients") +
  guides(fill=FALSE, color=FALSE) + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Treatment_time), label="p.format", method = "wilcox.test", paired = T, label.x= 1.42, size=2.8, label.y.npc = "top", vjust=-0.4, hjust=0.4)
ggsave(file="T0_vs_T24/Alfa diversity_treatment_saliva.png", width = 5,height =4.3, dpi=300)


######## Healthy vs HIV
suppressWarnings(rm(pAlpha, data_sub))
data_sub<-subset_samples(data, Condition %in% c("Healthy","HIV"))
pAlpha<-plot_richness(data_sub, measures=c("Shannon", "Observed"), x="Condition", color="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
New_data<-New_data[order(New_data$Condition, New_data$Sample),]
pAlpha$data<-New_data
pAlpha +
  geom_point(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha= 0.1) +
  theme_classic(base_size = 9.5) + 
  scale_color_manual(values=c("Healthy"="deepskyblue","HIV"="coral")) +
  geom_boxplot(size= 0.5, alpha = 0.1, color="black") +
  labs(x="Condition",
       title="A) \n\nAlpha diversity between healthy and HIV groups") +
  guides(fill=FALSE, color=FALSE) + 
  theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=9),
        title = element_text( size=9)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", paired = F, label.x= 1.41, size=2.7, label.y.npc = "top", vjust=-0.4, hjust=0.4)
ggsave(file="Healthy_vs_HIV/Alfa diversity_healthy_HIV_saliva.png", width = 5,height =4.3, dpi=300)


# the difference in richness is caused by really low abundant ASV (noise)?
data_sub.prop<-transform_sample_counts(data_sub, function(x) (x/sum(x))*100)
data_sub_healthy<-subset_samples(data_sub.prop, Condition=="Healthy")
data_sub_healthy<-filter_taxa(data_sub_healthy, function(x) mean(x)<1, prune = T)
data_sub_healthy<-filter_taxa(data_sub_healthy, function(x) min(x)>0, prune = T) # there are also saliva ASV codes
head(otu_table(data_sub_healthy))
N_low_healthy<-length(unique(tax_table(data_sub_healthy)[,"Genus"]))
N_low_healthy
data_sub_HIV<-subset_samples(data_sub.prop, Condition=="HIV")
data_sub_HIV<-filter_taxa(data_sub_HIV, function(x) mean(x)<1, prune = T)
data_sub_HIV<-filter_taxa(data_sub_HIV, function(x) min(x)>0, prune = T) # there are also saliva ASV codes
head(otu_table(data_sub_HIV))
tail(otu_table(data_sub_HIV))
N_low_HIV<-length(unique(tax_table(data_sub_HIV)[,"Genus"]))
N_low_HIV

con<-file("Healthy_vs_HIV/Noisy_ASV_check.txt")
sink(con, append = T)
cat("The difference in ASV richness (see alpha diversity) is caused by a different background noise between two different sequencing batches???\n\n",fill=T)
cat("ASV with mean lower than 0.1% in HIV saliva samples:",N_low_HIV,"\n",fill=T)
cat("ASV with mean lower than 0.1% in Healthy saliva samples:",N_low_healthy,"\n\n",fill=T)
cat("--> According to this arbitrary threshold and considering the max richness >500, the 'noisy' seems to be comparable between the two batches (both in saliva and in stool analysis!)",fill=T)
sink()
close(con)

rm(N_low_HIV,N_low_healthy, con, data_sub_healthy, data_sub_HIV)


### Responders
suppressWarnings(rm(pAlpha, data_sub))
data_sub<-data.responders
pAlpha<-plot_richness(data_sub, measures=c("Shannon", "Observed"), x="Responders", color="Responders")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data$Responders<-gsub("1","Yes",New_data$Responders)
New_data$Responders<-gsub("0","No",New_data$Responders)
pAlpha$data<-New_data
pAlpha + 
  geom_point(data=pAlpha$data, aes(x=Responders, y=value, color=NULL), alpha=0.1) +
  geom_boxplot(size= 0.5, alpha = 0.1, color="black") +
  theme_classic(base_size = 9.5) + 
  labs(x="Responders",
       title="A) \n\nAlpha diversity between responders and not responders") +
  scale_color_manual(values=c("No"="chartreuse4","Yes"="chartreuse")) +
  guides(fill=FALSE, color=FALSE) + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=0.9, size=9),
        title = element_text(size=8.9)) +
  stat_compare_means(aes(group = Responders), label="p.format", method = "wilcox.test", paired = F, label.x= 1.425, size=2.8,
                     label.y.npc = "top", vjust=-0.4, hjust=0.4)
ggsave(file="Responders_CD4_Analysis/Alfa diversity_Responders_saliva.png", width = 5,height =4.3, dpi=300)

suppressWarnings(rm(Obser_value,factor, data_sub, pAlpha, obs, New_data, ev, H, alphadt))


################ BETA DIVERSITY Hellinger (SALIVA) #######################

dir.create("T0_vs_T24/Beta_div")
dir.create("Healthy_vs_HIV/Beta_div")
dir.create("Responders_CD4_Analysis/Beta_div")

############ Treatment
suppressWarnings(rm(data_sub, data_sub.prop))
{data_sub<-subset_samples(data, Treatment_time %in% c("T0","T24"))
  data_sub.genus<- tax_glom(data_sub, taxrank = "Genus", NArm = F)
  data_sub.phy<- tax_glom(data_sub, taxrank = "Phylum", NArm = F)
  data_sub.phy.prop <- transform_sample_counts(data_sub.phy, function(otu) otu/sum(otu)*100)
  data_sub.genus.prop <- transform_sample_counts(data_sub.genus, function(otu) otu/sum(otu)*100)
  data_sub.prop <- transform_sample_counts(data_sub, function(otu) otu/sum(otu)*100)
}

{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phy.prop))
}

#### PERMANOVA
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Treatment_time, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV$aov.tab
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Treatment_time, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_g$aov.tab
perm_g_euclidean<-perm_g$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Treatment_time, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_p$aov.tab

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
write.csv2(beta, file="T0_vs_T24/Beta_diversity_permanova Hellinger.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,as(sample_data(data_sub),"data.frame")$Treatment_time)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="T0_vs_T24/General Beta_dispersion_permanova Hellinger.csv",quote=F,row.names = T)

######### PCoA Treatment
data.sqrt_prop<-transform_sample_counts(data_sub.prop, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
plot_ordination(data.sqrt_prop, ordBC, color = "Treatment_time") +
  scale_color_manual(values = c("T0" = "coral", "T24"="green")) +
  geom_line(aes(group=sample_data(data.sqrt_prop)$Sample_code), col="black", size=0.2) +
  geom_point(size=3) + theme_classic(base_size = 14) + 
  labs(title="PCoA with Hellinger distance  (saliva samples)",
       color="Treatment_time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2.5, show.legend = FALSE)
ggsave(file="T0_vs_T24/Beta_div/PCoA Beta diversity Hellinger.png", width = 7, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Treatment_time") +
  scale_color_manual(values = c("T0" = "coral", "T24"="green")) +
  geom_line(aes(group=sample_data(data.sqrt_prop)$Sample_code), col="black", size=0.09) +
  geom_point(size=3.15, alpha= 0.4) +
  theme_classic(base_size = 9.5) + 
  labs(title="B) \n\nPCoA with Hellinger distance on ASVs",
       color="Treatment time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), 
       #subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1]))
       subtitle = paste("PERMANOVA Pr(>F) =", 0.9789)) # to perfectly match the one already written in the draft
ggsave(file="T0_vs_T24/Beta_div/PCoA Beta diversity Hellinger_no_names.png", width = 5, height = 4.4, dpi=300)


########## Healthy vs HIV
suppressWarnings(rm(data_sub, data_sub.prop))
{data_sub<-subset_samples(data, Condition %in% c("HIV","Healthy"))
  data_sub.genus<- tax_glom(data_sub, taxrank = "Genus", NArm = F)
  data_sub.phy<- tax_glom(data_sub, taxrank = "Phylum", NArm = F)
  data_sub.phy.prop <- transform_sample_counts(data_sub.phy, function(otu) otu/sum(otu)*100)
  data_sub.genus.prop <- transform_sample_counts(data_sub.genus, function(otu) otu/sum(otu)*100)
  data_sub.prop <- transform_sample_counts(data_sub, function(otu) otu/sum(otu)*100)
}

{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phy.prop))
}

### PERMANOVA
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV$aov.tab
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_g$aov.tab
perm_g_euclidean<-perm_g$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Condition, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_p$aov.tab

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
write.csv2(beta, file="Healthy_vs_HIV/Beta_div/Beta_diversity_permanova Hellinger.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,as(sample_data(data_sub),"data.frame")$Condition)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Healthy_vs_HIV/Beta_div/General Beta_dispersion_permanova Hellinger.csv",quote=F,row.names = T)

####### PCoA Condition
data.sqrt_prop<-transform_sample_counts(data_sub.prop, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values = c("Healthy" = "deepskyblue", "HIV"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) +stat_ellipse() +
  labs(title="PCoA with Hellinger distance on ASVs",
       color="Condition", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2.5, show.legend = FALSE)
ggsave(file="Healthy_vs_HIV/Beta_div/PCoA Beta diversity Hellinger.png", width = 7, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values = c("Healthy" = "deepskyblue", "HIV"="coral")) +
  geom_point(size=3.1, alpha= 0.4) +
  theme_classic(base_size = 9.5) +stat_ellipse(size=0.10) +
  labs(title="B) \n\nPCoA with Hellinger distance on ASVs",
       color="Condition", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       #subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1]))
       subtitle = paste("PERMANOVA Pr(>F) =", 0.0105)) # to perfectly match the one already written in the draft
ggsave(file="Healthy_vs_HIV/Beta_div/PCoA Beta diversity Hellinger_no names.png", width = 5, height = 4.4, dpi=300)


######## Responders yes vs no
suppressWarnings(rm(data_sub, data_sub.prop))
{data_sub<-data.responders
  data_sub.genus<- tax_glom(data_sub, taxrank = "Genus", NArm = F)
  data_sub.phy<- tax_glom(data_sub, taxrank = "Phylum", NArm = F)
  data_sub.phy.prop <- transform_sample_counts(data_sub.phy, function(otu) otu/sum(otu)*100)
  data_sub.genus.prop <- transform_sample_counts(data_sub.genus, function(otu) otu/sum(otu)*100)
  data_sub.prop <- transform_sample_counts(data_sub, function(otu) otu/sum(otu)*100)
}

{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phy.prop))
}

#### PERMANOVA
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Responders, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV$aov.tab
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Responders, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_g$aov.tab
perm_g_euclidean<-perm_g$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Responders, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_p$aov.tab

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
write.csv2(beta, file="Responders_CD4_Analysis/Beta_div/Beta_diversity_permanova Hellinger.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,as(sample_data(data_sub),"data.frame")$Responders)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Responders_CD4_Analysis/Beta_div/General Beta_dispersion_permanova Hellinger.csv",quote=F,row.names = T)

### PCoA Responders
data.temp <- transform_sample_counts(data_sub.prop, function(x) x/sum(x))
data.sqrt_prop<-transform_sample_counts(data.temp, sqrt) # square root of proportion
sample_data(data.sqrt_prop)$Responders<-gsub("1","Responding",sample_data(data.sqrt_prop)$Responders)
sample_data(data.sqrt_prop)$Responders<-gsub("0","Not responding",sample_data(data.sqrt_prop)$Responders)
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
plot_ordination(data.sqrt_prop, ordBC, color = "Responders") +
  scale_color_manual(values = c("Not responding" = "chartreuse4", "Responding"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) +stat_ellipse() +
  labs(title="PCoA with Hellinger distance  (saliva samples)",
       color="Responders", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2.5, show.legend = FALSE)
ggsave(file="Responders_CD4_Analysis/Beta_div/PCoA Beta diversity Hellinger.png", width = 7, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Responders") +
  scale_color_manual(values = c("Not responding" = "chartreuse4", "Responding"="chartreuse")) +
  geom_point(size=3.15, alpha =0.4) + theme_classic(base_size = 9.5) +stat_ellipse(size= 0.10) +
  labs(title="B) \n\nPCoA with Hellinger distance on ASVs",
       color="Responders", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), 
       # subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1]))
       subtitle = paste("PERMANOVA Pr(>F) =",0.6454)) # to perfectly match the one already written in the draft
ggsave(file="Responders_CD4_Analysis/Beta_div/PCoA Beta diversity Hellinger_no names.png",
       width = 5, height = 4.4, dpi=300)


################# DA WITH DESEQ2 _ Treatment (SALIVA) ##################

dir.create("T0_vs_T24/DA_DESeq2")

suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Treatment_time %in% c("T0","T24"))
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
# system(" echo 'In case of pair analysis, the log2foldchange displayes the change from T24 to T0, see https://support.bioconductor.org/p/105981/ \n' > T0_vs_T24/DA_DESeq2/NB.txt ")  

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sample+Treatment_time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Treatment_time", "T24", "T0"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res<-res[res$baseMean > 100, ] # threshold arbitrary chosen based on low abundances in plots
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
    #write.csv2(r, file=paste0("T0_vs_T24/DA_DESeq2/DA_",t,"_ratio_T24_vs_T0.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="T0_vs_T24/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="T0_vs_T24/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

tabella_g<-Table_tot[Table_tot$Taxa=="Family",]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Treatment_time)

unique(tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Treatment_time)) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("T24"="coral","T0"="chartreuse")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = "Family")~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Treatment_time), size=1.8, alpha=0.9) +
  geom_point(aes(color=Treatment_time), size=3.8, alpha=0.5) +
  geom_line(aes(group=Sample_code), size=0.12) +
  theme(strip.text.x=element_text(size=13,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line")  ) + 
  theme(axis.text.y = element_text(size=12))+
  scale_x_discrete(labels=rep(unique(levels(tabella_g$Treatment_time))), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
# scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))
plot_1 +
  labs(title= "Differently abundant families in saliva between T0 and T24", y="Proportional Abundance %", fill="Treatment_time", x="")
ggsave(filename = "T0_vs_T24/DA_DESeq2/Plot_result_DeSeq2_Treatment_time.png", width = 8, height = 5, dpi=300)
dev.off()

# for further comparisons ...
Family.DESEQ2_T0_T24<-unique(tabella_g[tabella_g$Taxa=="Family","Bacteria"])
# NB: it is called "_g" (not f) but just because it's a typo error during the design of the script (the obj is correct!)

rm(plot_2, plot_1, head_plot)


############## DA WITH DESeq2_ Conditions (SALIVA) #############

dir.create("Healthy_vs_HIV/DA_DESeq2")

suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition %in% c("Healthy","HIV"))
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
  res<-results(DE, contrast= c("Condition", "HIV", "Healthy"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res<-res[res$baseMean > 100, ] # arbitrary threshold to avoid the most noisy result
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
    #write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_HIV_vs_Healthy.csv"), row.names = F, quote=F, na = "")
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

#View(Res_tot)
write.csv(Res_tot, file="Healthy_vs_HIV/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Healthy_vs_HIV/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Condition<-factor(Table_tot$Condition, levels=c("Healthy","HIV") )
Table_tot<-Table_tot[Table_tot$Taxa=="Family",] # the other results are redundants

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, color=Condition, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.6, color="black", alpha=0, size= 0.4) + 
  geom_point(position=position_dodge(width=0.6), size=1, alpha= 0.9) +
  geom_point(position=position_dodge(width=0.6), size=2.5, alpha= 0.6) +
  theme_classic(base_size = 11) +
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  scale_color_manual(values=c("Healthy"="deepskyblue","HIV"="coral")) +
  guides( color=guide_legend(nrow=1), fill="none" ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=13),
        panel.grid.major.x = element_line(linewidth = 0.5),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=11), 
        axis.text.y = element_text(size=10),
        plot.title= element_text(size=12),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1, 4, seq(2,max(Table_tot$Abundance),4))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", 
       fill="Condition", x="")
ggsave(filename = "Healthy_vs_HIV/DA_DESeq2/DA_Condition_every_result.png", width = 6.5, height = 4.5, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Healthy_vs_HIV/DA_DESeq2/NB.txt ")

# for comparisons ...
Family.DESEQ2_Healthy_HIV<-unique(Table_tot[Table_tot$Taxa=="Family","Bacteria"]) # no genera!



############## DA WITH DESeq2_ Responders (SALIVA) #############

dir.create("Responders_CD4_Analysis/DA_DESeq2")

suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data.responders) > 10, data.responders) 
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Responders)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Responders", "1", "0"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res<-res[res$baseMean > 100, ] # arbitrary threshold to avoid the most noisy result
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
    #write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_1_vs_0.csv"), row.names = F, quote=F, na = "")
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

#View(Res_tot)
write.csv(Res_tot, file="Responders_CD4_Analysis/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Responders_CD4_Analysis/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Responders<-gsub("0","Not responding",Table_tot$Responders)
Table_tot$Responders<-gsub("1","Responding",Table_tot$Responders)
Table_tot$Responders<-factor(Table_tot$Responders, levels=c("Not responding","Responding") )
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, color=Responders, fill=Responders)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.6, color="black", alpha=0, size= 0.35) + 
  geom_point(position=position_dodge(width=0.6), size=1, alpha= 0.9) +
  geom_point(position=position_dodge(width=0.6), size=2.5, alpha= 0.6) +
  theme_classic(base_size = 11) +
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  scale_color_manual(values=c("Not responding"="chartreuse4","Responding"="chartreuse2")) +
  guides( color=guide_legend(nrow=1), fill="none" ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=13),
        panel.grid.major.x = element_line(linewidth = 0.5),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=11), 
        axis.text.y = element_text(size=9.8),
        plot.title= element_text(size=12),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1, 4, seq(2,max(Table_tot$Abundance),4))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", 
       fill="Immunological responders   ", x="")
ggsave(filename = "Responders_CD4_Analysis/DA_DESeq2/DA_Responders_every_result.png", width = 6.5, height = 4.5, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Responders_CD4_Analysis/DA_DESeq2/NB.txt ")

# for comparisons ...
Phyla.DESEQ2_responders<-unique(Table_tot[Table_tot$Taxa=="Phylum","Bacteria"])




######################## SAVING DESEQ2 RESULTS ABUNDANCES (SALIVA) ######################

suppressWarnings(rm(temp))
temp<-otu_table(subset_taxa(data.responders.phy.prop, Phylum %in% Phyla.DESEQ2_responders))
row.names(temp)<-as.character(tax_table(subset_taxa(data.phy.prop, Phylum %in% Phyla.DESEQ2_responders))[, "Phylum"])
if(identical(colnames(temp), sample_data(data.responders)$Sample_ID)){
  colnames(temp)<-sample_data(data.responders)$Sample
}
DA_Responders_Saliva<-temp

suppressWarnings(rm(temp))
temp<-otu_table(subset_taxa(data.fam.prop, Family %in% Family.DESEQ2_T0_T24))
row.names(temp)<-as.character(tax_table(subset_taxa(data.fam.prop, Family %in% Family.DESEQ2_T0_T24))[, "Family"])
if(identical(colnames(temp),sample_data(data)$Sample_ID)){
  colnames(temp)<-sample_data(data)$Sample
}
DA_Treatment_Saliva<-temp


################ RANDOM FOREST (HIV vs HEALTHY SALIVA) ###################

library("randomForest")
dir.create("Healthy_vs_HIV/Random_Forest")

# on family because the DA results are on families only --> searching for common results
data.target<-subset_samples(data.fam.prop, Treatment_time!="T24") # Proportional Abundance %s, e.g. see the script used in PMID: 31293616
tax_table(data.target)<-as.matrix(Taxa.fam.update)


### filtering noises
# according to relative abundance, e.g. see PMID: 24193493
filter_to_apply_RF<-function(x) mean(x)>0.1
data.target.cleaned <- filter_taxa(data.target, filter_to_apply_RF, TRUE)
# according to prevalence
min_sample<-round(4*length(sample_names(data.target.cleaned))/100, digits = 0) # at least 4% of samples, see PMID: 35140735
if(min_sample == 0 ){ min_sample<-1 }
who<-as.data.frame(otu_table(data.target.cleaned))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.5, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>min_sample,] # not if found in less than X "points"  (--> those one to be removed)
who<-as.vector(tax_table(data.target.cleaned)[row.names(who),"Family"]) # to translate in ASV code (of genus collapsed object)
data.target.cleaned<-subset_taxa(data.target.cleaned, ! Family %in% who)


predictors <- t(otu_table(data.target.cleaned)) # NB: samples as rows

### Scale the data should not be necessary, see https://stackoverflow.com/questions/8961586/do-i-need-to-normalize-or-scale-data-for-randomforest-r-package
# predictors <- as.data.frame(scale(predictors)) # NB: scale acts on columns --> scaled the abundances of the variable through the samples

y <- as.factor(sample_data(data.target.cleaned)$Condition)
rf_data <- cbind.data.frame(y, predictors)

# NB: error if the colnames starts with digit, then...
colnames(rf_data) <- make.names(colnames(rf_data)) #(... if start with a number, adds an X at the start of the name)

set.seed(1)
pdf("Healthy_vs_HIV/Random_Forest/Tuning_RF___BEST_NUMBER_OF_VARIABLES_FOR_SPLIT.pdf")
tuneRF(predictors, y, ntreeTry = 1000) #(to find the optimal number of variables at each split)
dev.off()
best <- 3 # see tuneRF results

set.seed(1)
RF <- randomForest(y ~ . , data = rf_data, ntree = 1000, mtry= best)
RF

Selected<-importance(RF) # mean decrease giny impurity for each taxa
row.names(Selected)<-gsub("^X", "", row.names(Selected)) # removing the added X
row.names(Selected)<-as.character(tax_table(data.target.cleaned)[ row.names(Selected), "Family"]) # translating from code to tax names
Selected<-Selected[order(Selected, decreasing = T), ] # mean DECREASE giny impurity --> higher is better
t_names<-factor(names(Selected), levels = unique(names(Selected))) # NB: now Selected is a "labelled vector!

ggplot(mapping=aes(x=t_names, y=as.numeric(Selected))) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 35, size = 12, hjust = 1, vjust = 1)) +
  labs(x=paste("\nTaxa remaining in the feature table after the filters"),
       y="mean decrease of Gini impurity", 
       title = paste("Importance of the taxa according to Random Forest (HIV vs HEALTHY) \n computed using only families with average abundance over 0.1 \n and found at least in more than",min_sample,"sample"),
       caption="\n dotted lines are arbitrarly plot at 50% of max Mean Decrease Gini value") +
  geom_hline(yintercept = max(as.numeric(Selected)) /2 , colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"),
        axis.text.x = element_text(
          colour= ifelse(t_names %in% Family.DESEQ2_Healthy_HIV, "red","black")
          )
        )
ggsave(filename = "Healthy_vs_HIV/Random_Forest/Importances_of_selected_Taxa.png", width = 12, height = 8, dpi=300)


#### EXPORTING THE SETTING AND THE RESULTS OF THIS RANDOM FOREST
suppressWarnings(rm(con))
con<-file("Healthy_vs_HIV/Random_Forest/Settings_and_results_RF.txt")
sink(con)
RF
cat("\n\n\nPrediction of the model: \n")
RF$predicted
cat("\n\n\nPrediction errors of the model: \n")
cat(RF$confusion[2,1]+RF$confusion[1,2],"errors and",RF$confusion[1,1]+RF$confusion[2,2],"samples predicted correctly")
cat("\n\n\nNumber of variable for split:",best,"(selected through tuneRF function, see the related plot)")
cat("\n\n\nFilter applied to the dataset:\n")
print(filter_to_apply_RF)
cat("Prevalence: found at least in more than",min_sample,"samples")
cat("\n\n\n\nVariables importance according to this forest: \n")
Selected
sink()
close(con)

pdf("Healthy_vs_HIV/Random_Forest/RF___Errors_associated_to_tree_number.pdf")
plot(RF)
dev.off()

suppressWarnings(rm(data.target, data.target.cleaned, filter_to_apply_RF, con, t_names))

# NB: the package has to be removed from the environment after being used, because it over-writes some ggplot2 function causing errors!
detach("package:randomForest") 



########### CORRELATIONS HIV DATA: DA TAXA SALIVA vs DA IL and FFA (Saliva) #############

dir.create("DA_taxa_vs_DA_SCFA_and_IL_along_T0_T24")

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data.fam.prop, Condition!="Healthy")
sample_names(data.temp)<-paste(sample_data(data.temp)$Sample_code, sample_data(data.temp)$Treatment_time, sep="_")
prune.dat_target <- subset_taxa(data.temp, Family %in% Family.DESEQ2_T0_T24)
target <-as.data.frame(otu_table(prune.dat_target))
target_tax<-as.data.frame(tax_table(prune.dat_target))
identical(row.names(target_tax),row.names(target))
row.names(target)<-target_tax$Family
rm(target_tax)
target<-as.data.frame(t(target))

SCFA <- read.delim2("../../SCFA.csv")
IL <- read.delim2("../../Interleuch.csv")

# same order with the phyloseq data
row.names(SCFA)<-paste(SCFA$Sample_Code, SCFA$Treatment_time, sep="_")
SCFA<-SCFA[row.names(target), colnames(SCFA) %in% c("Butyric","Propionic")]
row.names(IL)<-paste(IL$Sample_Code, IL$Treatment_time, sep="_")
IL<-IL[row.names(target), colnames(IL) %in% c("IL8","IL10")]

# removing every column with all zero
IL<-IL[,!colSums(IL)==0]
SCFA<-SCFA[,!colSums(SCFA)==0]

# Computing EVERY correlation to start...
x<-cbind.data.frame(target,IL,SCFA)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}

# ... focusing on Bacteria vs SCFA
corr<-subset(correlation, Var1 %in% colnames(target))
corr<-subset(corr, Var2 %in% colnames(SCFA))
colnames(corr)<-c("DA_Family","SCFA","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$DA_Family, y = corr$SCFA, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "DA SCFA", x= "DA Family in saliva", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="DA_taxa_vs_DA_SCFA_and_IL_along_T0_T24/Family_vs_SCFA.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="DA_taxa_vs_DA_SCFA_and_IL_along_T0_T24/Family_vs_SCFA.csv", quote = F, na="", row.names = F)


# ... focusing on Bacteria vs IL
corr<-subset(correlation, Var1 %in% colnames(target))
corr<-subset(corr, Var2 %in% colnames(IL))
colnames(corr)<-c("DA_Family","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$DA_Family, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "DA IL", x= "DA Family in saliva", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="DA_taxa_vs_DA_SCFA_and_IL_along_T0_T24/Family_vs_IL.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="DA_taxa_vs_DA_SCFA_and_IL_along_T0_T24/Family_vs_IL.csv", quote = F, na="", row.names = F)


########### CORRELATIONS HIV DATA (T0) WITH IL AND SCFA IN T0 #############

dir.create("T0_vs_T24/Spearm correlations")
dir.create("T0_vs_T24/Spearm correlations/Samples_at_time_T0")

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data.genus.prop, Treatment_time=="T0") # only at T0
sample_names(data.temp)<-sample_data(data.temp)$Sample_code
top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.temp)
Top5<-as.data.frame(otu_table(prune.dat_top))
Top5_tax<-as.data.frame(tax_table(prune.dat_top))
identical(row.names(Top5_tax),row.names(Top5))
row.names(Top5)<-Top5_tax$Genus
rm(Top5_tax)
Top5<-as.data.frame(t(Top5))

SCFA <- read.delim2("../../SCFA.csv")
IL <- read.delim2("../../Interleuch.csv")

SCFA_T0<-subset(SCFA, Treatment_time=="T0")
row.names(SCFA_T0)<-SCFA_T0$Sample_Code
SCFA_T0<-SCFA_T0[row.names(Top5), colnames(SCFA_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]
IL_T0<-subset(IL, Treatment_time=="T0")
row.names(IL_T0)<-IL_T0$Sample_Code
IL_T0<-IL_T0[row.names(Top5),colnames(IL_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]
# removing every column with all zero
IL_T0<-IL_T0[,!colSums(IL_T0)==0]
SCFA_T0<-SCFA_T0[,!colSums(SCFA_T0)==0]

# Computing EVERY correlation to start...
x<-cbind.data.frame(Top5,IL_T0,SCFA_T0)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}

# ... focusing on Bacteria vs SCFA
corr<-subset(correlation, Var1 %in% colnames(Top5))
corr<-subset(corr, Var2 %in% colnames(SCFA_T0))
colnames(corr)<-c("Top_Genera","SCFA","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
corr$SCFA<-gsub("_"," ",corr$SCFA)
corr$SCFA<-gsub("X","",corr$SCFA)
ggplot(corr, aes(x = corr$Top_Genera, y = corr$SCFA, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle = -22, hjust = 0, size= 11),
        axis.text.y=element_text(size= 10.5)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =10, vjust=0.7) +
  labs(title = "A) \n", y= "SCFA", x= "5 most abundant genera", 
       caption= "\n adjusted p-value (Holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=16)) +
  theme(legend.text = element_text(size=11), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(0.9, 'cm'), legend.title = element_text(size=13))
ggsave(file="T0_vs_T24/Spearm correlations/Samples_at_time_T0/Top5genera_vs_SCFA_T0.png", dpi=300, width = 6, height = 6)
write.csv2(corr, file="T0_vs_T24/Spearm correlations/Samples_at_time_T0/Spearman Top5gen vs SCFA T0.csv", quote = F, na="", row.names = F)

# ... focusing on Bacteria vs IL
corr<-subset(correlation, Var1 %in% colnames(Top5))
corr<-subset(corr, Var2 %in% colnames(IL_T0))
colnames(corr)<-c("Top_Genera","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$Top_Genera, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "5 most abundant genera in saliva", 
       caption= "\n adjusted p-value (Holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="T0_vs_T24/Spearm correlations/Samples_at_time_T0/Top5genera_vs_IL_T0.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="T0_vs_T24/Spearm correlations/Samples_at_time_T0/Spearman Top5gen vs IL T0.csv", quote = F, na="", row.names = F)

# ... focusing on SCFA vs IL
corr<-subset(correlation, Var1 %in% colnames(SCFA_T0))
corr<-subset(corr, Var2 %in% colnames(IL_T0))
colnames(corr)<-c("SCFA","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$SCFA, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "SCFA", 
       caption= "\n adjusted p-value (Holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="T0_vs_T24/Spearm correlations/Samples_at_time_T0/SCFA_vs_IL_T0.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="T0_vs_T24/Spearm correlations/Samples_at_time_T0/Spearman SCFA vs IL T0.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr, Top5, data.temp)


########### CORRELATIONS HIV DATA (T24) WITH IL AND SCFA IN T24 #############

dir.create("T0_vs_T24/Spearm correlations")
dir.create("T0_vs_T24/Spearm correlations/Samples_at_time_T24")

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data.genus.prop, Treatment_time=="T24") # only at T24
sample_names(data.temp)<-sample_data(data.temp)$Sample_code
top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.temp)
Top5<-as.data.frame(otu_table(prune.dat_top))
Top5_tax<-as.data.frame(tax_table(prune.dat_top))
identical(row.names(Top5_tax),row.names(Top5))
row.names(Top5)<-Top5_tax$Genus
rm(Top5_tax)
Top5<-as.data.frame(t(Top5))

SCFA <- read.delim2("../../SCFA.csv")
IL <- read.delim2("../../Interleuch.csv")

SCFA_T24<-subset(SCFA, Treatment_time=="T24")
row.names(SCFA_T24)<-SCFA_T24$Sample_Code
SCFA_T24<-SCFA_T24[row.names(Top5), colnames(SCFA_T24)!=c("Sample_ID","Sample_Code","Treatment_time")]
IL_T24<-subset(IL, Treatment_time=="T24")
row.names(IL_T24)<-IL_T24$Sample_Code
IL_T24<-IL_T24[row.names(Top5),colnames(IL_T24)!=c("Sample_ID","Sample_Code","Treatment_time")]
# removing every column with all zero
IL_T24<-IL_T24[,!colSums(IL_T24)==0]
SCFA_T24<-SCFA_T24[,!colSums(SCFA_T24)==0]

# Computing EVERY correlation to start...
x<-cbind.data.frame(Top5,IL_T24,SCFA_T24)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}

# ... focusing on Bacteria vs SCFA
corr<-subset(correlation, Var1 %in% colnames(Top5))
corr<-subset(corr, Var2 %in% colnames(SCFA_T24))
colnames(corr)<-c("Top_Genera","SCFA","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
corr$SCFA<-gsub("_"," ",corr$SCFA)
corr$SCFA<-gsub("X","",corr$SCFA)
ggplot(corr, aes(x = corr$Top_Genera, y = corr$SCFA, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle = -22, hjust = 0, size= 11),
        axis.text.y=element_text(size= 10.5)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =10, vjust=0.7) +
  labs(title = "B) \n", y= "SCFA", x= "5 most abundant genera", 
       caption= "\n adjusted p-value (Holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=16)) +
  theme(legend.text = element_text(size=11), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(0.9, 'cm'), legend.title = element_text(size=13))
ggsave(file="T0_vs_T24/Spearm correlations/Samples_at_time_T24/Top5genera_vs_SCFA_T24.png", dpi=300, width = 6, height = 6)
write.csv2(corr, file="T0_vs_T24/Spearm correlations/Samples_at_time_T24/Spearman Top5gen vs SCFA T24.csv", quote = F, na="", row.names = F)

# ... focusing on Bacteria vs IL
corr<-subset(correlation, Var1 %in% colnames(Top5))
corr<-subset(corr, Var2 %in% colnames(IL_T24))
colnames(corr)<-c("Top_Genera","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$Top_Genera, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "5 most abundant genera in saliva", 
       caption= "\n adjusted p-value (Holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="T0_vs_T24/Spearm correlations/Samples_at_time_T24/Top5genera_vs_IL_T24.png", dpi=300, width = 6, height = 6)
write.csv2(corr, file="T0_vs_T24/Spearm correlations/Samples_at_time_T24/Spearman Top5gen vs IL T24.csv", quote = F, na="", row.names = F)

# ... focusing on SCFA vs IL
corr<-subset(correlation, Var1 %in% colnames(SCFA_T24))
corr<-subset(corr, Var2 %in% colnames(IL_T24))
colnames(corr)<-c("SCFA","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$SCFA, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "SCFA", 
       caption= "\n adjusted p-value (Holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="T0_vs_T24/Spearm correlations/Samples_at_time_T24/SCFA_vs_IL_T24.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="T0_vs_T24/Spearm correlations/Samples_at_time_T24/Spearman SCFA vs IL T24.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr, Top5, data.temp)


############################ IL AND SCFA BETWEEN RESPONDERS #######################

{ if( grepl("_analysis",getwd()) ){setwd("..")}
  if( ! grepl("Results",getwd()) ){setwd("Results")} # it sets wd to "Results"
}

dir.create("Correlations_in_Responders_or_Not_Responders")

SCFA <- read.delim2("../SCFA.csv")
IL <- read.delim2("../Interleuch.csv")

### at T0

Responders<-as.data.frame(sample_data(data))
Responders<-Responders[Responders$Responders=="1", "Sample_code"]
# View(sample_data(data))
SCFA_T0<-subset(SCFA, Treatment_time=="T0")
SCFA_T0<-subset(SCFA_T0, Sample_Code %in% Responders$Sample_code)
row.names(SCFA_T0)<-SCFA_T0$Sample_Code
SCFA_T0<-SCFA_T0[Responders$Sample_code, colnames(SCFA_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

IL_T0<-subset(IL, Treatment_time=="T0")
IL_T0<-subset(IL_T0, Sample_Code %in% Responders$Sample_code)
row.names(IL_T0)<-IL_T0$Sample_Code
IL_T0<-IL_T0[Responders$Sample_code, colnames(IL_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

# removing every column with all zero
IL_T0<-IL_T0[,!colSums(IL_T0)==0]
SCFA_T0<-SCFA_T0[,!colSums(SCFA_T0)==0]


# Computing EVERY correlation to start...
x<-cbind.data.frame(IL_T0,SCFA_T0)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}
corr<-subset(correlation, Var1 %in% colnames(SCFA_T0))
corr<-subset(corr, Var2 %in% colnames(IL_T0))
colnames(corr)<-c("SCFA","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$SCFA, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "SCFA", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T0_in_RESPONDERS.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T0_in_RESPONDERS.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr)

### at T24

Responders<-as.data.frame(sample_data(data))
Responders<-Responders[Responders$Responders=="1", "Sample_code"]
# View(sample_data(data))
SCFA_T0<-subset(SCFA, Treatment_time=="T24")
SCFA_T0<-subset(SCFA_T0, Sample_Code %in% Responders$Sample_code)
row.names(SCFA_T0)<-SCFA_T0$Sample_Code
SCFA_T0<-SCFA_T0[Responders$Sample_code, colnames(SCFA_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

IL_T0<-subset(IL, Treatment_time=="T24")
IL_T0<-subset(IL_T0, Sample_Code %in% Responders$Sample_code)
row.names(IL_T0)<-IL_T0$Sample_Code
IL_T0<-IL_T0[Responders$Sample_code, colnames(IL_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

# removing every column with all zero
IL_T0<-IL_T0[,!colSums(IL_T0)==0]
SCFA_T0<-SCFA_T0[,!colSums(SCFA_T0)==0]


# Computing EVERY correlation to start...
x<-cbind.data.frame(IL_T0,SCFA_T0)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}
corr<-subset(correlation, Var1 %in% colnames(SCFA_T0))
corr<-subset(corr, Var2 %in% colnames(IL_T0))
colnames(corr)<-c("SCFA","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$SCFA, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "SCFA", 
       caption= "\n adjusted p-value (holm) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T24_in_RESPONDERS.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T24_in_RESPONDERS.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr)

#### EXTRA RIGUARDO QUELLA CORRELAZIONE 1
Sper<-cor.test(x[,"IL18"], x[,"Butyric"], method = "spearman")
Pears<- cor.test(x[,"IL18"], x[,"Butyric"], method = "pearson")
# x [ , c("IL18","Butyric")]
# rank(x$Butyric)
# rank(x$IL18)
png(filename = "Correlations_in_Responders_or_Not_Responders/About_that_perfect_correlation_at_T24_responders.png" )
par(mfrow=c(2,1))
plot(rank(x[,"Butyric"]),rank(x[,"IL18"]), xlab="Ranks Butyric", ylab="Ranks IL18")
abline(lm(rank(x[,"IL18"])~rank(x[,"Butyric"])),col="red")
title("With Spearman (p-value=0.001)")
plot(x[ ,"Butyric"],x[ ,"IL18"], xlab="Butyric", ylab="IL18")
abline(lm(x[ ,"IL18"]~x[ ,"Butyric"]),col="red")
title(paste0("With Pearson \n (p-value=",round(Pears$p.value,3)," and cor=",round(Pears$estimate,3) ,")"))
dev.off()


############################ IL AND SCFA BETWEEN NOT RESPONDERS #######################

{ if( grepl("_analysis",getwd()) ){setwd("..")}
  if( ! grepl("Results",getwd()) ){setwd("Results")} # it sets wd to "Results"
}

dir.create("Correlations_in_Responders_or_Not_Responders")

SCFA <- read.delim2("../SCFA.csv")
IL <- read.delim2("../Interleuch.csv")

### at T0

Responders<-as.data.frame(sample_data(data))
Responders<-Responders[Responders$Responders=="0", "Sample_code"]
# View(sample_data(data))
SCFA_T0<-subset(SCFA, Treatment_time=="T0")
SCFA_T0<-subset(SCFA_T0, Sample_Code %in% Responders$Sample_code)
row.names(SCFA_T0)<-SCFA_T0$Sample_Code
SCFA_T0<-SCFA_T0[Responders$Sample_code, colnames(SCFA_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

IL_T0<-subset(IL, Treatment_time=="T0")
IL_T0<-subset(IL_T0, Sample_Code %in% Responders$Sample_code)
row.names(IL_T0)<-IL_T0$Sample_Code
IL_T0<-IL_T0[Responders$Sample_code, colnames(IL_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

# removing every column with all zero
IL_T0<-IL_T0[,!colSums(IL_T0)==0]
SCFA_T0<-SCFA_T0[,!colSums(SCFA_T0)==0]


# Computing EVERY correlation to start...
x<-cbind.data.frame(IL_T0,SCFA_T0)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}
corr<-subset(correlation, Var1 %in% colnames(SCFA_T0))
corr<-subset(corr, Var2 %in% colnames(IL_T0))
colnames(corr)<-c("SCFA","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$SCFA, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "SCFA", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T0_in_NOT_RESPONDERS.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T0_in_NOT_RESPONDERS.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr)


### at T24

Responders<-as.data.frame(sample_data(data))
Responders<-Responders[Responders$Responders=="0", "Sample_code"]
# View(sample_data(data))
SCFA_T0<-subset(SCFA, Treatment_time=="T24")
SCFA_T0<-subset(SCFA_T0, Sample_Code %in% Responders$Sample_code)
row.names(SCFA_T0)<-SCFA_T0$Sample_Code
SCFA_T0<-SCFA_T0[Responders$Sample_code, colnames(SCFA_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

IL_T0<-subset(IL, Treatment_time=="T24")
IL_T0<-subset(IL_T0, Sample_Code %in% Responders$Sample_code)
row.names(IL_T0)<-IL_T0$Sample_Code
IL_T0<-IL_T0[Responders$Sample_code, colnames(IL_T0)!=c("Sample_ID","Sample_Code","Treatment_time")]

# removing every column with all zero
IL_T0<-IL_T0[,!colSums(IL_T0)==0]
SCFA_T0<-SCFA_T0[,!colSums(SCFA_T0)==0]


# Computing EVERY correlation to start...
x<-cbind.data.frame(IL_T0,SCFA_T0)
cor.test(x$IL18, x$Butyric)
#install.packages("Hmisc")
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Var1","Var2","Corr","pvalue")
}
corr<-subset(correlation, Var1 %in% colnames(SCFA_T0))
corr<-subset(corr, Var2 %in% colnames(IL_T0))
colnames(corr)<-c("SCFA","IL","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$SCFA, y = corr$IL, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "IL", x= "SCFA", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T24_in_NOT_RESPONDERS.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="Correlations_in_Responders_or_Not_Responders/SCFA_vs_IL_at_T24_in_NOT_RESPONDERS.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr)



################### ANALYSING LEFSE RESULTS SALIVA HEALTHY vs HIV ###########################

{ if( grepl("_analysis",getwd()) ){setwd("..")}
  if( ! grepl("Results",getwd()) ){setwd("Results")} # it sets wd to "Results"
}

dir.create("Saliva_bacteriota_analysis/Healthy_vs_HIV/PICRUST2_LEFSE_results")

suppressWarnings(rm(a, Significative_functions_LEFSE, Significative_functions_LEFSE_3))
a <- read.delim("../QIIME/PICRUST2_LEFSE/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("../QIIME/PICRUST2_LEFSE/Result_LEFSE_Healthy_vs_HIV_SALIVA.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="Saliva_bacteriota_analysis/Healthy_vs_HIV/PICRUST2_LEFSE_results/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)


# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Healthy"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Healthy"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), 
                                               fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), 
            hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Healthy",1,0), size=3) +
  scale_fill_manual(values=c("Healthy"="deepskyblue","HIV"="coral")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15))+
  theme(legend.position = "bottom")
ggsave(filename = "Saliva_bacteriota_analysis/Healthy_vs_HIV/PICRUST2_LEFSE_results/PICRUST2_LEFSE_Healthy_VS_HIV_METACYC_every_result.png", width = 10, height = 11, dpi = 300)


rm(Significative_functions_LEFSE, a)


################### ANALYSING LEFSE RESULTS SALIVA Not responders_vs_Responders ###########################

dir.create("Saliva_bacteriota_analysis/Responders_CD4_Analysis/PICRUST2_LEFSE_results")

suppressWarnings(rm(a, Significative_functions_LEFSE, Significative_functions_LEFSE_3))
a <- read.delim("../QIIME/PICRUST2_LEFSE/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("../QIIME/PICRUST2_LEFSE/Result_LEFSE_Responders_vs_Not_SALIVA.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="Saliva_bacteriota_analysis/Responders_CD4_Analysis/PICRUST2_LEFSE_results/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)


# modifing names for the plot
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Not responders"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Not responders"]

# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), 
                                               fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), 
            hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Not responders",1,0), size=3.2) +
  scale_fill_manual(values=c("Not responders"="chartreuse4","Responders"="chartreuse")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15))+
  theme(legend.position = "bottom")
ggsave(filename = "Saliva_bacteriota_analysis/Responders_CD4_Analysis/PICRUST2_LEFSE_results/PICRUST2_LEFSE_Responders_CD4_Analysis_METACYC_every_result.png", width = 10.5, height = 7, dpi = 300)


rm(Significative_functions_LEFSE, a)



###################### PCoA ON OTHER FACTORS (SALIVA) ################################

dir.create("Other_factors")

######## Oral infection
suppressWarnings(rm(data_sub, data_sub.prop))
{data_sub<-subset_samples(data, Treatment_time== "T0" & Condition=="HIV")
  data_sub.genus<- tax_glom(data_sub, taxrank = "Genus", NArm = F)
  data_sub.phy<- tax_glom(data_sub, taxrank = "Phylum", NArm = F)
  data_sub.phy.prop <- transform_sample_counts(data_sub.phy, function(otu) otu/sum(otu)*100)
  data_sub.genus.prop <- transform_sample_counts(data_sub.genus, function(otu) otu/sum(otu)*100)
  data_sub.prop <- transform_sample_counts(data_sub, function(otu) otu/sum(otu)*100)
}

{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phy.prop))
}

#### PERMANOVA
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Oral_infections, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV$aov.tab
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Oral_infections, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_g$aov.tab
perm_g_euclidean<-perm_g$aov.tab$`Pr(>F)` # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Oral_infections, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_p$aov.tab

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
write.csv2(beta, file="Other_factors/INFECTIONS_Beta_diversity_permanova_Hellinger.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,as(sample_data(data_sub),"data.frame")$Oral_infections)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Other_factors/INFECTIONS_General Beta_dispersion_permanova_Hellinger.csv",quote=F,row.names = T)


### PCoA INFECTIONS
data.temp <- transform_sample_counts(data_sub.prop, function(x) x/sum(x))
data.sqrt_prop<-transform_sample_counts(data.temp, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
plot_ordination(data.sqrt_prop, ordBC, color = "Oral_infections") +
  stat_ellipse(size= 0.15) +
  scale_color_manual(values = c("no" = "#52d0f8", "yes"="red")) +
  geom_point(size=2.8, alpha= 0.3) + theme_classic(base_size = 9) + 
  labs(title="PCoA with Hellinger distance (saliva samples)",
       color="Oral Infections", x=paste("PC1: ",eigval[1],"% variance"), 
       y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2, show.legend = FALSE)
ggsave(file="Other_factors/INFECTIONS_PCoA_Beta diversity_Hellinger_with_names.png", width = 4.5, height = 3.5, dpi=300)
plot_ordination(data.sqrt_prop, ordBC, color = "Oral_infections") +
  stat_ellipse(size= 0.15) +
  scale_color_manual(values = c("no" = "#52d0f8", "yes"="red")) +
  geom_point(size=3.15, alpha= 0.4) +
  theme_classic(base_size = 9.5) +
  labs(title="A) \n\nPCoA with Hellinger distance on ASVs",
       color="Oral infections", x=paste("PC1: ",eigval[1],"% variance"), 
       y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1]))
ggsave(file="Other_factors/INFECTIONS_PCoA_Beta diversity_Hellinger.png", width = 5, height = 4.5, dpi=300)


########## PCoA also on other factors

### caries
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~caries, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot
plot_ordination(data.sqrt_prop, ordBC, color = "caries") +
  scale_color_manual(values = c("no" = "#52d0f8", "yes"="red")) +
  stat_ellipse(size= 0.15) +
  geom_point(size=2.8, alpha= 0.3) + theme_classic(base_size = 9) + 
  labs(title="PCoA with Hellinger distance (saliva samples)",
       color="Caries", x=paste("PC1: ",eigval[1],"% variance"), 
       y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2, show.legend = FALSE)
ggsave(file="Other_factors/CARIES_PCoA_Beta diversity_Hellinger.png", width = 5, height =5, dpi=300)
### Smoker
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Smoker, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot
plot_ordination(data.sqrt_prop, ordBC, color = "Smoker") +
  stat_ellipse(size= 0.2) +
  scale_color_manual(values = c("no" = "lightblue", "yes"="grey2")) +
  geom_point(size=2.8, alpha= 0.3) + theme_classic(base_size = 9) + 
  labs(title="PCoA with Hellinger distance (saliva samples)",
       color="Smoker", x=paste("PC1: ",eigval[1],"% variance"), 
       y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2, show.legend = FALSE)
ggsave(file="Other_factors/SMOKER_PCoA_Beta diversity_Hellinger.png", width = 5, height = 5, dpi=300)
### Abitudinal_Mouthwash
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Abitudinal_Mouthwash, data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_ASV_euclidean<-perm_ASV$aov.tab$`Pr(>F)` # needed later for plot
plot_ordination(data.sqrt_prop, ordBC, color = "Abitudinal_Mouthwash") +
  stat_ellipse(size= 0.15) +
  scale_color_manual(values = c("no" = "#52d0f8", "yes"="red")) +
  geom_point(size=2.8, alpha= 0.3) + theme_classic(base_size = 9) + 
  labs(title="PCoA with Hellinger distance (saliva samples)",
       color="Mouthwash", x=paste("PC1: ",eigval[1],"% variance"), 
       y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_euclidean[1])) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=2, show.legend = FALSE)
ggsave(file="Other_factors/MOUNTHWASH_PCoA_Beta diversity_Hellinger.png", width = 5, height = 5, dpi=300)


#################### ALPHA DIVERSITY ON OTHER FACTORS (SALIVA) #################

######## ORAL INFECTIONS
suppressWarnings(rm(pAlpha, data_sub))
data_sub<-subset_samples(data, Condition == "HIV") # this also excludes T24
pAlpha<-plot_richness(data_sub, measures=c("Shannon", "Observed"), x="Oral_infections")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
New_data<-New_data[order(New_data$Oral_infections, New_data$Sample),]
pAlpha$data<-New_data
pAlpha +
  geom_point(data=pAlpha$data, aes(x=Oral_infections, y=value, color=NULL), alpha= 0.1) +
  theme_classic(base_size = 10) + 
  geom_boxplot(size= 0.5, alpha = 0.9) +
  labs(x="Oral_infections", title="Alpha diversity in saliva \nbetween HIV patients with and without oral infection") +
  guides(fill=FALSE, color=FALSE) + 
  theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=9),
        title = element_text( size=9)) +
  stat_compare_means(aes(group = Oral_infections), label="p.format", method = "wilcox.test", paired = F, label.x= 1.45, size=3.5, label.y.npc = "top", vjust=-0.4, hjust=0.4)
ggsave(file="Other_factors/Alfa diversity_INFECTIONS.png", width = 6,height =5, dpi=300)



############# DA WITH DESEQ2 (OTHER FACTORS : ORAL INFECTIONS, SALIVA) ################

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_sub<-subset_samples(data, Condition == "HIV") # this also excludes T24
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Oral_infections)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Oral_infections", "no", "yes"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res<-res[res$baseMean > 100, ] # arbitrary threshold to avoid the most noisy result
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_no_vs_yes.csv"), row.names = F, quote=F, na = "")
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

#View(Res_tot)
write.xlsx(Res_tot, file="Other_factors/DA_DESeq2_INFECTIONS_every_result.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Oral_infections)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("no"="blue","yes"="red")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(
    legend.margin=margin(-20, 0, 0, 0), 
    legend.position="bottom", 
    legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
    axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
    axis.text.y = element_text(size=12), plot.title= element_text(size=18),
    panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", 
       fill="Oral_infections", x="")
# ggsave(filename = "Other_factors/DA_DESeq2_INFECTIONS.png", width = 15, height = 8, dpi=300)
dev.off()

unique(Table_tot$Bacteria)
Redund<- c("Spirochaetaceae","Spirochaetia","Spirochaetales","Spirochaetota")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, color= Oral_infections, fill=Oral_infections)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.6, color="black", alpha=0, size= 0.4) + 
  geom_point(position=position_dodge(width=0.6), size=1, alpha= 0.9) +
  geom_point(position=position_dodge(width=0.6), size=2.5, alpha= 0.6) +
  theme_classic(base_size = 11) +
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  scale_color_manual(values=c("no"="#52d0f8","yes"="red")) +
  guides( color=guide_legend(nrow=1), fill="none" ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=13),
        panel.grid.major.x = element_line(linewidth = 0.5),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=11), 
        axis.text.y = element_text(size=10),
        plot.title= element_text(size=12),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1, 2, 3, 5, seq(5,max(Table_tot$Abundance+2),3))) +
  labs(title= "B) \n\nDifferently abundant Taxa (only HIV patients)", y="Proportional Abundance %", 
       color="Oral infections", x="")

DE_plot
ggsave(filename = "Other_factors/DA_Oral_infections_no_redundants.png", width = 5.7, height = 4.9, dpi=300)
dev.off()

rm(DE_plot, data_sub, data_pruned)



##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

setwd("..")

package<-sessionInfo()

con <- file("R_version_and_packages.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"

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
package$otherPkgs$pca3d[c(1,4)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$mixOmics[c(1,4)]

sink() # restore STR OUTPUT to R console
