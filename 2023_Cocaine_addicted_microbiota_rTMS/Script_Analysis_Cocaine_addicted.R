#################### PREPARING THE ENVIRONMENT #################

setwd("/home/Desktop")

library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2")
library("dendextend")
library("qiime2R")

Sample_data<- as.data.frame(readxl::read_excel("Sample_Data.xlsx"))

data.saliva<-qza_to_phyloseq(features="QIIME2/Saliva/table.qza", taxonomy="QIIME2/Saliva/Bayes_saliva_taxonomy.qza", tree = "QIIME2/Saliva/tree.qza")
# updating names and adding sample data
Updating<-sample_names(data.saliva)
{ Updating
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-[0-9][0-9]-","",Updating)
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-16S-","",Updating)
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9]-","",Updating)
Updating<-gsub("[A-Z][0-9][0-9]-[0-9][0-9]-","",Updating)
Updating<-gsub("-[A-Z][0-9]-[A-Z][0-9][0-9]","",Updating)
Updating
}
row.names(Sample_data)<-Sample_data$ID_FASTQ
Saliva<-Sample_data[Updating,]
sample_names(data.saliva)<-Updating
sample_data(data.saliva)<-Saliva
head(sample_data(data.saliva))
identical(length(sample_names(data.saliva)),length(Updating)) # TRUE --> no sample cutted out

data.Stool<-qza_to_phyloseq(features="QIIME2/Stool/table.qza", taxonomy ="QIIME2/Stool/Bayes_Stool_taxonomy.qza", tree = "QIIME2/Stool/Tree/rooted-tree.qza")
# updating names and adding sample data
Updating<-sample_names(data.Stool)
{ Updating
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-[0-9][0-9]-","",Updating)
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-[0-9]-","",Updating)
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-16S-","",Updating)
Updating<-gsub("ID[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9]-","",Updating)
Updating<-gsub("[A-Z][0-9][0-9]-[0-9][0-9]-","",Updating)
Updating<-gsub("-[A-Z][0-9]-[A-Z][0-9][0-9]","",Updating)
Updating
}
row.names(Sample_data)<-Sample_data$ID_FASTQ
Stool<-Sample_data[Updating,]
sample_names(data.Stool)<-Updating
sample_data(data.Stool)<-Stool
head(sample_data(data.Stool))
identical(length(sample_names(data.Stool)),length(Updating)) # TRUE --> no sample cutted out

# to merge the two objects I need to remove both trees
data.sal.temp<-phyloseq(otu_table(data.saliva),tax_table(data.saliva), sample_data(data.saliva))
data.Stool.temp<-phyloseq(otu_table(data.Stool),tax_table(data.Stool), sample_data(data.Stool))
data.total<-merge_phyloseq(data.Stool.temp,data.sal.temp)
identical(length(sample_names(data.total)),(length(sample_names(data.saliva))+length(sample_names(data.Stool)))) # TRUE
rm(Stool,Saliva,Updating,data.sal.temp,data.Stool.temp)

# save.image("data.RData")

dir.create("Results")
setwd("Results")

################# STARTING GENERAL ANALYSIS ###############

dir.create("General_analysis")
setwd("General_analysis")

##################### GENERAL PCoA BRAY CURTIS #####################

dir.create("plots")

data.prop <- transform_sample_counts(data.total, function(ASV) ASV/sum(ASV)*100)
{ data.sqrt_prop<-transform_sample_counts(data.prop, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues 
eigval<- round((eigval/sum(eigval))*100, 1) 
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type", shape = "Condition") +
  geom_point(size=3.5) + theme_classic() + geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n on proportional ASV after sqrt transformation", color="Condition", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_classic()
ggsave(file="plots/PCoA Beta diversity Bray Curtis.png", width = 8, height = 6, dpi=300)

# jaccard
{ DistBC = phyloseq::distance(data.sqrt_prop, method = "jaccard")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues 
eigval<- round((eigval/sum(eigval))*100, 1) 
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type", shape = "Condition") +
  geom_point(size=3.5) + theme_classic() + geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Jaccard distnce \n on proportional ASV after sqrt transformation", color="Condition", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) + theme_classic()
ggsave(file="plots/PCoA Beta diversity Jaccard.png", width = 8, height = 6, dpi=300)

# same sample at the same time between saliva and stool (only untreated and pre_rTMS)
{ p<-subset(Sample_data, Three_times=="No")
pre_rTMS<-subset(Sample_data, Three_times=="Yes" & Time=="pre_rTMS")
p<-rbind.data.frame(p, pre_rTMS, make.row.names = F)
s<-subset(p, Sample_type=="Stool")
f<-subset(p, Sample_type=="Saliva")
pair<- unique(s$Sample_Condition_Time[s$Sample_Condition_Time %in% f$Sample_Condition_Time]) # with both saliva and stool
# NB: healthy samples are automatically removed because not in pair
}

{ data.total.pair<-subset_samples(data.prop, Sample_Condition_Time %in% pair)
data.sqrt_prop<-transform_sample_counts(data.total.pair, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues 
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_type") +
  geom_point(size=4) + theme_classic() + geom_text(aes(label=sample_data(data.sqrt_prop)$Sample), color="black", size=1.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance on proportional ASV after sqrt transformation: \n only cocaine addict samples untreated or treated at pre_rTMS", caption = "the lines are connetting the same sample", color="Sampling type", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) +
  theme_classic() + geom_line(aes(group=Sample_Condition_Time), size=0.1, color="black")
ggsave(file="plots/PCoA Beta diversity_only cocain pairs_Bray Curtis.png", width = 8, height = 6, dpi=300)

rm(f, p, s, DistBC, eigval, pair, ordBC, pre_rTMS)

#################### GENERAL ANALYSIS COUNT EXPORT #####################

{ data.phy = tax_glom(data.total, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data.total, taxrank = "Class", NArm = F)
data.order = tax_glom(data.total, taxrank = "Order", NArm = F)
data.fam = tax_glom(data.total, taxrank = "Family", NArm = F)
data.genus = tax_glom(data.total, taxrank = "Genus", NArm = F)

data.prop <- transform_sample_counts(data.total, function(ASV) ASV/sum(ASV)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)

Taxa.genus<-as.data.frame(tax_table(data.genus))
Taxa.fam<-as.data.frame(tax_table(data.fam))
Taxa.phy<-as.data.frame(tax_table(data.phy))
Taxa.class<-as.data.frame(tax_table(data.class))
Taxa.order<-as.data.frame(tax_table(data.order))
}

{ dir.create("Raw_counts")
write.csv2(cbind(as(otu_table(data.total),"matrix"),as(tax_table(data.total),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)

options(scipen = 100)
dir.create("Relative_abundances")
write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Relative_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Relative_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Relative_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Relative_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Relative_abundances/counts_genus.csv",quote=F)
}

######################### CLOSING GENERAL PART #########################

to_remove<-ls()
rm(list=to_remove[! to_remove %in% c("data.total", "data.saliva", "data.Stool","Sample_data")])
setwd("~/Desktop/Results")

#################### PREPARATION OF SALIVA DATA #######################

dir.create("Saliva")
setwd("Saliva")

{data<-data.saliva
data.pair<-subset_samples(data, Three_times=="Yes") # only "treated" subjects sampled three times
data.pair.placebo<-subset_samples(data.pair, Treatment_type=="Placebo")
data.pair.treated<-subset_samples(data.pair, Treatment_type=="Treated")
}

# data total
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

# pair data total
{data.pair.phy = subset_samples(data.phy, Three_times=="Yes")
data.pair.class = subset_samples(data.class, Three_times=="Yes")
data.pair.order = subset_samples(data.order, Three_times=="Yes")
data.pair.fam = subset_samples(data.fam, Three_times=="Yes")
data.pair.genus = subset_samples(data.genus, Three_times=="Yes")
data.pair.prop <- transform_sample_counts(data.pair, function(ASV) ASV/sum(ASV)*100)
data.pair.phy.prop <- transform_sample_counts(data.pair.phy, function(ASV) ASV/sum(ASV)*100)
data.pair.class.prop <- transform_sample_counts(data.pair.class, function(ASV) ASV/sum(ASV)*100)
data.pair.order.prop <- transform_sample_counts(data.pair.order, function(ASV) ASV/sum(ASV)*100)
data.pair.fam.prop <- transform_sample_counts(data.pair.fam, function(ASV) ASV/sum(ASV)*100)
data.pair.genus.prop <- transform_sample_counts(data.pair.genus, function(ASV) ASV/sum(ASV)*100)
}

# pair data placebo
{data.pair.placebo.phy = subset_samples(data.pair.phy, Treatment_type=="Placebo")
data.pair.placebo.class = subset_samples(data.pair.class, Treatment_type=="Placebo")
data.pair.placebo.order = subset_samples(data.pair.order, Treatment_type=="Placebo")
data.pair.placebo.fam = subset_samples(data.pair.fam, Treatment_type=="Placebo")
data.pair.placebo.genus = subset_samples(data.pair.genus, Treatment_type=="Placebo")
data.pair.placebo.prop <- transform_sample_counts(data.pair.placebo, function(ASV) ASV/sum(ASV)*100)
data.pair.placebo.phy.prop <- transform_sample_counts(data.pair.placebo.phy, function(ASV) ASV/sum(ASV)*100)
data.pair.placebo.class.prop <- transform_sample_counts(data.pair.placebo.class, function(ASV) ASV/sum(ASV)*100)
data.pair.placebo.order.prop <- transform_sample_counts(data.pair.placebo.order, function(ASV) ASV/sum(ASV)*100)
data.pair.placebo.fam.prop <- transform_sample_counts(data.pair.placebo.fam, function(ASV) ASV/sum(ASV)*100)
data.pair.placebo.genus.prop <- transform_sample_counts(data.pair.placebo.genus, function(ASV) ASV/sum(ASV)*100)
}

# pair data treated
{data.pair.treated.phy = subset_samples(data.pair.phy, Treatment_type=="Treated")
data.pair.treated.class = subset_samples(data.pair.class, Treatment_type=="Treated")
data.pair.treated.order = subset_samples(data.pair.order, Treatment_type=="Treated")
data.pair.treated.fam = subset_samples(data.pair.fam, Treatment_type=="Treated")
data.pair.treated.genus = subset_samples(data.pair.genus, Treatment_type=="Treated")
data.pair.treated.prop <- transform_sample_counts(data.pair.treated, function(ASV) ASV/sum(ASV)*100)
data.pair.treated.phy.prop <- transform_sample_counts(data.pair.treated.phy, function(ASV) ASV/sum(ASV)*100)
data.pair.treated.class.prop <- transform_sample_counts(data.pair.treated.class, function(ASV) ASV/sum(ASV)*100)
data.pair.treated.order.prop <- transform_sample_counts(data.pair.treated.order, function(ASV) ASV/sum(ASV)*100)
data.pair.treated.fam.prop <- transform_sample_counts(data.pair.treated.fam, function(ASV) ASV/sum(ASV)*100)
data.pair.treated.genus.prop <- transform_sample_counts(data.pair.treated.genus, function(ASV) ASV/sum(ASV)*100)
}

# vocabulary ASV-taxonomic name
{Taxa.genus<-as.data.frame(tax_table(data.genus))
Taxa.fam<-as.data.frame(tax_table(data.fam))
Taxa.phy<-as.data.frame(tax_table(data.phy))
Taxa.class<-as.data.frame(tax_table(data.class))
Taxa.order<-as.data.frame(tax_table(data.order))
}

#################### TAXA NOT ASSIGNED IN SILVA #########################

# counting not assigned
{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Total","assigned","%","Taxa")
}
assigned
write.csv2(assigned,file="Percent_assigned_in_SILVA.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)

############ ADDING INFORMATIONS TO MISSING TAXA NAMES ##############

{taxa_temp<-Taxa.genus
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
Taxa.genus<-taxa_temp
}

rm(taxa_temp)

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
Taxa.fam<-taxa_temp
}

rm(taxa_temp, x)

######################## SALIVA COUNTS EXPORT ####################################

dir.create("Abundances")
setwd("Abundances")

{dir.create("Raw_counts")
write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)

options(scipen = 100) # disable scientific annotation
dir.create("Relative_abundances")
write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Relative_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Relative_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Relative_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Relative_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Relative_abundances/counts_genus.csv",quote=F)
}

setwd("..")

################### SALIVA RAREFACTION CURVE ######################

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

png(file="Sample_saturation.png",width=3000,height=2100, res=300)
r<-rarecurve(t(otu_table(data)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()

########################### SALIVA ABUNDANCE PLOT ##########################

dir.create("Abundances")
setwd("Abundances")

# TOP 5 Phylum con stacked bar plot
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.prop)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-subset(tabella, Time!= "during_rTMS")
tabella$Time<-gsub("pre_rTMS","pre_rTMS and untreated",tabella$Time)
tabella$Time<-factor(tabella$Time, levels = c("Healthy","pre_rTMS and untreated","post_rTMS"))
tabella$sample_Sample<-gsub("Unknown_","Healthy",tabella$sample_Sample)
}
ggplot(data=tabella, aes(x=Abundance, y=sample_Sample, fill=Phylum)) + facet_grid(rows = vars(Time),scales = "free_y",space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 17.5) +
  theme(panel.spacing = unit(0.8,"lines"), axis.text.y=element_text(angle=0, size = 9), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(y="Sample", x="Relative abundance", title = "Five most abundant phyla and other taxa in saliva samples")
ggsave(file="TOP5_phyla.png",width=9,height=12, dpi=300)
dev.off()

# TOP 5 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.genus.prop)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-subset(tabella, Time!= "during_rTMS")
tabella$Time<-gsub("pre_rTMS","pre_rTMS and untreated",tabella$Time)
tabella$Time<-factor(tabella$Time, levels = c("Healthy","pre_rTMS and untreated","post_rTMS"))
tabella$sample_Sample<-gsub("Unknown_","Healthy",tabella$sample_Sample)
}
ggplot(data=tabella, aes(x=Abundance, y=sample_Sample, fill=Genus)) + facet_grid(rows = vars(Time),scales = "free_y",space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 17.5) +
  theme(panel.spacing = unit(0.8,"lines"), axis.text.y=element_text(angle=0, size = 9), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(y="Sample", x="Relative abundance", title = "Five most abundant genera and other taxa in saliva samples")
ggsave(file="TOP5_genera.png",width=9,height=12,dpi=300)
dev.off()

library(ggh4x) # needed for facet_nested

# TOP 5 phyla in pairs (two times)
{top5 <- names(sort(taxa_sums(data.pair.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.pair.phy.prop)
others<-taxa_names(data.pair.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.pair.phy.prop)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-subset(tabella, Time!= "during_rTMS")
}
ggplot(data=tabella, aes(x=Abundance, y=Time, fill=Phylum)) + facet_nested(Treatment_type+sample_Sample ~ ., scales = "free", space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 15) +
  theme(panel.spacing = unit(0.4,"lines"), axis.text.y=element_text(angle=0, size = 10), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + 
  labs(y="Sample", x="Relative abundance", title = "Five most abundant phyla and other taxa in saliva of treated samples")
ggsave(file="TOP5_phyla_in pair two times.png",width=10,height=10, dpi=300)
dev.off()

# TOP 5 genera in pairs (three times)
{top <- names(sort(taxa_sums(data.pair.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.pair.genus.prop)
others<-taxa_names(data.pair.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.pair.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-subset(tabella, Time!= "during_rTMS")
}
ggplot(data=tabella, aes(x=Abundance, y=Time, fill=Genus)) + facet_nested(Treatment_type + sample_Sample ~ .,scales = "free_y",space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 15) +
  theme(panel.spacing = unit(0.4,"lines"), axis.text.y=element_text(angle=0, size = 10), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) +  labs(y="Sample", x="Relative abundance", title = "Five most abundant genera and other taxa in saliva of treated samples")
ggsave(file="TOP5_genera_in pairs two times.png",width=10,height=10,dpi=300)
dev.off()

setwd("..")

###################### SALIVA HIERARCHICAL CLUSTERING ###################

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))
c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
c<-as.dendrogram(c)
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("Treated","coral",tabella_colore$Colore) # orange
colors <- gsub("Healthy","chartreuse2",colors) #green
colors <- gsub("Addicted","red4",colors)
colors_to_use <- colors[order.dendrogram(c)]
labels_colors(c) <- colors_to_use
}
png(file="Hierarchical cluster EUCLIDEAN on normalized ASVs.png",width=3000,height=1800, res=300)
par(mar=c(5,2,3,0),cex=0.75, cex.main=1.75, cex.sub = 1.25, cex.axis= 1.25)
plot(c, main="Community structure computed with Euclidean distance \n on sqrt proportional ASVs", 
     sub="Healthy= green          Addicted (untreated and placebo) = red           Treated= orange")
dev.off()

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))
  c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))), method = "bray"))
  c<-as.dendrogram(c)
  tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
  colnames(tabella_colore)<-c("Gruppo","Colore")
  colors <- gsub("Treated","coral",tabella_colore$Colore) # orange
  colors <- gsub("Healthy","chartreuse2",colors) #green
  colors <- gsub("Addicted","red4",colors)
  colors_to_use <- colors[order.dendrogram(c)]
  labels_colors(c) <- colors_to_use
}
png(file="Hierarchical cluster BRAY on normalized ASVs.png",width=3000,height=1800, res=300)
par(mar=c(5,2,3,0),cex=0.75, cex.main=1.75, cex.sub = 1.25, cex.axis= 1.25)
plot(c, main="Community structure computed with Bray-Curtis distance \n on sqrt proportional ASVs", 
     sub="Healthy= green          Addicted (untreated and placebo) = red           Treated= orange")
dev.off()

##################### ALFA DIVERSITY IN SALIVA ############################

dir.create("Alpha_diversity")
setwd("Alpha_diversity")

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

# healthy vs cocain addicted (untreated)
{pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition")
pAlpha$data<-dplyr::filter(pAlpha$data, Condition!="Treated")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev)
New_data$Condition<-gsub("pre_rTMS","Addicted",New_data$Condition)
New_data$Condition<-factor(New_data$Condition, levels = c("Healthy","Addicted"))
pAlpha$data<-New_data
}
library(ggpubr)
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw(base_size = 12.2) + 
  labs(x="Condition", title="Alpha diversity in healthy and cocaine addicted (pre_rTMS and untreated) \n in saliva samples")+
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=10)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4)
ggsave(file="Alfa diversity Healthy vs Addicted pre_rTMS.png", width = 8,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor) # OK

#############       now pre_rTMS vs post_rTMS (only samples with all pre_rTMS and post_rTMS time)

####################################################### placebo

rm(data.temp)
data.temp<-subset_samples(data.pair.placebo, Time!="during_rTMS") # only pre vs post
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
New_data<-pAlpha$data
# selecting only samples in pair with *all* two times 
{ pre_rTMS<-subset(New_data, Time=="pre_rTMS")
post_rTMS<-subset(New_data, Time=="post_rTMS")
pre_rTMS<-subset(pre_rTMS, Sample %in% post_rTMS$Sample)
post_rTMS<-subset(post_rTMS, Sample %in% pre_rTMS$Sample)
pre_rTMS<-pre_rTMS[order(pre_rTMS$Sample),]
post_rTMS<-post_rTMS[order(post_rTMS$Sample),]
New_data<-rbind.data.frame(pre_rTMS, post_rTMS)
identical(New_data[New_data$Time=="pre_rTMS","Sample"],New_data[New_data$Time=="post_rTMS","Sample"] )
}
# TRUE --> same order --> updating
{pAlpha$data<-New_data
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev)
New_data$Time<-factor(New_data$Time, levels=c("pre_rTMS","post_rTMS"))
pAlpha$data<-New_data
}
library(ggpubr)
pAlpha + theme_bw(base_size = 12.3) + 
  labs(x="Times", title="Alpha diversity in pre_rTMS and post_rTMS of placebo group \n in saliva samples", caption = "the lines are connetting the same sample") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=10)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired=T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.3) +
  geom_line(aes(group=Sample), size=0.3, col="blue")
ggsave(file="Alfa diversity placebo pre_rTMS vs post_rTMS.png", width = 8,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Time
wilcox.test(Obser_value$value~factor, paired=T) # OK

rm(con, pre_rTMS, post_rTMS, Obser_value, New_data, alphadt, obs, H)

####################################################### treated

rm(data.temp)
data.temp<-subset_samples(data.pair.treated, Time!="during_rTMS") # only pre vs post
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
pAlpha<-plot_richness(data.pair.treated, measures=c("Shannon", "Observed"), x="Time")
New_data<-pAlpha$data
# selecting only samples in pair with *all* two times 
{ pre_rTMS<-subset(New_data, Time=="pre_rTMS")
  post_rTMS<-subset(New_data, Time=="post_rTMS")
  pre_rTMS<-subset(pre_rTMS, Sample %in% post_rTMS$Sample)
  post_rTMS<-subset(post_rTMS, Sample %in% pre_rTMS$Sample)
  pre_rTMS<-pre_rTMS[order(pre_rTMS$Sample),]
  post_rTMS<-post_rTMS[order(post_rTMS$Sample),]
  New_data<-rbind.data.frame(pre_rTMS, post_rTMS)
}
 identical(New_data[New_data$Time=="pre_rTMS","Sample"],New_data[New_data$Time=="post_rTMS","Sample"] )
# TRUE --> same order --> updating
{pAlpha$data<-New_data
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  New_data$Time<-factor(New_data$Time, levels=c("pre_rTMS","post_rTMS"))
  pAlpha$data<-New_data
}
library(ggpubr)
pAlpha + theme_bw(base_size = 12.3) + 
  labs(x="Times", title="Alpha diversity between pre_rTMS and post_rTMS of treated group \n in saliva samples", caption = "the lines are connetting the same sample") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=10)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired=T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.3) +
  geom_line(aes(group=Sample), size=0.3, col="blue")
ggsave(file="Alfa diversity treated pre_rTMS vs post_rTMS.png", width = 8,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Time
wilcox.test(Obser_value$value~factor, paired=T) # OK

rm(con, pre_rTMS, post_rTMS, Obser_value,  New_data, alphadt, obs, H)

setwd("..")

######################## BETA DIVERSITY BRAY CURTIS #######################

dir.create("Bray_Curtis_Beta_diversity")
setwd("Bray_Curtis_Beta_diversity")

# NB: only between pre and post
{ASV.prop<-as.data.frame(otu_table(data.prop))
ASV.prop<-ASV.prop[ , ! colnames(ASV.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
ASV.genus.prop<-ASV.genus.prop[ , ! colnames(ASV.genus.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
ASV.fam.prop<-ASV.fam.prop[ , ! colnames(ASV.fam.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
ASV.class.prop<-ASV.class.prop[, ! colnames(ASV.class.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
ASV.order.prop<-ASV.order.prop[, ! colnames(ASV.order.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
ASV.phy.prop<-ASV.phy.prop[, ! colnames(ASV.phy.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")
metadata<-subset(metadata, Time != "during_rTMS")
metadata$Treatment_time<-paste(metadata$Treatment_type,metadata$Time, sep="_")
library(vegan)

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_ASV

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_g

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_f

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_o

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_c

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_p

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

# pairwise beta diversity
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=metadata$Treatment_time, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV
pair_ASV$pairs<-gsub("Healthy_Healthy","Healthy", pair_ASV$pairs)
pair_ASV$pairs<-gsub("pre_rTMS","pre-rTMS", pair_ASV$pairs)
pair_ASV$pairs<-gsub("post_rTMS","post-rTMS", pair_ASV$pairs)

# exporting beta diversity
rm(con)
con<-file("Beta diversity general and pairwise between Healthy Untreated and T1_post_rTMS_Treated")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n", fill=TRUE)
pair_ASV
sink()
close(con)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
disper<-vegan::betadisper(BC.dist,metadata$Treatment_time)
disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
disp_ASV$tab
a<-as.data.frame(disp_ASV$pairwise$permuted)
colnames(a)<-c("permuted_p_value")
a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
a$Significative<-a$padj_BH
a$Significative[a$padj_BH < 0.05]<- "*"
a$Significative[a$Significative > 0.05]<- ""
a
row.names(a)<-gsub("Healthy_Healthy","Healthy", row.names(a))
row.names(a)<-gsub("pre_rTMS","pre-rTMS", row.names(a))
row.names(a)<-gsub("post_rTMS","post-rTMS", row.names(a))

#export dispersion
rm(con)
con<-file("Beta dispersion General and Pairwise between Healthy Untreated and T1_post_rTMS_Treated")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n")
a
sink()
close(con)

rm(disp_ASV,a)

########################### PCoA BRAY CURTIS #####################

# NB: only pre and post treatment

# in base ad ASV
data.prop.labels<-subset_samples(data.prop, Time!="during_rTMS")
{sample_data(data.prop.labels)$Time<-gsub("pre_rTMS","Addicted",sample_data(data.prop.labels)$Time)
sample_data(data.prop.labels)$True_treated_Sign<-sample_data(data.prop.labels)$Under_TRUE_treatment
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Healthy","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Untreated","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Treated","*",sample_data(data.prop.labels)$True_treated_Sign)
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.005))+
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral")) +
  labs(title="PCoA computed with Bray-Curtis distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign")
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_Conditions.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + 
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.005))+
  labs(title="PCoA computed with Bray-Curtis distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign") +
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_no_ellipses.png", width = 8, height = 6, dpi=300)
# without post_rTMS
data.temp<-subset_samples(data.sqrt_prop, Time!="post_rTMS")
data.temp<-transform_sample_counts(data.temp, function(x) x/sum(x))
DistBC = phyloseq::distance(data.temp, method = "bray")
ordBC = ordinate(data.temp, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
plot_ordination(data.temp, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  labs(title="PCoA computed with Bray-Curtis distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) +
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_solo healthy vs addicted.png", width = 8, height = 6, dpi=300)

setwd("..")

######################## BETA DIVERSITY WEIGHTED UNIFRAC #######################

rm(ASV.prop, ASV.genus.prop, perm_ASV, perm_g, perm_f)

dir.create("Weighted Unifrac Beta diversity")
setwd("Weighted Unifrac Beta diversity")

# NB: only between pre and post

#### PERMANOVA
{metadata<-as(sample_data(data.prop),"data.frame")
metadata<-subset(metadata, Time != "during_rTMS")
metadata$Treatment_time<-paste(metadata$Treatment_type,metadata$Time, sep="_")
library(vegan)
}

{rm(data.temp)
data.temp<-subset_samples(data.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac") # needing phyloseq to compute unifrac
perm_ASV<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
perm_ASV
}
{rm(data.temp)
data.temp<-subset_samples(data.genus.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_g<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
perm_g
}
{rm(data.temp)
data.temp<-subset_samples(data.fam.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_f<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
perm_f
}
{rm(data.temp)
data.temp<-subset_samples(data.order.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_o<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
perm_o
}
{rm(data.temp)
data.temp<-subset_samples(data.class.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_c<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
perm_c
}
{rm(data.temp)
data.temp<-subset_samples(data.phy.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_p<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
perm_p
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

# pairwise beta diversity
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# on genera because on ASV level there is to much noise
rm(data.temp)
data.temp<-subset_samples(data.genus.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
pair_g<- pairwise.adonis(DistBC,  factors=metadata$Treatment_time, p.adjust.m = "BH", sim.method="wunifrac", perm = 9999)
pair_g
pair_g$pairs<-gsub("Healthy_Healthy","Healthy", pair_g$pairs)
pair_g$pairs<-gsub("pre_rTMS","pre-rTMS", pair_g$pairs)
pair_g$pairs<-gsub("post_rTMS","post-rTMS", pair_g$pairs)

# exporting beta diversity
rm(con)
con<-file("Beta diversity general and pairwise between Healthy Addicted and Untreated VS post_rTMS_")
sink(con, append = TRUE)
cat("General beta diversity Weighted Unifrac on Genera \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Weighted Unifrac (correction Benjamini-Hochberg on p-value) on Genera \n", fill=TRUE)
pair_ASV
sink()

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
rm(data.temp)
data.temp<-subset_samples(data.genus.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
disper<-vegan::betadisper(DistBC,metadata$Treatment_time)
disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
disp_ASV$tab
a<-as.data.frame(disp_ASV$pairwise$permuted)
colnames(a)<-c("permuted_p_value")
a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
a$Significative<-a$padj_BH
a$Significative[a$padj_BH < 0.05]<- "*"
a$Significative[a$Significative > 0.05]<- ""
a
row.names(a)<-gsub("Healthy_Healthy","Healthy", row.names(a))
row.names(a)<-gsub("pre_rTMS","pre-rTMS", row.names(a))
row.names(a)<-gsub("post_rTMS","post-rTMS", row.names(a))

#export dispersion
rm(con)
con<-file("Beta dispersion General and Pairwise between Healthy Addicted and Untreated vs post_rTMS")
sink(con, append=TRUE)
cat("General beta dispersion Weighted Unifrac \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) on Genera \n \n")
a
sink()

rm(disp_ASV,a)

########################### PCoA WEIGHTED UNIFRAC #####################

# in base ad ASV
{ data.prop.labels<-subset_samples(data.prop, Time!="during_rTMS")
sample_data(data.prop.labels)$Time<-gsub("pre_rTMS","Addicted",sample_data(data.prop.labels)$Time)
sample_data(data.prop.labels)$True_treated_Sign<-sample_data(data.prop.labels)$Under_TRUE_treatment
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Healthy","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Untreated","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Treated","*",sample_data(data.prop.labels)$True_treated_Sign)
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "wunifrac")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.000005))+
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral")) +
  labs(title="PCoA computed with Weighted Unifrac distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign")
ggsave(file="PCoA_Beta_diversity_W_Unifrac_on_ASV.png", width = 8, height = 6, dpi=300)

# it's much better on genera (less noise)
{data.prop.labels<-subset_samples(data.genus.prop, Time!="during_rTMS")
  sample_data(data.prop.labels)$Time<-gsub("pre_rTMS","Addicted",sample_data(data.prop.labels)$Time)
  sample_data(data.prop.labels)$True_treated_Sign<-sample_data(data.prop.labels)$Under_TRUE_treatment
  sample_data(data.prop.labels)$True_treated_Sign<-gsub("Healthy","",sample_data(data.prop.labels)$True_treated_Sign)
  sample_data(data.prop.labels)$True_treated_Sign<-gsub("Untreated","",sample_data(data.prop.labels)$True_treated_Sign)
  sample_data(data.prop.labels)$True_treated_Sign<-gsub("Treated","*",sample_data(data.prop.labels)$True_treated_Sign)
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "wunifrac")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.0002))+
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral")) +
  labs(title="PCoA computed with Weighted Unifrac distance \n on sqrt proportional genera", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign")
ggsave(file="PCoA_Beta_diversity_W_Unifrac_sui generi.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + 
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.0002))+
  labs(title="PCoA computed with Weighted Unifrac distance \n on sqrt proportional genera", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign") +
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral"))
ggsave(file="PCoA_Beta_diversity_W_Unifrac_no_ellipses_generi.png", width = 8, height = 6, dpi=300)

setwd("..")

############## DA WITH DESEQ2 - Healthy vs Untreated SALIVA #################

dir.create("DESeq2")
setwd("DESeq2")

# Trimming under sum of 10, see DESeq2  tutorial
data_pruned<- prune_taxa(taxa_sums(data) > 10, data)
# following are needed in case of different ASV names due to pruning
{data_pruned.phy = tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
data_pruned.class = tax_glom(data_pruned, taxrank = "Class", NArm = F)
data_pruned.order = tax_glom(data_pruned, taxrank = "Order", NArm = F)
data_pruned.fam = tax_glom(data_pruned, taxrank = "Family", NArm = F)
data_pruned.genus = tax_glom(data_pruned, taxrank = "Genus", NArm = F)
data_pruned.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data_pruned.phy.prop <- transform_sample_counts(data_pruned.phy, function(ASV) ASV/sum(ASV)*100)
data_pruned.class.prop <- transform_sample_counts(data_pruned.class, function(ASV) ASV/sum(ASV)*100)
data_pruned.order.prop <- transform_sample_counts(data_pruned.order, function(ASV) ASV/sum(ASV)*100)
data_pruned.fam.prop <- transform_sample_counts(data_pruned.fam, function(ASV) ASV/sum(ASV)*100)
data_pruned.genus.prop <- transform_sample_counts(data_pruned.genus, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured taxa
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
  taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
tax_table(data_pruned.genus.prop)<-as.matrix(taxa_temp)
rm(taxa_temp)
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
tax_table(data_pruned.fam.prop)<-as.matrix(taxa_temp)
}
rm(taxa_temp)

# function to automatize the comparison between groups
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Condition)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Condition)
    }
    # DE<-DESeq(d, test="LRT",reduced= ~ 1)
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Condition", nume[i], deno[i]), test="Wald")
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        #
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("Healthy")        #("nume" e "deno" choose who confront)
deno<-c("Addicted")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100) # disable scientific annotation
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)  #(automatically computes all DESeq2 analyses)
{DE_results <- read.table("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
write.csv2(DE_results, file="DE results_Healthy vs Untreated.csv", quote=F, row.names = F)
unlink("DE_results.txt")
}

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
target<-prune_taxa(target, data_pruned.genus.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_g<-tabella
}
{target<-DE_results[DE_results$Rank=="Family","ASV"]
target<-prune_taxa(target, data_pruned.fam.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_f<-tabella
}
{target<-DE_results[DE_results$Rank=="Order","ASV"]
target<-prune_taxa(target, data_pruned.order.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_o<-tabella
}
{target<-DE_results[DE_results$Rank=="Class","ASV"]
target<-prune_taxa(target, data_pruned.class.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_c<-tabella
}
{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
target<-prune_taxa(target, data_pruned.phy.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_p<-tabella
}
# unique boxplot of DESeq2  results
{tabella_g$Taxa<-"Genera"
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
}

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_c,tabella_o,tabella_p)
tabella_tot<-tabella_tot[tabella_tot$Condition!="Treated",]
# tabella_tot<-tabella_tot[tabella_tot$Under_TRUE_treatment!="Treated",]
tabella_tot<-tabella_tot[tabella_tot$Time!="during_rTMS",]
tabella_tot$Condition<-factor(tabella_tot$Condition, levels=c("Healthy","Addicted"))

# there are reduntants ( same ASV, same abundance, same result! ) --> manteining only the lower level
Redund<-DE_results
Redund<-Redund[duplicated(round(Redund$log2FoldChange, digits = 0), fromLast=T), ]
Redund<-as.character(c(Redund[Redund$Rank=="Phylum","Phylum"],Redund[Redund$Rank=="Class","Class"],Redund[Redund$Rank=="Order","Order"],Redund[Redund$Rank=="Family","Family"]))

tabella_tot2<-subset(tabella_tot, ! Bacteria %in% Redund)

ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) + facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space = "free") +
  geom_boxplot(width=0.8) + theme_bw( base_size= 16) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values = c("Healthy"="Chartreuse","Addicted"="red4")) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = -40, vjust=0.8, hjust=0, size=12), axis.text.y = element_text(size=13)) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  theme(plot.title= element_text(size=18) ,
        legend.key.size=unit(0.7,"cm"), legend.text =element_text(size=18), legend.title = element_text(size=16)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant taxa between Healthy and Addicted (placebo and untreated) patients", y="Sqrt Proportional Abundance", fill="Condition", x="") +
  scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  theme(plot.margin=margin(1,1.5,1,0.3, unit = "cm"))
ggsave(filename = "Total_differences_DESeq2 Healthy vs Untreated_no redundant.png", width = 13, height = 8, dpi=300)

dev.off()

setwd("..")

############## DA WITH DESEQ2 - PLACEBO PAIRS pre_rTMS vs post_rTMS SALIVA #################

setwd("DESeq2")

# Trimming under sum of 10, see DESeq2  tutorial
{data_pruned<- prune_taxa(taxa_sums(data.pair.placebo) > 10, data.pair.placebo)
# following are needed in case of different ASV names due to pruning
data_pruned.phy = tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
data_pruned.class = tax_glom(data_pruned, taxrank = "Class", NArm = F)
data_pruned.order = tax_glom(data_pruned, taxrank = "Order", NArm = F)
data_pruned.fam = tax_glom(data_pruned, taxrank = "Family", NArm = F)
data_pruned.genus = tax_glom(data_pruned, taxrank = "Genus", NArm = F)
data_pruned.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data_pruned.phy.prop <- transform_sample_counts(data_pruned.phy, function(ASV) ASV/sum(ASV)*100)
data_pruned.class.prop <- transform_sample_counts(data_pruned.class, function(ASV) ASV/sum(ASV)*100)
data_pruned.order.prop <- transform_sample_counts(data_pruned.order, function(ASV) ASV/sum(ASV)*100)
data_pruned.fam.prop <- transform_sample_counts(data_pruned.fam, function(ASV) ASV/sum(ASV)*100)
data_pruned.genus.prop <- transform_sample_counts(data_pruned.genus, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured taxa
{
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
  taxa_temp$Genus[which(is.na(taxa_temp$Genus))[1]]<- paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"] )}
for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
  taxa_temp$Genus[which(is.na(taxa_temp$Genus))[1]]<- paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"] )}
for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
  taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
  taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  tax_table(data_pruned.genus.prop)<-as.matrix(taxa_temp)
}
rm(taxa_temp)
{
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
for( x in 1: length(which(is.na(taxa_temp$Family))) ) {
  taxa_temp$Family[which(is.na(taxa_temp$Family))[1]]<- paste("NA_ o",taxa_temp[which(is.na(taxa_temp$Family))[1],"Order"] )}
  tax_table(data_pruned.fam.prop)<-as.matrix(taxa_temp)
}
rm(taxa_temp)

# function to automatize the comparison between groups /// NB: ~Sample + Time
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Sample + Time)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Sample + Time)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Time", nume[i], deno[i]), test="Wald")
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

nume<-c("pre_rTMS")        #("nume" e "deno" choose who confront)
deno<-c("post_rTMS")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100) # disable scientific annotation
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
write.csv2(DE_results, file="DE results_Placebo pre_rTMS vs post_rTMS.csv", quote=F, row.names = F)
unlink("DE_results.txt")

rm(tabella_g,tabella_o,tabella_f,tabella_c,tabella_p, tabella)

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
target<-prune_taxa(target, data_pruned.genus.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_g<-tabella
}

{target<-DE_results[DE_results$Rank=="Family","ASV"]
target<-prune_taxa(target, data_pruned.fam.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_f<-tabella
}

# unique boxplot of DESeq2  results

{tabella_g$Taxa<-"Genera"   # cut off from kingdowm column to last not interesting taxonomic level
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_f$Taxa<-"Families"
tabella_f[,c("Phylum","Order","Class")]<-NULL
colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"
}

tabella_tot<-rbind.data.frame(tabella_g,tabella_f)
tabella_tot<-subset(tabella_tot, Time !="during_rTMS")
tabella_tot$Time<-factor(tabella_tot$Time, levels = c("pre_rTMS","post_rTMS"))

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Time)) + facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space = "free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 16 ) + theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values = c("pre_rTMS"="red4","post_rTMS"="coral")) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = -40, vjust=0.8, hjust=0, size=12), axis.text.y = element_text(size=13)) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  theme(plot.title= element_text(size=18) ,
        legend.key.size=unit(0.7,"cm"), legend.text =element_text(size=16), legend.title = element_text(size=15)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between pre and post rTMS in placebo group", y="Sqrt Proportional Abundance", fill="Time", x="") +
  scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  theme(plot.margin=margin(1,2,1,0.3, unit = "cm"))
ggsave(filename = "Total_differences_DESeq2 Placebo two_times.png", width = 12, height = 8, dpi=300)
dev.off()

setwd("..")

############## DA WITH DESEQ2 - treated PAIRS pre_rTMS vs post_rTMS SALIVA #################

setwd("DESeq2")

# Trimming under sum of 10, see DESeq2  tutorial
{data_pruned<- prune_taxa(taxa_sums(data.pair.treated) > 10, data.pair.treated)
# following are needed in case of different ASV names due to pruning
data_pruned.phy = tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
data_pruned.class = tax_glom(data_pruned, taxrank = "Class", NArm = F)
data_pruned.order = tax_glom(data_pruned, taxrank = "Order", NArm = F)
data_pruned.fam = tax_glom(data_pruned, taxrank = "Family", NArm = F)
data_pruned.genus = tax_glom(data_pruned, taxrank = "Genus", NArm = F)
data_pruned.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data_pruned.phy.prop <- transform_sample_counts(data_pruned.phy, function(ASV) ASV/sum(ASV)*100)
data_pruned.class.prop <- transform_sample_counts(data_pruned.class, function(ASV) ASV/sum(ASV)*100)
data_pruned.order.prop <- transform_sample_counts(data_pruned.order, function(ASV) ASV/sum(ASV)*100)
data_pruned.fam.prop <- transform_sample_counts(data_pruned.fam, function(ASV) ASV/sum(ASV)*100)
data_pruned.genus.prop <- transform_sample_counts(data_pruned.genus, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured taxa
{
  taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
  for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
    taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[which(is.na(taxa_temp$Genus))[1]]<- paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"] )}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[which(is.na(taxa_temp$Genus))[1]]<- paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"] )}
  for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  tax_table(data_pruned.genus.prop)<-as.matrix(taxa_temp)
}
rm(taxa_temp)
{
  taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
  for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
    taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  for( x in 1: length(which(is.na(taxa_temp$Family))) ) {
    taxa_temp$Family[which(is.na(taxa_temp$Family))[1]]<- paste("NA_ o",taxa_temp[which(is.na(taxa_temp$Family))[1],"Order"] )}
  tax_table(data_pruned.fam.prop)<-as.matrix(taxa_temp)
}
rm(taxa_temp)

# function to automatize the comparison between groups /// NB: ~Sample + Time
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Sample + Time)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Sample + Time)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Time", nume[i], deno[i]), test="Wald")
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        #			cat(paste0(rank,":",nume[i],"/",dena[i],"\n"),file=outfile,append=T)
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("pre_rTMS")        #("nume" e "deno" choose who confront)
deno<-c("post_rTMS")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100) # disable scientific annotation
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)             #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
write.csv2(DE_results, file="DE_results_TRUE treated two times.csv", quote=F, row.names = F)
unlink("DE_results.txt")

rm(tabella_g,tabella_o,tabella_f,tabella_c,tabella_p)

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
target<-prune_taxa(target, data_pruned.genus.prop)
tabella<-psmelt(target)
# tabella$Abundance<-sqrt(tabella$Abundance)
tabella_g<-tabella
}
{target<-DE_results[DE_results$Rank=="Family","ASV"]
target<-prune_taxa(target, data_pruned.fam.prop)
tabella<-psmelt(target)
# tabella$Abundance<-sqrt(tabella$Abundance)
tabella_f<-tabella
}
{target<-DE_results[DE_results$Rank=="Order","ASV"]
target<-prune_taxa(target, data_pruned.order.prop)
tabella<-psmelt(target)
# tabella$Abundance<-sqrt(tabella$Abundance)
tabella_o<-tabella
}
{target<-DE_results[DE_results$Rank=="Class","ASV"]
target<-prune_taxa(target, data_pruned.class.prop)
tabella<-psmelt(target)
# tabella$Abundance<-sqrt(tabella$Abundance)
tabella_c<-tabella
}
{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
target<-prune_taxa(target, data_pruned.phy.prop)
tabella<-psmelt(target)
# tabella$Abundance<-sqrt(tabella$Abundance)
tabella_p<-tabella
}

# unique boxplot of DESeq2  results

{tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"
} # all other results are redundants (same result, different level)

tabella_tot<-rbind.data.frame(tabella_g)
tabella_tot<-subset(tabella_tot, Time!="during_rTMS") # only pre vs post
tabella_tot<-subset(tabella_tot, Bacteria!="NA_ o NA") 
tabella_tot$Time<-factor(tabella_tot$Time, levels=c("pre_rTMS","post_rTMS"))

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Time)) + facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space = "free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 16 ) + theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values = c("pre_rTMS"="red4","post_rTMS"="coral")) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = -40, vjust=0.8, hjust=0, size=12), axis.text.y = element_text(size=13)) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  theme(plot.title= element_text(size=18) ,
        legend.key.size=unit(0.7,"cm"), legend.text =element_text(size=16), legend.title = element_text(size=15)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between pre and post rTMS in treated group", y="Sqrt Proportional Abundance", fill="Time", x="") +
  scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  theme(plot.margin=margin(1,3.25,1,0.3, unit = "cm"))
ggsave(filename = "Total_differences_DESeq2 treated Three_times.png", width = 15, height = 8, dpi=300)
dev.off()

setwd("..")

################## HANDLING PICRUST AND LEFSE OUTPUT ################

rm(a)
{a <- read.delim("path_abun_unstrat_descrip.tsv.gz") # output of PICRUST2
colnames(a)<-gsub("ID2117.[0-9][0-9][0-9].","",colnames(a))
colnames(a)<-gsub("ID2117.[0-9][0-9].","",colnames(a))
colnames(a)<-gsub(".A[0-9].[A-Z][0-9][0-9]","",colnames(a))
colnames(a)<-gsub("ID1642.16S.B[0-9][0-9].[0-9][0-9].","",colnames(a))
colnames(a)
Descriptions<-a[,1:2]
row.names(Descriptions)<-Descriptions[,1]
a[,2]<-NULL # no descriptions in LefSE
metadata<-as(sample_data(data),"data.frame")
identical(metadata$ID_FASTQ,colnames(a)[-1]) # TRUE
row.names(a)<-a[,1]
a[,1]<-NULL
a<-rbind.data.frame(metadata$Condition,a)
head(a)
row.names(a)[1]<-"Condition"
}
a<-a[,a[1,]!="Treated"] # filtered out "treated" sample to study healthy vs untreated only

write.table(a,file="Output_KO_ready_for_LefSe_online.tsv",row.names = T,quote = F, sep="\t")
# remember to add a rownames to ex colnames

# --> LEFSe
# <--
Significative_functions_LEFSE<- read.delim("Result.lefse_internal_res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
Significative_functions_LEFSE<-dplyr::left_join(Significative_functions_LEFSE,Descriptions, by=c("Pathway"="pathway"))
head(Significative_functions_LEFSE)

write.csv2(Significative_functions_LEFSE,file = "Significative_functions_LEFSE.csv",na="",quote = F, row.names = F)

######################### CLOSING SALIVA PART #########################

to_remove<-ls()
rm(list=to_remove[! to_remove %in% c("data.total", "data.saliva", "data.Stool","Sample_data")])
setwd("~/Desktop/Results")

#################### PREPARATION OF STOOL DATA #######################

dir.create("Stool")
setwd("Stool")

{data<-data.Stool
  data.pair<-subset_samples(data, Three_times=="Yes") # only "treated" subjects sampled three times
  data.pair.placebo<-subset_samples(data.pair, Treatment_type=="Placebo")
  data.pair.treated<-subset_samples(data.pair, Treatment_type=="Treated")
}

# data total
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
  data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

# pair data total
{data.pair.phy = subset_samples(data.phy, Three_times=="Yes")
  data.pair.class = subset_samples(data.class, Three_times=="Yes")
  data.pair.order = subset_samples(data.order, Three_times=="Yes")
  data.pair.fam = subset_samples(data.fam, Three_times=="Yes")
  data.pair.genus = subset_samples(data.genus, Three_times=="Yes")
  data.pair.prop <- transform_sample_counts(data.pair, function(ASV) ASV/sum(ASV)*100)
  data.pair.phy.prop <- transform_sample_counts(data.pair.phy, function(ASV) ASV/sum(ASV)*100)
  data.pair.class.prop <- transform_sample_counts(data.pair.class, function(ASV) ASV/sum(ASV)*100)
  data.pair.order.prop <- transform_sample_counts(data.pair.order, function(ASV) ASV/sum(ASV)*100)
  data.pair.fam.prop <- transform_sample_counts(data.pair.fam, function(ASV) ASV/sum(ASV)*100)
  data.pair.genus.prop <- transform_sample_counts(data.pair.genus, function(ASV) ASV/sum(ASV)*100)
}

# pair data placebo
{data.pair.placebo.phy = subset_samples(data.pair.phy, Treatment_type=="Placebo")
  data.pair.placebo.class = subset_samples(data.pair.class, Treatment_type=="Placebo")
  data.pair.placebo.order = subset_samples(data.pair.order, Treatment_type=="Placebo")
  data.pair.placebo.fam = subset_samples(data.pair.fam, Treatment_type=="Placebo")
  data.pair.placebo.genus = subset_samples(data.pair.genus, Treatment_type=="Placebo")
  data.pair.placebo.prop <- transform_sample_counts(data.pair.placebo, function(ASV) ASV/sum(ASV)*100)
  data.pair.placebo.phy.prop <- transform_sample_counts(data.pair.placebo.phy, function(ASV) ASV/sum(ASV)*100)
  data.pair.placebo.class.prop <- transform_sample_counts(data.pair.placebo.class, function(ASV) ASV/sum(ASV)*100)
  data.pair.placebo.order.prop <- transform_sample_counts(data.pair.placebo.order, function(ASV) ASV/sum(ASV)*100)
  data.pair.placebo.fam.prop <- transform_sample_counts(data.pair.placebo.fam, function(ASV) ASV/sum(ASV)*100)
  data.pair.placebo.genus.prop <- transform_sample_counts(data.pair.placebo.genus, function(ASV) ASV/sum(ASV)*100)
}

# pair data treated
{data.pair.treated.phy = subset_samples(data.pair.phy, Treatment_type=="Treated")
  data.pair.treated.class = subset_samples(data.pair.class, Treatment_type=="Treated")
  data.pair.treated.order = subset_samples(data.pair.order, Treatment_type=="Treated")
  data.pair.treated.fam = subset_samples(data.pair.fam, Treatment_type=="Treated")
  data.pair.treated.genus = subset_samples(data.pair.genus, Treatment_type=="Treated")
  data.pair.treated.prop <- transform_sample_counts(data.pair.treated, function(ASV) ASV/sum(ASV)*100)
  data.pair.treated.phy.prop <- transform_sample_counts(data.pair.treated.phy, function(ASV) ASV/sum(ASV)*100)
  data.pair.treated.class.prop <- transform_sample_counts(data.pair.treated.class, function(ASV) ASV/sum(ASV)*100)
  data.pair.treated.order.prop <- transform_sample_counts(data.pair.treated.order, function(ASV) ASV/sum(ASV)*100)
  data.pair.treated.fam.prop <- transform_sample_counts(data.pair.treated.fam, function(ASV) ASV/sum(ASV)*100)
  data.pair.treated.genus.prop <- transform_sample_counts(data.pair.treated.genus, function(ASV) ASV/sum(ASV)*100)
}

# vocabulary ASV-taxonomic name
{Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

#################### TAXA NOT ASSIGNED IN SILVA #########################

# counting not assigned
{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Totale","assigned","%","Taxa")
}
assigned
write.csv2(assigned,file="Percent_assigned_in_SILVA_Stool.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)

############ ADDING INFORMATIONS TO MISSING TAXA NAMES ##############

{taxa_temp<-Taxa.genus
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
Taxa.genus<-taxa_temp
}

rm(taxa_temp)

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
  Taxa.fam<-taxa_temp
}

rm(taxa_temp, x)

######################## STOOL COUNTS EXPORT ####################################

dir.create("Abundances")
setwd("Abundances")

{dir.create("Raw_counts")
  write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Raw_counts/counts_genus.csv",quote=F)
  
  options(scipen = 100) # disable scientific annotation
  dir.create("Relative_abundances")
  write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Relative_abundances/counts_genus.csv",quote=F)
}

setwd("..")

################### STOOL RAREFACTION CURVE ######################

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

png(file="Sample_saturation.png",width=3000,height=2100, res=300)
r<-rarecurve(t(otu_table(data)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()

########################### STOOL ABUNDANCE PLOT ##########################

dir.create("Abundances")
setwd("Abundances")

# TOP 5 Phylum con stacked bar plot
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-subset(tabella, Time!= "during_rTMS")
  tabella$Time<-gsub("pre_rTMS","pre_rTMS and untreated",tabella$Time)
  tabella$Time<-factor(tabella$Time, levels = c("Healthy","pre_rTMS and untreated","post_rTMS"))
  tabella$sample_Sample<-gsub("Unknown_","Healthy",tabella$sample_Sample)
}
ggplot(data=tabella, aes(x=Abundance, y=sample_Sample, fill=Phylum)) + facet_grid(rows = vars(Time),scales = "free_y",space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 17.1) +
  theme(panel.spacing = unit(0.8,"lines"), axis.text.y=element_text(angle=0, size = 9), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(y="Sample", x="Relative abundance", title = "Five most abundant phyla and other taxa in stool samples")
ggsave(file="TOP5_phyla.png",width=9,height=12, dpi=300)
dev.off()

# TOP 5 Generi
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-subset(tabella, Time!= "during_rTMS")
  tabella$Time<-gsub("pre_rTMS","pre_rTMS and untreated",tabella$Time)
  tabella$Time<-factor(tabella$Time, levels = c("Healthy","pre_rTMS and untreated","post_rTMS"))
  tabella$sample_Sample<-gsub("Unknown_","Healthy",tabella$sample_Sample)
}
ggplot(data=tabella, aes(x=Abundance, y=sample_Sample, fill=Genus)) + facet_grid(rows = vars(Time),scales = "free_y",space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 17.1) +
  theme(panel.spacing = unit(0.8,"lines"), axis.text.y=element_text(angle=0, size = 9), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(y="Sample", x="Relative abundance", title = "Five most abundant genera and other taxa in stool samples")
ggsave(file="TOP5_genera.png",width=9,height=12,dpi=300)
dev.off()

library(ggh4x) # needed for facet_nested

# TOP 5 phyla in pairs (two times)
{top5 <- names(sort(taxa_sums(data.pair.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.pair.phy.prop)
  others<-taxa_names(data.pair.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.pair.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-subset(tabella, Time!= "during_rTMS")
}
ggplot(data=tabella, aes(x=Abundance, y=Time, fill=Phylum)) + facet_nested(Treatment_type+sample_Sample ~ ., scales = "free", space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 15) +
  theme(panel.spacing = unit(0.4,"lines"), axis.text.y=element_text(angle=0, size = 10), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + 
  labs(y="Sample", x="Relative abundance", title = "Five most abundant phyla and other taxa in stool of treated samples")
ggsave(file="TOP5_phyla_in pair two times.png",width=10,height=10, dpi=300)
dev.off()

# TOP 5 genera in pairs (three times)
{top <- names(sort(taxa_sums(data.pair.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.pair.genus.prop)
  others<-taxa_names(data.pair.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.pair.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-subset(tabella, Time!= "during_rTMS")
}
ggplot(data=tabella, aes(x=Abundance, y=Time, fill=Genus)) + facet_nested(Treatment_type + sample_Sample ~ .,scales = "free_y",space="free") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size = 15) +
  theme(panel.spacing = unit(0.4,"lines"), axis.text.y=element_text(angle=0, size = 10), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) +  labs(y="Sample", x="Relative abundance", title = "Five most abundant genera and other taxa in stool of treated samples")
ggsave(file="TOP5_genera_in pairs two times.png",width=10,height=10,dpi=300)
dev.off()

setwd("..")

###################### STOOL HIERARCHICAL CLUSTERING ###################

library(dendextend)
{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))
  c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
  c<-as.dendrogram(c)
  tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
  colnames(tabella_colore)<-c("Gruppo","Colore")
  colors <- gsub("Treated","coral",tabella_colore$Colore) # orange
  colors <- gsub("Healthy","chartreuse2",colors) #green
  colors <- gsub("Addicted","red4",colors)
  colors_to_use <- colors[order.dendrogram(c)]
  labels_colors(c) <- colors_to_use
}
png(file="Hierarchical cluster EUCLIDEAN on normalized ASVs.png",width=3000,height=1800, res=300)
par(mar=c(5,2,3,0),cex=0.75, cex.main=1.75, cex.sub = 1.25, cex.axis= 1.25)
plot(c, main="Community structure computed with Euclidean distance \n on sqrt proportional ASVs", 
     sub="Healthy= green          Addicted (untreated and placebo) = red           Treated= orange")
dev.off()

library(dendextend)
{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))
  c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))), method = "bray"))
  c<-as.dendrogram(c)
  tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
  colnames(tabella_colore)<-c("Gruppo","Colore")
  colors <- gsub("Treated","coral",tabella_colore$Colore) # orange
  colors <- gsub("Healthy","chartreuse2",colors) #green
  colors <- gsub("Addicted","red4",colors)
  colors_to_use <- colors[order.dendrogram(c)]
  labels_colors(c) <- colors_to_use
}
png(file="Hierarchical cluster BRAY on normalized ASVs.png",width=3000,height=1800, res=300)
par(mar=c(5,2,3,0),cex=0.75, cex.main=1.75, cex.sub = 1.25, cex.axis= 1.25)
plot(c, main="Community structure computed with Bray-Curtis distance \n on sqrt proportional ASVs", 
     sub="Healthy= green          Addicted (untreated and placebo) = red           Treated= orange")
dev.off()

##################### ALFA DIVERSITY IN STOOL ############################

dir.create("Alpha diversity")
setwd("Alpha diversity")

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

# healthy vs cocain addicted (untreated)
{pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition")
  pAlpha$data<-dplyr::filter(pAlpha$data, Condition!="Treated")
  pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evanness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  New_data$Condition<-gsub("pre_rTMS","Addicted",New_data$Condition)
  New_data$Condition<-factor(New_data$Condition, levels = c("Healthy","Addicted"))
  pAlpha$data<-New_data
}
library(ggpubr)
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw(base_size = 12.2) + 
  labs(x="Condition", title="Alpha diversity in healthy and cocaine addicted (pre_rTMS and untreated) \n in stool samples")+
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=0.6, hjust=1, size=10)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4)
ggsave(file="Alfa diversity Healthy vs Addicted pre_rTMS.png", width = 8,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor) # OK

#############       now pre_rTMS vs post_rTMS (only samples with all pre_rTMS and post_rTMS time)

####################################################### placebo

rm(data.temp)
data.temp<-subset_samples(data.pair.placebo, Time!="during_rTMS") # only pre vs post
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
New_data<-pAlpha$data
# selecting only samples in pair with *all* two times 
{ pre_rTMS<-subset(New_data, Time=="pre_rTMS")
  post_rTMS<-subset(New_data, Time=="post_rTMS")
  pre_rTMS<-subset(pre_rTMS, Sample %in% post_rTMS$Sample)
  post_rTMS<-subset(post_rTMS, Sample %in% pre_rTMS$Sample)
  pre_rTMS<-pre_rTMS[order(pre_rTMS$Sample),]
  post_rTMS<-post_rTMS[order(post_rTMS$Sample),]
  New_data<-rbind.data.frame(pre_rTMS, post_rTMS)
  identical(New_data[New_data$Time=="pre_rTMS","Sample"],New_data[New_data$Time=="post_rTMS","Sample"] )
}
# TRUE --> same order --> updating
{pAlpha$data<-New_data
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  New_data$Time<-factor(New_data$Time, levels=c("pre_rTMS","post_rTMS"))
  pAlpha$data<-New_data
}
library(ggpubr)
pAlpha + theme_bw(base_size = 12.3) + 
  labs(x="Times", title="Alpha diversity in pre_rTMS and post_rTMS of placebo group \n in stool samples", caption = "the lines are connetting the same sample") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=10)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired=T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.3) +
  geom_line(aes(group=Sample), size=0.3, col="blue")
ggsave(file="Alfa diversity placebo pre_rTMS vs post_rTMS.png", width = 8,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Time
wilcox.test(Obser_value$value~factor, paired=T) # OK

rm(con, pre_rTMS, post_rTMS, Obser_value, New_data, alphadt, obs, H)

####################################################### treated

rm(data.temp)
data.temp<-subset_samples(data.pair.treated, Time!="during_rTMS") # only pre vs post
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
pAlpha<-plot_richness(data.pair.treated, measures=c("Shannon", "Observed"), x="Time")
New_data<-pAlpha$data
# selecting only samples in pair with *all* two times 
{ pre_rTMS<-subset(New_data, Time=="pre_rTMS")
  post_rTMS<-subset(New_data, Time=="post_rTMS")
  pre_rTMS<-subset(pre_rTMS, Sample %in% post_rTMS$Sample)
  post_rTMS<-subset(post_rTMS, Sample %in% pre_rTMS$Sample)
  pre_rTMS<-pre_rTMS[order(pre_rTMS$Sample),]
  post_rTMS<-post_rTMS[order(post_rTMS$Sample),]
  New_data<-rbind.data.frame(pre_rTMS, post_rTMS)
}
identical(New_data[New_data$Time=="pre_rTMS","Sample"],New_data[New_data$Time=="post_rTMS","Sample"] )
# TRUE --> same order --> updating
{pAlpha$data<-New_data
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  New_data$Time<-factor(New_data$Time, levels=c("pre_rTMS","post_rTMS"))
  pAlpha$data<-New_data
}
library(ggpubr)
pAlpha + theme_bw(base_size = 12.3) + 
  labs(x="Times", title="Alpha diversity between pre_rTMS and post_rTMS of treated group \n in stool samples", caption = "the lines are connetting the same sample") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=10)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired=T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.3) +
  geom_line(aes(group=Sample), size=0.3, col="blue")
ggsave(file="Alfa diversity treated pre_rTMS vs post_rTMS.png", width = 8,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Time
wilcox.test(Obser_value$value~factor, paired=T) # OK

rm(con, pre_rTMS, post_rTMS, Obser_value,  New_data, alphadt, obs, H)

setwd("..")

######################## BETA DIVERSITY BRAY CURTIS #######################

dir.create("Bray_Curtis_Beta_diversity")
setwd("Bray_Curtis_Beta_diversity")

# NB: only between pre and post
{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.prop<-ASV.prop[ , ! colnames(ASV.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.genus.prop<-ASV.genus.prop[ , ! colnames(ASV.genus.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.fam.prop<-ASV.fam.prop[ , ! colnames(ASV.fam.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.class.prop<-ASV.class.prop[, ! colnames(ASV.class.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.order.prop<-ASV.order.prop[, ! colnames(ASV.order.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
  ASV.phy.prop<-ASV.phy.prop[, ! colnames(ASV.phy.prop) %in% Sample_data[Sample_data$Time=="during_rTMS","ID_FASTQ"] ]
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")
metadata<-subset(metadata, Time != "during_rTMS")
metadata$Treatment_time<-paste(metadata$Treatment_type,metadata$Time, sep="_")
library(vegan)

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_ASV

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_g<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_g

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_f

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_o

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_c

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Treatment_time, data=metadata, permutations = 9999, method="bray")
perm_p

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

# pairwise beta diversity
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=metadata$Treatment_time, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV
pair_ASV$pairs<-gsub("Healthy_Healthy","Healthy", pair_ASV$pairs)
pair_ASV$pairs<-gsub("pre_rTMS","pre-rTMS", pair_ASV$pairs)
pair_ASV$pairs<-gsub("post_rTMS","post-rTMS", pair_ASV$pairs)

# exporting beta diversity
rm(con)
con<-file("Beta diversity general and pairwise between Healthy Untreated and T1_post_rTMS_Treated")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n", fill=TRUE)
pair_ASV
sink()
close()

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
disper<-vegan::betadisper(BC.dist,metadata$Treatment_time)
disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
disp_ASV$tab
a<-as.data.frame(disp_ASV$pairwise$permuted)
colnames(a)<-c("permuted_p_value")
a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
a$Significative<-a$padj_BH
a$Significative[a$padj_BH < 0.05]<- "*"
a$Significative[a$Significative > 0.05]<- ""
a
row.names(a)<-gsub("Healthy_Healthy","Healthy", row.names(a))
row.names(a)<-gsub("pre_rTMS","pre-rTMS", row.names(a))
row.names(a)<-gsub("post_rTMS","post-rTMS", row.names(a))

#export dispersion
rm(con)
con<-file("Beta dispersion General and Pairwise between Healthy Untreated and T1_post_rTMS_Treated")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n")
a
sink()
close()

rm(disp_ASV,a)

########################### PCoA BRAY CURTIS #####################

# NB: only pre and post treatment

# in base ad ASV
data.prop.labels<-subset_samples(data.prop, Time!="during_rTMS")
{sample_data(data.prop.labels)$Time<-gsub("pre_rTMS","Addicted",sample_data(data.prop.labels)$Time)
sample_data(data.prop.labels)$True_treated_Sign<-sample_data(data.prop.labels)$Under_TRUE_treatment
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Healthy","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Untreated","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Treated","*",sample_data(data.prop.labels)$True_treated_Sign)
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.005))+
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral")) +
  labs(title="PCoA computed with Bray-Curtis distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign")
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_Conditions.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + 
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.005))+
  labs(title="PCoA computed with Bray-Curtis distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign") +
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_no_ellipses.png", width = 8, height = 6, dpi=300)
# without post_rTMS
{data.temp<-subset_samples(data.sqrt_prop, Time!="post_rTMS")
data.temp<-transform_sample_counts(data.temp, function(x) x/sum(x))
DistBC = phyloseq::distance(data.temp, method = "bray")
ordBC = ordinate(data.temp, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.temp, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  labs(title="PCoA computed with Bray-Curtis distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance")) +
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4"))
ggsave(file="PCoA_Beta_diversity_Bray_Curtis_solo healthy vs addicted.png", width = 8, height = 6, dpi=300)
setwd("..")

######################## BETA DIVERSITY WEIGHTED UNIFRAC #######################

rm(ASV.prop, ASV.genus.prop, perm_ASV, perm_g, perm_f)

dir.create("Weighted Unifrac Beta diversity")
setwd("Weighted Unifrac Beta diversity")

# NB: only between pre and post

#### PERMANOVA
{metadata<-as(sample_data(data.prop),"data.frame")
  metadata<-subset(metadata, Time != "during_rTMS")
  metadata$Treatment_time<-paste(metadata$Treatment_type,metadata$Time, sep="_")
  library(vegan)
}

{rm(data.temp)
  data.temp<-subset_samples(data.prop, Time != "during_rTMS")
  data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
  DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac") # needing phyloseq to compute unifrac
  perm_ASV<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
  perm_ASV
}
{rm(data.temp)
  data.temp<-subset_samples(data.genus.prop, Time != "during_rTMS")
  data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
  DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
  perm_g<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
  perm_g
}
{rm(data.temp)
  data.temp<-subset_samples(data.fam.prop, Time != "during_rTMS")
  data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
  DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
  perm_f<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
  perm_f
}
{rm(data.temp)
  data.temp<-subset_samples(data.order.prop, Time != "during_rTMS")
  data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
  DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
  perm_o<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
  perm_o
}
{rm(data.temp)
  data.temp<-subset_samples(data.class.prop, Time != "during_rTMS")
  data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
  DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
  perm_c<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
  perm_c
}
{rm(data.temp)
  data.temp<-subset_samples(data.phy.prop, Time != "during_rTMS")
  data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
  DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
  perm_p<- vegan::adonis(DistBC ~Treatment_time, data=metadata, permutations = 9999)
  perm_p
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

# pairwise beta diversity
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# on genera because on ASV level there is to much noise
rm(data.temp)
data.temp<-subset_samples(data.genus.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
pair_g<- pairwise.adonis(DistBC,  factors=metadata$Treatment_time, p.adjust.m = "BH", sim.method="wunifrac", perm = 9999)
pair_g
pair_g$pairs<-gsub("Healthy_Healthy","Healthy", pair_g$pairs)
pair_g$pairs<-gsub("pre_rTMS","pre-rTMS", pair_g$pairs)
pair_g$pairs<-gsub("post_rTMS","post-rTMS", pair_g$pairs)

# exporting beta diversity
rm(con)
con<-file("Beta diversity general and pairwise between Healthy Addicted and Untreated VS post_rTMS_")
sink(con, append = TRUE)
cat("General beta diversity Weighted Unifrac on Genera \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Weighted Unifrac (correction Benjamini-Hochberg on p-value) on Genera \n", fill=TRUE)
pair_ASV
sink()
close()

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
rm(data.temp)
data.temp<-subset_samples(data.genus.prop, Time != "during_rTMS")
data.sqrt_prop_perm<-transform_sample_counts(data.temp, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
disper<-vegan::betadisper(DistBC,metadata$Treatment_time)
disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
disp_ASV$tab
a<-as.data.frame(disp_ASV$pairwise$permuted)
colnames(a)<-c("permuted_p_value")
a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
a$Significative<-a$padj_BH
a$Significative[a$padj_BH < 0.05]<- "*"
a$Significative[a$Significative > 0.05]<- ""
a
row.names(a)<-gsub("Healthy_Healthy","Healthy", row.names(a))
row.names(a)<-gsub("pre_rTMS","pre-rTMS", row.names(a))
row.names(a)<-gsub("post_rTMS","post-rTMS", row.names(a))

#export dispersion
rm(con)
con<-file("Beta dispersion General and Pairwise between Healthy Addicted and Untreated vs post_rTMS")
sink(con, append=TRUE)
cat("General beta dispersion Weighted Unifrac \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) on Genera \n \n")
a
sink()
close()

rm(disp_ASV,a)

########################### PCoA WEIGHTED UNIFRAC #####################

# in base ad ASV
{ data.prop.labels<-subset_samples(data.prop, Time!="during_rTMS")
sample_data(data.prop.labels)$Time<-gsub("pre_rTMS","Addicted",sample_data(data.prop.labels)$Time)
sample_data(data.prop.labels)$True_treated_Sign<-sample_data(data.prop.labels)$Under_TRUE_treatment
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Healthy","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Untreated","",sample_data(data.prop.labels)$True_treated_Sign)
sample_data(data.prop.labels)$True_treated_Sign<-gsub("Treated","*",sample_data(data.prop.labels)$True_treated_Sign)
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "wunifrac")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.0001))+
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral")) +
  labs(title="PCoA computed with Weighted Unifrac distance \n on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign")
ggsave(file="PCoA_Beta_diversity_W_Unifrac_on_ASV.png", width = 8, height = 6, dpi=300)

# it's much better on genera (less noise)
{data.prop.labels<-subset_samples(data.genus.prop, Time!="during_rTMS")
  sample_data(data.prop.labels)$Time<-gsub("pre_rTMS","Addicted",sample_data(data.prop.labels)$Time)
  sample_data(data.prop.labels)$True_treated_Sign<-sample_data(data.prop.labels)$Under_TRUE_treatment
  sample_data(data.prop.labels)$True_treated_Sign<-gsub("Healthy","",sample_data(data.prop.labels)$True_treated_Sign)
  sample_data(data.prop.labels)$True_treated_Sign<-gsub("Untreated","",sample_data(data.prop.labels)$True_treated_Sign)
  sample_data(data.prop.labels)$True_treated_Sign<-gsub("Treated","*",sample_data(data.prop.labels)$True_treated_Sign)
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "wunifrac")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + stat_ellipse() +
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.003))+
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral")) +
  labs(title="PCoA computed with Weighted Unifrac distance \n on sqrt proportional genera", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign")
ggsave(file="PCoA_Beta_diversity_W_Unifrac_sui generi.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_point(size=3) + theme_classic() + 
  geom_text(aes(label=sample_data(data.prop.labels)$True_treated_Sign), size=5, col="white", position = position_nudge(0,-0.001))+
  labs(title="PCoA computed with Weighted Unifrac distance \n on sqrt proportional genera", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"),
       caption="the patients under true treatment are marked with the * sign") +
  scale_color_manual(values = c("Healthy"="chartreuse2", "Addicted"="red4", "post_rTMS"="coral"))
ggsave(file="PCoA_Beta_diversity_W_Unifrac_no_ellipses_generi.png", width = 8, height = 6, dpi=300)

setwd("..")

############## DA WITH DESEQ2 - Healthy vs Untreated Stool #################

dir.create("DESeq2")
setwd("DESeq2")

# Trimming under sum of 10, see DESeq2  tutorial
{data_pruned<- prune_taxa(taxa_sums(data) > 10, data)
sample_data(data_pruned)$Condition<-gsub("pre-rTMS","Addicted",sample_data(data_pruned)$Condition)
# following are needed in case of different ASV names due to pruning
data_pruned.phy = tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
data_pruned.class = tax_glom(data_pruned, taxrank = "Class", NArm = F)
data_pruned.order = tax_glom(data_pruned, taxrank = "Order", NArm = F)
data_pruned.fam = tax_glom(data_pruned, taxrank = "Family", NArm = F)
data_pruned.genus = tax_glom(data_pruned, taxrank = "Genus", NArm = F)
data_pruned.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data_pruned.phy.prop <- transform_sample_counts(data_pruned.phy, function(ASV) ASV/sum(ASV)*100)
data_pruned.class.prop <- transform_sample_counts(data_pruned.class, function(ASV) ASV/sum(ASV)*100)
data_pruned.order.prop <- transform_sample_counts(data_pruned.order, function(ASV) ASV/sum(ASV)*100)
data_pruned.fam.prop <- transform_sample_counts(data_pruned.fam, function(ASV) ASV/sum(ASV)*100)
data_pruned.genus.prop <- transform_sample_counts(data_pruned.genus, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured taxa
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
      d<-phyloseq_to_deseq2(data_pruned, ~Condition)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Condition)
    }
    # DE<-DESeq(d, test="LRT",reduced= ~ 1)
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Condition", nume[i], deno[i]), test="Wald")
      #res<-results(DE, contrast=c("Condition", nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      #res<-lfcShrink(DE, type="normal", contrast=c("Condition",nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        #			cat(paste0(rank,":",nume[i],"/",dena[i],"\n"),file=outfile,append=T)
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("Healthy")        #("nume" e "deno" choose who confront)
deno<-c("Addicted")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100) # disable scientific annotation
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)    #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
write.csv2(DE_results, file="DE results_Healthy vs Untreated.csv", quote=F, row.names = F)
unlink("DE_results.txt")

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
target<-prune_taxa(target, data_pruned.genus.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_g<-tabella
}
{target<-DE_results[DE_results$Rank=="Family","ASV"]
target<-prune_taxa(target, data_pruned.fam.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_f<-tabella
}
{target<-DE_results[DE_results$Rank=="Order","ASV"]
target<-prune_taxa(target, data_pruned.order.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_o<-tabella
}
{target<-DE_results[DE_results$Rank=="Class","ASV"]
target<-prune_taxa(target, data_pruned.class.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_c<-tabella
}
{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
target<-prune_taxa(target, data_pruned.phy.prop)
tabella<-psmelt(target)
tabella$Abundance<-sqrt(tabella$Abundance)
tabella_p<-tabella
}
# unique boxplot of DESeq2  results
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
tabella_tot<-tabella_tot[tabella_tot$Condition!="post_rTMS",]
tabella_tot<-tabella_tot[tabella_tot$Under_TRUE_treatment!="Treated",]
tabella_tot$Condition<-factor(tabella_tot$Condition, levels=c("Healthy","Addicted"))

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space = "free") +
  geom_boxplot(width=0.8) + theme_bw( ) +
  scale_fill_manual(values = c("Healthy"="Chartreuse","Addicted"="red4","post_rTMS"="coral")) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,axis.text.x = element_text(angle = 90, vjust=0.3, hjust=1, size=10), axis.text.y = element_text(size=10)) + 
  scale_x_discrete(expand=c(-0.2, 1.5)) + theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between Healthy subjects and Untreated addicted", y="Sqrt Proportional Abundance", fill="Condition", x="") + scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))
ggsave(filename = "Total_differences_DESeq2 Healthy vs Untreated.png", width = 15, height = 8, dpi=300)
dev.off()

# if there are reduntants ( same ASV, same abundance, same result! ), manteining only the lower level
Redund<-DE_results
Redund<-Redund[duplicated(round(Redund$log2FoldChange, digits = 0), fromLast=T), ]
Redund<-as.character(c(Redund[Redund$Rank=="Phylum","Phylum"],Redund[Redund$Rank=="Class","Class"],Redund[Redund$Rank=="Order","Order"],Redund[Redund$Rank=="Family","Family"]))

tabella_tot2<-subset(tabella_tot, ! Bacteria %in% Redund)
tabella_tot2<-subset(tabella_tot2, ! Bacteria %in% "Bacteroidaceae") # redundant too

tabella_tot2[tabella_tot2$Bacteria=="Incertae_Sedis","Bacteria"]<-"Incertae_Sedis _ p Firmicutes"
# the following are almost absent!
tabella_tot2<-subset(tabella_tot2, ! Bacteria %in% c("Pseudomonadales","Hafniaceae","Hafnia-Obesumbacterium")) 
# name too long, shortening it
tabella_tot2[tabella_tot2$Bacteria=="uncultured_ f Desulfovibrionaceae", "Bacteria"] <- "uncultured_ f Desulfovib."
tabella_tot2[tabella_tot2$Bacteria=="uncultured_ f Christensenellaceae", "Bacteria"] <- "uncultured_ f Christensen."

ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) +
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space = "free") +
  geom_boxplot(width=0.8) + theme_bw( base_size= 16) + 
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  scale_fill_manual(values = c("Healthy"="Chartreuse","Addicted"="red4")) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = -40, vjust=0.8, hjust=0, size=12), axis.text.y = element_text(size=13)) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  theme(plot.title= element_text(size=18) ,
        legend.key.size=unit(0.7,"cm"), legend.text =element_text(size=18), legend.title = element_text(size=16)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant taxa between Healthy and Addicted (placebo and untreated) patients", y="Sqrt Proportional Abundance", fill="Condition", x="") +
  scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  theme(plot.margin=margin(1,3.25,1,0.3, unit = "cm"))
ggsave(filename = "Total_differences_DESeq2 Healthy vs Untreated_no redundant.png", width = 15.5, height = 8, dpi=300)

dev.off()

setwd("..")

############## DA WITH DESEQ2 - PLACEBO PAIRS pre_rTMS vs post_rTMS Stool #################

setwd("DESeq2")

# Trimming under sum of 10, see DESeq2  tutorial
data_pruned<- prune_taxa(taxa_sums(data.pair.placebo) > 10, data.pair.placebo)

# function to automatize the comparison between groups /// NB: ~Sample + Time
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Sample + Time)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Time)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Time", nume[i], deno[i]), test="Wald")
      #res<-results(DE, contrast=c("Time", nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      #res<-lfcShrink(DE, type="normal", contrast=c("Time",nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind("num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,], "matrix")))
        #			cat(paste0(rank,":",nume[i],"/",dena[i],"\n"),file=outfile,append=T)
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("pre_rTMS")        #("nume" e "deno" choose who confront)
deno<-c("post_rTMS")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100) # disable scientific annotation
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)             #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
write.csv2(DE_results, file="DE results_Placebo three times.csv", quote=F, row.names = F)
unlink("DE_results.txt")
# nothing! 

setwd("..")

############## DA WITH DESEQ2 - treated PAIRS pre_rTMS vs post_rTMS Stool #################

setwd("DESeq2")

# Trimming under sum of 10, see DESeq2  tutorial
data_pruned<- prune_taxa(taxa_sums(data.pair.treated) > 10, data.pair.treated)

# function to automatize the comparison between groups /// NB: ~Sample + Time
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Sample + Time)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Time)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Time", nume[i], deno[i]), test="Wald")
      #res<-results(DE, contrast=c("Time", nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      #res<-lfcShrink(DE, type="normal", contrast=c("Time",nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind("num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,], "matrix")))
        #			cat(paste0(rank,":",nume[i],"/",dena[i],"\n"),file=outfile,append=T)
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("pre_rTMS")        #("nume" e "deno" choose who confront)
deno<-c("post_rTMS")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100) # disable scientific annotation
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)             #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-1]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus")
write.csv2(DE_results, file="DE_results_true treated three times.csv", quote=F, row.names = F)
unlink("DE_results.txt")

# nothing again!

setwd("..")

################# PLS-DA and sPLS-DA on Stool ###########################

setwd("/home/Scrivania/Results/Stool")
dir.create("PLS-DA")
setwd("PLS-DA")

# BiocManager::install('mixOmics')
library("mixOmics")
data.no.treated<-subset_samples(data.genus.prop, Condition != "Treated")
metadata<-as(sample_data(data.no.treated),"data.frame")

X <- as.matrix(t(otu_table(data.no.treated)))
Y <- as.factor(metadata$Condition)

stool.splsda <- splsda(X, Y, ncomp = 10)

perf.splsda <- perf(stool.splsda, validation = "Mfold", # to assest the optimal number of comp
                    folds = 10, nrepeat = 4*(dim(X)[1]), cpus = 6,
                    progressBar = FALSE, auc = TRUE)
comp<-perf.splsda$choice.ncomp
comp #4

png(filename = "Error for each component PLS-DA", width = 2500, height = 1800, res=300)
plot(perf.splsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
title(main="Classification rate error \n for each component used to compute PLS-DA")
dev.off()

stool.splsda <- splsda(X, Y, ncomp = comp[1]) 

n= 2 # choose the component to use depending on clustering efficiency in plot
png(filename = paste("PLS-DA stool_comp 1 and", n), width = 2500, height = 1500, res=300)
plotIndiv(stool.splsda, comp = c(1,n),
          group = metadata$Condition, ind.names = TRUE, ellipse = TRUE, legend = TRUE, 
          title = paste0("PLSDA between Healthy and Untreated addicted \n on scaled proportional genus data \n (computed with ", comp[1], " components, plotted on comp 1 and ", n,")"))
dev.off()

################ sPLS-DA for variable selection

# cutting out same sample for test (otherwise clear overfitting during the test)
set.seed(1)
train <- sample(1:nrow(X), 50)
train
test <- setdiff(1:nrow(X), train)

X.test <- X[test,]
X <- X[train, ]
Y.test <- Y[test]
Y<- Y[train]
metadata<-metadata[row.names(X),]

tune.splsda.stool <- tune.splsda(X, Y, ncomp = comp[1], validation = 'Mfold',
                                 folds = 10, nrepeat = 4*(dim(X)[1]), cpus = 6, # use repeated cross-validation
                                 dist = 'max.dist',
                                 test.keepX =  seq(10,round(dim(X)[2]/5),10), # testing from 10 to 1/5 MAX number of variables
                                 measure = "BER") # use balanced error rate of dist measure

optimal.comp <- tune.splsda.stool$choice.ncomp$ncomp
optimal.comp # 2
optimal.keepX <- tune.splsda.stool$choice.keepX[1:optimal.comp]
optimal.keepX # 60, 10

final.splsda<-splsda(X,Y, ncomp=optimal.comp, keepX = optimal.keepX)
n=2 # choose which component plot depending on plot
png(filename = paste("sparse PLS-DA stool_comp 1 and",n), width = 2500, height = 1500, res=300)
plotIndiv(final.splsda, comp = c(1,n), group = metadata$Condition, ind.names = TRUE, legend = TRUE, ellipse = TRUE,
          title = paste0("sparse PLSDA between Healthy and Untreated addicted \n on scaled proportional genus data \n (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"))
dev.off()

################ plotting loadings

loadings<-as.data.frame(final.splsda$loadings$X)
identical(row.names(Taxa.genus[row.names(loadings),]),row.names(loadings))
loadings$Genus<-Taxa.genus[row.names(loadings),"Genus"]

a<-loadings[loadings$comp1!=0,]
a$comp1<-as.numeric(a$comp1)
a<-a[order(abs(a$comp1)),]
a$Genus <- factor(a$Genus, levels = a$Genus) # otherwise it would be re-ordered in plot
ggplot(mapping=aes(x=a$Genus,y=a$comp1)) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, size = 11, hjust = 1, vjust = 1)) +
  labs(x="Selected genera on component 1", 
       y="loadings", title = "Loadings of selected genera for component 1 in sPLSDA", 
       caption="\n dotted lines are plot at 50% of both min and max loadings value") +
  geom_hline(yintercept = max(a$comp1)/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  geom_hline(yintercept = min(a$comp1)/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = "Loadings of choosen genera for sPLSDA comp 1.png", width = 14, height = 8, dpi=300)

b<-loadings[loadings[[n]]!=0,]
tax_b<-Taxa.genus
tax_b$ASV<- row.names(tax_b)
tax_b<-subset(tax_b, ASV %in% row.names(b))
tax_b<-tax_b[row.names(b),]
identical(tax_b$ASV,row.names(b))
b$Genus<-tax_b$Genus

b[,n]<-as.numeric(b[,n])
b<-b[order( abs(b[,n])) ,]
b$Genus <- factor(b$Genus, levels = b$Genus)
ggplot(mapping=aes(x=b$Genus,y=b[,n])) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, size = 11, hjust = 1, vjust = 1)) +
  labs(x=paste("Selected genera on component",n), 
       y="loadings", title = paste("Loadings of selected genera for component",n,"in sPLSDA"), 
       caption="\n dotted lines are plot at 50% of both min and max loadings value") +
  geom_hline(yintercept = max(b[,n])/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  geom_hline(yintercept = min(b[,n])/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = paste("Loadings of choosen genera for sPLSDA comp",n,".png"), width = 14, height = 8, dpi=300)

############### testing the sPLS-DA model

predict.splsda <- predict(final.splsda, newdata =as.matrix(X.test), dist = "all")

# evaluating the prediction accuracy
predict.comp <- predict.splsda$class$centroids.dist[,n]
predict.comp
table(factor(predict.comp, levels = levels(Y)), as.factor(Y.test))
# now evaluating the prediction accuracy using ONLY the first component
predict.comp1 <- predict.splsda$class$centroids.dist[,1]
predict.comp1
table(factor(predict.comp1, levels = levels(Y)), as.factor(Y.test))

##### exporting impostations and values of sPLSDA

con<-file("Results_and_settings_sPLSDA")
sink(con)
cat("Number of training samples and their genera", fill = TRUE)
cat(dim(X), fill = TRUE)
cat("\n Samples casually selected as test (and then discarded from training data set)", fill=TRUE)
cat(row.names(X.test), fill = TRUE)
cat("\n number of components originally selected by perf function for normal PLSDA", fill=TRUE)
cat(comp[1], fill=TRUE)
cat("\n number of components selected by tune.splsda function for sPLSDA", fill=TRUE)
cat(optimal.comp, fill = TRUE)
cat("\n number of genera selected by tune.splsda function for each component", fill=TRUE)
cat(optimal.keepX, fill=TRUE)
cat("\n Possible number of genera that could be selected by function tune.splsda", fill=TRUE)
cat(seq(10,round(dim(X)[2]/5),10), fill = TRUE)
cat("\n \n \n ### testing the prediction efficiency of sPLSDA through confusion matrix using ONLY the first component \n", fill=TRUE)
table(factor(predict.comp1, levels = levels(Y)), Y.test)
cat("\n", length(which(predict.comp1!="")),"on", length(predict.comp1), "samples predicted" )
cat("\n \n ### testing the prediction efficiency of sPLSDA through confusion matrix component 1 and",n,"\n", fill=TRUE)
table(factor(predict.comp, levels = levels(Y)), Y.test)
cat("\n", length(which(predict.comp!="")),"on", length(predict.comp), "samples predicted" )
cat("\n \n type of distance used for predictions: centroid distance")
sink()

rm(a,b,con,tax_a,tax_b,loadings,metadata,predict.comp,predict.comp1,train,X,Y,X.test,Y.test,final.splsda,tune.splsda.stool,predict.splsda,comp,n,test,optimal.comp,optimal.keepX)

setwd("..")

################## HANDLING PICRUST AND LEFSE OUTPUT ################

a <- read.delim("path_abun_unstrat_descrip.tsv.gz") # PICRUTS2 output
colnames(a)<-gsub("ID2117.[0-9][0-9][0-9].","",colnames(a))
colnames(a)<-gsub("ID2117.[0-9][0-9].","",colnames(a))
colnames(a)<-gsub("ID2117.[0-9].","",colnames(a))
colnames(a)<-gsub(".A[0-9].[A-Z][0-9][0-9]","",colnames(a))
colnames(a)<-gsub("ID1548.16S.[A-Z][0-9][0-9].[0-9][0-9].","",colnames(a))
colnames(a)<-gsub("ID1642.16S.[A-Z][0-9][0-9].[0-9][0-9].","",colnames(a))
colnames(a)
Descriptions<-a[,1:2]
row.names(Descriptions)<-Descriptions[,1]
a[,2]<-NULL # no descriptions in LefSE
metadata<-as(sample_data(data),"data.frame")
identical(metadata$ID_FASTQ,colnames(a)[-1]) # TRUE
row.names(a)<-a[,1]
a[,1]<-NULL
a<-rbind.data.frame(metadata$Condition,a)
head(a)
row.names(a)[1]<-"Condition"
a<-a[,a[1,]!="Treated"] # filtered out "treated" sample to study healthy vs untreated only

write.table(a,file="Output_KO_ready_for_LefSe_online.tsv",row.names = T,quote = F, sep="\t")
# remember to add a rownames to ex colnames

# --> LEFSe
# <--
Significative_functions_LEFSE<- read.delim("Result.lefse_internal_res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)
# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
Significative_functions_LEFSE<-dplyr::left_join(Significative_functions_LEFSE,Descriptions, by=c("Pathway"="pathway"))
head(Significative_functions_LEFSE)

write.csv2(Significative_functions_LEFSE,file = "Significative_functions_LEFSE.csv",na="",quote = F, row.names = F)

##################### R and PACKAGES VERSION #########################

setwd("/home/Desktop")

package<-sessionInfo()

con <- file("R_and_packages.txt")
sink(con, append = TRUE) 
sink(con, append=TRUE, type="message")

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

sink()
sink(type="message") # restore STR ERROR to R console
close(con)