##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
library("ggplot2")
library("vegan")
library("ggpubr")
library("ggh4x")
library("egg")
library("dendextend")
library("DESeq2")
library("xlsx")  
library("dplyr")
library("qiime2R")
library("Hmisc")
}

{dir.create("Data_check")
dir.create("Data_check/PCoA_test")
dir.create("Results")
dir.create("Results/Abundances")
dir.create("Results/Abundances/Ratio_Firmi_Bacteroi_Healthy_vs_AMD")
dir.create("Results/Abundances/Ratio_Firmi_Bacteroi_T0_vs_T1")
dir.create("Results/Beta_div")
dir.create("Results/Beta_div/Healthy_vs_AMD")
dir.create("Results/Beta_div/T0_vs_T1")
dir.create("Results/DA_DESeq2/")
dir.create("Results/Tests_and_Correlations")
}

options(scipen = 100) # disable scientific annotation


####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample

Metadata<-as.data.frame(read.table(file="Metadata.tsv", header = T))

row.names(Metadata)<-Metadata$FASTQ_Code # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,sample)

sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels=c("Healthy","AMD","Post_treatment")) # decides the order in plots
sample_names(data)<-sample_data(data)$Sample_name
head(sample_data(data))
row.names(Metadata)<-Metadata$Sample_name

data_all<-data  # preparing for the subsets
rm(data)


########## IMPORTING ALSO THE CONTAMINATED ONE (FROM HOST DNA)

data_contam<-qza_to_phyloseq(features="QIIME/table.qza", tree = "QIIME/eucar_contaminants/rooted-tree_with_contam.qza")
# changing names (already done before, check up there!)
if(identical(original_names,sample_names(data_contam))){
  sample_names(data_contam)<-sample_names(data_all)
  }
if(identical(sample_names(data_contam), row.names(Metadata))) {
  sample_data(data_contam)<-Metadata
  cat("Everything alright there!")} else { cat("\n Sample names are not matching, check them \n\n") }


# save.image("data_after_import.RData")


#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
unfiltered_data<-data_all
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under the average of 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_average_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured"] # to avoid the removal of other uncultured genera
filtered<-filtered[! is.na(filtered)] # to avoid the removal of other NA genera
data_all<-subset_taxa(data_all, ! Genus %in% filtered)
data_all<-subset_taxa(data_all, ! is.na(Phylum)) # manually checked, there is a really low abundance of NA phylum is some sample, potential off-target (Bayesian Classification)
data_all<-prune_taxa(taxa_sums(data_all) > 0, data_all) 
rm(filtered, data.genus.temp)


proof1<- "Marker of the filtering, it is required for the script"


################# CHECKING THE COMPOSITION AFTER FILTERING ####################

data.unf.prop<-transform_sample_counts(data_contam, function(x) (x/sum(x))*100)


if(! "proof1" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}
data.prop<-transform_sample_counts(data_all, function(x) (x/sum(x))*100)


################# BRAY HELLINGER

##### unfiltered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.unf.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


################# sqrt BRAY PROPORTIONAL AB

##### unfiltered sqrt BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.unf.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "bray")
  DistBC = sqrt(DistBC)
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.prop.labels, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA sqrt Bray-Curtis (on proportional ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.prop.labels))
data.prop.labels<-data.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "bray")
  DistBC = sqrt(DistBC)
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_sqrt_BRAY_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


################# EUCLIDEAN HELLINGER

##### unfiltered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
data.prop.labels<-data.unf.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))



################# wUNIFRAC PROP

##### unfiltered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
data.prop.labels<-data.unf.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.prop.labels, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Condition") +
  scale_color_manual(values=c("AMD"="red3","Healthy"="chartreuse","Post_treatment"="coral")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_name), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


############################ RAREFACTION ANALYSIS ################################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as generated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of the original samples (and in the same order!!!)
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    check<-"red3"
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

png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(otu_table(data_all)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)



############# PREPARATION OF THE DATA FOR A GENERAL VIEW #######################

{data.phy = tax_glom(data_all, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data_all, taxrank = "Class", NArm = F)
data.order = tax_glom(data_all, taxrank = "Order", NArm = F)
data.fam = tax_glom(data_all, taxrank = "Family", NArm = F)
data.genus = tax_glom(data_all, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data_all, function(ASV) ASV/sum(ASV)*100)
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
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_AFTER_filters.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)


######################## GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(unfiltered_data)==1, unfiltered_data)
length(which(taxa_sums(unfiltered_data)==1)) # if zero there are no singletons
{n<-as.data.frame(otu_table(filter))
  N<-as.data.frame(otu_table(unfiltered_data))
  G<-1-(colSums(n)/colSums(N))
}
con<-file("Data_check/Percentuale_di_singletons_Good's_coverage.txt")
sink(con)
cat("GOOD'S COVERAGE ESTIMATOR \n", fill=TRUE)
cat("1-(n/N) for each sample, where n is number of singletons and N is Total ASV \n \n", fill=TRUE)
G
sink()
close(con)

rm(con, filter, G, n, N)


########################### COUNTS EXPORT ##########################################

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data_all),"matrix"),as(tax_table(data_all),"matrix")),file="Results/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
  }

options(scipen = 100)
dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT (GENERAL VIEW) ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one

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
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Condition<-gsub("_"," ",tabella$Condition)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =13) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/TOP_5_phyla.png",width=10,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Abundances/TOP5_phyla_average_abundance_EVERY_SAMPLE.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)

# TOP 5 Genera
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
tabella$Condition<-gsub("_"," ",tabella$Condition)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =13) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Abundances/TOP_5_genera.png",width=10,height=5,dpi=300)
dev.off()

rm(top, tabella, tabella_top)

# TOP 8 Genera
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
tabella$Condition<-gsub("_"," ",tabella$Condition)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =13) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.35, "cm"),legend.text = element_text ( size = 11 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Abundances/TOP_8_genera.png",width=10,height=5.1,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


#################### \\\\\\\\ STARTING THE ANALYSIS OF HEALTHY VS AMD \\\\\\\\ ######################

# Healthy controls VS AMDs (T0, both patients with paired T1 and not)

if(! "proof1" %in% ls()){
  cat("\nWait! Did you performed the filtering steps???\n\n")
  Sys.sleep(2)
}

remove<-ls()[!ls() %in% c("proof1","data_all","data_contam","Metadata")] # to reset the environment
rm(list = remove)

data<-subset_samples(data_all, Condition != "Post_treatment") # no T1, they are the same sample of some T0
data<-prune_taxa(taxa_sums(data)>0, data)

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


#################### RATIO FIRMICUTES/BACTEROIDES (HEALTHY VS AMD) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

ratio_fb<-ratio_fb[sample_data(data)$Sample_name, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb,sample_data(data)[,c("Sample_ID","Condition")])

# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Condition, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Condition=="AMD","Ratios_Mean"]<-Ratios_Mean["AMD"]
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_Mean"]<-Ratios_Mean["Healthy"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Condition, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Condition=="AMD","Ratios_st_err"]<-Ratios_st["AMD"]/sqrt(length(which(ratio_fb$Condition=="AMD")))
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_st_err"]<-Ratios_st["Healthy"]/sqrt(length(which(ratio_fb$Condition=="Healthy")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/Abundances/Ratio_Firmi_Bacteroi_Healthy_vs_AMD/Firmi_Bacter_Ratio.csv", ratio_fb)

ratio_fb$Condition<-factor(ratio_fb$Condition, levels = c("AMD","Healthy"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-wilcox.test(Ratio ~ Condition, paired=F, data=ratio_fb)
p_val<-round(res$p.value, digits = 4)
p_val
# exporting the results
con <- file("Results/Abundances/Ratio_Firmi_Bacteroi_Healthy_vs_AMD/Mann_Whi_Wilcoxon.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Condition:   V=",res$statistic, "p-value=", p_val, "\n",fill=T)
cat("\n Shapiro Wilk test p-value: ",distr$p.value)
sink()
close(con)


######## BAR PLOT
ratio_fb_plot<-ratio_fb
ratio_fb_plot[ratio_fb_plot$Ratio>100, "Ratio"]<-100    # almost no Bacteroidetes --> ratio too high --> setting a limit
ggplot(data=ratio_fb_plot, aes(x=Sample_ID, fill=Condition, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Healthy"="chartreuse",
                               "AMD"="coral")) +
  facet_grid2(.~Condition, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Condition"))+
  scale_y_sqrt(breaks=c(1,5,seq(10,100,10)) ) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Mann Whitney p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi_Healthy_vs_AMD/Ratio_barplot.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means)
set.seed(2)
ggplot(data=ratio_fb_plot, aes(x=Condition, color=Condition,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("AMD"="coral","Healthy"="chartreuse")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(1,5,seq(10,100,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= ifelse(Ratios_Mean + Ratios_st_err > 100, # standard deviation, the if else is needed to avoid geom bar below zero
                                 100, Ratios_Mean + Ratios_st_err),
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0,
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", subtitle = paste("Mann Whitney p-value:",p_val),
       caption="The upper limit of the Ratio is manually setted to 100\n because the absence of Bacteroidetes causes a really high ratio")
ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi_Healthy_vs_AMD/Ratios_jitter_with_mean_and_STerr.png",width=4,height=3.5, dpi=300) 
dev.off()

suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY (Healthy vs AMD) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition", color="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating
New_data<-rbind.data.frame(obs,H,ev)
pAlpha$data<-New_data
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha= 0) +
  scale_color_manual(values = c("Healthy"="chartreuse", "AMD"="coral")) +
  theme_classic2() +
  labs(x="Condition") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=11),
        panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 0.8,
                     size=3.2, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="Results/Alfa_diversity_HEALTHY_vs_AMD_with_Mann_W.png", width = 6,height =5, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


######################## BETA DIVERSITY  (Healthy VS AMD) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[1]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[1]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Beta_div/Healthy_vs_AMD/Beta_diversity_permanova_Helling.csv",quote=F,row.names = T)


### PLUS: checking it with Bray too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="bray")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Beta_div/Healthy_vs_AMD/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Condition)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Beta_div/Healthy_vs_AMD/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS (Healthy vs AMD) #####################

# on ASV
data.prop.labels<-data.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","AMD"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/Healthy_vs_AMD/PCoA_Beta_diversity_Hellinger_on_ASV.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","AMD"="coral")) +
  geom_point(size=3, alpha=0.5) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  labs(color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H),
       caption = "NB: significative but there is really low the variance along the axes")
ggsave(file="Results/Beta_div/Healthy_vs_AMD/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 8, height = 6, dpi=300)

rm(data.sqrt_prop, eigval, DistBC, ordBC,data.prop.labels)


# again but on genera
{data.prop.labels<-data.genus.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","AMD"="coral")) +
  geom_point(size=3, alpha=0.5) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  labs(color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Beta_div/Healthy_vs_AMD/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 6, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



################### DA WITH DESEQ2 (Healthy vs AMD) #####################

if(! "proof1" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Healthy", "AMD"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res<-res[res$baseMean > 50, ] # arbitrary threshold to avoid the most noisy result
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Healthy_vs_AMD.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2_Healthy_vs_AMD.csv", row.names = F)
write.xlsx(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2_Healthy_vs_AMD.xlsx", showNA = F, col.names = T)

Table_tot$Bacteria<-gsub("[","",Table_tot$Bacteria,fixed=T)
Table_tot$Bacteria<-gsub("]","",Table_tot$Bacteria,fixed=T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
  geom_boxplot(width=0.9, size= 0.25, alpha= 0.21, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.6),
             aes(color=Condition), size= 0.5, alpha= 0.95) +
  scale_fill_manual(values=c("AMD"="coral","Healthy"="chartreuse")) +
  scale_color_manual(values=c("AMD"="coral2","Healthy"="chartreuse3")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme_classic2(base_size = 11.2) + 
  theme(strip.text.x=element_text(size=11.4,colour="black"),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 28, vjust=1, hjust=1, size=9.5), 
        axis.text.y = element_text(size=9),
        plot.margin = margin(2,2,5,2),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray")
  ) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,2.5,5,7.5, seq(5,max(Table_tot$Abundance+3),5))) +
  labs(y="Percent abundance", 
       fill="Condition", x="")
ggsave(filename = "Results/DA_DESeq2/DA_Healthy_vs_AMD_every_result.png", width = 10.2, height = 6.8, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/DA_DESeq2/NB.txt ")



############################ PICRUST2 - LEFSE (HEALTHY vs AMD) ##################################

dir.create("Results/PICRUST2_LEFSE")

suppressWarnings(rm(a))
a <- read.delim("QIIME/PICRUST2_LEFSE/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("QIIME/PICRUST2_LEFSE/Result_LEFSE.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="Results/PICRUST2_LEFSE_results/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)


# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="AMD"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="AMD"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Healthy",1,0), size=3.2) +
  scale_fill_manual(values=c("AMD"="coral", "Healthy"="deepskyblue")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 18),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_every_result.png", width = 11.81, height = 9.5, dpi = 300)

system('touch Results/PICRUST2_LEFSE_results/no_result_over_3')

rm(Significative_functions_LEFSE, a)


#################### \\\\\\\\ STARTING THE ANALYSIS OF TREATMENT \\\\\\\\ ######################

# T0 vs T1, only paired samples

if(! "proof1" %in% ls()){
  cat("\nWait! Did you performed the filtering steps???\n\n")
Sys.sleep(2)
}

remove<-ls()[!ls() %in% c("proof1","data_all","data_contam","Metadata")] # to reset the environment
rm(list = remove)

data<-subset_samples(data_all, Pair == "Yes") # only samples with also T1
data<-prune_taxa(taxa_sums(data)>0, data)

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


###################### ABUNDANCES BAR PLOT (TREATMENT, IN PAIR) ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("red3","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","red3","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one

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
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Phylum)) + theme_classic(base_size =14) + 
  facet_grid2(Sample_ID+Treatment_time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/TOP5_phyla_PAIR_ALONG_TIME.png",width=7,height=11, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Abundances/TOP5_phyla_average_abundance_ONLY_PAIRED_SAMPLES.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))


# TOP 5 Genera
suppressWarnings(rm(top, others, tabella))
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
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(Sample_ID+Treatment_time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_5) +
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Abundances/TOP5_genera_PAIR_ALONG_TIME.png",width=7,height=11, dpi=300) 
dev.off()

# TOP 8 Genera
suppressWarnings(rm(top, others, tabella))
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
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(Sample_ID+Treatment_time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_8) +
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) +
  labs(y="Patients", x="Relative abundance", title = "Eight most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Abundances/TOP8_genera__PAIR_ALONG_TIME.png",width=7,height=11, dpi=300) 
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))


# #################### RATIO FIRMICUTES/BACTEROIDES (TREATMENT, IN PAIR) ###################
# 
# suppressWarnings(rm(data_fb, ratio_fb))
# 
# data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
# F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
# B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
# ratio_fb<-otu_table(data_fb)
# ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
# row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
# ratio_fb<-t(ratio_fb)
# 
# ratio_fb<-ratio_fb[sample_data(data)$Sample_name, ] # same order
# identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
# ratio_fb<-cbind.data.frame(ratio_fb,sample_data(data)[,c("Sample_ID","Treatment_time")])
# 
# # mean
# Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Treatment_time, mean)
# ratio_fb$Ratios_Mean<-rep("temp")
# ratio_fb[ratio_fb$Treatment_time=="T1","Ratios_Mean"]<-Ratios_Mean["T1"]
# ratio_fb[ratio_fb$Treatment_time=="T0","Ratios_Mean"]<-Ratios_Mean["T0"]
# ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# # st err
# Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Treatment_time, sd)
# ratio_fb$Ratios_st_err<-rep("temp")
# ratio_fb[ratio_fb$Treatment_time=="T1","Ratios_st_err"]<-Ratios_st["T1"]/sqrt(length(which(ratio_fb$Treatment_time=="T1")))
# ratio_fb[ratio_fb$Treatment_time=="T0","Ratios_st_err"]<-Ratios_st["T0"]/sqrt(length(which(ratio_fb$Treatment_time=="T0")))
# ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)
# 
# head(ratio_fb, n=2)
# write.csv2(file="Results/Abundances/Ratio_Firmi_Bacteroi_T0_vs_T1/Firmi_Bacter_Ratio.csv", ratio_fb)
# 
# ratio_fb$Treatment_time<-factor(ratio_fb$Treatment_time, levels = c("T1","T0"))
# 
# 
# ####### STATISTICAL TEST
# 
# distr<-shapiro.test(ratio_fb[,"Ratio"])
# distr
# hist(ratio_fb[,"Ratio"])
# 
# res<-wilcox.test(Ratio~Treatment_time, paired=TRUE, data=ratio_fb)
# p_val<-round(res$p.value, digits = 4)
# p_val
# # exporting the results
# con <- file("Results/Abundances/Ratio_Firmi_Bacteroi_T0_vs_T1/Wilcoxon_paired.txt") # create this file and connect it to R through "con" object 
# sink(con, append = TRUE) # redirect STR ERROR to "con"
# cat("Wilcoxon test (Paired)\n", fill=T)
# cat("Ratio~Treatment_time:   V=",res$statistic, "p-value=", p_val, "\n",fill=T)
# cat("\n Shapiro Wilk test p-value: ",distr$p.value)
# sink()
# close(con)
# 
# 
# ######## BAR PLOT
# ratio_fb_plot<-ratio_fb
# ratio_fb_plot[ratio_fb_plot$Ratio>100, "Ratio"]<-100    # almost no Bacteroidetes --> ratio too high --> setting a limit
# ggplot(data=ratio_fb_plot, aes(x=Sample_ID, fill=Treatment_time, y=Ratio)) +
#   theme_bw(base_size =12) + 
#   scale_fill_manual(values = c("T0"="coral",
#                                "T1"="red3")) +
#   facet_grid2(.~Treatment_time, scales = "free_x", 
#               space="free", strip = strip_nested(size="constant"))+
#   theme(panel.spacing.x = unit(2,"pt"))+
#   scale_x_discrete (expand = c(0.01,0) ) +
#   geom_line(aes(y= ratio_fb$Ratios_Mean, group="Treatment_time"))+
#   scale_y_sqrt(breaks=c(1,5,seq(10,100,10)) ) +
#   geom_bar(stat="identity", position="stack", width = 0.8) +
#   guides(fill="none") +
#   theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
#   labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio",
#        subtitle = paste("Mann Whitney p-value:",p_val),
#        caption = "the line is the average ratio of the group")
# ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi_T0_vs_T1/Ratio_barplot.png",width=8,height=4, dpi=300) 
# dev.off()
# 
# 
# ####### JITTER PLOT (with means)
# set.seed(2)
# ggplot(data=ratio_fb_plot, aes(x=Treatment_time, color=Treatment_time,y=Ratio)) +
#   theme_bw(base_size =9) +
#   scale_color_manual(values = c("T1"="red3","T0"="coral")) +
#   facet_grid(.~Treatment_time, scales = "free_x", space="free")+
#   theme(panel.spacing.x = unit(20,"pt"))+
#   scale_x_discrete (expand = c(0.2,0.3) ) +
#   scale_y_sqrt(breaks=c(1,5,seq(10,100,10)) ) +
#   geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
#                     ymax = Ratios_Mean),
#                 size=0.35,width = 0.35, color= "black") +
#   geom_errorbar(aes(ymax= ifelse(Ratios_Mean + Ratios_st_err > 100, # standard deviation, the if else is needed to avoid geom bar below zero
#                                  100, Ratios_Mean + Ratios_st_err),
#                     ymin= ifelse(Ratios_Mean - Ratios_st_err < 0,
#                                  0, Ratios_Mean - Ratios_st_err)),
#                 width=.1, color="black", size= 0.15) +
#   geom_jitter(width = 0.25, size=0.5, show.legend = F) +
#   guides(color="none") +
#   theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
#         axis.text.y=element_text(size=6)) +
#   labs(x="", y="Firmicutes/Bacteroidetes Ratio", subtitle = paste("Mann Whitney p-value:",p_val),
#        caption="The upper limit of the Ratio is manually setted to 100\n because the absence of Bacteroidetes causes a really high ratio")
# ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi_T0_vs_T1/Ratios_jitter_with_mean_and_STerr.png",width=4,height=3.5, dpi=300) 
# dev.off()
# 
# suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY (TREATMENT, IN PAIR) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Treatment_time")
pAlpha
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Treatment_time, New_data$Sample_ID),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha + theme_bw() + 
  geom_line(aes(group = pAlpha$data$Sample_ID),col="grey",size=0.15) +
  labs(x="Treatment_time", title="Alpha diversity between T1 and T0 patients") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=1, size=11)) +
  stat_compare_means(aes(group = Treatment_time), label="p.format", method = "wilcox.test", paired = T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4)
ggsave(file="Results/Alfa_diversity_T0_vs_T1_TREATMENT.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
Obser_value<-Obser_value[order(Obser_value$Treatment_time, Obser_value$Sample_ID),]
factor<-Obser_value$Treatment_time
factor
Obser_value$Sample_ID # to re-check if they have the same order
wilcox.test(Obser_value$value~factor, paired=T)

rm(pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY (TREATMENT, IN PAIR) #######################

suppressWarnings(rm(ASV.prop))

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

###### PERMANOVA

metadata<-as(sample_data(data.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Beta_div/T0_vs_T1/Beta_divers_permanova_Helling.csv",quote=F,row.names = T)


### PLUS: checking it with Bray too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Sample_ID+Treatment_time, data=metadata, permutations = 9999, method="bray")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Beta_div/T0_vs_T1/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Treatment_time)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Beta_div/T0_vs_T1/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)


######################## PCoA BETA DIV (TREATMENT, IN PAIR) #########################

# on ASV
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Treatment_time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  scale_color_manual(values=c("T0"="red3","T1"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", color="Treatment_time", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), 
       subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H),
       caption="lines connect paired sample")
ggsave(file="Results/Beta_div/T0_vs_T1/PCoA_Beta_diversity_Hellinger_on_ASV.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Treatment_time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  scale_color_manual(values=c("T0"="red3","T1"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", color="Treatment_time", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired sample")
ggsave(file="Results/Beta_div/T0_vs_T1/PCoA_Beta_div_Hellinger_ASV_no_ellipse.png", width = 8, height = 6, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Treatment_time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  scale_color_manual(values=c("T0"="red3","T1"="coral")) +
  geom_point(size=2.5) + theme_classic(base_size = 14) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", color="Treatment_time", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =", perm_ASV_H),
       caption="lines connect paired sample")
ggsave(file="Results/Beta_div/T0_vs_T1/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 8, height = 6, dpi=300)


# again but on genera
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Treatment_time") +
  scale_color_manual(values=c("T0"="red3","T1"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  # geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Treatment_time", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), 
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H),
       caption="lines connect paired sample")
ggsave(file="Results/Beta_div/T0_vs_T1/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 8, height = 6, dpi=300)


############## DA WITH DESEQ2 (TREATMENT, IN PAIR) #############

if(! "proof1" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sample_ID+Treatment_time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Treatment_time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res<-res[res$baseMean > 50, ] # arbitrary threshold to avoid the most noisy result
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
    write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_T0_vs_T1.csv"), row.names = F, quote=F, na = "")
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

system(" echo 'Any Differential abundance regarding the treatment (t0 vs t1)' > Results/DA_DESeq2/Anything_about_treatment.txt ")



################### ANALYSIS ON CLINICAL VALUES AT EACH TREATMENT TIME ########################

Values<-read.table(file="SCFA_and_clinical_values.csv", sep=",", dec=".", header = T)
head(Values, n=3)

FFA <- Values[ , ! colnames(Values) %in% c("Visus","OCT_change")]
Visus<- Values[ , colnames(Values) %in% c("Treatment_time","Visus","OCT_change")] # NB: same original table --> SAME ORDER already

Signif_FFA <- FFA[ , c("Treatment_time", "IsoHexanoic", "PhenylAcetic", "PhenylPropionic")] # significant according to the analyses of another team
Signif_FFA_T0 <- Signif_FFA[Signif_FFA$Treatment_time=="T0", ! colnames(Signif_FFA) %in% "Treatment_time"] 
Signif_FFA_T1 <- Signif_FFA[Signif_FFA$Treatment_time=="T1",  ! colnames(Signif_FFA) %in% "Treatment_time"] 

# FFA_T0 <- FFA[FFA$Treatment_time=="T0", ]
# FFA_T1 <- FFA[FFA$Treatment_time=="T1", ]

Visus_T0 <- Visus[Visus$Treatment_time=="T0", ! colnames(Visus) %in% c("Sample_ID","Treatment_time") ]
Visus_T1 <- Visus[Visus$Treatment_time=="T1", ! colnames(Visus) %in% c("Sample_ID","Treatment_time") ]


# testing the distributions
png(filename = "Results/Tests_and_Correlations/Testing_some_distribution.png")
par(mfrow=c(2,2))
hist(Values[Values$Treatment_time=="T0", "OCT_change"], main="Only the OCT change is gaussian", xlab = "OCT change at T0", sub=paste("Gaussian probability:",round(shapiro.test(Values[Values$Treatment_time=="T0", "OCT_change"])$p.value,3) ) )
hist(Values[Values$Treatment_time=="T0", "Visus"], main="then is not possible to use Pearson", xlab = "Visus at T0", sub=paste("Gaussian probability:",round(shapiro.test(Values[Values$Treatment_time=="T0", "Visus"])$p.value,3) ) )
hist(Values[Values$Treatment_time=="T1", "IsoHexanoic"], main="", xlab = "IsoHexanoic at T1", sub=paste("Gaussian probability:",round(shapiro.test(Values[Values$Treatment_time=="T1", "IsoHexanoic"])$p.value,3) ) )
hist(Values[Values$Treatment_time=="T1", "Visus"], main="" , xlab = "Visus at T1", sub=paste("Gaussian probability:",round(shapiro.test(Values[Values$Treatment_time=="T1", "Visus"])$p.value,3) ) )
dev.off()


######## CORRELATIONS 

# Computing FFAs vs Visus at T0
x<- cbind(Signif_FFA_T0, Visus_T0)
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("FFA","Visus","Corr","pvalue")
}
corr<-subset(correlation, FFA %in% colnames(Signif_FFA_T0))
corr<-subset(corr, Visus %in% colnames(Visus_T0))
colnames(corr)<-c("FFA","Visus","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$FFA, y = corr$Visus, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =8) +
  labs(title = "Pre treatment", y= "Visus", x= "FFA", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Tests_and_Correlations/FFA_vs_Visus_only_T0.png", dpi=300, width = 7, height = 7)

write.csv2(corr, file="Results/Tests_and_Correlations/FFA_vs_Visus_only_T0.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr)


# Computing FFAs vs Visus at T1
x<- cbind(Signif_FFA_T1, Visus_T1)
{x<-as.matrix(x)
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("FFA","Visus","Corr","pvalue")
}
corr<-subset(correlation, FFA %in% colnames(Signif_FFA_T1))
corr<-subset(corr, Visus %in% colnames(Visus_T1))
colnames(corr)<-c("FFA","Visus","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "BH")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$FFA, y = corr$Visus, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =8) +
  labs(title = "", y= "Visus", x= "FFA", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Tests_and_Correlations/FFA_vs_Visus_only_T1.png", dpi=300, width = 7, height = 7)
### again but plotting the p-value, not the p-adjusted
corr$Sign<-corr$pvalue
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$FFA, y = corr$Visus, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =8) +
  labs(title = "Post treatment", y= "Visus", x= "FFA", 
       caption= "\n p-values (not adjusted) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Tests_and_Correlations/FFA_vs_Visus_only_T1_NOT_ADJUSTED.png", dpi=300, width = 7, height = 7)

write.csv2(corr, file="Results/Tests_and_Correlations/FFA_vs_Visus_only_T1.csv", quote = F, na="", row.names = F)

rm(correlation, correlation_corr, correlation_pvalue, corr)


######## PCoA

dist_matrix <- vegdist(FFA[ , ! colnames(FFA) %in% c("Sample_ID","Treatment_time")], method = "bray")
dist_matrix

png("Results/Tests_and_Correlations/PCoA_on_every_FFA_Bray-Curtis.png", res=300,width = 1800, height = 1800)
pcoaVS <- ecodist::pco(dist_matrix, negvals = "zero")
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], xlab = "",
       ylab = "", axes = TRUE, type="n", sub = "T0 = black\nT1 = red")
text(pcoaVS$vectors[,1], pcoaVS$vectors[,2], labels(FFA$Sample_ID), cex = 0.9, xpd = TRUE,
     col = ifelse(FFA$Treatment_time=="T0",1,2) )
dev.off()



##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results/R_version_and_packages.txt")
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
cat("\n", "\n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
