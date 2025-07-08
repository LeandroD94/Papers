##################### PREPARING THE ENVIRONMENT #################
{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggvenn")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("DESeq2")
  library("vegan")
  # utilities
  library("xlsx")  
  library("dplyr")
  library("qiime2R")
}

setwd ('../../../SARA/Breast_cancer_H_vs_Tumor_Fazia')
{dir.create("Data_check")
 dir.create("Data_check/PCoA_test")
 dir.create("Results/")
 dir.create("Results/General_Abundances")

for(x in c("H_vs_Cancer/",
           # "Age_group_ONLY_cancer_tissue/",
           # "Age_group_ONLY_Healthy_tissue/",
           "City_ONLY_cancer_tissue/",
           "City_ONLY_Healthy_tissue/"
           )) {
  dir.create(paste0("Results/",x))
  dir.create(paste0("Results/",x,"Abundances"))
  dir.create(paste0("Results/",x,"Beta_div"))
  dir.create(paste0("Results/",x,"DA_DESeq2"))
  }

dir.create("Results/Radiotherapy_ONLY_cancer_tissue/")
dir.create("Results/Chemotherapy_ONLY_cancer_tissue/")
dir.create("Results/PCoA_EXTRA_FACTORS_only_cancer_tissues")

}
rm(x)

options(scipen = 100) # disable scientific annotation


####################### IMPORTING DATA #####################
# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample

Metadata <- read.csv("Metadata.csv", header = T)
row.names(Metadata)<-Metadata$FASTQ_Code # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,sample)

sample_data(data)$Tissue<-factor(sample_data(data)$Tissue, levels=c("Healthy","Cancer")) # decides the order in plots
sample_names(data)<-sample_data(data)$FASTQ_ID
head(sample_data(data))

# Estrai la matrice taxonomica
tax_table <- tax_table(data)

# Estrai la matrice di counts (abundance)
otu_table <- otu_table(data)

# Controlla se taxa senza assegnazione a livello Genus sono NA o ""
# Qui assumiamo che siano NA, altrimenti sostituisci con ""
unassigned_taxa <- is.na(tax_table[, "Genus"]) | tax_table[, "Genus"] == ""

# Calcola il totale di reads in ogni campione
total_reads_per_sample <- sample_sums(data)

# Calcola la somma delle reads dei taxa non assegnati a livello Genus per ogni campione
unassigned_reads_per_sample <- colSums(otu_table[unassigned_taxa, ])

# Calcola la percentuale per campione
percent_unassigned <- (unassigned_reads_per_sample / total_reads_per_sample) * 100
library(writexl)

# Creiamo un data frame con le percentuali per campione
df_percent_unassigned <- data.frame(
  SampleID = names(percent_unassigned),
  PercentUnassigned = percent_unassigned
)

# Salva il data frame come file xlsx
write_xlsx(df_percent_unassigned, path = "percent_unassigned_reads_per_sample.xlsx")

# Visualizza la percentuale media su tutti i campioni
mean_percent_unassigned <- mean(percent_unassigned)

mean_percent_unassigned

########## IMPORTING ALSO THE CONTAMINATED ONE (FROM HOST DNA)

data_contam<-qza_to_phyloseq(features="QIIME/table.qza", tree = "QIIME/eucar_contaminants/rooted-tree_with_contam.qza")
# changing names (already done before, check up there!)
if(identical(original_names,sample_names(data_contam))){
  sample_names(data_contam)<-sample_names(data)
}
if(identical(sample_names(data_contam), Metadata$FASTQ_ID)) {
  sample_data(data_contam)<-sample_data(data)
  cat("Everything alright there!")} else { cat("\n Sample names are not matching, check them \n\n") }


#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)

###### cutting under 0.01% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.01, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_001_MEAN_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.01, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured" & !is.na(filtered)] # to avoid the removal of other uncultured genera
data<-subset_taxa(data, ! Genus %in% filtered)
data<-prune_taxa(taxa_sums(data) > 0, data) 
rm(filtered, data.genus.temp)


######## Checking unassigned in prokaryote kingdom
{Unass<-tax_glom(data, taxrank = "Kingdom", NArm = F) # or domain
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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Unassigned_domain_checking.csv", row.names = T, quote = F)
####### now filtering out every Unassigned ASV
data<- subset_taxa(data, ! Kingdom %in% c("Unassigned","d__Eukaryota") )
data<- subset_taxa(data, ! is.na(Phylum)  ) # with this classifier, a NA at phylum level is probably an eucaryotic contaminat
head(tax_table(data))
write.csv2(tax_table(data), file="Data_check/Every_filtered_ASV_and_taxonomic_assignment.csv", row.names = T)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


#### Applying a prevalence filter: at least 1 read at least in 2 samples in the whole dataset
data.temp<-tax_glom(data, taxrank = "Genus", NArm = F)
who<-as.data.frame(otu_table(data.temp) )
who<-apply(who, MARGIN=2, function(x) ifelse(x>1, 1, 0)) # if at least one read in the sample then it gets "a point"
who<-who[!rowSums(who)>1,] # more than 1 "points" --> at least in 2 samples
who<-as.vector(tax_table(data.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
who<-who[!is.na(who) & who!="uncultured"] # otherwise EVERY other NA and uncultured would be also eliminated
who # --> ok!
write.csv(who, "Data_check/Filtered_genera_that_were_ONLY_in_ONE_sample.csv", row.names = F)
data<-subset_taxa(data, ! Genus %in% who)
rm(who, data.temp)

# removing the Chloroplasts and Mithocondria
data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mithocondria"))


# View(tax_table(data))
proof1<- "Marker of the filtering, it is required for the script"


################# CHECKING THE COMPOSITION AFTER FILTERING ####################

data.unf.prop<-transform_sample_counts(data_contam, function(x) (x/sum(x))*100)


if(! "proof1" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}
data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100)


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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Tissue", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Tissue", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3200, height = 1800, res=300)
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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Tissue", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Tissue", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


################# wUNIFRAC PROP

##### unfiltered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
data.prop.labels<-data.unf.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.prop.labels, ordBC, color = "Tissue") +
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data", 
       color="Tissue", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Tissue") +
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Tissue", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels,data.prop,data.unf.prop))


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

png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=0.5)
dev.off()
rm(r)

# the sample DM70 is not saturated
who_unsaturated<-Metadata[Metadata$FASTQ_ID=="DM70", "ID"] # patient 8


#################### BACKUP OF THE COMPLETE DATASET BEFORE THE SUBSETTINGS ##################

if( length(sample_names(data)) >=31 ){
  complete_Metadata<-Metadata
  complete_data<-data
  rm(data,Metadata)
}


suppressWarnings(rm(original_names,evalslopes))
# save.image("data.RData")
### NB: this backup still has every sample


################## REMOVING AND ANNOTATING ANOMALOUS SAMPLES ######################

complete_data <- subset_samples(complete_data, FASTQ_ID!= "DM70") # removing the unsaturated sample

# same sample has not been sequenced successfully, then its pair is also removable IN THE PAIR analysis
pairs<- complete_Metadata$ID[duplicated(complete_Metadata$ID)] # only the complete pairs

proof2<- "The unsaturated sample has been removed"

####################### ***** PREPARATION OF THE DATA (GENERAL CHECK) #######################

# this dataset differs from the PAIRED ANALYSIS one because it still has the unpaired samples
data<-complete_data
Metadata<-complete_Metadata

sample_data(data)



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


#################### % ASSIGNED IN SILVA (GENERAL CHECK) #########################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Total","Assigned","%","Taxa")
assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)


######################## GOOD'S COVERAGE ESTIMATOR (GENERAL CHECK) #########################

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


########################### COUNTS EXPORT (GENERAL CHECK) ##########################################


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
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


####################### ***** PREPARATION OF THE DATA (H_vs_Cancer) #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  cat("\n\nDid you perform the decontamination and the removal of unsatureted samples already???\n\n")
  Sys.sleep(2)
}

# to reset the environment
remove<-ls()[! ls() %in% c("complete_Metadata", "complete_data", "unfiltered_data", "data_contam", "pairs", "who_unsaturated", "proof1","proof2")]
rm(list=remove)

  
data<-subset_samples(complete_data, ID %in% pairs)
data<-subset_samples(data, ! ID %in% who_unsaturated)
Metadata<-complete_Metadata[complete_Metadata$ID %in% pairs, ]
Metadata<-Metadata[! Metadata$ID %in% who_unsaturated, ]

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


###################### ABUNDANCES BAR PLOT  (H_vs_Cancer) ##########################

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
tabella$ID<-gsub("Patient_"," ",tabella$ID)
tabella$Tissue<-gsub("Healthy","H",tabella$Tissue)
tabella$Tissue<-gsub("Cancer","C",tabella$Tissue)
tabella$Tissue<-factor(tabella$Tissue, levels=c("H","C"))
ggplot(data=tabella, aes(x=Abundance, y=Work_ID, fill=Phylum)) + theme_classic(base_size =14) + 
  facet_grid2(ID+Tissue~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP5_phyla_abundances.pdf",width=8.5,height=11, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/H_vs_Cancer/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))


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
tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$ID<-gsub("Patient_"," ",tabella$ID)
tabella$Tissue<-gsub("Healthy","H",tabella$Tissue)
tabella$Tissue<-gsub("Cancer","C",tabella$Tissue)
tabella$Tissue<-factor(tabella$Tissue, levels=c("H","C"))
ggplot(data=tabella, aes(x=Abundance, y=Work_ID, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(ID+Tissue~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_8) +
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
  labs(y="Patients", x="Relative abundance", title = "Eight most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="TOP8_genera_abundances.png",width=7,height=8, dpi=300) 
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))


########################## ALFA DIVERSITY  (H_vs_Cancer) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Tissue", color="Tissue")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$ID, obs$ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Tissue, New_data$ID),]
pAlpha$data<-New_data
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
# with points only
pAlpha + theme_classic2() +
  scale_color_manual(values = c("Healthy"="chartreuse", "Cancer"="coral")) +
  geom_line(aes(group = pAlpha$data$ID),col="grey",size=0.15) +
  labs(x="Tissue", title="Alpha diversity between Cancer and Healthy patients") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Tissue), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/H_vs_Cancer/Alfa_diversity_GENUS_paired_wilcoxon.png", width = 6,height =5.6, dpi=300)
# with ID
pAlpha + theme_classic2() +
  scale_color_manual(values = c("Healthy"="chartreuse", "Cancer"="coral")) +
  geom_line(aes(group = pAlpha$data$ID),col="grey",size=0.15) +
  geom_text(aes(label=ID), size= 2.5, color= "grey") +
  labs(x="Tissue", title="Alpha diversity between Cancer and Healthy patients") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Tissue), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/H_vs_Cancer/Alfa_diversity_GENUS_paired_wilcoxon_with_ID.png", width = 6,height =5.6, dpi=300)


suppressWarnings(rm(pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor))


##################### BETA DIVERSITY  (H_vs_Cancer) #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~ID + Tissue, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~ID + Tissue, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~ID + Tissue, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~ID + Tissue, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~ID + Tissue, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~ID + Tissue, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/H_vs_Cancer/Beta_div/Beta_divers_permanova_Helling.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Tissue)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/H_vs_Cancer/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)


######################## PCoA BETA DIV  (H_vs_Cancer) #########################

# on genera
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  geom_line(aes(group=ID),col="grey", size=0.15)+
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  geom_point(size=3.5, alpha=0.3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Tissue", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), 
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H),
       caption="lines connect paired samples")
ggsave(file="Results/H_vs_Cancer/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 7, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  geom_line(aes(group=ID),col="grey", size=0.15)+
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  geom_point(size=3.5, alpha=0.3) + theme_classic(base_size = 13) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=2.3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Tissue", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples")
ggsave(file="Results/H_vs_Cancer/Beta_div/PCoA_Beta_div_Hellinger_genera_no_ellipse.png", width = 6.5, height = 5, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Tissue") +
  geom_line(aes(group=ID),col="grey", size=0.15)+
  scale_color_manual(values=c("Cancer"="coral","Healthy"="chartreuse")) +
  geom_point(size=3.5, alpha=0.3) + theme_classic(base_size = 13) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Tissue", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples")
ggsave(file="PCoA_Beta_diversity_Hellinger_on_genera_points_HC_cancer.pdf", width = 8.5, height = 11, dpi=300)


#################### VEEN DIAGRAM (H_vs_Cancer Tissue) ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.1% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Healthy<-subset_samples(data.venn, Tissue=="Healthy")
Healthy<-as.character(tax_table(prune_taxa(taxa_sums(Healthy)>0, Healthy))[,"Genus"])

Cancer<-subset_samples(data.venn, Tissue=="Cancer")
Cancer<-as.character(tax_table(prune_taxa(taxa_sums(Cancer)>0, Cancer))[,"Genus"])


ONLY_IN_Cancer<- Cancer[! Cancer %in% Healthy]
ONLY_IN_Cancer<- paste(ONLY_IN_Cancer, collapse = ", ")
head(ONLY_IN_Cancer)

ONLY_IN_Healthy<- Healthy[! Healthy %in% Cancer]
ONLY_IN_Healthy<- paste(ONLY_IN_Healthy, collapse = ", ")
head(ONLY_IN_Healthy)


x<-list(Healthy=Healthy,Cancer=Cancer)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("chartreuse","coral")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = " Venn Diagram of genera in Healthy and Cancer tissue \n (only genera with minimal abundance > 0.1% \n at least in three samples)\n")
ggsave(filename = "Venn_Diagramm.pdf", width = 8.5, height = 11, dpi=300, bg = "white")
dev.off()


suppressWarnings(rm(ONLY_IN_Cancer, ONLY_IN_Healthy, x, con, Cancer, Healthy, data.venn, who))


############## DIFFERENTIAL ABUNDANCES WITH DESEQ2  (H_vs_Cancer) #################

if(! "unfiltered_data" %in% ls() | ! "proof2" %in% ls()){
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
  
  # if the paired samples are seen as numbers (it depends from on levels name)
  sample_data(d)[["ID"]]<-as.character(sample_data(d)[["ID"]])
  
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
  DEseq_data<-phyloseq_to_deseq2(d, ~ID+Tissue)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Tissue", "Healthy", "Cancer"))
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
    write.csv2(r, file=paste0("Results/H_vs_Cancer/DA_DESeq2/DA_",t,"_ratio_Healthy_vs_Cancer.csv"), row.names = F, quote=F, na = "")
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

# View(Res_tot)
system(" echo 'Any significative result here!' >> Results/H_vs_Cancer/DA_DESeq2/No_result.txt ")  


################### ***** PREPARATION OF THE DATA (Radiotherapy _ Cancer Tissue) #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  cat("\n\nDid you perform the decontamination and the removal of unsatureted samples already???\n\n")
  Sys.sleep(2)
}

# to reset the environment
remove<-ls()[! ls() %in% c("complete_Metadata", "complete_data", "unfiltered_data", "data_contam", "pairs", "who_unsaturated", "proof1","proof2")]
rm(list=remove)


data <- subset_samples(complete_data, Tissue=="Cancer") # Only cancer Samples


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


###################### ABUNDANCES BAR PLOT  (Radiotherapy _ Cancer Tissue) ##########################

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
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Radiotherapy),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (only Cancer Tissues)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP_5_phyla_2.pdf",width=8.5,height=11, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Radiotherapy_ONLY_cancer_tissue/TOP5_phyla_average_abundance_CANCER.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
  tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Radiotherapy),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (only Cancer Tissues)", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Radiotherapy_ONLY_cancer_tissue/TOP_8_genera.png",width=7,height=5,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


########################### PCoA  (Radiotherapy _ Cancer Tissue) ########################

# on genera
{data.prop.labels<-data.genus.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Radiotherapy") +
  scale_color_manual(values=c("no"="royalblue3","yes"="deepskyblue")) +
  geom_point(size=3.5, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       color="Radiotherapy", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Radiotherapy_ONLY_cancer_tissue/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 6.5, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### ***** PREPARATION OF THE DATA (Chemotherapy _ Cancer Tissue) #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  cat("\n\nDid you perform the decontamination and the removal of unsatureted samples already???\n\n")
  Sys.sleep(2)
}

# to reset the environment
remove<-ls()[! ls() %in% c("complete_Metadata", "complete_data", "unfiltered_data", "data_contam", "pairs", "who_unsaturated", "proof1","proof2")]
rm(list=remove)


data <- subset_samples(complete_data, Tissue=="Cancer") # Only cancer Samples

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


###################### ABUNDANCES BAR PLOT  (Chemotherapy _ Cancer Tissue) ##########################

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
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Chemotherapy),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (only Cancer Tissues)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP_5_phyla_onlycancer.pdf",width=8.5,height=11, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Chemotherapy_ONLY_cancer_tissue/TOP5_phyla_average_abundance_CANCER.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
  tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Chemotherapy),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (only Cancer Tissues)", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Chemotherapy_ONLY_cancer_TOP_8_genera.png",width=7,height=5,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


########################### PCoA  (Chemotherapy _ Cancer Tissue) ########################

# on genera
{data.prop.labels<-data.genus.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Chemotherapy") +
  scale_color_manual(values=c("no"="royalblue4","yes"="deepskyblue")) +
  geom_point(size=3.5, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       color="Chemotherapy", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Chemotherapy_ONLY_cancer_tissue/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 6.5, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



################### ***** PREPARATION OF THE DATA (City  Cancer Tissue) #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  cat("\n\nDid you perform the decontamination and the removal of unsatureted samples already???\n\n")
  Sys.sleep(2)
}

# to reset the environment
remove<-ls()[! ls() %in% c("complete_Metadata", "complete_data", "unfiltered_data", "data_contam", "pairs", "who_unsaturated", "proof1","proof2")]
rm(list=remove)

data <- subset_samples(complete_data, Tissue=="Cancer") # Only cancer Samples

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


###################### ABUNDANCES BAR PLOT  (City  Cancer Tissue) ##########################

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
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(City),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (only Cancer Tissues)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP_5_phyla_urban_rural.pdf",width=8.5,height=11, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/City_ONLY_cancer_tissue/Abundances/TOP5_phyla_average_abundance_CANCER.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
  tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(City),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) +
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (only Cancer Tissues)", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="TOP_8_genera_urban_rural.pdf",width=8.5,height=11,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


########################## ALFA DIVERSITY  (City Cancer Tissue) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="City", color="City")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$ID, obs$ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Tissue, New_data$ID),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
  pAlpha$data$City<-gsub("more_than_","> ",pAlpha$data$City, fixed=T)
  pAlpha$data$City<-gsub("less_than_","< ",pAlpha$data$City, fixed=T)
}
# with points only
pAlpha + theme_classic2() +
  geom_boxplot(aes(group=City), size =0.8, color = "black", alpha= 0) +
  scale_color_manual(values = c("Rural"="coral2", "Urban"="royalblue")) +
  labs(x="Tissue", title="Alpha diversity between cities (Cancer Tissues)") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = City), label="p.format", method = "wilcox.test", paired = F,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/City_ONLY_cancer_tissue/Alfa_diversity_GENUS.png", width = 6,height =5.6, dpi=300)



##################### BETA DIVERSITY  (City Cancer Tissue) #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[1]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[1]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1]

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/City_ONLY_cancer_tissue/Beta_div/Beta_div_permanova_Hellinger_Cancer.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
set.seed(1)
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$City)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/City_ONLY_cancer_tissue/Beta_div/Beta_dispersion_Hellinger_on_genera.csv",quote=F,row.names = T)


###### FURTHER TEST : it is still significative without that strange sample???
data.temp<-subset_samples(data.genus.prop, ID!="Patient_13")
ASV.genus.prop2<-as.data.frame(otu_table(data.temp))
set.seed(1)
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop2)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,sample_data(data.temp)$City)
disp_ASV2<-vegan::permutest(disper, permutations=9999)
disp_ASV2



suppressWarnings(rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV))


########################### PCoA  (City  Cancer Tissue) ########################

# on genera
{data.prop.labels<-data.genus.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.6, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/City_ONLY_cancer_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_PATIENT.png", width = 6.5, height = 5, dpi=300)
# only points
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.2, alpha=1, shape=16) + geom_point(size=4.5, alpha=0.65) +
  theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F):",perm_g_H,"\nBeta dispersion Pr(>F):",disp_ASV$tab$`Pr(>F)`[1]))
ggsave(file="PCoA_Beta_diversity_Hellinger_on_genera_POINT.pdf", width = 8.5, height = 11, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.2, alpha=1, shape=16) + geom_point(size=4.5, alpha=0.65) +
  theme_classic(base_size = 13) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F):",perm_g_H,"\nBeta dispersion Pr(>F):",disp_ASV$tab$`Pr(>F)`[1]))
ggsave(file="Results/City_ONLY_cancer_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_NO_ELLIPSES.png", width = 6.5, height = 5, dpi=300)

#### checking without that sample
{data.prop.labels<-subset_samples(data.genus.prop, ID != "Patient_13")
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.6, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="This is a further check: the dispersion is still significative\n even after the removal of that strange sample --> good!",
       subtitle = paste("Betadispersion Pr(>F):", disp_ASV2$tab$`Pr(>F)`),
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/City_ONLY_cancer_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_JUST_A_CHECK.png", width = 6.5, height = 5, dpi=300)


suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


#################### VEEN DIAGRAM (City Cancer Tissue) ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.1% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Rural<-subset_samples(data.venn, City=="Rural")
Rural<-as.character(tax_table(prune_taxa(taxa_sums(Rural)>0, Rural))[,"Genus"])

Urban<-subset_samples(data.venn, City=="Urban")
Urban<-as.character(tax_table(prune_taxa(taxa_sums(Urban)>0, Urban))[,"Genus"])


ONLY_IN_Urban<- Urban[! Urban %in% Rural]
ONLY_IN_Urban<- paste(ONLY_IN_Urban, collapse = ", ")
head(ONLY_IN_Urban)

ONLY_IN_Rural<- Rural[! Rural %in% Urban]
ONLY_IN_Rural<- paste(ONLY_IN_Rural, collapse = ", ")
head(ONLY_IN_Rural)


x<-list(Rural=Rural,Urban=Urban)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("coral2","royalblue")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = " Venn Diagram of genera in Cancer Tissue \n (only genera with minimal abundance > 0.1% \n at least in three samples)\n")
ggsave(filename = "Venn_Diagramm_urban_rural.pdf", width = 8.5, height = 11, dpi=300, bg = "white")
dev.off()

write.csv2(ONLY_IN_Rural, file="Results/City_ONLY_cancer_tissue/Beta_div/Exclusive_genera_list.txt", quote = F, row.names = F)

suppressWarnings(rm(ONLY_IN_Urban, ONLY_IN_Rural, x, con, Urban, Rural, data.venn, who))



################### DA WITH DESEQ2 (City Cancer Tissue) #####################

if(! "proof1" %in% ls() | ! "proof2" %in% ls() ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~City)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("City", "Rural", "Urban"))
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
write.csv(Res_tot, file="Results/City_ONLY_cancer_tissue/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/City_ONLY_cancer_tissue/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

# boxplot
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=City)) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + 
  theme_bw(base_size = 12) +
  theme(strip.text.x=element_text(size=13,colour="black")) +
  scale_fill_manual(values=c("Rural"="coral2","Urban"="royalblue")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom",
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=9.5),
        axis.text.y = element_text(size=9.5), plot.title= element_text(size=16),
        panel.grid.minor.y= element_blank() ) +
  scale_x_discrete(expand=c(-0.2, 1)) +
  # scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance)+2,3))) +
  scale_y_sqrt(breaks=c(0, 0.1, 0.5, 1,1.5,2,2.5,3,3.5,4)) +
  labs(title= "Differently abundant Taxa (Cancer Tissue)", y="Percent Abundance",
       fill="City group", x="")
ggsave(filename = "DA_City_Cancer_Tissue_BOXPLOT.png", width = 7, height = 5, dpi=300)
dev.off()

# boxplot AND points
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, color=City, fill=City)) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8, size=0.4, outlier.size = 1, color="black",alpha=0.3) +
  geom_point(position = position_dodge(width=0.8),size=2.4, alpha=0.3, color="black", aes(group=City)) +
  geom_point(position = position_dodge(width=0.8),size=2, alpha=1, aes(group=City)) +
  theme_bw(base_size = 12) +
  theme(strip.text.x=element_text(size=13,colour="black")) +
  scale_fill_manual(values=c("Rural"="coral2","Urban"="royalblue")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom",
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=9.5),
        axis.text.y = element_text(size=9.5), plot.title= element_text(size=16),
        panel.grid.minor.y= element_blank() ) +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0, 0.1, 0.5, 1,1.5,2,2.5,3,3.5,4)) +
  labs(title= "Differently abundant Taxa (Cancer Tissue)", y="Percent Abundance",
       fill="City group", x="") +
  guides(fill="none")
ggsave(filename = "DA_City_Cancer_Tissue_BOXPLOT2.pdf", width = 8.5, height = 11, dpi=300)
dev.off()

# IDs
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, color=City)) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_text(position = position_dodge(width=0.75),size=2, alpha=1, aes(group=City, label=ID), color="black") +
  geom_point(position = position_dodge(width=0.75),size=3.5, alpha=0.3, aes(group=City)) +
  theme_bw(base_size = 12) +
  theme(strip.text.x=element_text(size=13,colour="black")) +
  scale_color_manual(values=c("Rural"="coral2","Urban"="royalblue")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom",
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=9.5),
        axis.text.y = element_text(size=9), plot.title= element_text(size=16),
        panel.grid.minor.y= element_blank() ) +
  scale_x_discrete(expand=c(-0.2, 1)) +
  # scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance)+2,3))) +
  scale_y_sqrt(breaks=c(0, 0.1, 0.5, 1,1.5,2,2.5,3,3.5,4)) +
  labs(title= "Differently abundant Taxa (Cancer Tissue)", y="Percent Abundance",
       fill="City group", x="")
ggsave(filename = "Results/City_ONLY_cancer_tissue/DA_DESeq2/DA_City_Cancer_Tissue_ID.png", width = 7, height = 5, dpi=300)
dev.off()


system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/City_ONLY_cancer_tissue/DA_DESeq2/NB.txt ")


# for further comparisons ...
Genera.DESEQ2_City_Cancer<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##################### PCoA on EXTRA FACTORS (Only Cancer Tissue) ###########################

data <- subset_samples(complete_data, Tissue=="Cancer") # Only cancer Samples
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)

{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
# colnames(sample_data(data.sqrt_prop))

plot_ordination(data.sqrt_prop, ordBC, color = "Family_background") +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 13) + 
  stat_ellipse(size=0.2) +
  theme(plot.caption =element_text(size=8)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption ="NB: the difference seems to be caused only by that unique sample")
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_cancer_tissues/CITY.png", width = 6, height = 4.5, dpi=300)


plot_ordination(data.sqrt_prop, ordBC, color = "Cancer_Type") +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_cancer_tissues/Cancer_TYPE.png", width = 6, height = 4.5, dpi=300)

sample_data(data.sqrt_prop)$Cancer_grade<-gsub(" ","",sample_data(data.sqrt_prop)$Cancer_grade) # to remove any invisible space
plot_ordination(data.sqrt_prop, ordBC, color = "Cancer_grade") +
  scale_color_manual(values=c("I"="orange","II"="red","III"="red4")) +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_cancer_tissues/Cancer_Grade.png", width = 6, height = 4.5, dpi=300)

plot_ordination(data.sqrt_prop, ordBC, color = "Anticonceptional_pill") +
  scale_color_manual(values=c("no"="royalblue3","yes"="red2")) +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 12) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="Anticonceptional pill", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_cancer_tissues/ANTICONCEPTIONAL_PILL.png", width = 6, height = 4.5, dpi=300)


suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### ***** PREPARATION OF THE DATA (City  Healthy Tissue) #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  cat("\n\nDid you perform the decontamination and the removal of unsatureted samples already???\n\n")
  Sys.sleep(2)
}

# to reset the environment
remove<-ls()[! ls() %in% c("complete_Metadata", "complete_data", "unfiltered_data", "data_contam", "pairs", "who_unsaturated", "proof1","proof2")]
rm(list=remove)

data <- subset_samples(complete_data, Tissue=="Healthy") # Only Healthy Samples

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


###################### ABUNDANCES BAR PLOT  (City  Healthy Tissue) ##########################

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
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(City),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (only Healthy Tissues)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP_5_phyla_HC_urban_rural.pdf",width=8.5,height=11, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/City_ONLY_Healthy_tissue/Abundances/TOP5_phyla_average_abundance_Healthy.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
  tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$ID<-gsub("Patient_"," ",tabella$ID)
ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(City),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) +
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (only Healthy Tissues)", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="TOP_8_genera_HC_urban_rural.pdf",width=8.5,height=11,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


########################## ALFA DIVERSITY  (City Healthy Tissue) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="City", color="City")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$ID, obs$ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Tissue, New_data$ID),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
  pAlpha$data$City<-gsub("more_than_","> ",pAlpha$data$City, fixed=T)
  pAlpha$data$City<-gsub("less_than_","< ",pAlpha$data$City, fixed=T)
}
# with points only
pAlpha + theme_classic2() +
  geom_boxplot(aes(group=City), size =0.8, color = "black", alpha= 0) +
  scale_color_manual(values = c("Rural"="coral2", "Urban"="royalblue")) +
  labs(x="Tissue", title="Alpha diversity between cities (Healthy Tissues)") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = City), label="p.format", method = "wilcox.test", paired = F,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/City_ONLY_Healthy_tissue/Alfa_diversity_GENUS.png", width = 6,height =5.6, dpi=300)



##################### BETA DIVERSITY  (City Healthy Tissue) #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[1]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[1]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1]

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~City, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/City_ONLY_Healthy_tissue/Beta_div/Beta_div_permanova_Hellinger_Healthy.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
set.seed(1)
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$City)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/City_ONLY_Healthy_tissue/Beta_div/Beta_dispersion_Hellinger_on_genera.csv",quote=F,row.names = T)

suppressWarnings(rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV))


########################### PCoA  (City  Healthy Tissue) ########################

# on genera
{data.prop.labels<-data.genus.prop
data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.6, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/City_ONLY_Healthy_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_PATIENT.png", width = 6.5, height = 5, dpi=300)
# only points
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.2, alpha=1, shape=16) + geom_point(size=4.5, alpha=0.65) +
  theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F):",perm_g_H,"\nBeta dispersion Pr(>F):",disp_ASV$tab$`Pr(>F)`[1]))
ggsave(file="Results/City_ONLY_Healthy_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_POINT.png", width = 6.5, height = 5, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "City") +
  scale_color_manual(values=c("Urban"="royalblue","Rural"="coral2")) +
  geom_point(size=3.2, alpha=1, shape=16) + geom_point(size=4.5, alpha=0.65) +
  theme_classic(base_size = 13) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="City", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F):",perm_g_H,"\nBeta dispersion Pr(>F):",disp_ASV$tab$`Pr(>F)`[1]))
ggsave(file="Results/City_ONLY_Healthy_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_NO_ELLIPSES.png", width = 6.5, height = 5, dpi=300)


suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


#################### VEEN DIAGRAM (City Healthy Tissue) ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.1% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Rural<-subset_samples(data.venn, City=="Rural")
Rural<-as.character(tax_table(prune_taxa(taxa_sums(Rural)>0, Rural))[,"Genus"])

Urban<-subset_samples(data.venn, City=="Urban")
Urban<-as.character(tax_table(prune_taxa(taxa_sums(Urban)>0, Urban))[,"Genus"])


ONLY_IN_Urban<- Urban[! Urban %in% Rural]
ONLY_IN_Urban<- paste(ONLY_IN_Urban, collapse = ", ")
head(ONLY_IN_Urban)

ONLY_IN_Rural<- Rural[! Rural %in% Urban]
ONLY_IN_Rural<- paste(ONLY_IN_Rural, collapse = ", ")
head(ONLY_IN_Rural)


x<-list(Rural=Rural,Urban=Urban)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("coral2","royalblue")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = " Venn Diagram of genera in Healthy Tissue \n (only genera with minimal abundance > 0.1% \n at least in three samples)\n")
ggsave(filename = "Results/City_ONLY_Healthy_tissue/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()


suppressWarnings(rm(ONLY_IN_Urban, ONLY_IN_Rural, x, con, Urban, Rural, data.venn, who))



################### DA WITH DESEQ2 (City Healthy Tissue) #####################

if(! "proof1" %in% ls() | ! "proof2" %in% ls() ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~City)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("City", "Rural", "Urban"))
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
# write.csv(Res_tot, file="Results/City_ONLY_Healthy_tissue/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)

system(" echo 'No results here!' > Results/City_ONLY_Healthy_tissue/DA_DESeq2/Anything_there.txt ")
# 
# 
# # for further comparisons ...
# Genera.DESEQ2_City_Healthy<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##################### PCoA on EXTRA FACTORS (Only Healthy Tissue) ###########################

data <- subset_samples(complete_data, Tissue=="Healthy") # Only Healthy Samples
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)

{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
# colnames(sample_data(data.sqrt_prop))

plot_ordination(data.sqrt_prop, ordBC, color = "Family_background") +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_Healthy_tissues/CITY.png", width = 6, height = 4.5, dpi=300)


plot_ordination(data.sqrt_prop, ordBC, color = "Cancer_Type") +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_Healthy_tissues/Healthy_TYPE.png", width = 6, height = 4.5, dpi=300)

sample_data(data.sqrt_prop)$Cancer_grade<-gsub(" ","",sample_data(data.sqrt_prop)$Cancer_grade) # to remove any invisible space
plot_ordination(data.sqrt_prop, ordBC, color = "Cancer_grade") +
  scale_color_manual(values=c("I"="orange","II"="red","III"="red4")) +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_Healthy_tissues/Healthy_Grade.png", width = 6, height = 4.5, dpi=300)

plot_ordination(data.sqrt_prop, ordBC, color = "Anticonceptional_pill") +
  scale_color_manual(values=c("no"="royalblue3","yes"="red2")) +
  geom_point(size=3, alpha=1) + theme_classic(base_size = 12) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="Anticonceptional pill", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_EXTRA_FACTORS_only_Healthy_tissues/ANTICONCEPTIONAL_PILL.png", width = 6, height = 4.5, dpi=300)


suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


# ################### ***** PREPARATION OF THE DATA (Age_group  Cancer Tissue) #######################
# 
# if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
#   cat("\n\nDid you perform the decontamination and the removal of unsatureted samples already???\n\n")
#   Sys.sleep(2)
# }
# 
# # to reset the environment
# remove<-ls()[! ls() %in% c("complete_Metadata", "complete_data", "unfiltered_data", "data_contam", "pairs", "who_unsaturated", "proof1","proof2")]
# rm(list=remove)
# 
# data <- subset_samples(complete_data, Tissue=="Cancer") # Only cancer Samples
# 
# {data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
#   data.class = tax_glom(data, taxrank = "Class", NArm = F)
#   data.order = tax_glom(data, taxrank = "Order", NArm = F)
#   data.fam = tax_glom(data, taxrank = "Family", NArm = F)
#   data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
# }
# 
# { data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
#   data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
#   data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
#   data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
#   data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
#   data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
# }
# 
# { Taxa.genus<-as.data.frame(tax_table(data.genus))
#   Taxa.fam<-as.data.frame(tax_table(data.fam))
#   Taxa.phy<-as.data.frame(tax_table(data.phy))
#   Taxa.class<-as.data.frame(tax_table(data.class))
#   Taxa.order<-as.data.frame(tax_table(data.order))
# }
# 
# # adding informations to missing names
# taxa_temp<-Taxa.genus
# {for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
#   taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
#   for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
#     taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
#   for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
#     taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
#   for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
#     taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
#   for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
#     taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
#   Taxa.genus.update<-taxa_temp
# }
# 
# rm(taxa_temp)
# 
# 
# ###################### ABUNDANCES BAR PLOT  (Age_group  Cancer Tissue) ##########################
# 
# # choosing colors  (see grDevices::colors() )
# fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
# fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
# 
# # TOP 5 Phyla
# suppressWarnings(rm(top, others, tabella))
# {top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
#   prune.dat_top <- prune_taxa(top,data.phy.prop)
#   others<-taxa_names(data.phy.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.phy.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Phylum<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
# }
# tabella$ID<-gsub("Patient_"," ",tabella$ID)
# ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Phylum)) +
#   facet_grid(cols= vars(Age_group),scales = "free_x", space = "free_x") + 
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_5) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
#   labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (only Cancer Tissues)", caption = " 'Others' includes every phylum below rank 5 ")
# ggsave(file="Results/Age_group_ONLY_cancer_tissue/Abundances/TOP_5_phyla.png",width=7,height=5, dpi=300)
# dev.off()
# 
# # means of TOP5 phyla
# write.xlsx(file = "Results/Age_group_ONLY_cancer_tissue/Abundances/TOP5_phyla_average_abundance_CANCER.xlsx", row.names = F,
#            cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))
# 
# rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)
# 
# 
# # TOP 8 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# tabella$ID<-gsub("Patient_"," ",tabella$ID)
# ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Age_group),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") +
#   theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_8) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
#   labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (only Cancer Tissues)", caption = " 'Others' includes every genus below rank 8 ")
# ggsave(file="Results/Age_group_ONLY_cancer_tissue/Abundances/TOP_8_genera.png",width=7,height=5,dpi=300)
# dev.off()
# 
# suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))
# 
# 
# ########################## ALFA DIVERSITY  (Age_group Cancer Tissue) ############################
# 
# # no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# # moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )
# 
# pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Age_group", color="Age_group")
# pAlpha
# # plot_richness( ) compute diversity like estimate_diversity( )
# H<-dplyr::filter(pAlpha$data, variable=="Shannon")
# obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# # adding evenness
# {identical(H$ID, obs$ID) # TRUE
#   ev<-H
#   ev$value<-(H$value)/log((obs$value))
#   ev$variable<-rep("Evenness")
#   # updating and ordering samples for pairwise wilcoxon
#   New_data<-rbind.data.frame(obs,H,ev)
#   head(New_data)
#   New_data<-New_data[order(New_data$Tissue, New_data$ID),]
#   pAlpha$data<-New_data
#   pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
#   pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
#   pAlpha$data$Age_group<-gsub("more_than_","> ",pAlpha$data$Age_group, fixed=T)
#   pAlpha$data$Age_group<-gsub("less_than_","< ",pAlpha$data$Age_group, fixed=T)
# }
# # with points only
# pAlpha + theme_classic2() +
#   geom_boxplot(aes(group=Age_group), size =0.8, color = "black", alpha= 0) +
#   scale_color_manual(values = c("< 50"="coral2", "> 50"="royalblue")) +
#   labs(x="Tissue", title="Alpha diversity between age groups (Cancer Tissue)") +
#   guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
#   stat_compare_means(aes(group = Age_group), label="p.format", method = "wilcox.test", paired = F,
#                      label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
#   theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
#         axis.title.x = element_text(vjust = -1))
# ggsave(file="Results/Age_group_ONLY_cancer_tissue/Alfa_diversity_GENUS.png", width = 6,height =5.6, dpi=300)
# 
# 
# 
# ##################### BETA DIVERSITY  (Age_group Cancer Tissue) #######################
# 
# suppressWarnings(rm(ASV.prop))
# 
# {ASV.prop<-as.data.frame(otu_table(data.prop))
#   ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
#   ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
#   ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
#   ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
#   ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
# }
# 
# #### PERMANOVA
# metadata<-as(sample_data(data.prop),"data.frame")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# perm_ASV$aov.tab$`Pr(>F)`[1]
# perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
# perm_g<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# perm_g$aov.tab$`Pr(>F)`[1]
# perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1] 
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
# perm_f<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
# perm_o<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
# perm_c<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
# perm_p<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
# row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
# beta
# write.csv2(beta, file="Results/Age_group_ONLY_cancer_tissue/Beta_div/Beta_div_permanova_Hellinger_Cancer.csv",quote=F,row.names = T)
# 
# 
# # Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# # on ASV
# BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
# disper<-vegan::betadisper(BC.dist,metadata$Age_group)
# disp_ASV<-vegan::permutest(disper, permutations=9999)
# disp_ASV
# write.csv2(disp_ASV$tab, file="Results/Age_group_ONLY_cancer_tissue/Beta_div/Beta_dispersion_Hellinger_on_genera.csv",quote=F,row.names = T)
# 
# suppressWarnings(rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV))
# 
# 
# ########################### PCoA  (Age_group  Cancer Tissue) ########################
# 
# # on genera
# {data.prop.labels<-data.genus.prop
# data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
# DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
# ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
# eigval<-ordBC$values$Eigenvalues
# eigval<- round((eigval/sum(eigval))*100, 1)
# }
# sample_data(data.sqrt_prop)$Age_group<-gsub("more_than_","> ",sample_data(data.sqrt_prop)$Age_group, fixed=T)
# sample_data(data.sqrt_prop)$Age_group<-gsub("less_than_","< ",sample_data(data.sqrt_prop)$Age_group, fixed=T)
# plot_ordination(data.sqrt_prop, ordBC, color = "Age_group") +
#   scale_color_manual(values=c("> 50"="royalblue","< 50"="coral2")) +
#   geom_point(size=3.6, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
#        color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/Age_group_ONLY_cancer_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_PATIENT.png", width = 6.5, height = 5, dpi=300)
# # only points
# plot_ordination(data.sqrt_prop, ordBC, color = "Age_group") +
#   scale_color_manual(values=c("> 50"="royalblue","< 50"="coral2")) +
#   geom_point(size=3.5, alpha=1, shape=16) + geom_point(size=4.2, alpha=0.65) +
#   theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
#        color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        subtitle = "PERMANOVA Pr(>F):",perm_g_H)
#   ggsave(file="Results/Age_group_ONLY_cancer_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_POINT.png", width = 6.5, height = 5, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Age_group") +
#   scale_color_manual(values=c("> 50"="royalblue","< 50"="coral2")) +
#   geom_point(size=3.5, alpha=1, shape=16) + geom_point(size=4.2, alpha=0.65) +
#   theme_classic(base_size = 13) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
#        color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        subtitle = "PERMANOVA Pr(>F):",perm_g_H)
# ggsave(file="Results/Age_group_ONLY_cancer_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_NO_ELLIPSES.png", width = 6, height = 5, dpi=300)
# 
# 
# suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))
# 
# 
# #################### VEEN DIAGRAM (Age_group Cancer Tissue) ##########################
# 
# data.genus.temp<-data.genus.prop
# tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
# data.venn<-data.genus.temp
# 
# ### abundance AND prevalence filter (0.5% abundance at least in 3 sample)
# who<-as.data.frame(otu_table(data.venn))
# who<-apply(who, MARGIN=2, function(x) ifelse(x>0.5, 1, 0)) # if more than 0.5 --> "a point"
# who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
# who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
# data.venn<-subset_taxa(data.venn, ! Genus %in% who)
# 
# 
# less_than_50<-subset_samples(data.venn, Age_group=="less_than_50")
# less_than_50<-as.character(tax_table(prune_taxa(taxa_sums(less_than_50)>0, less_than_50))[,"Genus"])
# 
# more_than_50<-subset_samples(data.venn, Age_group=="more_than_50")
# more_than_50<-as.character(tax_table(prune_taxa(taxa_sums(more_than_50)>0, more_than_50))[,"Genus"])
# 
# 
# ONLY_IN_more_than_50<- more_than_50[! more_than_50 %in% less_than_50]
# ONLY_IN_more_than_50<- paste(ONLY_IN_more_than_50, collapse = ", ")
# head(ONLY_IN_more_than_50)
# 
# ONLY_IN_less_than_50<- less_than_50[! less_than_50 %in% more_than_50]
# ONLY_IN_less_than_50<- paste(ONLY_IN_less_than_50, collapse = ", ")
# head(ONLY_IN_less_than_50)
# 
# 
# x<-list(less_than_50=less_than_50,more_than_50=more_than_50)
# ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
#        fill_color = c("chartreuse","coral")) +
#   theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
#   labs(title = " Venn Diagram of genera in Cancer Tissue \n (only genera with minimal abundance > 0.5% \n at least in three samples)\n")
# ggsave(filename = "Results/Age_group_ONLY_cancer_tissue/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
# dev.off()
# 
# 
# suppressWarnings(rm(ONLY_IN_more_than_50, ONLY_IN_less_than_50, x, con, more_than_50, less_than_50, data.venn, who))
# 
# 
# 
# ################### DA WITH DESEQ2 (Age_group Cancer Tissue) #####################
# 
# if(! "proof1" %in% ls() ){
#   cat("\nDid you perform the filtering step yet?\n", fill = T)
#   Sys.sleep(2)
# }
# 
# 
# ##### STARTING THE DIFFERENTIAL ANALYSIS
# suppressWarnings(rm(data_pruned, data.genus_pruned))
# data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# # Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
# 
# Table_tot<-NULL
# Res_tot<-NULL
# 
# for( t in c("Genus","Family","Class","Order","Phylum") ){
#   cat("\nWorking on",t,"level...\n")
#   suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
#   d <- tax_glom(data_pruned, taxrank = t, NArm = F)
#   d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
#   
#   if(t=="Genus"){ # updating missing names (NA and uncultured) but only for genus level
#     taxa_temp<-as.data.frame(tax_table(d))
#     for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
#       taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
#     for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
#       taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
#     for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
#       taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
#     for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
#       taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
#     for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
#       taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
#     tax_table(d)<-as.matrix(taxa_temp)
#     tax_table(d.prop)<-as.matrix(taxa_temp)
#     rm(taxa_temp) }
#   
#   ### starting the analysis
#   DEseq_data<-phyloseq_to_deseq2(d, ~Age_group)
#   DE<-DESeq(DEseq_data)
#   res<-results(DE, contrast= c("Age_group", "less_than_50", "more_than_50"))
#   res = res[order(res$padj, na.last=NA), ]
#   res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
#   res<-res[res$baseMean > 50, ] # arbitrary threshold to avoid the most noisy result
#   res
#   if(length(res$log2FoldChange)>0){ # if there are results...
#     cat(paste(length(res$log2FoldChange),"results for the",t,"level\n"))
#     Sys.sleep(1)
#     r<-as.data.frame(res)
#     r$ASV<-row.names(r)
#     Taxa.d<-as.data.frame(tax_table(d))
#     Taxa.d$ASV<-row.names(Taxa.d)
#     r<-dplyr::left_join(r, Taxa.d, by="ASV")
#     r$Kingdom<-NULL
#     r$Species<-NULL
#     assign(paste(t,"results",sep="_"), r)
#     r_level<-r
#     r_level[, "Taxon"]<- rep(t)
#     Res_tot<-rbind.data.frame(Res_tot,r_level)
#     ### single box plots
#     target<-r[[t]]
#     colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
#     target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
#     Table_DE<-psmelt(target)
#     colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
#     Table_DE$ASV<-NULL
#     # Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
#     assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
#     ### appending to unique box plot
#     index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t)
#     index<- index[-length(index)] # removing the last index, regarding the taxa of interest
#     Table_DE[,index]<-NULL
#     Table_DE$Taxa<-t
#     colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
#     Table_tot<-rbind.data.frame(Table_tot, Table_DE)
#   } else {
#     cat("Any results for the",t,"level\n")
#     
#     Sys.sleep(1)
#   }
# }
# 
# View(Res_tot)
# write.csv(Res_tot, file="Results/Age_group_ONLY_cancer_tissue/DA_DESeq2/Every_result_DESeq2_.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/Age_group_ONLY_cancer_tissue/DA_DESeq2/Every_result_DESeq2_.xlsx", showNA = F, col.names = T)
# 
# Table_tot$Age_group<-gsub("less_than_50","< 50", Table_tot$Age_group, fixed=T)
# Table_tot$Age_group<-gsub("more_than_50","> 50", Table_tot$Age_group, fixed=T)
# ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Age_group)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   scale_fill_manual(values=c("< 50"="chartreuse","> 50"="coral")) +
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-28, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
#         axis.text.x = element_text(angle = 28, vjust=1, hjust=1, size=11), 
#         axis.text.y = element_text(size=9.5), plot.title= element_text(size=16),
#         panel.grid.minor.y= element_blank() ) +   
#   scale_x_discrete(expand=c(-0.2, 1)) +
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,22,2), 24.5, seq(27,max(Table_tot$Abundance)+3,3))) +
#   labs(title= "Differently abundant Taxa (Cancer Tissue)", y="Percent Abundance", 
#        fill="Age group", x="")
# ggsave(filename = "Results/Age_group_ONLY_cancer_tissue/DA_DESeq2/DA_Age_group_Cancer_Tissue.png", width = 11, height = 6.5, dpi=300)
# dev.off()
# 
# system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Age_group_ONLY_cancer_tissue/DA_DESeq2/NB.txt ")
# 
# 
# # for further comparisons ...
# Genera.DESEQ2_Age_group_Cancer<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])
# 
# 


# ################### ***** PREPARATION OF THE DATA (Age_group  Healthy Tissue) #######################
# 
# data <- subset_samples(complete_data, Tissue=="Healthy") # Only healthy Samples
# 
# {data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
#   data.class = tax_glom(data, taxrank = "Class", NArm = F)
#   data.order = tax_glom(data, taxrank = "Order", NArm = F)
#   data.fam = tax_glom(data, taxrank = "Family", NArm = F)
#   data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
# }
# 
# { data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
#   data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
#   data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
#   data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
#   data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
#   data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
# }
# 
# { Taxa.genus<-as.data.frame(tax_table(data.genus))
#   Taxa.fam<-as.data.frame(tax_table(data.fam))
#   Taxa.phy<-as.data.frame(tax_table(data.phy))
#   Taxa.class<-as.data.frame(tax_table(data.class))
#   Taxa.order<-as.data.frame(tax_table(data.order))
# }
# 
# # adding informations to missing names
# taxa_temp<-Taxa.genus
# {for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
#   taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
#   for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
#     taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
#   for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
#     taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
#   for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
#     taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
#   for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
#     taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
#   Taxa.genus.update<-taxa_temp
# }
# 
# rm(taxa_temp)
# 
# 
# ###################### ABUNDANCES BAR PLOT  (Age_group  Healthy Tissue) ##########################
# 
# # choosing colors  (see grDevices::colors() )
# fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
# fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
# 
# # TOP 5 Phyla
# suppressWarnings(rm(top, others, tabella))
# {top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
#   prune.dat_top <- prune_taxa(top,data.phy.prop)
#   others<-taxa_names(data.phy.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.phy.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Phylum<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
# }
# tabella$ID<-gsub("Patient_"," ",tabella$ID)
# ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Phylum)) +
#   facet_grid(cols= vars(Age_group),scales = "free_x", space = "free_x") + 
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_5) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
#   labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (only Healthy Tissues)", caption = " 'Others' includes every phylum below rank 5 ")
# ggsave(file="Results/Age_group_ONLY_Healthy_tissue/Abundances/TOP_5_phyla.png",width=7,height=5, dpi=300)
# dev.off()
# 
# # means of TOP5 phyla
# write.xlsx(file = "Results/Age_group_ONLY_Healthy_tissue/Abundances/TOP5_phyla_average_abundance_HEALTHY.xlsx", row.names = F,
#            cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))
# 
# rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)
# 
# 
# # TOP 8 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.",tabella$Genus)
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# tabella$ID<-gsub("Patient_"," ",tabella$ID)
# ggplot(data=tabella, aes(x=ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Age_group),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") +
#   theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_8) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
#   labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (only Healthy Tissues)", caption = " 'Others' includes every genus below rank 8 ")
# ggsave(file="Results/Age_group_ONLY_Healthy_tissue/Abundances/TOP_8_genera.png",width=7,height=5,dpi=300)
# dev.off()
# 
# suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))
# 
# 
# ########################## ALFA DIVERSITY  (Age_group Healthy Tissue) ############################
# 
# # no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# # moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )
# 
# pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Age_group", color="Age_group")
# pAlpha
# # plot_richness( ) compute diversity like estimate_diversity( )
# H<-dplyr::filter(pAlpha$data, variable=="Shannon")
# obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# # adding evenness
# {identical(H$ID, obs$ID) # TRUE
#   ev<-H
#   ev$value<-(H$value)/log((obs$value))
#   ev$variable<-rep("Evenness")
#   # updating and ordering samples for pairwise wilcoxon
#   New_data<-rbind.data.frame(obs,H,ev)
#   head(New_data)
#   New_data<-New_data[order(New_data$Tissue, New_data$ID),]
#   pAlpha$data<-New_data
#   pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
#   pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
#   pAlpha$data$Age_group<-gsub("more_than_","> ",pAlpha$data$Age_group, fixed=T)
#   pAlpha$data$Age_group<-gsub("less_than_","< ",pAlpha$data$Age_group, fixed=T)
# }
# # with points only
# pAlpha + theme_classic2() +
#   geom_boxplot(aes(group=Age_group), size =0.8, color = "black", alpha= 0) +
#   scale_color_manual(values = c("< 50"="coral2", "> 50"="royalblue")) +
#   labs(x="Tissue", title="Alpha diversity between age groups (Healthy Tissue)") +
#   guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
#   stat_compare_means(aes(group = Age_group), label="p.format", method = "wilcox.test", paired = F,
#                      label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
#   theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
#         axis.title.x = element_text(vjust = -1))
# ggsave(file="Results/Age_group_ONLY_Healthy_tissue/Alfa_diversity_GENUS.png", width = 6,height =5.6, dpi=300)
# 
# 
# ##################### BETA DIVERSITY  (Age_group Healthy Tissue) #######################
# 
# suppressWarnings(rm(ASV.prop))
# 
# {ASV.prop<-as.data.frame(otu_table(data.prop))
#   ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
#   ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
#   ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
#   ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
#   ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
# }
# 
# #### PERMANOVA
# metadata<-as(sample_data(data.prop),"data.frame")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# perm_ASV$aov.tab$`Pr(>F)`[1]
# perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
# perm_g<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# perm_g$aov.tab$`Pr(>F)`[1]
# perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1] 
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
# perm_f<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
# perm_o<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
# perm_c<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
# perm_p<- vegan::adonis(sample_OTU ~Age_group, data=metadata, permutations = 9999, method="euclidean")
# 
# beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
# row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
# beta
# write.csv2(beta, file="Results/Age_group_ONLY_Healthy_tissue/Beta_div/Beta_div_permanova_Hellinger_Healthy.csv",quote=F,row.names = T)
# 
# 
# # Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# # on ASV
# BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
# disper<-vegan::betadisper(BC.dist,metadata$Age_group)
# disp_ASV<-vegan::permutest(disper, permutations=9999)
# disp_ASV
# write.csv2(disp_ASV$tab, file="Results/Age_group_ONLY_Healthy_tissue/Beta_div/Beta_dispersion_Hellinger_on_genera.csv",quote=F,row.names = T)
# 
# suppressWarnings(rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV))
# 
# 
# ########################### PCoA  (Age_group  Healthy Tissue) ########################
# 
# # on genera
# {data.prop.labels<-data.genus.prop
# data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
# DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
# ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
# eigval<-ordBC$values$Eigenvalues
# eigval<- round((eigval/sum(eigval))*100, 1)
# }
# sample_data(data.sqrt_prop)$Age_group<-gsub("more_than_","> ",sample_data(data.sqrt_prop)$Age_group, fixed=T)
# sample_data(data.sqrt_prop)$Age_group<-gsub("less_than_","< ",sample_data(data.sqrt_prop)$Age_group, fixed=T)
# plot_ordination(data.sqrt_prop, ordBC, color = "Age_group") +
#   scale_color_manual(values=c("> 50"="royalblue","< 50"="coral2")) +
#   geom_point(size=3.6, alpha=1) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_data(data.sqrt_prop)$ID), color="black", size=1.8, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
#        color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/Age_group_ONLY_Healthy_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_PATIENT.png", width = 6.5, height = 5, dpi=300)
# # only points
# plot_ordination(data.sqrt_prop, ordBC, color = "Age_group") +
#   scale_color_manual(values=c("> 50"="royalblue","< 50"="coral2")) +
#   geom_point(size=3.5, alpha=1, shape=16) + geom_point(size=4.2, alpha=0.65) +
#   theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
#        color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        subtitle = "PERMANOVA Pr(>F):",perm_g_H)
# ggsave(file="Results/Age_group_ONLY_Healthy_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_POINT.png", width = 6.5, height = 5, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Age_group") +
#   scale_color_manual(values=c("> 50"="royalblue","< 50"="coral2")) +
#   geom_point(size=3.5, alpha=1, shape=16) + geom_point(size=4.2, alpha=0.65) +
#   theme_classic(base_size = 13) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
#        color="Age", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        subtitle = "PERMANOVA Pr(>F):",perm_g_H)
# ggsave(file="Results/Age_group_ONLY_Healthy_tissue/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_NO_ELLIPSES.png", width = 6, height = 5, dpi=300)
# 
# 
# suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))
# 
# 
# #################### VEEN DIAGRAM (Age_group Healthy Tissue) ##########################
# 
# data.genus.temp<-data.genus.prop
# tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
# data.venn<-data.genus.temp
# 
# ### abundance AND prevalence filter (0.5% abundance at least in 3 sample)
# who<-as.data.frame(otu_table(data.venn))
# who<-apply(who, MARGIN=2, function(x) ifelse(x>0.5, 1, 0)) # if more than 0.5 --> "a point"
# who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
# who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
# data.venn<-subset_taxa(data.venn, ! Genus %in% who)
# 
# 
# less_than_50<-subset_samples(data.venn, Age_group=="less_than_50")
# less_than_50<-as.character(tax_table(prune_taxa(taxa_sums(less_than_50)>0, less_than_50))[,"Genus"])
# 
# more_than_50<-subset_samples(data.venn, Age_group=="more_than_50")
# more_than_50<-as.character(tax_table(prune_taxa(taxa_sums(more_than_50)>0, more_than_50))[,"Genus"])
# 
# 
# ONLY_IN_more_than_50<- more_than_50[! more_than_50 %in% less_than_50]
# ONLY_IN_more_than_50<- paste(ONLY_IN_more_than_50, collapse = ", ")
# head(ONLY_IN_more_than_50)
# 
# ONLY_IN_less_than_50<- less_than_50[! less_than_50 %in% more_than_50]
# ONLY_IN_less_than_50<- paste(ONLY_IN_less_than_50, collapse = ", ")
# head(ONLY_IN_less_than_50)
# 
# 
# x<-list(less_than_50=less_than_50,more_than_50=more_than_50)
# ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
#        fill_color = c("chartreuse","coral")) +
#   theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
#   labs(title = " Venn Diagram of genera in Healthy Tissue \n (only genera with minimal abundance > 0.5% \n at least in three samples)\n")
# ggsave(filename = "Results/Age_group_ONLY_Healthy_tissue/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
# dev.off()
# 
# 
# suppressWarnings(rm(ONLY_IN_more_than_50, ONLY_IN_less_than_50, x, con, more_than_50, less_than_50, data.venn, who))
# 
# 
# 
# ################### DA WITH DESEQ2 (Age_group Healthy Tissue) #####################
# 
# if(! "proof1" %in% ls() ){
#   cat("\nDid you perform the filtering step yet?\n", fill = T)
#   Sys.sleep(2)
# }
# 
# 
# ##### STARTING THE DIFFERENTIAL ANALYSIS
# suppressWarnings(rm(data_pruned, data.genus_pruned))
# data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# # Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
# 
# Table_tot<-NULL
# Res_tot<-NULL
# 
# for( t in c("Genus","Family","Class","Order","Phylum") ){
#   cat("\nWorking on",t,"level...\n")
#   suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
#   d <- tax_glom(data_pruned, taxrank = t, NArm = F)
#   d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
#   
#   if(t=="Genus"){ # updating missing names (NA and uncultured) but only for genus level
#     taxa_temp<-as.data.frame(tax_table(d))
#     for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
#       taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
#     for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
#       taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
#     for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
#       taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
#     for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
#       taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
#     for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
#       taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
#     tax_table(d)<-as.matrix(taxa_temp)
#     tax_table(d.prop)<-as.matrix(taxa_temp)
#     rm(taxa_temp) }
#   
#   ### starting the analysis
#   DEseq_data<-phyloseq_to_deseq2(d, ~Age_group)
#   DE<-DESeq(DEseq_data)
#   res<-results(DE, contrast= c("Age_group", "less_than_50", "more_than_50"))
#   res = res[order(res$padj, na.last=NA), ]
#   res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
#   res<-res[res$baseMean > 50, ] # arbitrary threshold to avoid the most noisy result
#   res
#   if(length(res$log2FoldChange)>0){ # if there are results...
#     cat(paste(length(res$log2FoldChange),"results for the",t,"level\n"))
#     Sys.sleep(1)
#     r<-as.data.frame(res)
#     r$ASV<-row.names(r)
#     Taxa.d<-as.data.frame(tax_table(d))
#     Taxa.d$ASV<-row.names(Taxa.d)
#     r<-dplyr::left_join(r, Taxa.d, by="ASV")
#     r$Kingdom<-NULL
#     r$Species<-NULL
#     assign(paste(t,"results",sep="_"), r)
#     r_level<-r
#     r_level[, "Taxon"]<- rep(t)
#     Res_tot<-rbind.data.frame(Res_tot,r_level)
#     ### single box plots
#     target<-r[[t]]
#     colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
#     target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
#     Table_DE<-psmelt(target)
#     colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
#     Table_DE$ASV<-NULL
#     # Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
#     assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
#     ### appending to unique box plot
#     index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t)
#     index<- index[-length(index)] # removing the last index, regarding the taxa of interest
#     Table_DE[,index]<-NULL
#     Table_DE$Taxa<-t
#     colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
#     Table_tot<-rbind.data.frame(Table_tot, Table_DE)
#   } else {
#     cat("Any results for the",t,"level\n")
#     
#     Sys.sleep(1)
#   }
# }
# 
# View(Res_tot)
# write.csv(Res_tot, file="Results/Age_group_ONLY_Healthy_tissue/DA_DESeq2/Every_result_DESeq2_.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/Age_group_ONLY_Healthy_tissue/DA_DESeq2/Every_result_DESeq2_.xlsx", showNA = F, col.names = T)
# 
# Table_tot$Age_group<-gsub("less_than_50","< 50", Table_tot$Age_group, fixed=T)
# Table_tot$Age_group<-gsub("more_than_50","> 50", Table_tot$Age_group, fixed=T)
# ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Age_group)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   scale_fill_manual(values=c("< 50"="chartreuse","> 50"="coral")) +
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-28, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
#         axis.text.x = element_text(angle = 28, vjust=1, hjust=1, size=11), 
#         axis.text.y = element_text(size=9.5), plot.title= element_text(size=16),
#         panel.grid.minor.y= element_blank() ) +   
#   scale_x_discrete(expand=c(-0.2, 1)) +
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance)+2,2))) +
#   labs(title= "Differently abundant Taxa (Healthy Tissue)", y="Percent Abundance", 
#        fill="Age group", x="")
# ggsave(filename = "Results/Age_group_ONLY_Healthy_tissue/DA_DESeq2/DA_Age_group_Healthy_Tissue.png", width = 9.5, height = 6.5, dpi=300)
# dev.off()
# 
# system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Age_group_ONLY_Healthy_tissue/DA_DESeq2/NB.txt ")
# system(" echo 'The Rhyzobiales are plant related bacteria --> they are contaminants!' > Results/Age_group_ONLY_Healthy_tissue/DA_DESeq2/BEWARE_OF_RHYZOBIALES.txt ")
# 
# 
# # for further comparisons ...
# Genera.DESEQ2_Age_group_HEALTHY<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])
# 

############ PICURST ############ 
#healthy vs cancer
a <- read.delim("QIIME/exported-feature-table/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
# a <- read.delim("picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz")
colnames(a)[1:4]
if (!require('stringr')) {install.packages('stringr')}
temp<-substring(colnames(a)[stringr::str_starts(colnames(a),pattern = "X")],first = 10)
temp
colnames(a)[stringr::str_starts(colnames(a),pattern = "X")]<-temp

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]

Metadata <- as.data.frame(read.csv("Metadata.csv")) # import metadata file to regroups
head(Metadata)
Metadata$FASTQ_ID<-as.character(Metadata$FASTQ_ID)
rownames(Metadata)<-Metadata$FASTQ_ID

colnames(a)
extract_PN <- function(name) {
  sub(".*(DM\\d+).*", "\\1", name)
}

new_colnames <- sapply(colnames(a), extract_PN)
colnames(a) <- new_colnames
Metadata<-Metadata[colnames(a),]
if(identical(Metadata$FASTQ_ID,colnames(a))){a<-rbind.data.frame(Metadata$Tissue,a)} else {cat ("\nsomething gone wrong...\n\n")}
colnames(a)
a<-cbind.data.frame(c("Tissue",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("Tissue", Descriptions$pathway)']<-"Sample"
head(a,n=3)

write.table(a,file="Output_for_LefSe_C_H.tsv",row.names = F,quote = F, sep="\t")

#urban vs rural
a <- read.delim("QIIME/exported-feature-table/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
# a <- read.delim("picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz")
colnames(a)[1:4]
if (!require('stringr')) {install.packages('stringr')}
temp<-substring(colnames(a)[stringr::str_starts(colnames(a),pattern = "X")],first = 10)
temp
colnames(a)[stringr::str_starts(colnames(a),pattern = "X")]<-temp

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]

Metadata <- as.data.frame(read.csv("Metadata.csv")) # import metadata file to regroups
head(Metadata)
Metadata$FASTQ_ID<-as.character(Metadata$FASTQ_ID)
rownames(Metadata)<-Metadata$FASTQ_ID

colnames(a)
extract_PN <- function(name) {
  sub(".*(DM\\d+).*", "\\1", name)
}

new_colnames <- sapply(colnames(a), extract_PN)
colnames(a) <- new_colnames
Metadata<-Metadata[colnames(a),]
if(identical(Metadata$FASTQ_ID,colnames(a))){a<-rbind.data.frame(Metadata$City,a)} else {cat ("\nsomething gone wrong...\n\n")}
colnames(a)
a<-cbind.data.frame(c("City",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("City", Descriptions$pathway)']<-"Sample"
head(a,n=3)

write.table(a,file="Output_for_LefSe_city.tsv",row.names = F,quote = F, sep="\t")

############  Picrust PLOT ############ 
Significative_functions_LEFSE <- read.delim("Result_LEFSE.res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]

for(x in 1:length(Significative_functions_LEFSE$Pathway)) {
  if(stringr::str_detect(Significative_functions_LEFSE$Pathway[x],"^f_[0-9]")) {Significative_functions_LEFSE$Pathway[x]<-gsub("f_","",Significative_functions_LEFSE$Pathway[x]) }
}

Significative_functions_LEFSE$Pathway<-gsub("_", "-", Significative_functions_LEFSE$Pathway)
head(Descriptions)
vett<-Significative_functions_LEFSE$Pathway
prova<-subset(Descriptions, pathway %in% vett)

Significative_functions_LEFSE$Pathway<-prova$description

colnames(Significative_functions_LEFSE)[colnames(Significative_functions_LEFSE)=="Pathway"]<-"MetaCyc_ID"
Significative_functions_LEFSE$MetaCyc_ID[is.na(Significative_functions_LEFSE$MetaCyc_ID)] # there has not to be NA here

# modifing names for the plot
Significative_functions_LEFSE$MetaCyc_ID<-gsub("&beta;-","", Significative_functions_LEFSE$MetaCyc_ID, fixed = T)
Significative_functions_LEFSE$MetaCyc_ID<-gsub("_"," ", Significative_functions_LEFSE$MetaCyc_ID, fixed = T)
{Significative_functions_LEFSE$MetaCyc_ID<-paste("",Significative_functions_LEFSE$MetaCyc_ID,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$MetaCyc_ID<-factor(Significative_functions_LEFSE$MetaCyc_ID, levels = Significative_functions_LEFSE$MetaCyc_ID) # to prevent alphabetical re-sorting
}
# plotting results only over 3
Significative_functions_LEFSE<- Significative_functions_LEFSE[Significative_functions_LEFSE$logLDA_score>2, ]
write.xlsx(Significative_functions_LEFSE,file = "Significative_functions_LEFSE.xlsx")
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Healthy"] <- - Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Cancer"]

ggplot(data=Significative_functions_LEFSE, aes(y=MetaCyc_ID, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=MetaCyc_ID), hjust = ifelse(Significative_functions_LEFSE$logLDA_score>0,1,0),size=4)+
  scale_fill_manual(values = c("Cancer"="lightblue","Healthy"="coral"),labels=c("Cancer","Healthy")) +
  theme_classic(base_size = 10) +
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.4)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 7)) +
  scale_x_continuous(breaks = c(1,2,3))+
  theme(legend.position = "bottom")

ggsave(filename = "PICRUST2_LEFSE.png", width = 15, height = 12, dpi = 300)

#urban rural
Significative_functions_LEFSE <- read.delim("Result_LEFSE_city.res", header=FALSE)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
head(Significative_functions_LEFSE)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]

for(x in 1:length(Significative_functions_LEFSE$Pathway)) {
  if(stringr::str_detect(Significative_functions_LEFSE$Pathway[x],"^f_[0-9]")) {Significative_functions_LEFSE$Pathway[x]<-gsub("f_","",Significative_functions_LEFSE$Pathway[x]) }
}

Significative_functions_LEFSE$Pathway<-gsub("_", "-", Significative_functions_LEFSE$Pathway)
head(Descriptions)
vett<-Significative_functions_LEFSE$Pathway
prova<-subset(Descriptions, pathway %in% vett)

Significative_functions_LEFSE$Pathway<-prova$description

colnames(Significative_functions_LEFSE)[colnames(Significative_functions_LEFSE)=="Pathway"]<-"MetaCyc_ID"
Significative_functions_LEFSE$MetaCyc_ID[is.na(Significative_functions_LEFSE$MetaCyc_ID)] # there has not to be NA here

# modifing names for the plot
Significative_functions_LEFSE$MetaCyc_ID<-gsub("&beta;-","", Significative_functions_LEFSE$MetaCyc_ID, fixed = T)
Significative_functions_LEFSE$MetaCyc_ID<-gsub("_"," ", Significative_functions_LEFSE$MetaCyc_ID, fixed = T)
{Significative_functions_LEFSE$MetaCyc_ID<-paste("",Significative_functions_LEFSE$MetaCyc_ID,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$MetaCyc_ID<-factor(Significative_functions_LEFSE$MetaCyc_ID, levels = Significative_functions_LEFSE$MetaCyc_ID) # to prevent alphabetical re-sorting
}
# plotting results only over 2
Significative_functions_LEFSE<- Significative_functions_LEFSE[Significative_functions_LEFSE$logLDA_score>2, ]
write.xlsx(Significative_functions_LEFSE,file = "Significative_functions_LEFSE_city.xlsx")
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Urban"] <- - Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Rural"]

ggplot(data = Significative_functions_LEFSE, 
       aes(y = MetaCyc_ID, x = as.numeric(logLDA_score), fill = Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  labs(x = "log LDA score", y = "", fill = "") +
  geom_text(aes(x = 0, label = MetaCyc_ID), 
            hjust = ifelse(Significative_functions_LEFSE$logLDA_score > 0, 1, 0), size = 4) +
  scale_fill_manual(values = c("Urban" = "lightblue", "Rural" = "coral"),
                    labels = c("Rural", "Urban")) +
  theme_classic(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),     # Rimuove solo le griglie verticali
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "bottom"
  ) +
  scale_x_continuous(breaks = c(1, 2, 3))


ggsave(filename = "PICRUST2_LEFSE_city_2.png", width = 15, height = 12, dpi = 300)


###### HEATMAP AND CLUSTER ON PHYLA ###### 
library(pheatmap)

otu_table(data.phy.prop)
sample_data(data.phy.prop)

ps_rel <- data.phy.prop
top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.phy.prop)

ps_phylum_top <- prune.dat_top
otu_table(prune.dat_top)
mat <- as.matrix(otu_table(ps_phylum_top))

tax <- tax_table(ps_phylum_top)
rownames(mat) <- tax[rownames(mat), "Phylum"]

wd_ids <- sample_data(ps_phylum_top)$Work_ID
colnames(mat) <- wd_ids

png("heatmap_phylum_workID.png", width = 3000, height = 2000, res = 300)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 8)      

dev.off() 


###### HEATMAP AND CLUSTER ON GENERA ###### 
otu_table(data.phy.prop)
sample_data(data.phy.prop)

ps_rel <- data.phy.prop
top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
prune.dat_top <- prune_taxa(top,data.genus.prop)

ps_phylum_top <- prune.dat_top
otu_table(prune.dat_top)
mat <- as.matrix(otu_table(ps_phylum_top))

tax <- tax_table(ps_phylum_top)
rownames(mat) <- tax[rownames(mat), "Genus"]

wd_ids <- sample_data(ps_phylum_top)$Work_ID
colnames(mat) <- wd_ids

png("heatmap_genus_workID.png", width = 3000, height = 2000, res = 300)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 8)       

dev.off() 


###### HEATMAP AND CLUSTER ON TOP 5 PHYLA - CANCER ###### 
library(pheatmap)
data <- subset_samples(complete_data, Tissue=="Cancer")
data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)

sd <- sample_data(data.phy.prop)
sd$ID <- gsub("Patient_", "", sd$ID)
sample_data(data.phy.prop) <- sd

otu_table(data.phy.prop)
sample_data(data.phy.prop)

ps_rel <- data.phy.prop
top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.phy.prop)

ps_phylum_top <- prune.dat_top
otu_table(prune.dat_top)
mat <- as(otu_table(ps_phylum_top), "matrix")
summary(as.vector(mat))
head(mat[,1:5])

tax <- tax_table(ps_phylum_top)
rownames(mat) <- tax[rownames(mat), "Phylum"]

wd_ids <- sample_data(ps_phylum_top)$ID
colnames(mat) <- wd_ids

annot_col <- data.frame(
  City = sample_data(ps_phylum_top)$City  
)

rownames(annot_col) <- sample_data(ps_phylum_top)$ID
annot_col$City <- gsub("Urban", "U", annot_col$City)
annot_col$City <- gsub("Rural", "R", annot_col$City)

graphics.off()
png("heatmap_phylum_cancer_city.png", width = 3000, height = 2000, res = 300)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annot_col,
         fontsize_col = 8)

dev.off() 


###### HEATMAP AND CLUSTER ON TOP 5 PHYLA - HEALTHY ###### 
data <- subset_samples(complete_data, Tissue=="Healthy")
data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)

sd <- sample_data(data.phy.prop)
sd$ID <- gsub("Patient_", "", sd$ID)
sample_data(data.phy.prop) <- sd

otu_table(data.phy.prop)
sample_data(data.phy.prop)

ps_rel <- data.phy.prop
top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.phy.prop)

ps_phylum_top <- prune.dat_top
otu_table(prune.dat_top)
mat <- as.matrix(otu_table(ps_phylum_top))

tax <- tax_table(ps_phylum_top)
rownames(mat) <- tax[rownames(mat), "Phylum"]

wd_ids <- sample_data(ps_phylum_top)$ID
colnames(mat) <- wd_ids

annot_col <- data.frame(
  City = sample_data(ps_phylum_top)$City 
)

rownames(annot_col) <- sample_data(ps_phylum_top)$ID
annot_col$City <- gsub("Urban", "U", annot_col$City)
annot_col$City <- gsub("Rural", "R", annot_col$City)

png("heatmap_phylum_healthy_city.png", width = 3000, height = 2000, res = 300)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annot_col, 
         fontsize_col = 7)

dev.off() 

###### HEATMAP AND CLUSTER ON TOP 8 GENUS - HEALTHY ###### 
data <- subset_samples(complete_data, Tissue=="Healthy")
data.phy = tax_glom(data, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)

sd <- sample_data(data.phy.prop)
sd$ID <- gsub("Patient_", "", sd$ID)
sample_data(data.phy.prop) <- sd

otu_table(data.phy.prop)
sample_data(data.phy.prop)

ps_rel <- data.phy.prop
top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:8]
prune.dat_top <- prune_taxa(top,data.phy.prop)

ps_phylum_top <- prune.dat_top
otu_table(prune.dat_top)
mat <- as(otu_table(ps_phylum_top), "matrix")

tax <- tax_table(ps_phylum_top)
rownames(mat) <- tax[rownames(mat), "Genus"]

wd_ids <- sample_data(ps_phylum_top)$ID
colnames(mat) <- wd_ids


annot_col <- data.frame(
  City = sample_data(ps_phylum_top)$City  
)


rownames(annot_col) <- sample_data(ps_phylum_top)$ID
annot_col$City <- gsub("Urban", "U", annot_col$City)
annot_col$City <- gsub("Rural", "R", annot_col$City)


png("heatmap_genus_healthy_city.png", width = 3000, height = 2000, res = 300)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annot_col,  
         fontsize_col = 8)

dev.off() 

###### HEATMAP WITH CLUSTER ON TOP 8 GENUS - CANCER ###### 
data <- subset_samples(complete_data, Tissue=="Cancer")
data.phy = tax_glom(data, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)

sd <- sample_data(data.phy.prop)
sd$ID <- gsub("Patient_", "", sd$ID)
sample_data(data.phy.prop) <- sd

otu_table(data.phy.prop)
sample_data(data.phy.prop)

ps_rel <- data.phy.prop
top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:8]
prune.dat_top <- prune_taxa(top,data.phy.prop)

ps_phylum_top <- prune.dat_top
otu_table(prune.dat_top)

mat <- as.matrix(otu_table(ps_phylum_top))

tax <- tax_table(ps_phylum_top)
rownames(mat) <- tax[rownames(mat), "Genus"]

wd_ids <- sample_data(ps_phylum_top)$ID
colnames(mat) <- wd_ids

annot_col <- data.frame(
  City = sample_data(ps_phylum_top)$City 
)

rownames(annot_col) <- sample_data(ps_phylum_top)$ID
annot_col$City <- gsub("Urban", "U", annot_col$City)
annot_col$City <- gsub("Rural", "R", annot_col$City)

png("heatmap_genus_Cancer_city.png", width = 3000, height = 2000, res = 300)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annot_col,  
         fontsize_col = 8)

dev.off() 

##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results/R_version_and_packages.txt") # create this file and connect it to R through "con" object 
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
print("egg")
packageVersion("egg")
cat("\n", "\n", fill=TRUE)
print("ggh4x")
packageVersion("ggh4x")
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink() # restore STR OUTPUT to R console
close(con)
suppressWarnings(rm(con))