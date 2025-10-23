# Bar plot e PCoA solo AGS diretti/conservati nel tempo



##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("DESeq2")
  library("ALDEx2")
  library("vegan")
  library("ecodist")
  # utilities
  library("reshape")
  library("xlsx")  
  library("qiime2R")
}

{dir.create("Data_check")
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Abundances/Abund_over_time")
  dir.create("Results/Beta_diversity_PCoA")
  dir.create("Results/Differently_abundant_bact")
}
options(scipen = 100) # disable scientific annotation

### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group



# !!! NB: the order of the row in the Metadata is employed also to set the order of certain elements in the plots (e.g. arrows directions)!
metadata <- as.data.frame(read.table("Metadata_aero_conserv_sequencing.txt", header = T, sep="\t"))



####################### IMPORTING EACH DATASET #####################

### IMPORTING RAW 16S seq

# devtools::install_github("jbisanz/qiime2R")
data1<-qza_to_phyloseq(features="QIIME/table_batch1.qza", taxonomy="QIIME/taxonomy_SILVA_batch1.qza")
original_names1<-sample_names(data1)
sample_names(data1)<-gsub(".*F","",sample_names(data1))

data2<-qza_to_phyloseq(features="QIIME/table_batch2.qza", taxonomy="QIIME/taxonomy_SILVA_batch2.qza")
original_names2<-sample_names(data2)
sample_names(data2)<-gsub(".*F","",sample_names(data2))

data3<-qza_to_phyloseq(features="QIIME/table_batch3.qza", taxonomy="QIIME/taxonomy_SILVA_batch3.qza")
original_names3<-sample_names(data3)
sample_names(data3)<-gsub(".*F","",sample_names(data3))

# using the sample data to subset the third batch (it includes both AERO and PN2 samples)
these_are_flocs <- metadata$FASTQ_ID [metadata$Seq_batch%in%c("5",5) & metadata$Biomass=="Flocs"]
data_AERO3 <- prune_samples( sample_names(data3)%in% these_are_flocs, data3)

# removing the replicates from data2 (not ecologically saturated)
# replicated <- metadata$FASTQ_ID [grepl("_unsat", metadata$Sample_name)]
# data2 <- prune_samples( ! sample_names(data2)%in% replicated, data2)



### GENERATING THE MAIN RAW OBJS

data_from_AS<-merge_phyloseq(data1,data2,data_AERO3)
data_from_PN2 <- prune_samples( ! sample_names(data3)%in% these_are_flocs, data3)


#importing also AGS IDRO and PN    (from other works)
load("Data_from_IDRO.RData")
load("Data_from_PN.RData")
sample_names(data_from_PN) <- sample_data(data_from_PN)$FASTQ_ID
data_from_IDRO<- samples_from_IDRO # renaming
rm(samples_from_IDRO)

sample_names(data_from_IDRO) <- sample_data(data_from_IDRO)$FASTQ_ID
data_from_PN <- phyloseq(otu_table(data_from_PN),tax_table(data_from_PN))  # removing the sample data (otherwise the "data" samples will be removed, as they still lack the sample data!)
data_from_IDRO <- phyloseq(otu_table(data_from_IDRO),tax_table(data_from_IDRO))




########### FILTERING NOISES FROM DATA SETS (EACH REACTOR ON ITS OWN) ##############

list_data <- list("data_from_IDRO"=data_from_IDRO,   # NB: the name in the list MUST be the same name of the obj to overwrite the latter
                  "data_from_PN"=data_from_PN,
                  "data_from_PN2"=data_from_PN2,
                  "data_from_AS"=data_from_AS
                  )
for(y in 1:length(list_data)){
  cat("Filtering ",names(list_data)[y],"...")
  data_y <- list_data[[y]]
  data.genus.temp<-tax_glom(data_y, taxrank = "Genus", NArm = F)
  data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
  target_filter<- filter_taxa(data.genus.temp, function(x) mean(x) < 0.005, TRUE)  # to remove noises/contaminants, see PMID:23202435 and DOI: 10.1128/mSystems.00290-19

  filtered<- taxa_names( target_filter )    # collecting the genera names
  # updating taxonomy
  taxa_temp<-as.data.frame(tax_table(target_filter))
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
    genera_to_filter<-taxa_temp[["Genus"]]  # every NA or uncultured will now have a more detailed path in the genus name
  }
  data_y<-subset_taxa(data_y, ! Genus %in% genera_to_filter)
  ASVs_code_to_maintain<- taxa_names(data_y)  # NB: taxa_names here means ASVs codes
  data_y<-prune_taxa(ASVs_code_to_maintain, data_y)  # working on ASVs searched with specific names between different levels of taxonomy (raw vs genus) 
  
               # name in env                   # modified obj
  assign(   names(list_data)[y]  , prune_taxa(taxa_sums(data_y)>0,data_y)    )  # overwriting the obj
}

rm(x,y,list_data,target_filter,taxa_temp,data_y,ASVs_code_to_maintain)




############# \\\\\\ CREATING AN UNIQUE PHY OBJ WITH SAMPLE DATA \\\\\\ ##############

### Merging
original_length <- length(c(sample_names(data_from_IDRO),sample_names(data_from_PN),sample_names(data_from_AS), sample_names(data_from_PN2) ))
sample <- c(sample_names(data_from_AS), sample_names(data_from_IDRO), sample_names(data_from_PN), sample_names(data_from_PN2) )
data <- merge_phyloseq(data_from_PN, data_from_IDRO, data_from_AS, data_from_PN2 )

row.names(metadata)<-metadata$FASTQ_ID # column with FASTQ/SAMPLE name

# head(Metadata)
#metadata$FASTQ_ID [ ! metadata$FASTQ_ID %in% sample]

metadata<-metadata[sample, ]
identical(as.numeric(length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])),as.numeric(original_length))


# ordering the factors
metadata$Sample_name<- factor(metadata$Sample_name, levels = unique(metadata$Sample_name) )
metadata$Sample_code<- factor(metadata$Sample_code, levels = unique(metadata$Sample_code) )
metadata$Reactor_ID<- factor(metadata$Reactor_ID, levels = unique(metadata$Reactor_ID) )
metadata$Collection_data <- factor(metadata$Collection_data, levels = unique(metadata$Collection_data) )
metadata$Glycerol<- factor(metadata$Glycerol, levels = c(" ","noGly","Gly") )
metadata$Storage <- factor(metadata$Storage , levels= c("Direct","4C","-20C","-20C_noGly","-80C","-80C_noGly","No_storage" ) )
metadata$Time_before_extr <- factor(metadata$Time_before_extr , levels= c("None","1day","1week","1month","2month","3month","6month","9month","1year","3year") )
metadata$Time_before_extr2 <- factor(metadata$Time_before_extr2 , levels= c(" ","1d","1w","1m","2m","3m","6m","9m","1y","3y") )
metadata$Season <- factor(metadata$Season , levels= c("Autumn","Winter","Spring","Summer") )

sample_data(data)<-metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(sample)))  # TRUE --> no samples cut away

# restoring the original row order of the metadata (it sets also the order in the plots!)
row.names(metadata)<-metadata$Sample_name

rm(original_length)



### removing the sample A21 (is the only one without direct extraction, it was paired with another RNAseq analysis )
data <- subset_samples(data, Sample_name!="A21")
metadata <- metadata[metadata$Sample_name!="A21", ]




############## FILTERING ANOMALOUS CLASSIFICATIONS ###################

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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Bacteria_Archaea_proportion_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)

# removing every Unassigned Phylum (even the phylum was not in SILVA? Probable contaminant!)
data<- subset_taxa(data, ! is.na(Phylum) )


# removing chloroplast and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplasts\n\n")
  data<-subset_taxa(data, Genus != "Chloroplast")
}
if( "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Mitochondria\n\n")
  data<-subset_taxa(data, Genus != "Mitochondria")
}



proof1<-"This is the proof of having performed the filtering steps"




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
      differ<- length(v)-t
      if(differ<0){differ<-0} # to avoid the error with negative values
      slope<-mean(diff(v[ (differ):length(v) ])/dx)
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


if("proof2" %in% ls()){
  stop("\nLooks like you already performed this filter, you won't see the unsaturated samples!\n")
} else {
  data.genus.temp <- tax_glom(data, taxrank = "Genus", NArm = F)
  png(file="Data_check/Rarefaction_curve.png",width=2450,height=1800, res=300)
  r<-rarecurve(t(as(otu_table(data.genus.temp),"matrix")), step=100,label=F, ylab = "Genera", xlab= "Reads amount")
  evalslopes(r,
             paste0("       ",sample_data(data.genus.temp)$Sample_name),  # adding a space before the name, in order to emulate a hjust
             lim=0.001,cex=0.85)
  dev.off()
}
rm(r)

unsaturated <- c("A43") #  AS 08/23 at 6m
{
  proof2<-"Rarefaction curves have been evaluated"
  data<- subset_samples(data, !Sample_name%in%unsaturated)
  metadata <- metadata[!metadata$Sample_name%in%unsaturated, ]
}



#################### \\\\\\ PREPARING THE MAIN OBJECTS \\\\\ #######################

if( ! "proof2" %in% ls() ){
  stop("\nDid you perform the filtering steps yet?\n")
}


### everything ready, let's start with the analysis
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
#  data.class = tax_glom(data, taxrank = "Class", NArm = F)
#  data.order = tax_glom(data, taxrank = "Order", NArm = F)
#  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
#  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
#  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
#  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
#  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
#  Taxa.class<-as.data.frame(tax_table(data.class))
#  Taxa.order<-as.data.frame(tax_table(data.order))
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

suppressWarnings( rm(sample, taxa_temp, data.genus.temp, these_are_flocs, replicated, x, evalslopes))


save.image("data_prepared_for_analysis.RData")





########################### % ASSIGNED IN SILVA ##################################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
# b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
# c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
#  d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,e)
colnames(assigned)<-c("Total","Assigned","%","Taxa")
assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
suppressWarnings( rm(a,b,c,d,e,assigned) )




########################### COUNTS EXPORT ##########################################

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
    #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
}




##################### \\\\\\ R AND PACKAGES VERSION \\\\\\ #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results/R_version_and_packages.txt")
sink(con, append = TRUE)

cat(package$R.version$version.string)
cat("   running on", package$running)
cat("\n", "\n", fill=TRUE)
package$otherPkgs$phyloseq[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$ggplot2[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$vegan[c(1,3)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$DESeq2[1:2]
cat("\n", "\n", fill=TRUE)
print("ggh4x")
packageVersion("ggh4x")
cat("\n", "\n", fill=TRUE)
print("vegan")
packageVersion("vegan")
cat("\n", "\n", fill=TRUE)
print("microbiome")
packageVersion("microbiome")
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))

