##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  # analysis packages
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
}
options(scipen = 100) # disable scientific annotation


### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","indianred","springgreen1","wheat","black","coral1","yellow","darkmagenta","cyan","gray","firebrick3", "royalblue3","gold3","darkgreen","violet","chocolate3", "deepskyblue1","red","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza")

# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

# NB: the order of the row in the Metadata is employed also to set the order of the samples in the plots
metadata <- as.data.frame(read.table("metadata_PNFlorence.txt", header = T, sep="\t"))
row.names(metadata)<-metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])
original_order<-metadata$Sample_name # to maintain this sample order in plots

metadata<-metadata[sample, ]
identical(as.numeric(length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])),as.numeric(original_length))

# metadata$Sample_name<- factor(metadata$Sample_name, levels = unique(metadata$Sample_name) )
metadata$Exp_day_factor<- factor(metadata$Exp_day, levels = unique(metadata$Exp_day[order(as.numeric(metadata$Exp_day))]) )
metadata<- metadata[ order(as.numeric(metadata$Exp_day)) , ]
metadata$Sample_name <- factor(metadata$Sample_name, levels = unique(metadata$Sample_name) )
metadata$Date <- factor(metadata$Date, levels = unique(metadata$Date) )


sample_data(data)<-metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

# restoring the original row order of the metadata (it sets also the order in the plots!)
row.names(metadata)<-metadata$Sample_name
metadata<-metadata[original_order, ]

rm(original_length)




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
# write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.001% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
target_filter<- filter_taxa(data.genus.temp, function(x) max(x) <= 0.001, TRUE)
filtered<-taxa_names( target_filter )
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0001_max_cutoff.csv")


# if the ASVs are filtered at ASV level then only the representative of that genera will be filtered ...
# otherwise, if the genera name are filtered then also every eventual NA and uncultured will be filtered also in other families
# then, a more specific taxonomy table (with a partial path included in the unidentified genus) will be now build

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

unf_data2<-unfiltered_data
taxa_temp<-as.data.frame(tax_table(unf_data2))
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
  tax_table(unf_data2)<-as.matrix(taxa_temp)
}

unf_data2<-subset_taxa(unf_data2, ! Genus %in% genera_to_filter)
ASVs_code_to_maintain<- taxa_names(unf_data2)  # NB: taxa_names here means ASVs codes
data<-prune_taxa(ASVs_code_to_maintain, data)  # working on ASVs searched with specific names between different levels of taxonomy (raw vs genus) 
data<-prune_taxa(taxa_sums(data) > 0, data) 
rm(filtered, data.genus.temp, unf_data2, genera_to_filter, target_filter)



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
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria\n\n")
  data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mitochondria"))
}



proof1<-"This is the proof of having performed the filtering steps"




############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################
 
# if(! "proof1" %in% ls()){
#   cat("\n Wait! Did you perform the filtering step??? \n\n")
#   Sys.sleep(2)
# }
# 
# Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
# Original_number$sample.ID<-gsub("^.*F","",Original_number$sample.ID)
# row.names(Original_number)<-Original_number$sample.ID
# Original_number<- Original_number[as.character(sample_data(data)$FASTQ_ID), ]   # the files come from the whole sequencing batch
# Original_number<-sum( c(Original_number$forward.sequence.count, Original_number$reverse.sequence.count))
# Remaining_number<-sum(otu_table(data))
# 
# 
# con<-file("Data_check/DETAILS_ABOUT_ORIGINAL_AND_PROCESSED_READS_NUMBER.txt")
# sink(con, append=TRUE)
# cat("Raw reads abundance (sequenced, considering the pairs as one)", fill=TRUE)
# cat(Original_number/2)
# cat("Reads abundance after every filter", fill=TRUE)
# cat(Remaining_number)
# cat("Percentual of remaining vs original", fill=TRUE)
# cat(paste0(round((Remaining_number/(Original_number/2))*100,digits=2)),"%")
# sink()
# close(con)
# rm(con)
# 
# 
# ### barplot of reads
# Original_read_number<-as.data.frame(read.table(file="QIIME/Original_number_of_reads_for_sample.tsv", header = T, sep="\t", row.names = 1))
# row.names(Original_read_number)<-gsub("^.*F","",row.names(Original_read_number))
# Original_read_number <- Original_read_number[as.character(sample_data(data)$FASTQ_ID), ]   # the file come from the whole sequencing batch
# DADA2_read_number<-as.data.frame(read.table(file="QIIME/DADA2_Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
# row.names(DADA2_read_number)<-gsub("^.*F","",row.names(DADA2_read_number))
# DADA2_read_number <- DADA2_read_number[as.character(sample_data(data)$FASTQ_ID), ]   # the file come from the whole sequencing batch
# after_filter_number<-as.data.frame(colSums(otu_table(data)))
# # updating names from FASTQ codes to sample names
# row.names(Original_read_number)<-sample_data(data)$Sample_name # same order, already modified
# row.names(DADA2_read_number)<-sample_data(data)$Sample_name # same order, already modified
# row.names(after_filter_number)<-sample_data(data)$Sample_name # same order, already modified
# # creating the table
# table<-cbind.data.frame(Original=Original_read_number$forward.sequence.count, # merged --> just one column, not both of them
#                         After_quality_filter=DADA2_read_number$non.chimeric,
#                         After_contam_filter=after_filter_number[,1])
# table$FASTQ<-as.character(sample_data(data)$FASTQ_ID)
# # re-naming and adding groups
# if(identical( table$FASTQ , as.character(sample_data(data)$FASTQ_ID) ) ){ # same order
#   table$Samples<-as.factor(sample_data(data)$Sample_name)
#   table$Exp_day_factor<-sample_data(data)$Exp_day_factor
#   # table$factor <- sample_data(data)$Type
#   # table$factor <- as.factor(table$factor)
#   cat("\n\nOK!\n\n")
# }
# # plotting
# ggplot(aes(x=Exp_day_factor, y=Original), data=table) +
#   geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
#   geom_bar( aes(y= After_quality_filter, fill="Read quality filters"), alpha= 0.8,
#             width = 0.35, stat = "identity") +
#   geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
#             width= 0.45, stat = "identity") +
#   theme_classic( base_size = 8.5) +
#   theme(axis.text.x = element_text(size=8, angle = 40, vjust=1, hjust=1),
#         axis.text.y = element_text(size=5.5),
#         axis.title.y = element_text(size=7),
#         strip.text = element_text(size=9.5 ),
#         panel.grid.major.y = element_line(linewidth =0.1, color="grey"),
#         legend.position = "bottom",
#         plot.margin = margin(2, 1, 2, 1),
#         legend.margin = margin(0, 30, 0, 0),
#         legend.text = element_text(size=9),
#         legend.title = element_text(size=9.5),
#         legend.key.height = unit(0.4,"cm")
#   ) +
#   scale_fill_manual(name='',
#                     breaks=c('Original number', 'Read quality filters', 'Relative abundance filters'),
#                     values=c('Original number'='green3', 'Read quality filters'='coral', 'Relative abundance filters'='red3')) +
#   # scale_y_continuous(breaks = c(10000, seq(0, max(table$Original+10000), 10000))) +
#   labs(y="Reads abundance", x="Experiment day")
# ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 4.5, height = 3.8, dpi=300)
# # saving also the table itself
# write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)
# 
# 
# suppressWarnings( rm(table, table2, Original_read_number, DADA2_read_number, after_filter_number) )
 



########################### % ASSIGNED IN SILVA ##################################

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}
{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

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
    if(length(x) < 2) {
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


png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data.genus),"matrix")), step=100,label=F, ylab = "Genera", xlab= "Reads amount")
evalslopes(r,paste0("Day",sample_data(data.genus)$Exp_day),lim=0.001,cex=1)
dev.off()
rm(r)


proof2<-"No unsaturated sample to remove"




#################### \\\\\\ STARTING THE ANALYSIS \\\\\ #######################


if(! "proof1" %in% ls() ){
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

rm(taxa_temp)


# save.image("data_prepared_for_analysis.RData")




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

write.csv2(cbind(as(otu_table(unfiltered_data),"matrix"),as(tax_table(unfiltered_data),"matrix")),file="Data_check/ASVs_TAXA_counts_before_filters.csv",quote=F)




###################### ABUNDANCES BAR PLOT ##########################

### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.phy.prop)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-gsub ("Candidatus_","Ca. ", tabella$Phylum)
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}

ggplot(data=tabella, aes(x=Exp_day_factor, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_10) +
  theme_classic(base_size =10.5) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 9.5
                                ),
  axis.text.y=element_text(size=8),
  axis.title.x = element_text(vjust=5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.44, "cm"),
  legend.text = element_text ( size = 10.5 ),
  legend.position="bottom",
  legend.margin = margin(-5 ,15, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="\nExperiment day", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results/Abundances/TOP_phyla_EVERY_SAMPLE.png",width=4.5,height=4, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_phyla__EVERY_SAMPLE_averages.xlsx", row.names = F, to_save)

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)



### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tax_table(prune.data.others)[,"Genus"]<- "Others"
  prune.data.others <- tax_glom(prune.data.others, taxrank = "Genus")  # to avoid graphical glitches
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  #tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Rhizobi.","Rhizobiaceae", tabella$Genus, fixed=T)
  tabella$Genus<-gsub ("SM1A02","Planctomyc. SM1A02", tabella$Genus)
  tabella$Genus<-gsub ("A4b","Chloroflexi A4b", tabella$Genus)
  tabella$Genus<-gsub ("AKYH767","Sphingob. AKYH767", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
fill_color_modified<-c("brown4","orange","aquamarine4","darkmagenta","violet","red3","navajowhite2","black","cyan","slateblue1","gray85","gold3", "deeppink","darkgreen","blue2","springgreen1", "orangered","chartreuse3", "royalblue1","darkslategray3")
ggplot(data=tabella, aes(x=Exp_day_factor, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_modified) +
  # scale_fill_manual(values=fill_color_19) +
  theme_classic(base_size =10.5) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 size= 9
  ),
  axis.text.y=element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.175, "cm"),
  legend.text = element_text ( size = 7.7 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-13 ,40, 1 ,5),
  plot.margin = margin(3,1,2,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="\nExperiment day", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="")
ggsave(file="Results/Abundances/TOP_genera_EVERY_SAMPLE_Experiment_Day.png",width=4.5,height=4.5, dpi=300)
# again, but with dates
ggplot(data=tabella, aes(x=Date, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  scale_fill_manual(values=fill_color_modified) +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=25,
                                 vjust=1,
                                 hjust=1,
                                 size= 9.8
  ),
  axis.text.y=element_text(size=8),
  axis.title.x = element_text(vjust=5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.28, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.text = element_text ( size = 8.1 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-7 ,38, 1 ,5),
  plot.margin = margin(4,1,3,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_EVERY_SAMPLE_Dates.png",width=5.1,height=4.5, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame( "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Overall average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2),
                            "Average in Inoculum"= round( as.numeric(apply( otu_table(subset_samples(prune.dat_top, Sample_name=="PN1" )) ,1,mean)) ,2),
                            "Average on Last Day"= round( as.numeric(apply( otu_table(subset_samples(prune.dat_top, Sample_name=="PN9" )) ,1,mean)) ,2)
)
                           
to_save<-to_save[order(to_save$`Overall average`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_genera_EVERY_SAMPLE_averages.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




###################### ABUNDANCES OF ARCHAEA ##########################

### TOP Genera Archaea
suppressWarnings(rm(top, others, tabella, unass_data))
data.temp <- subset_taxa(data.genus.prop, Kingdom=="d__Archaea")
{top <- names(sort(taxa_sums(data.temp), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.temp)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
custom_colors<-c("chartreuse","black","yellow", 
                 "darkblue","red","darkslategray3")
ggplot(data=tabella, aes(x=Exp_day_factor, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  # facet_nested(~ Group, scales = "free_x", space="free") +
  theme_bw(base_size =10) +
  # scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 size= 9.5
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=10),
  axis.title.x = element_text(vjust=0.85),
  axis.title =element_text(size=10),
  # strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin(-1 ,35, 1 ,5),
  plot.margin = margin(2,1,1,1)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Experiment Day", y="Percentual abundance of clades (Archaea)",
       #title = "Most abundant identified phyla (every sample)",
       fill="")
ggsave(file="Results/Abundances/Most_Abundant_Archaea.png",width=4.5,height=4, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




################# TRENDS OF AOB AND NOB BACTERIA ####################

# AOBs derive from from MIDAS, also in other articles I can't find other names...
AOB_list <- c("Nitrosomonas", "Nitrosococcus", "Nitrosomonadales",   # synonymous (see MIDAS) or higher levels
              "Nitrosospira",  "Nitrosolobus", "Nitrosovibrio",   # synonymous, see MIDAS
              "966-1",     # https://www.sciencedirect.com/science/article/pii/S0045653523014649
              "Ellin6067" )   # https://www.sciencedirect.com/science/article/pii/S0960852423008970

# adding also AOA, from PMID: 24559743
AOA_list <- c("Nitrososphaera", "Nitrososphaerota", 
              "Nitrosocaldus", 
              "Nitrosotalea",
              "Nitrocosmicus")  
AOA_list <- c(AOA_list,
              "Nitrosarchaeum", # DOI: 10.1099/ijsem.0.002926
              "Nitrosacidococcus", # https://doi.org/10.1016/j.wroa.2022.100157
              "Nitrosocosmicus",  # https://www.sciencedirect.com/science/article/pii/S0038071720300225
              "Nitrosocaldus", "Nitrosothermus" , "Nitrosocaldaceae", # https://doi.org/10.3389/fmicb.2020.608832
              "Nitrosoabyssus", "Nitrosymbiomonas", # doi: 10.1186/s40168-025-02146-2
              "Nitrosopelagicus", # doi.org/10.1002/9781118960608.gbm01969
              "Nitrosopolaris", # https://doi.org/10.1093/femsmc/xtac019
              "Nitrosopumilus", # https://doi.org/10.1002/9781118960608.gbm01290
              "Nitrosopumilales", "Nitrosopumilaceae",  # https://doi.org/10.1002/9781118960608.obm00122
              "Nitrososphaeraceae"   # yes, also as "Genus", ref from https://www.sciencedirect.com/science/article/pii/S0038071720300225
)

# All of Nitrosomonadaceae are AOB ... REF: https://link.springer.com/rwe/10.1007/978-3-642-30197-1_372
Nitrosomonadaceae <- Taxa.genus[ Taxa.genus$Family=="Nitrosomonadaceae" , "Genus"]
Nitrosomonadaceae <- Nitrosomonadaceae[!is.na(Nitrosomonadaceae)]
Nitrosomonadaceae <- Nitrosomonadaceae[Nitrosomonadaceae!="uncultured"]

# NOBs
NOB_list <- c("Nitrospira","Nitrospirae", "Nitrospirae", "Nitrospirales","Nitrospiraceae", "Nitrospirota",  # alcuni da MIDAS, ma tutti elencati in DOI: 10.1016/S0076-6879(11)86005-2
              # few Nitrospira names are quite similar, these are synonims or, at least, they still are of Nitrospirales order which is still known for AO activity, see https://www.nature.com/articles/s41396-023-01518-6
              "Nitromaritima",
              "Nitrospina","Nitrospinae","Nitrospinaceae","Nitrospinota",
              "Nitrobacter",
              "Nitrococcus",
              "Nitronereus" , # DOI: 10.1038/s41396-023-01518-6
              "Nitrohelix", "Nitronauta", "Nitrocaldera", "Nitrotheca", # doi: 10.1021/acs.est.3c00636
              "Nitrolancea", #doi: 10.1099/ijs.0.062232-0.
              "Nitrotoga")


### plot 1
{ top_data <- data.genus.prop
  vector_without_Ca <- gsub("Candidatus_","", tax_table( top_data )[ ,"Genus"])  # No "Candidatus" but same order, to a easier use of the logical below
  prune.dat_top <- subset_taxa(top_data , vector_without_Ca %in% c(AOB_list,AOA_list,NOB_list,Nitrosomonadaceae) )
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
    tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  } else {
    stop("Something went wrong when updating the taxonomy")
  }
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  # prune.dat_top <- tax_glom(prune.dat_top, taxrank = "Genus", NArm=T ) 
  tabella_top<-psmelt(prune.dat_top)
  tabella<-tabella_top[order(sort(tabella_top$Abundance, decreasing = T)), ]
  tabella$Genus2 <- gsub("Candidatus_","", tabella$Genus, fixed = T)
  tabella$Genus[tabella$Genus2 %in% c(AOB_list,Nitrosomonadaceae)] <- paste(tabella$Genus[tabella$Genus2 %in% c(AOB_list,Nitrosomonadaceae)], "(AOB)" )
  tabella$Genus[tabella$Genus2 %in% c(AOA_list)] <- paste(tabella$Genus[tabella$Genus2 %in% c(AOA_list)], "(AOA)" )
  tabella$Genus[tabella$Genus2 %in% c(NOB_list)] <- paste(tabella$Genus[tabella$Genus2 %in% c(NOB_list)], "(NOB)  " )
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
fill_color_AOB_NOB<-c("blue1","slateblue","lightblue3","firebrick","orangered1","royalblue1")
tabella$Exp_day <- factor (as.character(tabella$Exp_day) , levels = unique(metadata$Exp_day) )

ggplot(data=tabella, aes(x=Exp_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid( ~ Year, space= "free_x", scales = "free_x") +
  scale_fill_manual(values=fill_color_AOB_NOB) +
  theme_classic(base_size =10.5) +
  scale_y_continuous(expand=c(0,0.1) ) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 9.5
  ),
  axis.text.y=element_text(size=7),
  strip.text.x = element_text(size=9.5, angle=0),
  panel.grid.major.y = element_line(colour="darkgray", linewidth=0.15),
  panel.grid.minor.y = element_line(colour="gray", linewidth=0.09),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin(0 ,20, 1 ,5),
  plot.margin = margin(2,1,2,1)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment Day", y="Percentual abundance of clades (16S)",
       fill="")
ggsave(file="Results/Abundances/AOB_NOB_AlongTime.png",width=4.5,height=4, dpi=300)

dev.off()



#### Again, but proportions among themselves...
{ top_data <- data.genus.prop
  vector_without_Ca <- gsub("Candidatus_","", tax_table( top_data )[ ,"Genus"])  # No "Candidatus" but same order, to a easier use of the logical below
  prune.dat_top <- subset_taxa(top_data , vector_without_Ca %in% c(AOB_list,AOA_list,NOB_list,Nitrosomonadaceae) )
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
    tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  } else {
    stop("Something went wrong when updating the taxonomy")
  }
  prune.dat_top <- transform_sample_counts( prune.dat_top, function(x) x/sum(x)*100 )
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella_top<-psmelt(prune.dat_top)
  tabella<-tabella_top[order(sort(tabella_top$Abundance, decreasing = T)), ]
  tabella$Genus2 <- gsub("Candidatus_","", tabella$Genus, fixed = T)
  tabella$Genus[! tabella$Genus2 %in% c("Nitrospira","Nitrosomonas","Nitrotoga")] <- "Other nitrifiers"
  tabella$Genus[tabella$Genus2 %in% c("Nitrospira","Nitrotoga")] <- paste(tabella$Genus[tabella$Genus2 %in% c("Nitrospira","Nitrotoga")], "(NOB)" )
  tabella$Genus[tabella$Genus2 %in% "Nitrosomonas"] <- paste(tabella$Genus[tabella$Genus2 %in% "Nitrosomonas"], "(AOB)" )
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c("Ca. Nitrotoga (NOB)","Nitrospira (NOB)", "Nitrosomonas (AOB)", "Other nitrifiers") )
}
fill_color_AOB_NOB<-c("orangered1","firebrick4","royalblue1","navyblue")
tabella$Exp_day <- factor (as.character(tabella$Exp_day) , levels = unique(metadata$Exp_day) )
ggplot(data=tabella, aes(x=Exp_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid( ~ Year, space= "free_x", scales = "free_x") +
  scale_fill_manual(values=fill_color_AOB_NOB) +
  theme_classic(base_size =10.5) +
  scale_y_continuous(expand=c(0,0.1) ) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 9.5
  ),
  axis.text.y=element_blank(),
  strip.text.x = element_text(size=9.5, angle=0),
  panel.grid.major.y = element_line(colour="darkgray", linewidth=0.15),
  panel.grid.minor.y = element_line(colour="gray", linewidth=0.09),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 10.5 ),
  legend.position="bottom",
  legend.margin = margin(0 ,20, 1 ,5),
  plot.margin = margin(3,1,2,1)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment Day", y="Proportion of nitrifiers",
       fill="")
ggsave(file="Results/Abundances/AOB_NOB_ProportionsAmongThemselves.png",width=4.5,height=4, dpi=300)




######################## ALPHA DIVERSITY ALONG TIME ##############################

data.genus.temp<-data.genus
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Sample_name
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<- mix_alpha$Shannon /log(mix_alpha$Observed)
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$Exp_day <- metadata[mix_alpha$Sample_name, "Exp_day" ]  # NB: the metadata has to be ordered according to time!
mix_alpha$Date <- metadata[mix_alpha$Sample_name, "Date" ]
mix_alpha$Date <- factor(mix_alpha$Date , levels=unique(metadata[mix_alpha$Sample_name, "Date" ]) ) # NB: the metadata has to be ordered according to time!
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
days <- unique(metadata$Exp_day)[!unique(metadata$Exp_day)%in%c(0,437)]
mix_alpha$Time_axisX<- factor(as.character(mix_alpha$Exp_day), levels =c (0, 
                                                                          days[order(days)],
                                                                          437) )  # continuous value to print a shade of color
mix_alpha <- mix_alpha[ order(mix_alpha$Exp_day) ,]
# plot
ggplot(data=mix_alpha, aes(y=value, x=Time_axisX, color=Exp_day)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(0,50,100,150,184),
                        midpoint = 50) +
  geom_point(size =1.8) +
  geom_point(size =3, alpha=0.4) +
  geom_path(aes(group=""), color="darkgray" , linewidth = 0.18, alpha=0.9,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic2(base_size = 11.2) +
  labs(y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  theme(axis.text.x = element_text(size=9.5, angle=35, hjust=1, vjust=1),
        axis.text.y= element_text(angle=0, size=7),
        #axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  ) +
  labs(x="Experiment day")
ggsave(file="Results/Alfa_diversity_at_GENUS_level_with_Days.png", width = 4.65,height =4, dpi=300)
# again, but with dates
# ggplot(data=mix_alpha, aes(y=value, x=Date, color=Exp_day)) +
#   facet_grid( variable ~ . , scales = "free", space = "free_x") +
#   scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
#                         breaks=c(0,50,100,150,184),
#                         midpoint = 50) +
#   geom_point(size =1.8) +
#   geom_point(size =3, alpha=0.4) +
#   geom_path(aes(group=""), color="darkgray" , linewidth = 0.18, alpha=0.9,
#             arrow=arrow(length =unit(0.22,"cm"), type = "closed")
#   ) +
#   theme_classic2(base_size = 11.2) +
#   labs(y="Alpha Diversity Measure") +
#   guides(fill="none", color="none", shape="none") +
#   scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
#   theme(axis.text.x = element_text(size=9.2, angle=35, hjust=1, vjust=1),
#         axis.text.y= element_text(angle=0, size=7),
#         axis.title.x = element_blank(),
#         panel.grid.major.y = element_line(linewidth=0.4),
#         panel.grid.minor.y = element_line(linewidth=0.25),
#         plot.margin = margin(5,1,5,1)
#   )
# ggsave(file="Results/Alfa_diversity_at_GENUS_level_with_DATES.png", width = 5,height =4.2, dpi=300)




########################### PCoA BETA DIV ##########################

# on genera (every sample)
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
new_label <- paste0(sample_data(data.sqrt_prop)$Exp_day,rep("th") )
new_label[new_label=="0th"] <- ""
color_time <- c("PN1"="cyan",
                "PN2"="deepskyblue", "PN3"="deepskyblue3",
                "PN4"="royalblue","PN5"="royalblue3",
                "PN6"="royalblue3","PN7"="blue1",
                "PN8"="blue2","PN9"="blue3"
                )
plot_ordination(data.sqrt_prop, ordBC, color="Sample_name") + 
  scale_color_manual(values= color_time  )+
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(shape="none", color="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=  new_label ) ,
            color="grey30", 
            vjust=-0.65,
            size=2.8, show.legend = FALSE) +
  labs(
    # title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="")
ggsave(file="Results/Beta_div_Hellinger_on_EVERYSAMPLE_genera_WITH_Day.png", width = 3.8, height = 3.5, dpi=300)
# again, with line
meta_temp <- sample_data(data.sqrt_prop)
meta_temp$temp <- as.factor("rep")
sample_data(data.sqrt_prop) <- meta_temp[ order(meta_temp$Exp_day), ]
plot_data <- plot_ordination(data.sqrt_prop, ordBC, color="Sample_name") + 
  scale_color_manual(values= color_time  )+
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(shape="none", color="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=  new_label ) ,
            color="grey30", 
            vjust=-0.65,
            size=2.8, show.legend = FALSE) +
  labs( 
    #title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="")
plot_data$data <- plot_data$data[ order(plot_data$data$Exp_day) , ]   # to re-order the data (sample data in phyloseq not sortable) and then the geom path order
plot_data +
  geom_path(aes(group=temp),  col= "darkgray", linewidth = 0.24,
                      arrow=arrow(length =unit(0.25,"cm"), type = "closed")
            )
ggsave(file="Results/Beta_div_Hellinger_on_EVERYSAMPLE_genera_path.png", width = 3.8, height = 3.5, dpi=300)



# # on genera (first samples)
# data.prop.labels<-subset_samples(data.genus.prop, Exp_day < 100 )
# {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# new_label <- paste0(sample_data(data.sqrt_prop)$Exp_day,rep("th") )
# new_label[new_label=="0th"] <- ""
# plot_ordination(data.sqrt_prop, ordBC, color="Sample_name") + 
#   scale_color_manual(values= c("CA1"="cyan",
#                                "CA2"="deepskyblue",
#                                "CA3"="deepskyblue3",
#                                "CA4"="blue",
#                                "CA5"="blue4"
#   )
#   )+
#   geom_point(size=3.5, alpha=0.5,  show.legend = F) +
#   guides(shape="none", color="none") +
#   theme_classic(base_size = 10) +
#   theme(legend.text =element_text(size=9)) +
#   geom_text(aes(label=  new_label ) ,
#             color="grey30", 
#             vjust=-0.55,
#             size=2.8, show.legend = FALSE) +
#   labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
#        x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        color="")
# ggsave(file="Results/Beta_div_Hellinger_on_FIRSTfiveSAMPLES_genera_WITH_Day.png", width = 4, height = 3.5, dpi=300)
# # again, with line
# meta_temp <- sample_data(data.sqrt_prop)
# meta_temp$temp <- as.factor("rep")
# sample_data(data.sqrt_prop) <- meta_temp[ order(meta_temp$Exp_day), ]
# plot_data <- plot_ordination(data.sqrt_prop, ordBC, color="Sample_name") + 
#   scale_color_manual(values= c("CA1"="deepskyblue",
#                                "CA2"="deepskyblue2",
#                                "CA3"="deepskyblue3",
#                                "CA4"="deepskyblue4",
#                                "CA5"="blue",
#                                "CA6"="blue2",
#                                "CA7"="blue3",
#                                "CA8"="blue4"
#   )
#   )+
#   geom_point(size=3.5, alpha=0.5,  show.legend = F) +
#   guides(shape="none", color="none") +
#   theme_classic(base_size = 10) +
#   theme(legend.text =element_text(size=9)) +
#   labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
#        x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        color="")
# plot_data$data <- plot_data$data[ order(plot_data$data$Exp_day) , ]   # to re-order the data (sample data in phyloseq not sortable) and then the geom path order
# plot_data +
#   geom_path(aes(group=temp),  col= "darkgray", linewidth = 0.3,
#             arrow=arrow(length =unit(0.25,"cm"), type = "closed")
#   )
# ggsave(file="Results/Beta_div_Hellinger_on_FIRSTfiveSAMPLES_genera_path.png", width = 4, height = 3.5, dpi=300)




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

