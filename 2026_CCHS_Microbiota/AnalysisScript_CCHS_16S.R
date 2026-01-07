##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
  library("ggplot2")
  library("ggpubr")
  library("ggh4x") #
  library("ggvenn") #
  library("vegan")
  library("ecodist")
  library("DESeq2")
  library("mixOmics")
  library("qiime2R")
  library("Hmisc")
}

{
  dir.create("Data_check")
  dir.create("Data_check/PCoA_test")
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Abundances/Ratio_Firmi_Bacteroi")
  dir.create("Results/Beta_div")
  dir.create("Results/DA_DESeq2")
  dir.create("Results/PLS_DA")
}

options(scipen = 100) # disable scientific annotation



### CircoTax function from https://github.com/matteoramazzotti/CircoTax/blob/main/CircoTax.R
source("CircoTax.R")




####################### PREPARING THE 16S DATA #####################

### IMPORTING RAW DATA

# if from QIIME
# devtools::install_github("jbisanz/qiime2R")
data1<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza")
original_names1<-sample_names(data1)
sample_names(data1)<-gsub(".*F","",sample_names(data1))

data2<-qza_to_phyloseq(features="QIIME_batch2/table.qza", taxonomy="QIIME_batch2/taxonomy_SILVA.qza")
original_names2<-sample_names(data2)
sample_names(data2)<-gsub(".*F","",sample_names(data2))



### FILTERING NOISES FROM 16S' 1st and 2nd BATCH INDEPENDENTLY

if(! "proof1" %in% ls()){
  unfiltered_data1<-data1
  unfiltered_data2<-data2
}

# 1st batch
suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data1, taxrank = "Genus", NArm = F)
#write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)

# cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_cutoff_batch1.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured" & !is.na(filtered)] # to avoid the removal of other uncultured genera
data1<-subset_taxa(data1, ! Genus %in% filtered)
data1<-prune_taxa(taxa_sums(data1) > 0, data1) 
rm(filtered, data.genus.temp)

# 2nd batch
suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data2, taxrank = "Genus", NArm = F)
#write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)

# cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_cutoff_batch1.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured" & !is.na(filtered)] # to avoid the removal of other uncultured genera
data2<-subset_taxa(data2, ! Genus %in% filtered)
data2<-prune_taxa(taxa_sums(data2) > 0, data2) 
rm(filtered, data.genus.temp)


proof1<- "Marker of the filtering, it is required for the script"



### JOINING

data <- merge_phyloseq(data1, data2)



### REMOVING CHLOROPLASTS AND MITOCHONDRIA
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria\n\n")
  data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mitochondria"))
}



### CHECKING UNASSIGNED IN PROKARYOTE DOMAIN/KINGDOM

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
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)




##################### IMPORTING THE METADATA #########################

Metadata <- as.data.frame(read.delim(file="Metadata_CCHS.txt", sep="\t", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample_names(data), ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names1)+length(original_names2)))

rm(original_length)

sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels=c("Control","Control_parent","CCHS")) # decides the order in plots
sample_data(data)$Condition2<-factor(sample_data(data)$Condition2, levels=c("Control","CCHS")) 
sample_data(data)$Age_groups<-factor(sample_data(data)$Age_groups, levels=c("<25",">25"))
sample_data(data)$Age_3groups<-factor(sample_data(data)$Age_3groups, levels=c("<20","20-35",">35"))
# head(sample_data(data))


unsubsetted<-data
data <- subset_samples(data, ! Sample_name %in% c("C1_bis","C17G","H2_bis","H5","HP4_bis","HP9G" , "C5","C10","C11" ))
Metadata <- Metadata[ ! Metadata$Sample_name %in% c("C1_bis","C17G","H2_bis","H5","HP4_bis","HP9G" , "C5","C10","C11" ) , ]




#################### CHECKING SAMPLES BALANCE BETWEEN GROUPS ####################

Meta_no_replicates <- sample_data(data)
ages_C <- rbind.data.frame( round(mean( Meta_no_replicates$Age[Meta_no_replicates$Condition2=="CCHS"] ),0) ,
                            min (Meta_no_replicates$Age[Meta_no_replicates$Condition2=="CCHS"] ),
                            max (Meta_no_replicates$Age[Meta_no_replicates$Condition2=="CCHS"] )
                            )
ages_H <- rbind.data.frame( round(mean( Meta_no_replicates$Age[Meta_no_replicates$Condition2=="Control"] ),0) ,
                            min (Meta_no_replicates$Age[Meta_no_replicates$Condition2=="Control"] ),
                            max (Meta_no_replicates$Age[Meta_no_replicates$Condition2=="Control"] )
                            )

counting <- as.matrix(table( Meta_no_replicates$Age_3groups, Meta_no_replicates$Condition2 ))
counting <- rbind( counting , as.matrix(cbind(ages_C,ages_H)) )
row.names(counting)<- c("<20","20-35",">35","mean","min","max" )

write.csv2(counting , 
           file="Data_check/Contingency_table_COND_AGE.csv",
           quote=F)




############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

# if(! "proof1" %in% ls()){
#   cat("\n Wait! Did you perform the filtering step??? \n\n")
#   Sys.sleep(2)
# }
# 
# Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
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


# ### barplot of reads
# Original_read_number<-as.data.frame(read.table(file="QIIME/Original_number_of_reads_for_sample.tsv", header = T, sep="\t", row.names = 1))
# DADA2_read_number<-as.data.frame(read.table(file="QIIME/DADA2_Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
# after_filter_number<-as.data.frame(colSums(otu_table(data)))
# # updating names from FASTQ codes to sample names
# Original_read_number<-Original_read_number[original_names, ] # this vector came from the import section (ordered as "sample", which contains the modified names)
# DADA2_read_number<-DADA2_read_number[original_names, ]
# row.names(Original_read_number)<-sample # same order, already modified
# row.names(DADA2_read_number)<-sample # same order, already modified
# # same order of rows
# Original_read_number<-Original_read_number[row.names(DADA2_read_number), ]
# this_order<-codes[row.names(DADA2_read_number), "Sample"]
# after_filter_number<-after_filter_number[this_order, ]
# # creating the table
# table<-cbind.data.frame(Original=Original_read_number$forward.sequence.count, # merged --> just one column, not both of them
#                         After_quality_filter=DADA2_read_number$non.chimeric,
#                         After_contam_filter=after_filter_number)
# # adding factors to re-group the table
# table$factor<-as.data.frame(sample_data(data))[this_order, ]$Condition2
# table$Samples<-row.names(Original_read_number)
# # plotting
# ggplot(aes(x=Samples, y=Original), data=table) +
#   facet_grid( ~ factor, space="free_x", scales = "free_x") +
#   geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
#   geom_bar( aes(y= After_quality_filter, fill="Quality read filters"), alpha= 0.8,
#             width = 0.6, stat = "identity") +
#   geom_bar( aes(y= After_contam_filter, fill="Contaminant filters"),
#             width= 0.4, stat = "identity") +
#   theme_classic( base_size = 10.2) +
#   theme(axis.text.x = element_text(size=6, angle = 90, vjust=0.5),
#         axis.text.y = element_text(size=5),
#         panel.grid.major.y = element_line(size=0.1, color="grey"),
#         legend.position = "bottom",
#         legend.margin = margin(0, 0, 0, 0),
#         legend.text = element_text(size=9.2),
#         legend.title = element_text(size=9.8),
#         legend.key.height = unit(0.4,"cm")
#   ) +
#   scale_fill_manual(name='Number of reads:  ',
#                     breaks=c('Original number', 'Host read filters', 'Quality read filters', 'Contaminant filters'),
#                     values=c('Original number'='green3', 'Host read filters'='blue','Quality read filters'='coral', 'Contaminant filters'='red3')) +
#   scale_y_continuous(breaks = c(10000, seq(0, max(table$Original)+10000, 15000))) +
#   labs(y="Reads abundance", x="FASTQ")
# ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 7.2, height = 3.8, dpi=300)
# # saving also the table itself
# write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)
# 
# 
# rm(table, Original_read_number, DADA2_read_number, after_filter_number)
 



################# CHECKING THE COMPOSITION AFTER FILTERING ####################
# 
# data.unf.prop<- merge_phyloseq(unfiltered_data1, unfiltered_data2)
# data.unf.prop<-transform_sample_counts(data.unf.prop, function(x) (x/sum(x))*100) # if BOWTIE STEP IS NOT PERFORMED
# 
# 
# if(! "proof1" %in% ls()){
#   cat("\n Wait! Did you perform the filtering step??? \n\n")
#   Sys.sleep(2)
# }
# data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100)
# 
# 
# ################# BRAY HELLINGER
# 
# ##### unfiltered BRAY
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
# data.prop.labels<-data.unf.prop
# {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
#   scale_color_manual(values=c("CCHS"="coral","Control"="deepskyblue2","Control_parent"="deepskyblue")) +
#   guides(color="none") +
#   geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(linewidth=0.2) + 
#   geom_text(aes(label= Sample_name), 
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data", 
#        color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# ##### filtered BRAY
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
# data.prop.labels<-data.prop
# {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
#   scale_color_manual(values=c("CCHS"="coral","Control"="deepskyblue2","Control_parent"="deepskyblue")) +
#   geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(linewidth=0.2) + 
#   geom_text(aes(label= Sample_name), 
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="\n \n filtered data", 
#        color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3200, height = 1800, res=300)
# ggarrange(p1,p2, nrow = 1)
# dev.off()
# 
# suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))
# 
# 
# 
# 
# ################# EUCLIDEAN HELLINGER
# 
# ##### unfiltered EUCLIDEAN
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
# data.prop.labels<-data.unf.prop
# {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
#   scale_color_manual(values=c("CCHS"="coral","Control"="deepskyblue2","Control_parent"="deepskyblue")) +
#   guides(color="none") +
#   geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(linewidth=0.2) + 
#   geom_text(aes(label= Sample_name), 
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
#        color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# ##### filtered EUCLIDEAN
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
# data.prop.labels<-data.prop
# {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
#   scale_color_manual(values=c("CCHS"="coral","Control"="deepskyblue2","Control_parent"="deepskyblue")) +
#   geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(linewidth=0.2) + 
#   geom_text(aes(label= Sample_name), 
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="\n \n filtered data", 
#        color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3200, height = 1800, res=300)
# ggarrange(p1,p2, nrow = 1)
# dev.off()
# 
# suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))




################# CHECKING THE VARIABILITY AMONG REPLICATES #######################

data.genus_unsub = tax_glom(unsubsetted, taxrank = "Genus", NArm = F)
data.genus_unsub.prop <- transform_sample_counts(data.genus_unsub, function(ASV) ASV/sum(ASV)*100)

data.prop.labels<-data.genus_unsub.prop
sample_data(data.prop.labels)$ID<-gsub("_.*","",sample_data(data.prop.labels)$Sample_Family)       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition2") +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral")) + # ,"Control_parent"="deepskyblue")) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10) +
  geom_line(aes(group=ID),linewidth=0.35, color="grey25") +
  stat_ellipse(linewidth=0.1, aes(color=Condition2, shape=NULL))  +
  theme( legend.text = element_text(size=10) ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_Family), color="black", size=1.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (with replicates)",
       color="Condition", # shape="Hirschsprung", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Data_check/PCoA_with_replicates.png", width = 5, height = 4.2, dpi=300)

# PCoA BETWEEN BATCHES
plot_ordination(data.sqrt_prop, ordBC, color = "Batch", shape="Condition2") +
  scale_color_manual(values=c("1o"="green","2o"="red")) + # ,"Control_parent"="deepskyblue")) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 8.5) +
  geom_line(aes(group=ID),linewidth=0.35, color="grey25") +
  theme( legend.text = element_text(size=8) ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_Family), color="black",
            size=1.95, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (with lines connect replicates)",
       color="Batch" ,
       shape="Condition", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Data_check/PCoA_with_between_batches.png", width = 5, height = 4.2, dpi=300)




####################### PREPARATION OF THE 16S OBJECT TO ANALYSE #######################

if(! "proof1" %in% ls() ){
  stop("Filter the data before generating the main objects!")
}

sample_data(data)[["Sample_Family"]] <- gsub("_bis","", sample_data(data)[["Sample_Family"]] )   # the info about the replicates is removed in the main dataset

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
#data.class = tax_glom(data, taxrank = "Class", NArm = F)
#data.order = tax_glom(data, taxrank = "Order", NArm = F)
#data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  #data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  #data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  #data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  #Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  #Taxa.class<-as.data.frame(tax_table(data.class))
  #Taxa.order<-as.data.frame(tax_table(data.order))
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


# save.image(file="data_prepared_after_filters_CCHS.RData")




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

png(file="Data_check/Rarefaction_curve.png",width=2000,height=1600, res=300)
r<-rarecurve(t(as(otu_table(data.genus),"matrix")), ylab = "Genera" , xlab="Filtered reads amount", step=100,label=F)
evalslopes(r,sample_data(data.genus)$Sample_Family,lim=0.0001,cex=1)
dev.off()

rm(r)




#################### % ASSIGNED IN SILVA #########################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
#b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
#c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
#d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
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
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  # write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  # write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}




######################## ABUNDANCES BAR PLOT ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
# fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
# fill_color_16<-c("darkblue","brown4","wheat","coral","yellow3","darkmagenta", "deepskyblue3","firebrick4","gray60","pink3","darkgreen","gold","violet", "red","chartreuse3","darkslategray3")


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
  tabella$Condition2<-factor(tabella$Condition2, levels = c("CCHS","Control"))
}
ggplot(data=tabella, aes(x=Sample_Family, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Condition2),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size=9) + 
  scale_fill_manual(values=fill_color_5) +
  scale_y_continuous(expand = c(0,1)) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size = 6.5), 
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(-10,5,11,0),
        strip.text = element_text ( size = 10.5 ),
        legend.text = element_text ( size = 11 )
        ) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Percentual abundance",
       # title = "Five most abundant phyla",
       caption = " 'Others' includes every phylum below rank 5 ",
       fill="")
ggsave(file="Results/Abundances/TOP_5_phyla.png",width=4.5,height=3.35, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.csv2(file = "Results/Abundances/TOP5_phyla_average_abundance.csv", row.names = F, to_save )

rm(top, prune.dat_top,tabella, tabella_top)



# TOP Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:16]
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
  tabella$Genus <- gsub ("Prevotellaceae_NK3B31_group","Prevotellaceae_NK3B31",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
  tabella$Condition2<-factor(tabella$Condition2, levels = c("CCHS","Control"))
}
fill_color_16<-c(  "darkblue", "lightblue3", "brown","wheat","coral","yellow4",
                  "darkmagenta", "deepskyblue3","orange2","gray75",
                  "pink3","darkgreen","gold","violet", "red","chartreuse3","darkslategray3")

genera_plot <- ggplot(data=tabella, aes(x=Sample_Family , y=Abundance, fill=Genus)) + 
  facet_grid(cols= vars(Condition2),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + 
  theme_classic(base_size=9) + 
  scale_fill_manual(values=fill_color_16) +
  scale_y_continuous(expand = c(0,1)) +
  labs(x="", y="Percentual abundance",
       caption = " 'Others' includes every genus below rank 16 ", fill="")
# genera_plot + 
#   theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size = 6.5), 
#         legend.key.height  = unit(0.25, "cm"),
#         legend.key.width  = unit(0.15, "cm"),
#         legend.margin = margin(-15,30,0,0),
#         strip.text = element_text ( size = 10.5 ),
#         legend.text = element_text ( size = 8.15 ),
#         legend.position="bottom"
#         ) +
#   guides(fill=guide_legend(nrow=4) )
# ggsave(file="Results/Abundances/TOP_Genera_1.png",width=4.5,height=3.35,dpi=300)
### again, but this time on the legend will be on the right side
genera_plot + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size = 5.8), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=8.5),
        axis.text.y = element_text(size=6),
        legend.key.height  = unit(0.382, "cm"),
        legend.key.width  = unit(0.15, "cm"),
        legend.margin = margin(-12.2,-2,0,-4),
        plot.margin  = margin(1,1,1,1),
        strip.text = element_text ( size = 9.8 ),
        legend.text = element_text ( size = 7.2 ),
        legend.position="right"
  ) +
  guides(fill=guide_legend(nrow=17) )
ggsave(file="Results/Abundances/TOP_Genera_2.png",width=4.8,height=3.2,dpi=300)

dev.off()

# means of TOP Genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.csv2(file = "Results/Abundances/TOP_genera_average_abundance.csv", row.names = F, to_save )
suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))




#################### RATIO FIRMICUTES/BACTEROIDES ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
sample_names(data_fb) <- sample_data(data_fb)$Sample_Family
data_fb <- subset_samples(data_fb, Sample_name !="C1") # this is THE outlier sample
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-cbind.data.frame(otu_table(data_fb),tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.numeric(ratio_fb[F_index,]) / as.numeric(ratio_fb[B_index,]) )
# NB: both the numerator and denom are multiplied by 10 to compute the ratio to avoid the issues arising from fraction values (e.g. 10 / 0.08 , which improperly increase the final product)
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

ratio_fb<-ratio_fb[sample_data(data_fb)$Sample_Family, ] # same order
ratio_fb<-cbind.data.frame(ratio_fb,sample_data(data_fb)[,c("Sample_Family","Condition2")])

# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Condition2, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Condition2=="CCHS","Ratios_Mean"]<-Ratios_Mean["CCHS"]
ratio_fb[ratio_fb$Condition2=="Control","Ratios_Mean"]<-Ratios_Mean["Control"]
#ratio_fb[ratio_fb$Condition=="Control_parent","Ratios_Mean"]<-Ratios_Mean["Control_parent"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Condition2, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Condition2=="CCHS","Ratios_st_err"]<-Ratios_st["CCHS"]/sqrt(length(which(ratio_fb$Condition2=="CCHS")))
ratio_fb[ratio_fb$Condition2=="Control","Ratios_st_err"]<-Ratios_st["Control"]/sqrt(length(which(ratio_fb$Condition2=="Control")))
#ratio_fb[ratio_fb$Condition=="Control_parent","Ratios_st_err"]<-Ratios_st["Control_parent"]/sqrt(length(which(ratio_fb$Condition=="Control_parent")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/Abundances/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio.csv", ratio_fb)

# ratio_fb$Condition<-factor(ratio_fb$Condition, levels = c("CCHS","Control","Control_parent"))
ratio_fb$Condition2<-factor(ratio_fb$Condition2, levels = c("CCHS","Control"))



####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

# ratio_fb<-ratio_fb[order(ratio_fb$Sample_number,ratio_fb$Time),] # if paired
res<-wilcox.test(Ratio ~ Condition2, data=ratio_fb, pair=F)
p_val<-round(res$p.value, digits = 2)
p_val    # 0.38 with the outlier sample, 0.22 without

# DT<-FSA::dunnTest(Ratio ~ Condition, data=ratio_fb, method="bh")
# DT

# # exporting the results
# con <- file("Results/Abundances/Ratio_Firmi_Bacteroi/Statistical_tests_FB_ratios.txt") # create this file and connect it to R through "con" object 
# sink(con, append = TRUE) # redirect STR ERROR to "con"
# cat("Kruskal Wallis test \n", fill=T)
# cat("Ratio~Condition :   V=",res$statistic, "p-value=", p_val, "\n",fill=T) # da adattare al kruskal
# cat("\n\n Dunnet test (p-value adjusted through BH)\n", fill=T)
# print(DT$res)
# cat("\n\n\n Shapiro Wilk test p-value: ",distr$p.value)
# sink()
# close(con)



######## BAR PLOT
ggplot(data=ratio_fb, aes(x=Sample_Family, fill=Condition2, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Control"="deepskyblue2",
                               "CCHS"="coral"
                               # "Control_parent"="deepskyblue"
                               )
                    ) +
  facet_grid2(.~Condition2, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= Ratios_Mean, group="Condition"))+
  scale_y_sqrt(breaks=c(1,5,seq(10,150,10)) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Wilcoxon Wallis p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi/Ratio_barplot.png",width=6,height=5, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))




########################## ALPHA DIVERSITY ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Condition2", color="Condition2")
pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
{ H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
  obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
  # adding evenness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Condition2, New_data$Family_group),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition2, y=value, color=NULL), alpha=0) +
  theme_classic2(base_size = 10.5) + 
  scale_color_manual(values = c("Control"="deepskyblue2", "CCHS"="tomato")) + #, "Control_parent"="deepskyblue")) +
  labs(x=""
       # title="Alpha diversity in Control and CCHS"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=8)) +
  stat_compare_means(aes(group = Condition2), label="p.format", method = "wilcox.test",
                     label.x= 0.95, size=2.5, label.y.npc = "top", vjust=-0.5, hjust=-0.35)
ggsave(file="Results/Alfa_diversity_GENUS_CONDITION_Wilcox.png", width = 4.5,height =3.5, dpi=300)
# # again,  ID and lines
# pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition2, y=value, color=NULL), alpha=0) +
#   theme_classic2() + 
#   scale_color_manual(values = c("Control"="deepskyblue2", "CCHS"="coral")) + #, "Control_parent"="deepskyblue")) +
#   labs(x="", title="Alpha diversity in Control and CCHS",
#        caption="Lines connect the subject of the same family" ) +
#   guides(fill="none", color="none") +
#   geom_text(aes(label=pAlpha$data$Sample_Family), color="black", size=2, show.legend = FALSE) +
#   theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
#   geom_line(aes(group=Family_group),linewidth=0.2, color="grey") +
#   stat_compare_means(aes(group = Condition2), label="p.format", method = "wilcox.test", label.x= 0.95, size=3, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
# ggsave(file="Results/Alfa_diversity_GENUS_CONDITION_Wilcox_Lines.png", width = 5.2,height =4.5, dpi=300)
# again, age
# again,  ID and lines
pAlpha2<-pAlpha
pAlpha2$layers[[1]] <- NULL # removing geom point layer
pAlpha2$data$Age_Condition <- paste0( pAlpha2$data$Age_3groups , pAlpha2$data$Condition2 )
pAlpha2$data$ID_special <- pAlpha2$data$Sample_Family
pAlpha2$data$ID_special[! pAlpha2$data$ID_special %in% c("M1","L1","T3","P1","F1","W1")] <- ""
pAlpha2$data$ID_special_LEFT <- pAlpha2$data$ID_special
pAlpha2$data$ID_special_LEFT[ pAlpha2$data$Condition2!="CCHS" ] <- ""
pAlpha2$data$ID_special_RIGHT <- pAlpha2$data$ID_special
pAlpha2$data$ID_special_RIGHT[ pAlpha2$data$Condition2=="CCHS" ] <- ""

pAlpha2 +
  geom_boxplot(data=pAlpha2$data, aes(x=Condition2, y=value, color=NULL, fill=Condition2),
               alpha=0.08, size = 0.12 ) +
  geom_point( aes(color=Age_Condition , shape=HIRSCHSPRUNG) , size= 1.35, show.legend = F) +
  theme_classic2(base_size = 9.8) + 
  scale_color_manual(values = c("<20Control"="cadetblue2", "20-35Control"="deepskyblue3", ">35Control"="blue",
                                "<20CCHS"="indianred1", "20-35CCHS"="red1", ">35CCHS"="darkred"
                                )) +
  scale_fill_manual(values = c("Control"="cadetblue2","CCHS"="tomato2") ) +
  labs(x="" #,
       #  title="Alpha diversity in Control and CCHS" ,
       # caption="The ages are displayed on each point"
       ) +
  guides(fill="none", color="none", shape="none") +
  geom_text(aes(label=pAlpha2$data$ID_special_LEFT), color="gray25", size=2.4, hjust=-0.3, show.legend = FALSE) +
  geom_text(aes(label=pAlpha2$data$ID_special_RIGHT), color="gray25", size=2.4, hjust=1.3, show.legend = FALSE) +
  theme(axis.text.x= element_text(angle=22, vjust=1, hjust=1, size=8)) +
  stat_compare_means(aes(group = Condition2), label="p.format", method = "wilcox.test",
                     label.x= 0.85, size=2.45, label.y.npc = "top", vjust=-0.38, hjust=-0.35)
ggsave(file="Results/Alfa_diversity_GENUS_CONDITION_Wilcox_colorAGE.png", width = 3.95,height =3.2, dpi=300)
# again, but p-value computes the ages
pAlpha2$mapping$x<-NULL
pAlpha2 + geom_boxplot(data=pAlpha2$data, aes(x=Age_groups, y=value, color=NULL), alpha=0) +
  geom_point( aes(x=Age_groups, y=value, color=Condition2)) +
  theme_classic2( ) + 
  scale_color_manual(values = c("Control"="deepskyblue2", "CCHS"="coral")) +
  labs(x="", title="Alpha diversity in Control and  (p-value on ages)" ,
       caption="The ages are displayed on each point") +
  guides(fill="none", color="none") +
  geom_text(aes(x=Age_groups,label=pAlpha$data$Age), color="black", size=2.2, show.legend = FALSE) +
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(x=Age_groups,  group = Age_groups), label="p.format", method = "wilcox.test", label.x= 0.95, size=3, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="Results/Alfa_diversity_GENUS_AGE_Wilcox_colorCONDITION.png", width = 5.2,height =4.5, dpi=300)

suppressWarnings ( rm(con, pAlpha, alphadt,H, ev, obs, a, New_data) )




######################## TESTING THE CO-FACTORS ON BETA DIV #######################

suppressWarnings(rm(ASV.genus.prop))

metadata<-as(sample_data(data.genus.prop),"data.frame")

ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t

results_here <- NULL

contrast_lists <- c( "Condition2",
                     "Age_groups",
                     "Age_3groups",
                     "Family_group",
                     "Sex",
                     "Batch",
                     "Age_groups+Family_group",
                     "Age_groups+HIRSCHSPRUNG",
                     "Sex+Age_groups",
                     "Condition2+HIRSCHSPRUNG",
                     "Condition2+Family_group",
                     "Condition2+Age",
                     "Condition2+Age_groups",
                     "Age_3groups+Condition2",
                     "Condition2+Sex" ,
                     "Condition2+Batch",
                     "Condition2+Batch+Age",
                     "Condition2+Batch+Age+Sex",
                     "Condition2+Sex+Age_groups",
                     "Condition2+Sex+Age_groups+Family_group"
                     )

for (this_one in contrast_lists){
  set.seed(1994)
  perm_g<- vegan::adonis2( formula(paste( "sample_OTU ~", this_one)) , data=metadata, permutations = 9999, method="euclidean")
  which_rows <- row.names(perm_g) %in% colnames(Metadata)
  new_results_row<- cbind.data.frame( "Factor"= row.names(perm_g[which_rows, ]) ,
                                      "P-value" = as.data.frame(perm_g[which_rows, ]) [, "Pr(>F)" ] ,
                                      "Design" = rep(paste(this_one))
  )
  results_here <- rbind.data.frame( results_here , new_results_row , rep("") )
}

# IF THE ANOMALOUS SAMPLES ARE REMOVED...
set.seed(1994)
ASV.genus.prop2<-as.data.frame( otu_table(subset_samples(data.genus.prop , !Sample_Family%in%c("T3","P1","M1") ) ) )
meta2 <- as(sample_data(subset_samples(data.genus.prop,!Sample_Family%in%c("T3","P1","M1"))), "data.frame")
sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop2))) # samples has to be on rows --> t

perm_g<- vegan::adonis2( sample_OTU ~ Condition2 + Age , data=meta2, permutations = 9999, method="euclidean")
which_rows <- row.names(perm_g) %in% colnames(meta2)
new_results_row<- cbind.data.frame( "Factor"= row.names(perm_g[which_rows, ]) ,
                                    "P-value" = as.data.frame(perm_g[which_rows, ]) [, "Pr(>F)" ] ,
                                    "Design" = "Condition2+Age (NO OUTL nor YOUNG)"
                                    )
results_here <- rbind.data.frame( results_here , new_results_row , rep("") )

perm_diversity_res <- results_here
perm_diversity_res$Sign <- ""
perm_diversity_res$Sign[perm_diversity_res$`P-value`<=0.05 &perm_diversity_res$`P-value`>0 ] <- "*"

write.csv2(file= "Results/Beta_div/Co_factors_influence_on_genera_div.csv", x= perm_diversity_res,  row.names = F, quote = F)



# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
{BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  # Condition
  disper<-vegan::betadisper(BC.dist,metadata$Condition2)
  disp_genus<-vegan::permutest(disper,  permutations=9999)
  # Age
  disper<-vegan::betadisper(BC.dist,metadata$Age_3groups)
  disp_age<- vegan::permutest(disper, permutations=9999)
  a<-rbind.data.frame( disp_genus$tab$`Pr(>F)`[1] , disp_age$tab$`Pr(>F)`[1] )
  # Age no Outliers, Prob and <12y
  BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop2)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,meta2$Age_3groups)
  disp_no_outliers<- vegan::permutest(disper, permutations=9999)
  a<-rbind.data.frame( disp_genus$tab$`Pr(>F)`[1] , disp_age$tab$`Pr(>F)`[1] , disp_no_outliers$tab$`Pr(>F)`[1] )
}
colnames(a)<-c("permuted_p_value")
row.names(a)<-c("Condition","Age","Age no <12y or outliers")
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/Beta_div/Beta_disp_General_and_Pairwise_between_Conditions.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger (computed on genera)")
cat("\n", fill=TRUE)
a
sink()
close(con)

rm(disp_genus, disp_age, a, con)




####################### PCoA BETA DIV #########################

### on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Age_Condition<- paste0( sample_data(data.prop.labels)$Age_3groups , sample_data(data.prop.labels)$Condition2 )
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition2", shape="HIRSCHSPRUNG") +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral")) + # ,"Control_parent"="deepskyblue")) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10) +
  stat_ellipse(linewidth=0.15, aes(color=Condition2, shape=NULL))  +
  theme( legend.text = element_text(size=10)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_Family), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)",
       # subtitle= paste0("PERMANOVA Pr(>F) Control - CCHS: ", perm_g_H["Condition2", "Pr(>F)"] ),
       color="Condition", shape="Hirschsprung", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_1.png", width = 5, height = 4.2, dpi=300)
# only certain ID
these_IDs <- sample_data(data.sqrt_prop)$Sample_Family # to obtain the correct order
these_IDs [ ! these_IDs %in% c("M1","F1","L1","T3","P1","W1" ) ] <- ""   # everything else is just blank ...
# these_IDs[ these_IDs=="U2" ]<- "\nU2"  # to avoid sovrapositions

plot_ordination(data.sqrt_prop, ordBC, shape="HIRSCHSPRUNG") +
  geom_point(aes(color=Age_Condition), size=3.2, alpha= 0.55, show.legend = F) +
  scale_color_manual(values = c("<20Control"="cadetblue1", "20-35Control"="deepskyblue3", ">35Control"="darkblue",
                                "<20CCHS"="indianred1", "20-35CCHS"="red1", ">35CCHS"="darkred",
                                "Control"="royalblue", "CCHS"="red2"
                                )
                     ) +
  theme_classic(base_size = 9) +
  stat_ellipse(linewidth=0.12, aes(color=Condition2, shape=NULL))  +
  theme( legend.text = element_text(size=7),
         legend.margin = margin(1,-3,0,-5)
         ) +
  geom_text(aes(label= these_IDs ), lineheight=0.2, color="grey10", vjust=1.32 , size=3, show.legend = FALSE) +
  labs(color="Condition", shape="Hirschsprung", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")
       )
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_1_few_ID.png", width = 4, height = 3.45, dpi=300)
# CCHS + HISRCHS group
sample_data(data.sqrt_prop)$Combination <- paste0(sample_data(data.sqrt_prop)$Condition2, sample_data(data.sqrt_prop)$HIRSCHSPRUNG)
sample_data(data.sqrt_prop)$Combination <- gsub("No","", sample_data(data.sqrt_prop)$Combination)
sample_data(data.sqrt_prop)$Combination <- gsub("CCHSYes","CCHS+HIRSCHSPRUNG", sample_data(data.sqrt_prop)$Combination)
plot_ordination(data.sqrt_prop, ordBC, color = "Combination") +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral", "CCHS+HIRSCHSPRUNG"="red3")) + # ,"Control_parent"="deepskyblue")) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10) +
  stat_ellipse(linewidth=0.15)  +
  theme( legend.text = element_text(size=7.5)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_Family), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", 
       #subtitle= paste0("PERMANOVA Pr(>F) Control - CCHS: ", perm_g_H["Condition2", "Pr(>F)"] ),
       color="Condition", shape="Hirschsprung",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_2.png", width = 5, height = 4.2, dpi=300)
# without IDs, with family lines
plot_ordination(data.sqrt_prop, ordBC, color = "Condition2", shape="HIRSCHSPRUNG") +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral")) + # ,"Control_parent"="deepskyblue")) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10) +
  stat_ellipse(linewidth=0.15)  +
  geom_line(aes(group=Family_group),linewidth=0.15, color="grey") +
  theme( legend.text = element_text(size=10)) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)",
       #subtitle= paste0("PERMANOVA Pr(>F) Control - CCHS: ", perm_g_H["Condition2", "Pr(>F)"] ),
       color="Condition", shape="Hirschsprung", 
       caption= "Lines connects subjects of the same family",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_3.png", width = 5, height = 4.2, dpi=300)
# Digestive_issues  + Age
plot_ordination(data.sqrt_prop, ordBC, color = "Condition2", shape="Digestive_issues") +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral")) + # ,"Control_parent"="deepskyblue")) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10) +
  stat_ellipse(linewidth=0.15)  +
  theme( legend.text = element_text(size=10)) +
  # geom_text(aes(label=sample_data(data.sqrt_prop)$Age), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)",
       caption="The number on the points indicates the sample age",
       color="Condition", shape="Digestive_issues", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_4.png", width = 5, height = 4.2, dpi=300)
# without IDs, with family lines (again, no second group)
plot_ordination(data.sqrt_prop, ordBC, color="Condition2", shape="HIRSCHSPRUNG") +
  geom_point(size=3, alpha= 0.6 ) +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral")) + # ,"Control_parent"="deepskyblue")) +
  theme_classic(base_size = 10) +
  stat_ellipse(linewidth=0.15, aes(color=Condition2, group=Condition2))  +
  geom_line(aes(group=Family_group),linewidth=0.15, color="grey") +
  theme( legend.text = element_text(size=10)) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)",
       # subtitle= paste0("PERMANOVA Pr(>F) Control - CCHS: ", perm_g_H["Condition2", "Pr(>F)"] ),
       color="Condition", shape="Hirschsprung", 
       caption= "Lines connects subjects of the same family",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_5.png", width = 5, height = 4.2, dpi=300)
#NB: the warning message is attended when building the last plot


# EXTRA: only AGE
age_white<- sample_data(data.sqrt_prop)$Age
age_white[! age_white>=35]<- ""
age_black<- sample_data(data.sqrt_prop)$Age
age_black[! age_black<35]<- ""
plot_ordination(data.sqrt_prop, ordBC, color="Age_3groups", shape="Condition2") +
  geom_point(size=3.8, alpha= 0.6 ) +
  scale_color_manual(values=c("<25"="deepskyblue",
                              ">25"="blue",
                              "<20"="cadetblue3",
                              "20-35"="deepskyblue1",
                              ">35"="blue"
                              )
                     ) +
  theme_classic(base_size = 9.5) +
  geom_text(aes(label=age_white), size=2.35, show.legend = FALSE , color="grey95") +
  geom_text(aes(label=age_black), size=2.35, show.legend = FALSE , color="black"  ) +
  stat_ellipse(linewidth=0.15, aes(color=Age_3groups, group=Age_3groups) )  +
  theme( legend.text = element_text(size=9)) +
  labs(
       color="Age groups", shape="Condition", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Data_check/PCoA_Beta_diversity_AgeGroups.png", width = 4, height = 3.45, dpi=300)
# Yes, warnings are attended ... but the plotting still works!




# #################### VEEN DIAGRAM ##########################
# 
# data.genus.temp<-data.genus.prop
# tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
# data.venn<-data.genus.temp
# 
# ### abundance AND prevalence filter (0.1% abundance at least in 3 sample)
# who<-as.data.frame(otu_table(data.venn))
# who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
# who<-who[!rowSums(who)>4,] # more than 4 "points" --> at least in 5 samples (one third of the CCHS group sample size)
# who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
# data.venn<-subset_taxa(data.venn, ! Genus %in% who)
# 
# 
# Control<-subset_samples(data.venn, Condition=="Healthy")
# Control<-as.character(tax_table(prune_taxa(taxa_sums(Control)>0, Control))[,"Genus"])
# 
# CCHS<-subset_samples(data.venn, Condition=="CCHS")
# CCHS<-as.character(tax_table(prune_taxa(taxa_sums(CCHS)>0, CCHS))[,"Genus"])
# 
# Control_parent<-subset_samples(data.venn, Condition=="Control_parent")
# Control_parent<-as.character(tax_table(prune_taxa(taxa_sums(Control_parent)>0, Control_parent))[,"Genus"])
# 
# 
# ONLY_IN_CCHS<- CCHS[! CCHS %in% Control_parent & ! CCHS %in% Control]
# ONLY_IN_CCHS<- paste(ONLY_IN_CCHS, collapse = ", ")
# head(ONLY_IN_CCHS)
# 
# ONLY_IN_Control<- Control[! Control %in% Control_parent & ! Control %in% CCHS]
# ONLY_IN_Control<- paste(ONLY_IN_Control, collapse = ", ")
# head(ONLY_IN_Control)
# 
# ONLY_IN_Control_parent<- Control_parent[! Control_parent %in% Control & ! Control_parent %in% CCHS]
# ONLY_IN_Control_parent<- paste(ONLY_IN_Control_parent, collapse = ", ")
# head(ONLY_IN_Control_parent)
# 
# # IN_Control_AND_Control_parent_BUT_NOT_CCHS <- Control[Control %in% Control_parent & ! Control %in% CCHS ] 
# # head(IN_Control_AND_Control_parent_BUT_NOT_CCHS)
# # 
# IN_Control_parent_AND_CCHS_BUT_NOT_Control <- Control_parent[Control_parent %in% CCHS & ! Control_parent %in% Control ] 
# head(IN_Control_parent_AND_CCHS_BUT_NOT_Control)
# 
# 
# con<-file("Results/Beta_div/Venn_Diagr_exclusive_bacteria.txt")
# sink(con, append=TRUE)
# cat("ONLY IN Control_parent", fill=TRUE)
# cat(ONLY_IN_Control_parent)
# cat("\n\nONLY IN Control", fill=TRUE)
# cat(ONLY_IN_Control)
# cat("\n\nONLY IN CCHS", fill=TRUE)
# cat(ONLY_IN_CCHS)
# cat("\n\nIn Control_parent and CCHS but not in Control ", fill=TRUE)
# cat(IN_Control_parent_AND_CCHS_BUT_NOT_Control)
# # cat("\n\nIn Control and Control_parent but not in CCHS ", fill=TRUE)
# # cat(IN_Control_AND_Control_parent_BUT_NOT_CCHS)
# sink()
# close(con)
# 
# 
# x<-list(Control=Control,CCHS=CCHS,Control_parent=Control_parent)
# ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
#        fill_color = c("chartreuse","coral","deepskyblue")) +
#   theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
#   labs(title = "Venn Diagram \n(only genera with minimal abundance > 0.1% \n at least in five samples)")
# ggsave(filename = "Results/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
# dev.off()
# 
# 
# suppressWarnings(rm(ONLY_IN_CCHS, ONLY_IN_Control, ONLY_IN_Control_parent, 
#                     x, con, CCHS, Control, Control_parent, data.venn, who))



# con<-file("Results/Beta_div/NO_exclusive_bacteria.txt")
# sink(con, append=TRUE)
# cat("Any bacteria is exclusive of any group (NB: considered only bacteria with minimal abundance >0.1% and present at least in two sample", fill=TRUE)
# sink()
# close(con)




############## DIFFERENTIAL ABUNDANCES WITH DESEQ2 (CONDITIONS) #############

if(! "proof1" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

#system(" echo 'In case of pair analysis, the log2foldchange displayes the change from Control to CCHS, see https://support.bioconductor.org/p/105981/ \n' > Results/DA_DESeq2/NB.txt ")  
#system(" echo 'Every result with a baseMean value lower than 50 has been arbitrarly filtered.' >> Results/DA_DESeq2/NB.txt ")  

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition2)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition2", "Control", "CCHS"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj <= 0.05) & abs(res$log2FoldChange)>0, ]
  res<-res[res$baseMean >= 25, ] # arbitrary threshold to avoid the most noisy result
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Control_vs_CCHS.csv"), row.names = F, quote=F, na = "")
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

#View(Res_tot)   # the higher level results are redundant...
write.csv(Res_tot, file="Results/DA_DESeq2/CONDITIONS_Every_result_DESeq2.csv", row.names = F)

Table_tot$Condition2<-factor(Table_tot$Condition2, levels = unique(Table_tot$Condition2))
Table_tot<-Table_tot[order(match(Table_tot$Condition2, levels(Table_tot$Condition2))), ]
# View(Table_tot)

# Table_tot$Bacteria<-gsub("Acidaminococcus","Acidaminoco.", Table_tot$Bacteria, fixed=T)
# the lines connecting relatives of the same groups are hereby removed to avoid graphical issues
Table_tot$Combination<- paste0(Table_tot$Condition,Table_tot$Family_group,Table_tot$Bacteria)
# adding random numbers instead of the proper family group name whereas the Condition group is the same
Table_tot[ duplicated(Table_tot$Combination), "Family_group"] <- 1:length(Table_tot[ duplicated(Table_tot$Combination), "Family_group"])


tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Condition2)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Condition2)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

# to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
levels(tabella_g$Xaxis)
levels(tabella_g$Condition2)

tabella_g$Age_Condition <- paste0( tabella_g$Age_3groups , tabella_g$Condition2 )
these_colors <- c("<20Control"="lightblue2", "20-35Control"="deepskyblue2", ">35Control"="darkblue",
                  "<20CCHS"="orange", "20-35CCHS"="coral", ">35CCHS"="red" 
                  )

plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Condition2, shape=HIRSCHSPRUNG)) +
  theme_classic(base_size = 8.5) +
  geom_boxplot( aes( x=Xaxis, shape="none", fill="Condition2", color="black" ) ,
                size = 0.45, alpha = 0 , outliers = F ) +
  #scale_color_manual(values = c("CCHS"="coral","Control"="deepskyblue2")) +
  scale_color_manual(values = these_colors ) +
  scale_shape_manual(values = c(16,3,17)) +  # the "3" is the "none" level of the boxplot
  # facet_wrap2(nrow=1, factor(Taxa,levels = "Genus")~Bacteria,
  #             labeller = labeller(group = label_wrap_gen(width = 34)),
  #             scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  facet_wrap2(nrow=2, .~Bacteria, scale= "free") +
  geom_point(aes(color=Age_Condition), size=0.8, alpha=0.95) +
  geom_point(aes(color=Age_Condition), size=2.2, alpha=0.45) +
  # geom_text(aes(label=Family_group), size=2.3, color="darkgray") +
  # geom_line(aes(group=Family_group), linewidth=0.15, color="grey") +
  theme(strip.text.x=element_text(size=8,colour="black"), 
        strip.switch.pad.wrap = unit(20,"line") ,
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=5.85, angle = 0) ,
        axis.text.x = element_text(size=6.25, angle = 0, color="gray50" ),
        #axis.text.x = element_text(size=6.25, angle = 20, vjust=1, hjust=1 , color="gray45" ),
        panel.grid.major.y= element_line(linewidth = 0.12, color="gray45"),
        panel.grid.minor.y= element_line(linewidth = 0.085, color="gray70" ),
        plot.margin =margin(1, 0.5, 1, 0.5), 
        legend.position = "none"
        ) +
  scale_x_discrete(labels=unique(levels(tabella_g$Condition2)), expand=c(0,0.55))
  
plot_1 +
  labs(y="Percent abundance", fill="Condition", x="",
       caption="CCHS patients with Hirschsprung disease are rapresented as triangles"
       )
ggsave(filename = "Results/DA_DESeq2/Condition_Plot_genera_DeSeq2.png", width = 4.5, height = 4, dpi=300)

# extra: with Codes
plot_1 +
  geom_text(aes(label=Sample_Family), color="black", size=2.5, show.legend = FALSE) +
  labs(y="Proportional Abundance", fill="Condition", x="",
       caption="CCHS patients with Hirschsprung disease are rapresented as triangles")
ggsave(filename = "Results/DA_DESeq2/Condition_Plot_genera_DeSeq2_with_sample_code.png", width = 4.5, height = 4, dpi=300)

### with certain samples only ...
# these_IDs <- tabella_g$Sample_Family # to obtain the correct order
# these_IDs [ ! these_IDs %in% c("M1","L1","T3","P1","W1","F1") ] <- ""   # everything else is just blank ...
these_IDs_left <- tabella_g$Sample_Family
these_IDs_left[ ! these_IDs_left %in% c("T3","W1") ] <- ""    # everything else is just blank ...
these_IDs_right <- tabella_g$Sample_Family
these_IDs_right[ ! these_IDs_right %in% c("M1","L1","P1","F1") ] <- ""

plot_1 +
  geom_text(aes(label=these_IDs_left), color="gray25", size=2.25, hjust=-0.25, show.legend = FALSE) +
  geom_text(aes(label=these_IDs_right), color="gray25", size=2.25, hjust=1.25, show.legend = FALSE) +
  labs(y="Percent Abundance", fill="Condition", x="",
       caption="CCHS patients with Hirschsprung disease are rpresented as triangles")
ggsave(filename = "Results/DA_DESeq2/Plot_genera_DeSeq2_with_specific_samples.png", width = 4.5, height = 4, dpi=300)

dev.off()



# ##### removing redundant results
# { 
#   # 1) if redundant then same abundances
#   auto_Max<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, max )
#   auto_Min<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, min )
#   auto_Median<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, median )
#   # reordering the names (lower ranks first --> the redundants are those at higher levels)
#   ordered_names <- c( unique(Table_tot[Table_tot$Taxa=="Genus", "Bacteria" ]),
#                       unique(Table_tot[Table_tot$Taxa=="Family", "Bacteria" ]),
#                       unique(Table_tot[Table_tot$Taxa=="Order", "Bacteria" ]),
#                       unique(Table_tot[Table_tot$Taxa=="Class", "Bacteria" ]),
#                       unique(Table_tot[Table_tot$Taxa=="Phylum", "Bacteria" ])
#   )
#   ordered_names[!is.na(ordered_names)]
#   auto_Max<-auto_Max[ ordered_names  ]
#   auto_Min<-auto_Min[ ordered_names  ]
#   auto_Median<-auto_Median[ ordered_names ] # the table it self follows the correct taxonomic order because it is built in this way
#   # same abundances
#   auto_redund<-names(auto_Max)[ duplicated(auto_Min) & duplicated(auto_Max) & duplicated(auto_Median)  ]
#   
#   # 2) if redundant then same ASV
#   auto_ASV<-Res_tot$ASV[duplicated(Res_tot$ASV)]
#   auto_ASV_Names<-unique(Table_tot$Bacteria[Table_tot$OTU %in% auto_ASV])
#   auto_redund<-auto_redund[auto_redund %in% auto_ASV_Names]
# }
# 
# Redund<-auto_redund
# 
##############  ... everything is a redundant result !!!


# tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)
# plot_2<-ggplot(tabella_2, aes(x= Xaxis, y=Abundance, fill=Condition2, shape=HIRSCHSPRUNG)) +
#   theme_classic(base_size = 11) +
#   scale_color_manual(values = c("CCHS"="coral","Control"="deepskyblue2")) +
#   facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
#   geom_point(aes(color=Condition2), size=1.8) +
#   geom_point(aes(color=Condition2), size=3, alpha=0.5) +
#   # geom_text(aes(label=Family_group), size=2.3, color="darkgray") +
#   geom_line(aes(group=Family_group), linewidth=0.15, color="grey") +
#   theme(strip.text.x=element_text(size=9.2,colour="black"),
#         strip.switch.pad.wrap = unit(10,"line"),
#         axis.text.y = element_text(size=8.5))+
#   scale_x_discrete(labels=rep(unique(levels(tabella_2$Condition2))), expand=c(0,0.5)) +
#   theme(plot.title= element_text(size=11)) +
#   theme(panel.grid.minor.y= element_blank()) +
#   guides( fill=guide_legend(nrow=1) ) +
#   #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_2$Abundance+2),2))) +
#   theme(legend.margin=margin(0, 0, 0, 0), legend.position="none")
# plot_2 +
#   labs(title= "Differently abundant families and orders", y="Proportional Abundance", fill="Condition2", x="",
#        caption="Lines connect subjects of the same family between the two groups\nCCHS patients with Hirschsprung disease are rapresented as triangles")
# ggsave(filename = "Results/DA_DESeq2/Plot_higher_taxLevel_DESeq2_Condition2_no redundants.png", width = 6, height = 4, dpi=300)
# dev.off()
# 
# 
# 
# 
# ### UNIQUE PLOT (GENERA + HIGHER)
# table_reunited<- rbind(tabella_g, tabella_2[ ! tabella_2$Bacteria %in% Redund, ] )
# ggplot(table_reunited, aes(x= Xaxis, y=Abundance, fill=Condition2, shape=HIRSCHSPRUNG)) +
#   theme_classic(base_size = 10) +
#   scale_color_manual(values = c("CCHS"="coral","Control"="deepskyblue2")) +
#   facet_wrap2(nrow=1,factor(Taxa,levels = c("Order","Family","Genus"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
#   geom_point(aes(color=Condition2), size=1.8) +
#   geom_point(aes(color=Condition2), size=3, alpha=0.5) +
#   # geom_text(aes(label=Family_group), size=2.3, color="darkgray") +
#   geom_line(aes(group=Family_group), linewidth=0.15, color="grey") +
#   theme(strip.text.x=element_text(size=7.2,colour="black", face="italic"), 
#         strip.switch.pad.wrap = unit(10,"line"),
#         axis.text.x = element_text(size=8.4, angle=30, hjust=1, vjust=1),
#         axis.text.y = element_text(size=8) 
#         )+
#   scale_x_discrete(labels=rep(unique(levels(tabella_2$Condition2))), expand=c(0,0.5)) +
#   theme(panel.grid.minor.y= element_blank()) +
#   guides( fill=guide_legend(nrow=1) ) +
#   #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_2$Abundance+2),2))) +
#   theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
#   labs(y="Proportional Abundance", fill="Condition2", x="" )
# ggsave(filename = "Results/DA_DESeq2/Plot_all_good_results_DESeq2_no_redundants.png", width = 7, height = 4, dpi=300)
# dev.off()
# 


# for comparisons...
Genera.DESEQ2<-unique(tabella_g[tabella_g$Taxa=="Genus","Bacteria"])
suppressWarnings ( rm(tabella_2, plot_2, plot_1) )



### CIRCOTAX PLOT
#DE_for_Circo<-Res_tot[ Res_tot$ASV %in% unique(table_reunited$OTU) & Res_tot$Taxon!="Family" , c("log2FoldChange","Phylum","Class","Order","Family","Genus")]
DE_for_Circo<-Res_tot[ Res_tot$Taxon=="Genus" , c("log2FoldChange","Phylum","Class","Order","Family","Genus")]
DE_for_Circo$Genus <- gsub("_"," ",DE_for_Circo$Genus, fixed=T)
DE_for_Circo$Genus <- gsub("-","",DE_for_Circo$Genus, fixed=T)
DE_for_Circo$Genus <- gsub("[","",DE_for_Circo$Genus, fixed=T)
DE_for_Circo$Genus <- gsub("]","",DE_for_Circo$Genus, fixed=T)
DE_for_Circo$Genus <- gsub(" group","",DE_for_Circo$Genus, fixed=T)

CircoTax(DE_for_Circo,title="",
         # ramp=c("orange","lightskyblue1","cadetblue1"),
         fc_col=1,
         sort="no",
         fill_text = "log2FC",
         size_taxon_circo = 3.5
) +
  scale_fill_gradient2( high = "cadetblue2" , mid="azure1", low = "tomato1" ) +
  scale_x_discrete( expand= c(0,0))
ggsave("Results/DA_DESeq2/Condition_CircoTax_Condition_without_age.png", dpi= 300 , width = 4.8, height = 4.8)


 
 
############## PLOT DIFF ABUND THROUGH LINEAR MODELS ##################
# 
# table_pvalue<-NULL
# 
# for(i in unique(tabella_g$Bacteria)){
#   target<- psmelt (subset_taxa(data.genus.prop,Genus==i))
#   
#   model_tested <- lm( target$Abundance ~ target$Condition2 + target$Age)  # for the plot
#   ranked_model_tested <- lm( rank(target$Abundance) ~ target$Condition2 + target$Age)  # for the statistics
#   coef_model <- coef( model_tested )  # base intercept, slope and covariate-intercept-diff
#   test <- summary( ranked_model_tested )
#   p_value_value <- test$coefficients[2,4]  # Condition CCHS' p-value , after age correction
#   p_value_value <- round(p_value_value, 3)
#   p_value_value[p_value_value==0]<-"0.001"
#   age_pvalue <- round(test$coefficients[3,4] , 3)
#   
#   table_pvalue<- rbind.data.frame(table_pvalue, 
#                                   cbind(i,coef_model[1],coef_model[2],coef_model[3],p_value_value,age_pvalue)
#   )
# }
# 
# colnames(table_pvalue)<-c("variable",
#                           "intercept","Age_intercept_diff","Cond_slope",
#                           "p_value_COND_after_AGE","p_value_AGE"
# )
# 
# every_target <- psmelt (subset_taxa(data.genus.prop,Genus%in%unique(tabella_g$Bacteria)))
# # adding model coef to each sign genus
# row.names(table_pvalue)<-table_pvalue$variable
# every_target <- cbind.data.frame( every_target , table_pvalue[every_target$Genus, c("intercept","Age_intercept_diff","Cond_slope")] )
# 
# # plotting
# ggplot(data=every_target, aes(y=Abundance, x=Age, fill=Condition2) ) +
#   facet_grid2( rows = "Genus", scales = "free_y") +
#   geom_abline(intercept= as.numeric(every_target[["intercept"]]) , 
#               slope= as.numeric(every_target[["Cond_slope"]]) , 
#               colour="deepskyblue2" , linetype= 2 
#   ) +
#   geom_abline(intercept= as.numeric(every_target[["intercept"]])+ as.numeric(every_target[["Age_intercept_diff"]]) ,
#               slope= as.numeric(every_target[["Cond_slope"]]),
#               colour="coral", linetype= 2
#   ) +
#   scale_fill_manual(values=c( "CCHS"="coral","Control"="deepskyblue2"  )) +
#   scale_color_manual(values=c( "CCHS"="coral","Control"="deepskyblue2" )) +
#   geom_point(aes(color=Condition2), size= 1.3, alpha= 0.5) +
#   geom_point(aes(color=Condition2), size= 0.45, alpha= 0.9) +
#   theme_classic2(base_size = 7) +
#   theme(strip.text.x=element_text(size=12,colour="black"),
#         axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=5),
#         axis.text.y = element_text(size=5),
#         plot.margin = margin(5,2,0,2),
#         plot.title= element_text(size=8, hjust = 0.5, vjust=1.8) ,
#         plot.subtitle = element_text(size=6, vjust=1) ,
#         legend.key.size=unit(0.8,"cm"),
#         legend.text=element_text(size=7),
#         panel.grid.major.y = element_line(linewidth =0.12, color="gray"),
#         panel.grid.minor.y = element_line(linewidth =0.03, color="gray"),
#         panel.grid.major.x = element_blank()
#   ) +
#   guides( color="none", fill="none" ) +
#   labs(y="Percent Abundance", x="Age") +
#   scale_x_continuous(breaks = seq(10, max(every_target$Age)+1, 5)) #+
# # scale_y_continuous(breaks = c( seq(0, 15, 5) , seq(15,30,5), seq(30,100, 10 ) ) )
# 
# ggsave(file=paste0("Results/DA_DESeq2/Condition_PLOT_with_LM_model.png"), width = 2.3, height = 5, dpi=300)
 
 


############## DIFFERENTIAL ABUNDANCES WITH DESEQ2 IF THE ANOMALOUS SAMPLES ARE REMOVED #############

if(! "proof1" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
data_pruned<- subset_samples(data_pruned, ! Sample_Family %in% c("L1","T3","P1","M1") )  # L1 is under probiotic, T3 and P1 are younger than 12 years old, M1 is that strange CCHS outlier...
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

#system(" echo 'In case of pair analysis, the log2foldchange displayes the change from Control to CCHS, see https://support.bioconductor.org/p/105981/ \n' > Results/DA_DESeq2/NB.txt ")  
#system(" echo 'Every result with a baseMean value lower than 50 has been arbitrarly filtered.' >> Results/DA_DESeq2/NB.txt ")  

Table_tot<-NULL
Res_tot<-NULL

for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # if the paired samples are seen as numbers (it depends from on levels name)
  sample_data(d)[["Family_group"]]<-as.character(sample_data(d)[["Family_group"]])
  
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition2)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition2", "Control", "CCHS"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>0, ]
  res<-res[res$baseMean > 25, ] # arbitrary threshold to avoid the most noisy result
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
    #write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Control_vs_CCHS.csv"), row.names = F, quote=F, na = "")
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
#write.csv(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)

Table_tot$Condition2<-factor(Table_tot$Condition2, levels = unique(Table_tot$Condition2))
Table_tot<-Table_tot[order(match(Table_tot$Condition2, levels(Table_tot$Condition2))), ]
# View(Table_tot)

Table_tot$Bacteria<-gsub("_group","", Table_tot$Bacteria)
# the lines connecting relatives of the same groups are hereby removed to avoid graphical issues
Table_tot$Combination<- paste0(Table_tot$Condition,Table_tot$Family_group,Table_tot$Bacteria)
# adding random numbers instead of the proper family group name whereas the Condition group is the same
Table_tot[ duplicated(Table_tot$Combination), "Family_group"] <- 1:length(Table_tot[ duplicated(Table_tot$Combination), "Family_group"])


tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Condition2)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Condition2)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

# to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
levels(tabella_g$Xaxis)
levels(tabella_g$Condition2)


unique(tabella_g$Bacteria)
#tabella_g$Bacteria<-gsub("","",tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Condition2))  + # , shape=HIRSCHSPRUNG)) +
  theme_classic(base_size = 9.5) +
  scale_color_manual(values = c("CCHS"="coral","Control"="deepskyblue2")) +
  facet_wrap2(nrow=2, ~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Condition2), size=1.8, alpha=1) +
  geom_point(aes(color=Condition2), size=3, alpha=0.5) +
  # geom_text(aes(label=Family_group), size=2.3, color="darkgray") +
  # geom_line(aes(group=Family_group), linewidth=0.15, color="grey") +
  theme(strip.text.x=element_text(size=8,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_g$Condition2)), expand=c(0,0.5)) +
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), 
        legend.position = "none")
plot_1 +
  geom_text(aes(label=Age), color="gray40", size=2.8, show.legend = FALSE) +
  labs(title= "NO Lactogest, no <12 y old, no CCHS outlier", y="Proportional Abundance",
       fill="Condition2", x="",
       #caption="Lines connect subjects of the same family between the two groups\nCCHS patients with Hirschsprung disease are rapresented as triangles"
       #caption="Lines connect subjects of the same family between the two groups",
       )
ggsave(filename = "Results/DA_DESeq2/NO_OUTLIERS_OR_YOUNG_CONDITION_genera_DeSeq2.png", width = 6.25, height = 4, dpi=300)

  


############## DIFFERENTIAL ABUNDANCES WITH DESEQ2 (if Age is taken into account with COND) #############

if(! "proof1" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data)
sample_data(data_pruned)$Age_3groups <- factor( sample_data(data_pruned)$Age_3groups , levels = c("<20","20-35",">35"))
levels(sample_data(data_pruned)$Age_3groups) <- gsub("<","lesserT",levels(sample_data(data_pruned)$Age_3groups))
levels(sample_data(data_pruned)$Age_3groups) <- gsub(">","higherT",levels(sample_data(data_pruned)$Age_3groups))
levels(sample_data(data_pruned)$Age_3groups) <- gsub("-","TO",levels(sample_data(data_pruned)$Age_3groups), fixed = T)

#system(" echo 'In case of pair analysis, the log2foldchange displayes the change from Control to CCHS, see https://support.bioconductor.org/p/105981/ \n' > Results/DA_DESeq2/NB.txt ")  
#system(" echo 'Every result with a baseMean value lower than 50 has been arbitrarly filtered.' >> Results/DA_DESeq2/NB.txt ")  


Table_tot<-NULL
Res_tot<-NULL

#for( t in c("Genus","Family","Class","Order","Phylum") ){
for( t in c("Genus","Family") ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition2 + Age_3groups)
  # DEseq_data<-phyloseq_to_deseq2(d, ~Age_groups + Condition2)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition2", "Control", "CCHS"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj <= 0.05) & abs(res$log2FoldChange)>0, ]
  res<-res[res$baseMean > 25, ] # arbitrary threshold to avoid the most noisy result
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
    write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Control_vs_CCHS.csv"), row.names = F, quote=F, na = "")
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
# print(Res_tot)

#View(Res_tot)  # the Fam result is a redundant result
write.csv(Res_tot, file="Results/DA_DESeq2/CONDITION_plus_AGE_DESeq2.csv", row.names = F)


Table_tot$Condition2<-factor(Table_tot$Condition2, levels = unique(Table_tot$Condition2))
Table_tot<-Table_tot[order(match(Table_tot$Condition2, levels(Table_tot$Condition2))), ]

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Condition2)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))

# to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
levels(tabella_g$Xaxis)
levels(tabella_g$Condition2)

tabella_g$Age_Condition <- paste0( tabella_g$Age_3groups , tabella_g$Condition2 )
these_colors <- c("lesserT20Control"="lightblue2", "20TO35Control"="deepskyblue2", "higherT35Control"="darkblue",
                  "lesserT20CCHS"="orange", "20TO35CCHS"="coral", "higherT35CCHS"="red" 
)

plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Condition2, shape=HIRSCHSPRUNG)) +
  theme_classic(base_size = 8.5) +
  geom_boxplot( aes( x=Xaxis, shape="none", fill="Condition2", color="black" ) ,
                size = 0.45, alpha = 0 , outliers = F ) +
  #scale_color_manual(values = c("CCHS"="coral","Control"="deepskyblue2")) +
  scale_color_manual(values = these_colors ) +
  scale_shape_manual(values = c(16,3,17)) +  # the "3" is the "none" level of the boxplot
  # facet_wrap2(nrow=1, factor(Taxa,levels = "Genus")~Bacteria,
  #             labeller = labeller(group = label_wrap_gen(width = 34)),
  #             scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  facet_wrap2(nrow=2, .~Bacteria, scale= "free") +
  geom_point(aes(color=Age_Condition), size=0.8, alpha=0.95) +
  geom_point(aes(color=Age_Condition), size=2.2, alpha=0.45) +
  # geom_text(aes(label=Family_group), size=2.3, color="darkgray") +
  # geom_line(aes(group=Family_group), linewidth=0.15, color="grey") +
  theme(strip.text.x=element_text(size=8,colour="black"), 
        strip.switch.pad.wrap = unit(20,"line") ,
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=5.85, angle = 0) ,
        axis.text.x = element_text(size=6.25, angle = 0, color="gray50" ),
        #axis.text.x = element_text(size=6.25, angle = 20, vjust=1, hjust=1 , color="gray45" ),
        panel.grid.major.y= element_line(linewidth = 0.12, color="gray45"),
        panel.grid.minor.y= element_line(linewidth = 0.085, color="gray70" ),
        plot.margin =margin(1, 0.5, 1, 0.5), 
        legend.position = "none"
  ) +
  scale_x_discrete(labels=unique(levels(tabella_g$Condition2)), expand=c(0,0.55))

plot_1 +
  labs(y="Percent abundance", fill="Condition", x="",
       caption="CCHS patients with Hirschsprung disease are rapresented as triangles"
  )
ggsave(filename = "Results/DA_DESeq2/CONDITION_plus_AGE_Plot_genera_DeSeq2_ADJUSTED.png", width = 4.5, height = 4, dpi=300)

# extra: with Codes
plot_1 +
  geom_text(aes(label=Sample_Family), color="black", size=2.5, show.legend = FALSE) +
  labs(y="Proportional Abundance", fill="Condition", x="",
       caption="CCHS patients with Hirschsprung disease are rapresented as triangles")
ggsave(filename = "Results/DA_DESeq2/CONDITION_plus_AGE_Plot_genera_DeSeq2_ADJUSTED_Sample_codes.png", width = 4.5, height = 4, dpi=300)





################# PLS-DA and sPLS-DA ###########################

data_filtered<-filter_taxa( data.genus.prop ,function(x) mean(x)> 0.1, prune = T)
data_filtered<-subset_samples(data_filtered, !  Sample_Family %in% c("L1","T3","P1","M1") )  
# probiotics and/or under 12 years old and/or that anomalous CCHS


#### selecting three sample for each group (as test dataset)
{set.seed(1994)
  Control_train<-sample(sample_data(data_filtered)[sample_data(data_filtered)[["Condition2"]]=="Control", ] [["Sample_name"]] )[1:3]
}
{set.seed(1994)
  CCHS_train<-sample(sample_data(data_filtered)[sample_data(data_filtered)[["Condition2"]]=="CCHS", ] [["Sample_name"]] )[1:3]
}
# --> using the technical replicates to further test the sPLS-DA results ...
# replicates<-c("N1_bis","E2_bis","U1_bis")
# data.genus_unsub = tax_glom(unsubsetted, taxrank = "Genus", NArm = F)
# data.genus_unsub.prop <- transform_sample_counts(data.genus_unsub, function(ASV) ASV/sum(ASV)*100)
# Further_train<-sample(sample_data(data.genus_unsub.prop)[sample_data(data.genus_unsub.prop)[["Sample_Family"]]%in%replicates, ] [["Sample_name"]] )


# subsetting test and train datasets
data.selected<-subset_samples(data_filtered, ! Sample_name %in% c(Control_train, CCHS_train) )
metadata_PLSDA<-as(sample_data(data.selected),"data.frame")
data.train<-subset_samples(data_filtered, Sample_name %in% c(Control_train, CCHS_train) )
metadata_train<-as(sample_data(data.train),"data.frame")

# starting PLSDA
X <- as.matrix(t(otu_table(data.selected)))
Y <- as.factor(metadata_PLSDA$Condition2)

my.splsda <- splsda(X, Y, ncomp = 4) # test from 1 to 4 components


gc()
set.seed(1994)
perf.splsda <- perf(my.splsda, validation = "Mfold", # to assest the optimal number of comp
                    folds = 5, nrepeat = 100, cpus = 7,
                    progressBar = TRUE, auc = TRUE)

perf.splsda$choice.ncomp
#comp<-perf.splsda$choice.ncomp[1] # centroid distance (1)
### 1 is enough, but al least 2 are required for the plot!
comp <- 2

png(filename = "Results/PLS_DA/Error for each component PLS-DA.png", width = 2200, height = 1800, res=300)
plot(perf.splsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
title(main="Classification rate error \n for each component used to compute PLS-DA")
dev.off()

my.splsda <- splsda(X, Y, ncomp = comp)

n= 2 # choose the component to use depending on clustering efficiency in plot
png(filename = paste("Results/PLS_DA/normal_PLS_DA_comp_1_and_", n,".png", sep=""), width = 2400, height = 1600, res=300)
plotIndiv(my.splsda, comp = c(1,n), col = c("deepskyblue","coral"), cex = 1.8,
          group = metadata_PLSDA$Condition2, ind.names = TRUE, ellipse = TRUE, legend = F, 
          title = paste0("PLSDA among Control and CCHS \n on scaled proportional genus data \n (computed with ", comp, " components, plotted on comp 1 and ", n,")"))
dev.off()




################ sPLS-DA for variable selection

possible.pool<- c(1,3,5,8,10,15,20) # for test.keepX
gc()
set.seed(1994)
my.tuned.splsda <- tune.splsda(X, Y, ncomp = comp, validation = 'Mfold',
                               folds = 5, nrepeat = 100, cpus = 8, # use repeated cross-validation
                               dist = 'centroids.dist', test.keepX =  possible.pool, 
                               measure = "BER",  # use balanced error rate of dist measure
                               progressBar = T)

my.tuned.splsda$choice.ncomp$ncomp # 1

png(filename = paste("Results/PLS_DA/Error_of_splsda_depending_from_computing",comp,"components.png", sep="_"), width = 1500, height = 2000, res=300)
plot(my.tuned.splsda, col = color.jet(comp))
dev.off()
#optimal.comp<-my.tuned.splsda$choice.ncomp$ncomp  # 1
optimal.comp<-2

optimal.keepX <- my.tuned.splsda$choice.keepX[1:optimal.comp] # 20 and 1
optimal.keepX

final.splsda<-splsda(X,Y, ncomp=optimal.comp, keepX = optimal.keepX)
final.splsda$prop_expl_var$X   # 0.108  0.048

n=2 # choose which component plot depending on plot
# png(filename = paste("Results/PLS_DA/sparse_PLS_DA_comp_1_and",n, sep="_"), width = 2500, height = 1500, res=300)
# plotIndiv(final.splsda, comp = c(1,n), group = metadata_PLSDA2$Condition, ind.names = TRUE, legend = TRUE, ellipse = TRUE,
#           title = paste0("sparse PLSDA between Control and CCHS \n on scaled proportional genus data \n (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"))
# dev.off()
ggplot(mapping = aes(x=final.splsda$variates$X[,"comp1"], y=final.splsda$variates$X[,"comp2"], color=final.splsda$Y)) +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral"))+ #,"Control_parent"="red3")) +
  geom_point(size=2, alpha= 0.9) +
  geom_point(size=3.1, alpha= 0.5) +
  theme_classic(base_size = 10) +
  theme(title=element_text(size=10),
        legend.title =element_text(size=10)) +
  stat_ellipse(linewidth=0.15)  +  
  # geom_text(mapping = aes(label=row.names(final.splsda$X)), color="black", size=2.3) +
  labs(title="sparse PLS-DA between Control and CCHS", 
       subtitle=paste0("on scaled proportional genus data \n (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"),
       caption = paste0("selected ",optimal.keepX[1]," genera on LC1 and ",optimal.keepX[2]," genera on LC2"),
       color="Condition",
       x=paste("LC1: ",round(final.splsda$prop_expl_var$X[1]*100,digits = 2),"% variation"),
       y=paste("LC2: ",round(final.splsda$prop_expl_var$X[2]*100,digits = 2),"% variation"))
ggsave(filename = paste("Results/PLS_DA/sparsePLS-DA on comp 1 and",n,".png"), width = 5, height = 4.2, dpi=300)
# again, with sample names
ggplot(mapping = aes(x=final.splsda$variates$X[,"comp1"], y=final.splsda$variates$X[,"comp2"], color=final.splsda$Y)) +
  scale_color_manual(values=c("Control"="deepskyblue2","CCHS"="coral"))+ #,"Control_parent"="red3")) +
  geom_point(size=2, alpha= 0.9) +
  geom_point(size=3.1, alpha= 0.5) +
  theme_classic(base_size = 10) +
  theme(title=element_text(size=10),
        legend.title =element_text(size=10)) +
  stat_ellipse(linewidth=0.15)  +  
  geom_text(mapping = aes(label=row.names(final.splsda$X)), color="black", size=2.3) +
  labs(title="sparse PLS-DA between Control and CCHS", 
       subtitle=paste0("on scaled proportional genus data \n (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"),
       caption = paste0("selected ",optimal.keepX[1]," genera on LC1 and ",optimal.keepX[n]," genera on LC2"),
       color="Condition",
       x=paste("LC1: ",round(final.splsda$prop_expl_var$X[1]*100,digits = 2),"% variation"),
       y=paste("LC2: ",round(final.splsda$prop_expl_var$X[n]*100,digits = 2),"% variation"))
ggsave(filename = paste("Results/PLS_DA/sparsePLS-DA on comp 1 and",n,"_SAMPLE_NAMES.png"), width = 5, height = 4.2, dpi=300)



################ plotting loadings

loadings<-as.data.frame(final.splsda$loadings$X)
identical(row.names(Taxa.genus.update[row.names(loadings),]),row.names(loadings))
loadings$Genus<-Taxa.genus.update[row.names(loadings),"Genus"]

a<-loadings[loadings$comp1!=0,]
{tax_a<-Taxa.genus.update
  tax_a$ASV<- row.names(tax_a)
  tax_a<-subset(tax_a, ASV %in% row.names(a))
  tax_a<-tax_a[row.names(a),]
  identical(tax_a$ASV,row.names(a))
  a$Genus<-tax_a$Genus
}
a$comp1<-as.numeric(a$comp1)
a<-a[order(abs(a$comp1)),]
a$Genus <- factor(a$Genus, levels = a$Genus) # otherwise it would be re-ordered in plot
levels(a$Genus)<-gsub("_group", "", levels(a$Genus), fixed=T )
ggplot(mapping=aes(x=a$Genus,y=a$comp1)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, size = 10.5, hjust = 1, vjust = 1)) +
  labs(x="", 
       y="loadings",
       title = "Loadings of selected genera for component 1 in sPLSDA") +
  geom_hline(yintercept = max(a$comp1)/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, linewidth= 1) +
  geom_hline(yintercept = min(a$comp1)/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.2), "cm"))
ggsave(filename = "Results/PLS_DA/Loadings of choosen genera for sPLSDA comp 1.png", width = 8, height = 6, dpi=300)

b<-loadings[loadings[[2]]!=0,]
{tax_b<-Taxa.genus.update
  tax_b$ASV<- row.names(tax_b)
  tax_b<-subset(tax_b, ASV %in% row.names(b))
  tax_b<-tax_b[row.names(b),]
  identical(tax_b$ASV,row.names(b))
  b$Genus<-tax_b$Genus
}
n=2
b[,n]<-as.numeric(b[,n])
b<-b[order( abs(b[,n])) ,]
b$Genus <- factor(b$Genus, levels = b$Genus)
ggplot(mapping=aes(x=b$Genus,y=b[,n])) + geom_bar(stat = "identity", width = 0.6) + 
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, size = 12.5, hjust = 1, vjust = 1)) +
  labs(x="", 
       y="loadings", title = paste("Loadings of selected genera for component",n,"in sPLSDA")
       ) +
  #geom_hline(yintercept = max(b[,n])/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  geom_hline(yintercept = min(b[,n])/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = paste("Results/PLS_DA/Loadings_of_choosen_genera_for_sPLSDA_comp",n,".png"), width = 8, height = 6, dpi=300)

# c<-loadings[loadings[[3]]!=0,]
# {tax_c<-Taxa.genus.update
#   tax_c$ASV<- row.names(tax_c)
#   tax_c<-subset(tax_c, ASV %in% row.names(c))
#   tax_c<-tax_c[row.names(c),]
#   identical(tax_c$ASV,row.names(c))
#   c$Genus<-tax_c$Genus
#   c$Genus<-gsub("[Eubacterium]_", "Eubacterium_", c$Genus, fixed=T )
# }
# c[,3]<-as.numeric(c[,3])
# c<-c[order( abs(c[,3])) ,]
# c$Genus <- factor(c$Genus, levels = c$Genus)
# ggplot(mapping=aes(x=c$Genus,y=c[,3])) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 28, size = 10.5, hjust = 1, vjust = 1)) +
#   labs(x="", 
#        y="loadings", title = paste("Loadings of selected genera for component",3,"in sPLSDA")
#        )+
#   geom_hline(yintercept = max(c[,3])/2, colour="red", linetype="longdash") +
#   geom_hline(yintercept = 0, size= 1) +
#   geom_hline(yintercept = min(c[,3])/2, colour="red", linetype="longdash") +
#   theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
# ggsave(filename = paste("Results/PLS_DA/Loadings_of_choosen_genera_for_sPLSDA_comp",3,".png"), width = 8, height = 6, dpi=300)



############### testing the sPLS-DA model

X.test <- as.matrix(t(otu_table(data.train)))
Y.test <- as.factor(metadata_train$Condition2)
predict.splsda <- predict(final.splsda, newdata =as.matrix(X.test), dist = "all")

# evaluating the prediction accuracy
predict.comp <- predict.splsda$class$centroids.dist[,optimal.comp]
predict.comp
table(factor(predict.comp, levels = levels(Y)), as.factor(Y.test))
# now evaluating the prediction accuracy using ONLY the first component
predict.comp1 <- predict.splsda$class$centroids.dist[,1]
predict.comp1
table(factor(predict.comp1, levels = levels(Y)), as.factor(Y.test))




##### exporting settings and values of sPLSDA

suppressWarnings(rm(con))
con<-file("Results/PLS_DA/Settings_and_results_sPLSDA.txt")
sink(con)
cat("Number of training samples and their original number of genera", fill = TRUE)
cat(dim(X), fill = TRUE)
cat("\n Samples casually selected as test (and then discarded from training data set)", fill=TRUE)
cat(row.names(X.test), fill = TRUE)
cat("\n number of components suggested by perf function for normal PLSDA (centroid dist)", fill=TRUE)
cat(perf.splsda$choice.ncomp[3],fill=TRUE)
cat("\n number of components chosen for normal PLSDA (centroid dist)", fill=TRUE)
cat(comp, fill=TRUE)
cat("\n number of components suggested by tune.splsda function for sPLSDA", fill=TRUE)
cat(my.tuned.splsda$choice.ncomp$ncomp, fill=TRUE)
cat("\n number of components chosen for sPLSDA", fill=TRUE)
cat(optimal.comp, fill = TRUE)
cat("\n number of genera selected by tune.splsda function for each chosen component", fill=TRUE)
cat(optimal.keepX, fill=TRUE)
cat("\n Possible number of genera that could be selected by function tune.splsda", fill=TRUE)
cat(possible.pool, fill = TRUE)
cat("\n \n \n ### testing the prediction efficiency of sPLSDA through confusion matrix using ONLY the first component \n", fill=TRUE)
table(factor(predict.comp1, levels = levels(Y)), Y.test)
cat("\n \n ### testing the prediction efficiency of sPLSDA through confusion matrix component 1 and",n,"\n", fill=TRUE)
table(factor(predict.comp, levels = levels(Y)), Y.test)
cat("\n \n type of distance used: centroid distance")
sink()
close(con)


suppressWarnings( rm(a,b,con,tax_a,tax_b, colors_train, loadings, metadata_train, predict.comp,predict.comp1,train,X,Y,X.test,Y.test,predict.splsda,comp,n,test,optimal.comp,optimal.keepX) )
gc()




##################### R and PACKAGES VERSION #########################

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
print("Dendextend")
packageVersion("dendextend")
cat("\n", "\n", fill=TRUE)
package$otherPkgs$mixOmics[c(1,4)]
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
