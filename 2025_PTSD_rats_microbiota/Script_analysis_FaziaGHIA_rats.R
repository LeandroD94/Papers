##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("DESeq2")
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
  # dir.create("Results/Alpha_div")
  dir.create("Results/PCoA_and_BetaDivers")
  dir.create("Results/Differential_analyses")
}
options(scipen = 100) # disable scientific annotation



### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue4", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_vsearch.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
# sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

# NB: the order of the row in the Metadata is employed also to set the order of the samples in the plots
metadata <- as.data.frame(read.table("metadata_GHIA_rats.txt", header = T, sep="\t"))
row.names(metadata)<-metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])
original_order<-metadata$Sample_name # to maintein this sample order in plots
metadata<-metadata[sample, ]
identical(as.numeric(length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

# updating the names
sample_names(data) <- gsub( "GHIA_", "", sample_names(data) )

rm(original_length)

data<-subset_samples(data, Mouse_ID!="Rat1") # mislabeled (see preliminary PCoA)



########################### % ASSIGNED IN SILVA ##################################

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}
{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
}

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
  e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
  assigned<-rbind.data.frame(a,e)
  colnames(assigned)<-c("Total","Assigned","%","Taxa")
  assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
rm(a,e,assigned)




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
target_filter<- filter_taxa(data.genus.temp, function(x) mean(x)<= 0.005, TRUE)
filtered<-taxa_names( target_filter )
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ) , 
           file="Data_check/Filtered_genus_under_0005_mean_cutoff.csv")


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


# save.image("Fazia_16S_after_filters_RATS.RData")




############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

if(! "proof1" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
row.names(Original_number)<-Original_number$sample.ID
Original_number<- Original_number[as.character(sample_data(data)$FASTQ_ID), ]   # the files come from the whole sequencing batch
Original_number<-sum( c(Original_number$forward.sequence.count, Original_number$reverse.sequence.count))
Remaining_number<-sum(otu_table(data))


con<-file("Data_check/DETAILS_ABOUT_ORIGINAL_AND_PROCESSED_READS_NUMBER.txt")
sink(con, append=TRUE)
cat("Raw reads abundance (sequenced, considering the pairs as one)", fill=TRUE)
cat(Original_number/2)
cat("Reads abundance after every filter", fill=TRUE)
cat(Remaining_number)
cat("Percentual of remaining vs original", fill=TRUE)
cat(paste0(round((Remaining_number/(Original_number/2))*100,digits=2)),"%")
sink()
close(con)
rm(con)


### barplot of reads
Original_read_number<-as.data.frame(read.table(file="QIIME/Original_number_of_reads_for_sample.tsv", header = T, sep="\t", row.names = 1))
Original_read_number <- Original_read_number[as.character(sample_data(data)$FASTQ_ID), ]   # the file come from the whole sequencing batch
DADA2_read_number<-as.data.frame(read.table(file="QIIME/DADA2_Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
DADA2_read_number <- DADA2_read_number[as.character(sample_data(data)$FASTQ_ID), ]   # the file come from the whole sequencing batch
after_filter_number<-as.data.frame(colSums(otu_table(data)))
# updating names from FASTQ codes to sample names
row.names(Original_read_number)<-sample_data(data)$Sample_name # same order, already modified
row.names(DADA2_read_number)<-sample_data(data)$Sample_name # same order, already modified
row.names(after_filter_number)<-sample_data(data)$Sample_name # same order, already modified
# creating the table
table<-cbind.data.frame(Original=Original_read_number$forward.sequence.count, # merged --> just one column, not both of them
                        After_quality_filter=DADA2_read_number$non.chimeric,
                        After_contam_filter=after_filter_number[,1])
table$FASTQ<-as.character(sample_data(data)$FASTQ_ID)
# re-naming and adding groups
if(identical( table$FASTQ , as.character(sample_data(data)$FASTQ_ID) ) ){ # same order
  #table$Samples<-as.factor( gsub("GHIA_","", sample_data(data)$FASTQ_ID) )
  table$Samples<-as.factor( sample_data(data)$Sample_name)
  table$factor <- sample_data(data)$Sample_Type
  # table$factor <- as.factor(table$factor)
  cat("\n\nOK!\n\n")
}
# plotting
ggplot(aes(x=Samples, y=Original), data=table) +
  facet_grid( ~ factor, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Read quality filters"), alpha= 0.8,
            width = 0.65, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
            width= 0.35, stat = "identity") +
  theme_classic( base_size = 7.75) +
  theme(axis.text.x = element_text(size=5, angle = 45, vjust=1, hjust=1),
        axis.text.y = element_text(size=5),
        axis.title.y = element_text(size=7),
        panel.grid.major.y = element_line(linewidth =0.1, color="grey"),
        legend.position = "bottom",
        plot.margin = margin(2, 1, 2, 1),
        legend.margin = margin(0, 30, 0, 0),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9.5),
        legend.key.height = unit(0.4,"cm")
  ) +
  scale_fill_manual(name='Number of reads:  ',
                    breaks=c('Original number', 'Read quality filters', 'Relative abundance filters'),
                    values=c('Original number'='green3', 'Read quality filters'='coral', 'Relative abundance filters'='red3')) +
  #scale_y_continuous(breaks = c(10000, seq(0, max(table$Original+1000), 2500))) +  # non funziona perchè manca il valore del vecchio AGS, per ora
  labs(y="Reads abundance", x="FASTQ")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 6.5, height = 4.2, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)
dev.off()

suppressWarnings( rm(table, table2, Original_read_number, DADA2_read_number, after_filter_number) )




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
  #legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data.genus),"matrix")), step=100,label=F, 
             ylab = "Genera", xlab= "Reads amount")
evalslopes(r,sample_names(data.genus),lim=0.001,cex=0.7)
dev.off()
rm(r)

# data.temp <-subset_samples(data, Sample_Type == "Colon" )
# png(file="Data_check/Rarefaction_curve_focus_on_COLON_samples.png",width=3000,height=2100, res=300)
# r<-rarecurve(t(as(otu_table(data.temp),"matrix")), step=100,label=F, ylab = "ASVs", xlab= "Reads amount")
# evalslopes(r,sample_names(data.temp),lim=0.001,cex=0.7)
# dev.off()
# rm(r, data.temp)

# NB: the unsaturated samples is labelled as "faeces" in the FASTQ code... but this is a library preparation error!

# {
# data<-subset_samples(data, Sample_name!="M_F_C_1_WT")
# metadata<-metadata[metadata$Sample_name!="M_F_C_1_WT", ]
# proof2<-"The unsaturated sample has been discarded"
# }




############### \\\\\\ PREPARING THE ORIGINAL (COMPLETE) OBJECTS WITH COLON AND FEACES \\\\\\ #######################

if( "proof1" %in% ls() & length(unique(sample_data(data)$Sample_Type)) >1 )  { # then the objects are still complete
  # NB: the filter step is performed on each object !
  data_complete<-prune_taxa(taxa_sums(data)>0,data)
  metadata_complete<-metadata
  rm(data, metadata)
} 

data_complete.prop <- transform_sample_counts(data_complete, function(ASV) ASV/sum(ASV)*100)
data.genus_complete<-tax_glom(data_complete, taxrank ="Genus", NArm=F)
data.genus_complete.prop<-transform_sample_counts(data.genus_complete, function(x) (x/sum(x))*100)
data.phy_complete<-tax_glom(data_complete, taxrank ="Phylum", NArm=F)
data.phy_complete.prop<-transform_sample_counts(data.phy_complete, function(x) (x/sum(x))*100)




################# GENERAL PCoA BEFORE SUBSETTING (EVERY SAMPLE) ########################

suppressWarnings(rm(data.prop_temp, data.sqrt_prop))

if( ! "proof1" %in% ls()  ){
  stop("\n Wait! Did you perform the filtering steps??? \n\n")
  Sys.sleep(2)
}

data.prop.labels<-data.genus_complete.prop

{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}

sample_names(data.sqrt_prop)<-as.factor(sample_names(data.sqrt_prop))

# Type 1
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type", shape= "Sample_Group") +
  # scale_color_manual(values= c() ) +
  geom_point(size=2.8) +
  geom_point(size=4.8, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_line(aes(group= Mouse_ID),col="grey", linewidth=0.15)+
  # geom_text(aes(label= sample_names(data.prop.labels)), 
  #           color="black", size=1.8,
  #           show.legend = FALSE) +
  labs(title="",
       color="Sample Type",
       shape="Sample Group",
       caption = "The lines connect the same samples from the same mouse",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_and_BetaDivers/EVERY_SAMPLES_PCoA_Helling_on_genera_LINES.png", width = 5.2, height = 4.5, dpi=300)

# Type 2
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type", shape= "Sample_Group") +
  # scale_color_manual(values= c() ) +
  geom_point(size=2.8) +
  geom_point(size=4.8, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_text(aes(label= sample_data(data.prop.labels)$Sample_name), 
            color="black", size=1.5,
            show.legend = FALSE) +
  labs(title="",
       shape="Sample Group",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/PCoA_and_BetaDivers/EVERY_SAMPLE_PCoA_Helling_on_genera_Names.png", width = 5.2, height = 4.5, dpi=300)




####################### GENERAL COUNTS EXPORT ##############################

#write.csv2(cbind(as(otu_table(data_complete),"matrix"),as(tax_table(data_complete),"matrix")),file="Results/Abundances/BOTH_COLON_AND_FEACES_filtered_ASV.csv",quote=F)
options(scipen = 100)
#write.csv2(cbind(as(otu_table(data.genus_complete.prop),"matrix"),as(tax_table(data.genus_complete),"matrix")),file="Results/Abundances/EVERY_SAMPLE_COUNTS_genera_proportions_BOTH_COLON_AND_FEACES.csv",quote=F)
write.xlsx(cbind(as(otu_table(data.genus_complete.prop),"matrix"),as(tax_table(data.genus_complete),"matrix")),file="Results/Abundances/EVERY_SAMPLE_COUNTS_genera_proportions_BOTH_COLON_AND_FEACES.xlsx",showNA = F, col.names = T)

# raw feature_table
write.csv2(cbind(as(otu_table(unfiltered_data),"matrix"),as(tax_table(unfiltered_data),"matrix")),file="Data_check/raw_feature_table_BOTH_COLON_AND_FEACES.csv",quote=F)




################## ABUNDANCES BAR PLOTS (EVERY SAMPLE) ##################

### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- data.phy_complete.prop
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:5]
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
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Mouse_ID, y=Abundance, fill=Phylum) ) +
  facet_nested( . ~ Sample_Type + Sample_Group , scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =11) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust = 1, size= 6.5),
        axis.text.y=element_text( size= 7 ), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.45, "cm"),
        legend.text = element_text ( size = 10 ),
        legend.position="bottom", 
        legend.margin = margin(-9,25,10,0)
  ) + 
  scale_x_discrete(expand = c(0,0.5) ) +
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="", y="Percent abundance",
       title = "Most abundant identified phyla",
       fill="",
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/EVERY_SAMPLE_TOP_Phyla.png",width=7.2,height=5,dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/COLON_TOP_phyla_averages_EverySample.xlsx", row.names = F, to_save)

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)



### TOP Genera
data.target<-data.genus_complete.prop
# adding informations to missing names
taxa_temp<- as.data.frame( tax_table(data.target) )
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
  Taxa.target.update<-taxa_temp
}
rm(taxa_temp)

{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:15]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.target.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

ggplot(data=tabella, aes(x=Mouse_ID, y=Abundance, fill=Genus)) +
  facet_nested( . ~ Sample_Type + Sample_Group , scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =11) +
  scale_fill_manual(values=fill_color_15) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust = 1, size= 6.5),
        axis.text.y=element_text( size= 7 ), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.21, "cm"),
        legend.text = element_text ( size = 9.2 ),
        legend.position="bottom", 
        legend.margin = margin(-10,10,0,0)
  ) + 
  scale_x_discrete(expand = c(0,0.5) ) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="", y="Percent abundance",
       title = "Most abundant genera",
       fill= "" ,
       caption = " 'Others' includes every genus below rank 15 ")
ggsave(file="Results/Abundances/EVERY_SAMPLE_TOP_Genera.png",width=7.2,height=5,dpi=300)
dev.off()


write.xlsx(file = "Results/Abundances/EVERY_SAMPLE_TOP_genera.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)




#################### \\\\\\ STARTING THE ANALYSIS OF COLON \\\\\ #######################

# dir.create("Results/Abundances/")
data<-subset_samples(data_complete, Sample_Type=="Colon")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa

### everything ready, let's start with the analysis
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  # data.class = tax_glom(data, taxrank = "Class", NArm = F)
  # data.order = tax_glom(data, taxrank = "Order", NArm = F)
  # data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  # data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  # data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  # data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  # Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  # Taxa.class<-as.data.frame(tax_table(data.class))
  # Taxa.order<-as.data.frame(tax_table(data.order))
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

metadata<-metadata_complete[ metadata_complete$Sample_Type=="Colon",  ]




########################### COLON COUNTS EXPORT ##########################################

# {
#   #write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/COLON_ASV_counts.csv",quote=F)
#   write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/COLON_counts_phylum.csv",quote=F)
#   write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/COLON_counts_genus.csv",quote=F)
# }
# 
# {
#   #write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/COLON_counts_phylum.csv",quote=F)
#   write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/COLON_counts_genus.csv",quote=F)
#   write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/COLON_counts_genus.xlsx",showNA = F, col.names = T)
#   write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/COLON_counts_phylum.xlsx",showNA = F, col.names = T)
# }




###################### COLON ABUNDANCES BAR PLOT ##########################


### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.phy.prop)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:5]
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
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Mouse_ID, y=Abundance, fill=Phylum) ) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_nested(~ Sample_Group , scales = "free_x", space = "free_x" ) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_5) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=38,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=4),
  axis.title =element_text(size=10),
  strip.text = element_text(size=7.5),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin( -1.5 ,5, 8 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant identified phyla (Colon)",
       fill="",
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/COLON_TOP_phyla.png",width=6.5,height=4.5, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/COLON_TOP_phyla_averages_EverySample.xlsx", row.names = F, to_save)

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)



### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:15]
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
  tabella$Genus<- gsub("_group","",tabella$Genus)
  tabella$Genus<- gsub("uncultured_ f Lachnospiraceae","uncultured_ f Lachnospir",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Mouse_ID, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_nested(~ Sample_Group , scales = "free_x", space = "free_x" ) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_15) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=38,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=4),
  axis.title =element_text(size=10),
  strip.text = element_text(size=7.5),
  legend.key.height = unit(0.38, "cm"),
  legend.key.width = unit(0.21, "cm"),
  legend.text = element_text ( size = 8.25 ),
  legend.position="bottom",
  legend.margin = margin(-11,30,1,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant identified genera (Colon)",
       fill="",
       caption = " 'Others' includes every genus below rank 15 ")
ggsave(file="Results/Abundances/COLON_TOP_genera.png",width=6.5,height=4.5, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/COLON_TOP_genera_Averages.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




######################## COLON ALPHA DIVERSITY ##############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

data.label <- data.genus
pAlpha<-suppressWarnings( plot_richness(data.label, measures=c("Shannon", "Observed"), x="Sample_Group", color="Sample_Group") )
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H <- pAlpha$data[pAlpha$data$variable=="Shannon", ]
obs <- pAlpha$data[pAlpha$data$variable=="Observed", ]
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
  geom_boxplot(data=pAlpha$data, aes(x=Sample_Group, y=value, color=NULL, fill=Sample_Group), alpha=0.15, linewidth=0.15 ) +
  theme_classic2(base_size = 10) +
  # geom_text(aes(label=Sample_name), color="black", size=2.05, show.legend = FALSE) +
  labs(title="Alpha diversity (Colon)") +
  guides(fill="none", color="none", shape="none") +
  stat_compare_means(aes(group = Sample_Group), label="p.format", method = "wilcox.test", 
                     label.x= 1.45, size=3.2, label.y.npc = "top", vjust=-0.42, hjust=0.4) +
  theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=6.5),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25)
  )
ggsave(file="Results/COLON_Alpha_diversity_GENUS.png", width = 5.2,height =4, dpi=300)


suppressWarnings(rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor) )




####################### COLON PCoA and BETA DIV ##########################

#### PERMANOVA subsets

# C vs Shocked
suppressWarnings(rm(ASV.genus.prop))
sqrt_prop_table<-as.data.frame(t(sqrt(as.data.frame(otu_table(data)))))
perm_g<- vegan::adonis2(sqrt_prop_table ~Sample_Group, 
                        data=as(sample_data(data),"data.frame"),
                        permutations = 9999, method="euclidean")
perm_g_CvsS<-perm_g$`Pr(>F)`[1]
perm_g_CvsS


##### PCoA on genera
data.prop.labels<-data.genus.prop

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# Control vs Shocked
plot_ordination(data.sqrt_prop, ordBC, color="Sample_Group") + # , color = "Sample_Type") +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  stat_ellipse( linewidth= 0.25 ) +
  # guides(color="none", shape="none") +
  theme_classic(base_size = 12.5) +
  theme(title=element_text(size=10)) +
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances (Colon)",
       subtitle = paste("Control vs Shocked Pr(>F) =",perm_g_CvsS),
       x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       color="Group")
ggsave(file="Results/PCoA_and_BetaDivers/COLON_Beta_divers_Hellinger_on_genera_Control_vs_Shocked.png", width = 4.6, height = 4, dpi=300)




###################### COLON DA WITH DESEQ2  ########################

if(! "unfiltered_data" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS OF WT vs KO (Controls)
suppressWarnings(rm(data_pruned, data_target, data.genus_pruned))
data_target<- data
data_pruned<- prune_taxa(taxa_sums(data_target)  > 10, data_target)
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
  DEseq_data<-phyloseq_to_deseq2(d, ~ Sample_Group)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sample_Group", "Control", "Shocked"))
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_BBBB_vs_AAAA.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Differential_analyses/THERE_ARE_NOT_Signif_Diff_Abund_for_Control_vs_Shocked_in_COLON.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/Differential_analyses/COLON_THERE_ARE_NOT_Signif_Diff_Abundances_Control_vs_Shocked.xlsx", showNA = F, col.names = T)




#################### \\\\\\ STARTING THE ANALYSIS OF FAECES \\\\\ #######################

# dir.create("Results/Abundances/")
data<-subset_samples(data_complete, Sample_Type=="Faeces")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa

### everything ready, let's start with the analysis
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  # data.class = tax_glom(data, taxrank = "Class", NArm = F)
  # data.order = tax_glom(data, taxrank = "Order", NArm = F)
  # data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  # data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  # data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  # data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  # Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  # Taxa.class<-as.data.frame(tax_table(data.class))
  # Taxa.order<-as.data.frame(tax_table(data.order))
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

metadata<-metadata_complete[ metadata_complete$Sample_Type=="Faeces",  ]




########################### FAECES COUNTS EXPORT ##########################################

# {
#   #write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/FAECES_ASV_counts.csv",quote=F)
#   write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/FAECES_counts_phylum.csv",quote=F)
#   write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/FAECES_counts_genus.csv",quote=F)
# }
# 
# {
#   #write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/FAECES_counts_phylum.csv",quote=F)
#   write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/FAECES_counts_genus.csv",quote=F)
#   write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/FAECES_counts_genus.xlsx",showNA = F, col.names = T)
#   write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/FAECES_counts_phylum.xlsx",showNA = F, col.names = T)
# }




###################### FAECES ABUNDANCES BAR PLOT ##########################


### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.phy.prop)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:5]
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
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Mouse_ID, y=Abundance, fill=Phylum) ) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_nested(~ Sample_Group , scales = "free_x", space = "free_x" ) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_5) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=38,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=4),
  axis.title =element_text(size=10),
  strip.text = element_text(size=7.5),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin( -2 ,5, 8 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant identified phyla (Faeces)",
       fill="",
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/FAECES_TOP_phyla.png",width=6.5,height=4.5, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/FAECES_TOP_phyla_averages_EverySample.xlsx", row.names = F, to_save)

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)



### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:15]
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
  tabella$Genus<- gsub("_group","",tabella$Genus)
  tabella$Genus<- gsub("uncultured_ f Lachnospiraceae","uncultured_ f Lachnospir",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Mouse_ID, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_nested(~ Sample_Group , scales = "free_x", space = "free_x" ) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_15) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=38,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=4),
  axis.title =element_text(size=10),
  strip.text = element_text(size=7.5),
  legend.key.height = unit(0.38, "cm"),
  legend.key.width = unit(0.21, "cm"),
  legend.text = element_text ( size = 8.25 ),
  legend.position="bottom",
  legend.margin = margin(-11,30,1,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant identified genera (Faeces)",
       fill="",
       caption = " 'Others' includes every genus below rank 15 ")
ggsave(file="Results/Abundances/FAECES_TOP_genera.png",width=6.5,height=4.5, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/FAECES_TOP_genera_Averages.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




######################## FAECES ALPHA DIVERSITY ##############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

data.label <- data.genus
pAlpha<-suppressWarnings( plot_richness(data.label, measures=c("Shannon", "Observed"), x="Sample_Group", color="Sample_Group") )
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H <- pAlpha$data[pAlpha$data$variable=="Shannon", ]
obs <- pAlpha$data[pAlpha$data$variable=="Observed", ]
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
  geom_boxplot(data=pAlpha$data, aes(x=Sample_Group, y=value, color=NULL, fill=Sample_Group), alpha=0.15, linewidth=0.15 ) +
  theme_classic2(base_size = 10) +
  # geom_text(aes(label=Sample_name), color="black", size=2.05, show.legend = FALSE) +
  labs(title="Alpha diversity (Faeces)") +
  guides(fill="none", color="none", shape="none") +
  stat_compare_means(aes(group = Sample_Group), label="p.format", method = "wilcox.test", 
                     label.x= 1.45, size=3.2, label.y.npc = "top", vjust=-0.42, hjust=0.4) +
  theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=6.5),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25)
  )
ggsave(file="Results/FAECES_Alpha_diversity_GENUS.png", width = 5.2,height =4, dpi=300)


suppressWarnings(rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor) )




####################### FAECES PCoA and BETA DIV ##########################

#### PERMANOVA subsets

# C vs Shocked
suppressWarnings(rm(ASV.genus.prop))
sqrt_prop_table<-as.data.frame(t(sqrt(as.data.frame(otu_table(data)))))
perm_g<- vegan::adonis2(sqrt_prop_table ~Sample_Group, 
                        data=as(sample_data(data),"data.frame"),
                        permutations = 9999, method="euclidean")
perm_g_CvsS<-perm_g$`Pr(>F)`[1]
perm_g_CvsS


##### PCoA on genera
data.prop.labels<-data.genus.prop

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# Control vs Shocked
plot_ordination(data.sqrt_prop, ordBC, color="Sample_Group") + # , color = "Sample_Type") +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  stat_ellipse( linewidth= 0.25 ) +
  # guides(color="none", shape="none") +
  theme_classic(base_size = 12.5) +
  theme(title=element_text(size=10)) +
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances (Faeces)",
       subtitle = paste("Control vs Shocked Pr(>F) =",perm_g_CvsS),
       x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       color="Group")
ggsave(file="Results/PCoA_and_BetaDivers/FAECES_Beta_divers_Hellinger_on_genera_Control_vs_Shocked.png", width = 4.6, height = 4, dpi=300)




###################### FAECES DA WITH DESEQ2  ########################

if(! "unfiltered_data" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS OF WT vs KO (Controls)
suppressWarnings(rm(data_pruned, data_target, data.genus_pruned))
data_target<- data
data_pruned<- prune_taxa(taxa_sums(data_target)  > 10, data_target)
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
  DEseq_data<-phyloseq_to_deseq2(d, ~ Sample_Group)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sample_Group", "Control", "Shocked"))
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_BBBB_vs_AAAA.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Differential_analyses/THERE_ARE_NOT_Signif_Diff_Abund_for_Control_vs_Shocked_in_FAECES.csv", row.names = F)




##################### R AND PACKAGES VERSION #########################

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
package$otherPkgs$DESeq2 [1:2]
cat("\n", "\n", fill=TRUE)
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))


