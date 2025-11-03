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
  library("reshape")
  library("qiime2R")
}


{dir.create("Data_check")
  dir.create("Results")
  dir.create("Results/Abundances")
  #dir.create("Results/Time_series_in_R1")
  dir.create("Results/DA_DESeq2")
}
options(scipen = 100) # disable scientific annotation


### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
fill_color_19_modified<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "wheat3","deepskyblue2","red","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza", tree = "QIIME/rooted-tree.qza")

# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update
Metadata <- as.data.frame(read.csv2("metadata.csv", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)]) #- 1 # there is a sample (F3) collected but NOT sequenced
original_order<-Metadata$Sample_name # to maintein this sample order in plots
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

# adding a space to some day (technical replicate) to differentiate the labels in the plot
Metadata$Experiment_day[Metadata$Sample_name=="P2"]<-paste0(Metadata$Experiment_day[Metadata$Sample_name=="P2"]," ")
Metadata$Experiment_day[Metadata$Sample_name=="P7"]<-paste0(Metadata$Experiment_day[Metadata$Sample_name=="P7"]," ")

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)

# data<-subset_samples(data, sample_names!="S3")
# Metadata<-Metadata[Metadata$Sample_name!="S3", ]

sample_names(data)<-sample_data(data)$Sample_name
sample_data(data)$Sample_name<-factor(sample_data(data)$Sample_name,
                                      levels = original_order)
sample_data(data)$Experiment_day<-factor(sample_data(data)$Experiment_day,
                                         levels = unique(sample_data(data)$Experiment_day))



table_PHA<-read.table(file="PHA_Accumulators_maggio2025.tsv", fill = T, header = T, sep="\t")





#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_max_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE)))[["Genus"]]
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

if(! "proof1" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
Original_number$sample.ID<-gsub("^.*F","",Original_number$sample.ID)
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
row.names(Original_read_number)<-gsub("^.*F","",row.names(Original_read_number))
Original_read_number <- Original_read_number[as.character(sample_data(data)$FASTQ_ID), ]   # the file come from the whole sequencing batch
DADA2_read_number<-as.data.frame(read.table(file="QIIME/DADA2_Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
row.names(DADA2_read_number)<-gsub("^.*F","",row.names(DADA2_read_number))
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
  table$Samples<-as.factor(sample_data(data)$Sample_name)
  table$Experiment_day<-sample_data(data)$Experiment_day
  table$factor <- sample_data(data)$Sample_Type
  # table$factor <- as.factor(table$factor)
  cat("\n\nOK!\n\n")
}
# plotting
ggplot(aes(x=Experiment_day, y=Original), data=table) +
  facet_grid( ~ factor, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Read quality filters"), alpha= 0.8,
            width = 0.65, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
            width= 0.35, stat = "identity") +
  theme_classic( base_size = 7.75) +
  theme(axis.text.x = element_text(size=6, angle = 40, vjust=1, hjust=1),
        axis.text.y = element_text(size=5),
        axis.title.y = element_text(size=7),
        panel.grid.major.y = element_line(linewidth =0.1, color="grey"),
        legend.position = "bottom",
        strip.text = element_text( size= 10 ),
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
  labs(y="Reads abundance", x="Experimental Day")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 5.5, height = 4, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)


suppressWarnings( rm(table, table2, Original_read_number, DADA2_read_number, after_filter_number) )




############################ RAREFACTION ANALYSIS ################################

if(! "proof2" %in% ls()  ){
  data_with_unsat<-data
  data_with_unsat.genus = tax_glom(data_with_unsat, taxrank = "Genus", NArm = F)
}


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


png(file="Data_check/Rarefaction_curve.png",width=2100,height=1650, res=300)
r<-rarecurve(t(as(otu_table(data_with_unsat.genus),"matrix")), step=100,label=F, ylab = "Genera", xlab= "Reads amount")
evalslopes(r,sample_names(data_with_unsat.genus),lim=0.001,cex=1.05)
dev.off()
rm(r)


proof2<-"The sample S13 is strongly unsaturated"
data<-subset_samples(data_with_unsat, ! Sample_name %in% c("S13","F13"))   # removing also the paired sample





########################### % ASSIGNED IN SILVA ##################################

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
# data.class = tax_glom(data, taxrank = "Class", NArm = F)
# data.order = tax_glom(data, taxrank = "Order", NArm = F)
# data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}
{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  # Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  # Taxa.class<-as.data.frame(tax_table(data.class))
  # Taxa.order<-as.data.frame(tax_table(data.order))
}

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
  # b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
  # c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
  # d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
  e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
  # assigned<-rbind.data.frame(a,b,c,d,e)
  assigned<-rbind.data.frame(a,e)
  colnames(assigned)<-c("Total","Assigned","%","Taxa")
  assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
suppressWarnings (rm(a,b,c,d,e,assigned))




#######################  !!!!!  STARTING THE ANALYSIS OF R1 vs R2 !!!!!   #############################

if("proof1" %in% ls() & "proof2" %in% ls() ){
to_remove<-ls()
to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete","data","Metadata",
                                        "colors","fill_color_5","fill_color_10", "fill_color_19", "fill_color_19_modified",
                                        "table_PHA", "unfiltered_data", "data_with_unsat")]
rm(list=to_remove)
}

if( ! "data_complete" %in% ls() & "proof1" %in% ls() & "proof2" %in% ls()  ){
  data_complete<-data
}
data_complete.phy = tax_glom(data_complete, taxrank = "Phylum", NArm = F)
data_complete.genus = tax_glom(data_complete, taxrank = "Genus", NArm = F)
data_complete.genus.prop <- transform_sample_counts(data_complete.genus, function(x) x/sum(x)*100)
data_complete.phy.prop <- transform_sample_counts(data_complete.phy, function(x) x/sum(x)*100)

data<-subset_samples(data_complete, Experiment_state %in% c("Started","Almost_stabilized") )  # for comparisons of F/F ratios
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
 data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
 data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
 data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
 data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
 Taxa.genus<-as.data.frame(tax_table(data.genus))
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
# also with the "complete" data set
taxa_temp<- as.data.frame(tax_table(data_complete.genus))
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
  Taxa.complete.genus.update<-taxa_temp
}

rm(taxa_temp)




########################### COUNTS EXPORT ##########################################

{
  write.csv2(cbind(as(otu_table(data_with_unsat),"matrix"),as(tax_table(data_with_unsat),"matrix")),file="Results/Abundances/Counts_raw_ASVs.csv",quote=F)
  # with unsat --> feauture table to upload in GEO
  write.csv2(cbind(as(otu_table(data_complete.phy.prop),"matrix"),as(tax_table(data_complete.phy.prop),"matrix")),file="Results/Abundances/Counts_phyla_prop.csv",quote=F)
  write.csv2(cbind(as(otu_table(data_complete.genus.prop),"matrix"),as(tax_table(data_complete.genus.prop),"matrix")),file="Results/Abundances/Counts_genera_prop.csv",quote=F)
}




######################## ABUNDANCES BAR PLOT ##########################


### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data_complete.phy.prop)
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
  tabella$Phylum<-gsub ("Candidatus_","Ca. ", tabella$Phylum)
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}

tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
# tabella$Sampling_date<-gsub("_22","",tabella$Sampling_date)
unique(tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
tabella$Sampling_date<-gsub("_","/",tabella$Sampling_date, fixed=T)
# order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
order_time<-tabella$Sampling_date[c(grep("/04",tabella$Sampling_date),grep("/05",tabella$Sampling_date),grep("/06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))

# ggplot(data=tabella, aes(x=Abundance, y=Experiment_day, fill=Phylum)) +
#   theme_bw(base_size =10) + 
#   #facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
#   facet_grid2( Experiment_day+Sample_Type  ~ .,
#                scales = "free", 
#                space="free",
#                switch = "x", # This places the grid at the bottom (as x axis)
#                # switch = "y",
#                strip = strip_nested(size="constant"))+
#   scale_y_discrete (expand = c(0,0) ) +
#   scale_x_continuous(expand = c(0,0)) +  #  bars closer to the axis
#   scale_fill_manual(values=fill_color_10) +
#   geom_bar(stat="identity", position="stack", width = 1) +
#   theme( strip.text.y = element_text(size=6.85, angle= 0),
#          panel.spacing.y = unit(1.5, "pt"),
#          panel.background = element_rect(fill="white"),
#         axis.title.x = element_text(size=10, vjust=-0.95),
#         axis.text.x = element_text(size=8.5),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y = element_blank(),
#         plot.title = element_blank(),
#         #axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
#         legend.key.height = unit(0.2, "cm"),
#         legend.key.width = unit(0.45, "cm"),
#         legend.text = element_text ( size = 12.5 ),
#         legend.position="bottom",
#         plot.margin = margin(2,0.5,2,-4.5),
#         legend.margin = margin(0,0,11,-9)
#         ) +
#   guides(fill=guide_legend(nrow=3)) +
#   labs(y="Experiment day / Reactor", x="Percent abundance",
#        fill="",
#        caption = " 'Others' includes every phylum below rank 10 "
#        )
# ggsave(file="Results/Abundances/TOP_phyla_abundances.png",height=6,width=4, dpi=300) 
# dev.off()

# means of TOP phyla
mean_every<-round( as.numeric(apply(otu_table(prune.dat_top),1,mean)) ,2)
mean_R1_after_differentiation<-round( as.numeric(apply(otu_table(subset_samples(prune.dat_top, Experiment_state!="Still_preparing" & Sample_Type=="R1")),1,mean)) ,2 )
mean_R2_after_differentiation<-round( as.numeric(apply(otu_table(subset_samples(prune.dat_top, Experiment_state!="Still_preparing" & Sample_Type=="R2")),1,mean)) ,2 )
to_save<- cbind.data.frame("Mean_every_sample"=mean_every,
                           "mean_R1_after_differentiation"=mean_R1_after_differentiation,
                           "mean_R2_after_differentiation"=mean_R2_after_differentiation,
                           "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$Mean_every_sample, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_phyla_average_abundance.xlsx", row.names = F, to_save )




# TOP Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data_complete.genus.prop), decreasing=TRUE))[1:19]
prune.dat_top <- prune_taxa(top,data_complete.genus.prop)
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_selected<-Taxa.complete.genus.update[row.names(tax_selected),]
tax_table(prune.dat_top)<-as.matrix(tax_selected)
others<-taxa_names(data_complete.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data_complete.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-gsub("uncultured_","uncult. ",tabella$Genus)
tabella$Genus<-gsub("cteraceae","c.",tabella$Genus)
tabella$Genus<-gsub("Pseudopedobacter","Pseudopedobac.",tabella$Genus)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
tabella$Sampling_date<-gsub("_","/",tabella$Sampling_date, fixed=T)
# order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
order_time<-tabella$Sampling_date[c(grep("/04",tabella$Sampling_date),grep("/05",tabella$Sampling_date),grep("/06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
fill_color_modified<-c("brown4","springgreen2","chocolate4","cyan2",
                       "gray15","coral3","yellow4","darkmagenta",
                       "pink3","orange","gray","blue","gold",
                       "darkgreen","violet", "deepskyblue3","wheat3",
                       "red","chartreuse3","darkslategray3")
tabella <- aggregate( Abundance  ~ Genus + Sample_Type + Experiment_day, tabella , FUN=sum ) 
                        
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Genus)) + 
  theme_classic(base_size =12) + 
  facet_grid2( Sample_Type ~ Experiment_day,
               scales = "free", 
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0,0) , breaks = c(20,40,60,80,100)) +  #  bars closer to the x axis
  scale_fill_manual(values=fill_color_modified) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme( strip.text = element_text(size=10),
         plot.caption = element_text(size=7.5),
         panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
         panel.spacing.x = unit(0.75,"pt"),
         strip.clip = "off",
         axis.title.x = element_text(size=11),
         axis.title.y = element_text(size=11, vjust=0),
         axis.text.y = element_text(size=7.25, angle = 0),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.line.x = element_blank(),
         plot.title = element_blank(),
         plot.margin = margin(5,1,2,0),
         legend.key.height  = unit(0.4, "cm"),
         legend.key.width = unit(0.2, "cm"),
         legend.text = element_text ( size = 9.68 ),
         legend.position="bottom",
         legend.margin = margin(0,0,1,-14.5)
         ) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percent abundance",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_abundances.png",height =6,width = 4.95, dpi=300) 
dev.off()


# means of TOP genera
mean_every<-round( as.numeric(apply(otu_table(prune.dat_top),1,mean)) ,2)
mean_R1_after_differentiation<-round( as.numeric(apply(otu_table(subset_samples(prune.dat_top, Experiment_state!="Still_preparing" & Sample_Type=="R1")),1,mean)) ,2 )
mean_R2_after_differentiation<-round( as.numeric(apply(otu_table(subset_samples(prune.dat_top, Experiment_state!="Still_preparing" & Sample_Type=="R2")),1,mean)) ,2 )
to_save<- cbind.data.frame("Mean_every_sample"=mean_every,
                           "mean_R1_after_differentiation"=mean_R1_after_differentiation,
                           "mean_R2_after_differentiation"=mean_R2_after_differentiation,
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$Mean_every_sample, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_genera_Average_abundances.xlsx", row.names = F, to_save)


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top, to_save))




####################### ABUNDANCES OF PHA BACTERIA _ R1 versus R2 ###################

lista_PHA<-table_PHA[,1]     # they derive from papers

prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))

tabella$Experiment_state<-gsub("_"," ",tabella$Experiment_state)
# factoring also the sample names to ensure a certain order in the plot (=time) ; in the sample_names they are already ordered
tabella$Sample<-factor(tabella$Sample, levels=sample_names(data))
tabella$Genus<-gsub("Candidatus_","",tabella$Genus)

# fill_color_14<-c("bisque2","deeppink","darkgreen","darkmagenta","cyan","yellow2","brown","firebrick3","wheat3","violet","darkslategray3","springgreen3","blue","deepskyblue2") # , "gray", "pink3","yellow4","red","springgreen3","coral") 
fill_color_PHA<-c("chocolate","grey92","lightblue4","gold",
                  "lightcoral","lightgreen","grey20","deeppink",
                  "cyan","coral3","wheat4","darkmagenta","brown4",
                  "orange","yellow3",  "bisque3","blue",
                  "green3","darkslategray3", "chartreuse3","yellow4","red",
                  "magenta","deepskyblue2","gray","violet") 

ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=Genus)) +
  theme_classic( base_size =12 ) +
  facet_grid2( Sample_Type ~ Experiment_day,
               scales = "free", 
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0,0) , breaks = c(20,40,60,80,100), lim=c(0,100)) +  #  bars closer to the x axis
  scale_fill_manual(values=fill_color_PHA) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme( strip.text = element_text(size=11),
         plot.caption = element_text(size=7.5),
         panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
         panel.spacing.x = unit(0.75,"pt"),
         strip.clip = "off",
         axis.title.x = element_text(size=11),
         axis.title.y = element_text(size=11, vjust=0),
         axis.text.y = element_text(size=7.25, angle = 0),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.line.x = element_blank(),
         plot.title = element_blank(),
         plot.margin = margin(5,1,2,0),
         legend.key.height  = unit(0.4, "cm"),
         legend.key.width = unit(0.22, "cm"),
         legend.text = element_text ( size = 9.25 ),
         legend.position="bottom",
         legend.margin = margin(0,0,1,-14.5)
  ) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = paste("PHA accumulating genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"))
ggsave(file="Results/Abundances/PHA_accumulators_along_the_time.png",width=5.8,height=4.9,dpi=300)
dev.off()




######################## ALPHA DIVERSITY ##############################

data.genus.temp<-data_complete.genus
sample_names(data.genus.temp)<-paste(sample_data(data.genus.temp)$Sample_Type , sample_data(data.genus.temp)$Experiment_day , sep="_")
sample_names(data.genus.temp)<-gsub("R1_","S",sample_names(data.genus.temp))
sample_names(data.genus.temp)<-gsub("R2_","F",sample_names(data.genus.temp))

# values computation
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name" )
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
# aesthetics
mix_alpha$ID <- gsub("[1-9].*","",mix_alpha$Sample_name) # to connect the reactor points (without other reactors)
mix_alpha$ID_LineA<-mix_alpha$ID
mix_alpha$ID_LineA[ mix_alpha$ID =="F" ] <- seq(1,length(mix_alpha$ID_LineA[ mix_alpha$ID =="F" ])) # different labels --> diff groups --> these points will not be connected by this line
mix_alpha$ID_LineB<-mix_alpha$ID
mix_alpha$ID_LineB[ mix_alpha$ID =="S" ] <- seq(1,length(mix_alpha$ID_LineB[ mix_alpha$ID =="S" ])) # different labels --> diff groups --> these points will not be connected by this line
mix_alpha$ID<-gsub("S","1/1 ", mix_alpha$ID ) # NB: the space is added on purpose
mix_alpha$ID<-gsub("F","2/1", mix_alpha$ID )
mix_alpha$Time_color<- as.numeric(gsub("[A-Z]","",mix_alpha$Sample_name))   # invisible x axis
mix_alpha$ID [ mix_alpha$Time_color < 25 ] <-  "1/1 "
mix_alpha$Time_axisX<- factor(as.character(mix_alpha$Time), levels =c (0, unique(mix_alpha$Time)) )  # continuous value to print a shade of color
# plot
ggplot(data=mix_alpha, aes(y=value, x=Time_axisX, color=ID)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  geom_point(size = 1.25 ) +
  geom_point(size =3, alpha=0.4) +
  geom_path(aes(group=ID_LineA), col= "gray75", size = 0.28,
            arrow=arrow(length =unit(0.23,"cm"), type = "closed") 
  ) +
  geom_path(aes(group=ID_LineB), col= "gray50", linewidth = 0.28,
            arrow=arrow(length =unit(0.23,"cm"), type = "closed") 
  ) +
  theme_classic2(base_size = 9.5) +
  # guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  scale_y_continuous(expand = c(0.1,0)) + # to maintain the points inside the plot
  theme(axis.text.x = element_text(angle=25, size=7.5, vjust=1, hjust=1),
        axis.text.y= element_text(angle=0, size=6.5),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1),
        legend.title = element_text(size=9.5, margin = margin(0,0,0,35) ) ,
        legend.spacing.x = unit(0.15,"cm"),
        legend.margin = margin(0,15,0,-25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.25,"cm"),
        legend.key.width  = unit(0.6,"cm")
  ) +
  labs(y="Alpha Diversity Measure", color="F/F ratio" )
ggsave(file="Results/Alfa_diversity_GENUS_along_the_time.png", width = 3.5,height =3.5, dpi=300)




##################### PCoA OF R1 and R2 (including the period of equal feast/famine) ##################

# on genera
data.prop.labels<-data_complete.genus.prop
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day ))
sample_data(data.prop.labels)$F_F_ratio<-as.character( sample_data(data.prop.labels)$Sample_Type )
sample_data(data.prop.labels)$F_F_ratio<-gsub("R1","1/1", sample_data(data.prop.labels)$F_F_ratio )
sample_data(data.prop.labels)$F_F_ratio<-gsub("R2","2/1", sample_data(data.prop.labels)$F_F_ratio )
sample_data(data.prop.labels)$F_F_ratio[sample_data(data.prop.labels)$Experiment_state=="Still_preparing"]<-"1/1"
sample_data(data.prop.labels)$Sample_Type_LineA<-sample_data(data.prop.labels)$Sample_Type
sample_data(data.prop.labels)$Sample_Type_LineA[ sample_data(data.prop.labels)$Sample_Type =="R2" ] <- seq(1,length(sample_data(data.prop.labels)$Sample_Type_LineA[ sample_data(data.prop.labels)$Sample_Type =="R2" ])) # different labels --> diff groups --> these points will not be connected by this line
sample_data(data.prop.labels)$Sample_Type_LineB<-sample_data(data.prop.labels)$Sample_Type
sample_data(data.prop.labels)$Sample_Type_LineB[ sample_data(data.prop.labels)$Sample_Type =="R1" ] <- seq(1,length(sample_data(data.prop.labels)$Sample_Type_LineB[ sample_data(data.prop.labels)$Sample_Type =="R1" ])) # different labels --> diff groups --> these points will not be connected by this line

sample_data(data.prop.labels)$Experiment_day2[ ! sample_data(data.prop.labels)$Experiment_day2 %in% c(1, 38, 45, 63, 95 ) ] <- ""
sample_data(data.prop.labels)$Experiment_day2[ sample_data(data.prop.labels)$Experiment_day2 =="1" ] <- "1  "
sample_data(data.prop.labels)$Experiment_day2[ sample_data(data.prop.labels)$Experiment_day2 =="63" & sample_data(data.prop.labels)$F_F_ratio=="2/1" ] <- "    63"
sample_data(data.prop.labels)$Experiment_day2[ sample_data(data.prop.labels)$Experiment_day2 =="95" & sample_data(data.prop.labels)$F_F_ratio=="2/1" ] <- "      95"

{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
### with the color shift
plot_ordination(data.sqrt_prop, ordBC, color = "F_F_ratio") +
  # scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(1,95), midpoint = 15,
  #                       breaks=c(1, 35, 65, 95)
  # ) +
  geom_point(size=3.5, alpha=0.4) +
  geom_point(size=1.5, alpha=1) +
  geom_path(aes(group=Sample_Type_LineA), col= "gray75", linewidth = 0.28,
            arrow=arrow(length =unit(0.12,"cm"), type = "closed") 
  ) +
  geom_path(aes(group=Sample_Type_LineB), col= "gray50", linewidth = 0.28,
            arrow=arrow(length =unit(0.12,"cm"), type = "closed") 
  ) +
  geom_text( aes( label= Experiment_day2) , color = "black",
             show.legend = F, size= 3.5 , vjust= 1,
             lineheight = 3 ) +
  theme_classic(base_size = 8.85) +
  theme(
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.title = element_text(size=11, margin = margin(0,0,0,35) ),
        legend.text = element_text(size=9), 
        legend.spacing.x = unit(0.25,"cm"),
        legend.margin = margin(0,15,0,-35),
        legend.direction = "horizontal",
        legend.position = "bottom",
  ) +
  #guides(shape="F/F ratio") +
  labs(
       color="F/F ratio",
       x=paste("PC1: ",eigval[1],"% variation"),      
       y=paste("PC2: ",eigval[2],"% variation")
  )
ggsave(file="Results/PCoA_Beta_div_Hellinger_on_Genera_shades_R1R2_INCLUDING_DAYS_SAME_RATIO_FF.png", width = 4.2, height = 4, dpi=300)




########### CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R1 (CENTERED LOG RATIO) ################
# 
#     ALREADY DONE IN THE FIRST PART OF THIS PROJECT, WHEN ONLY THE P SAMPLES HAS BEEN USED!
#
# # selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
# data.genus.temp<-data_complete.genus
# tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
# data.genus.temp <- subset_samples(data.genus.temp, Sample_Type == "R1")
# min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
# ### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
# who<-as.data.frame(otu_table(data.genus.temp))
# who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
# who<-who[!rowSums(who)>=min_prevalence,]
# who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
# 
# # switching to log ratio transformation
# data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# # head(otu_table(data.genus.temp))
# # head(otu_table(data.genus))
# data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
# 
# # extracting the abundances
# abundances<-otu_table(data.genus.temp)
# row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])
# 
# # colnames (abundances)  # already ordered according to the time
# 
# 
# Corr_results<-NULL
# for(i in 1:length(row.names(abundances))){
#   save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "kendall") # correlated with the flow of time (the samples have been ordered accordingly)
#   new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
#   Corr_results<-rbind.data.frame(Corr_results, new_row)
# }
# Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
# {
#   row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
#   Corr_results<-Corr_results[ , -1]
#   colnames(Corr_results)<-c("rho","pvalue")
# }
# # Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
# Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
# Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")
# 
# write.csv2(Corr_results, "Results/Time_series_in_R1/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)
# 
# temp_for_save<-Corr_results[Corr_results$sign=="*", ]
# temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
# significative_abundances_R1_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
# significative_abundances_R1_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )
# 
# 
# 
# con <- file("Results/Time_series_in_R1/Only_certain_genera_have_been_tested.txt")
# sink(con, append = TRUE)
# cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R1 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
# sink()
# close(con)
# 
# suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, min_prevalence, save, new_row))
# 
# 
# 
# 
# 
# ########### LINE PLOTS OF STRONGEST CORRELATIONS
# 
# # fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!
# 
# # THE 12 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
# corr_top<-significative_abundances_R1_pos[1:12] # the table is already ordered
# suppressWarnings(rm(top, others, tabella))
# 
# # resetting to proportions
# samples_selected<-sample_names(data.genus.temp)
# data.genus.temp <- subset_samples(data.genus.prop, Sample_name %in% samples_selected)
# { 
#   tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update) 
#   prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% corr_top)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   tabella<-psmelt(prune.dat_top)
#   tabella$Genus<-gsub("Subgroup_7","Holophagae_Subgroup7",tabella$Genus)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
# }
# tabella<-tabella[order(tabella$Sample_name) , ] # days
# ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
#   theme_classic(base_size =12) + 
#   geom_point(aes(color=Genus), size =2) +
#   geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
#   geom_line(aes(group=Genus), linewidth= 0.4, alpha= 0.4) +
#   scale_color_manual(values= fill_color_19) +
#   theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
#         axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
#         panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
#         panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
#         legend.key.size = unit(0.4, "cm"),
#         legend.text = element_text ( size = 8.5 ),
#         legend.position="bottom",
#         legend.margin = margin(1,36,0,0)) +
#   guides(color=guide_legend(nrow=3)) + 
#   #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
#   scale_y_continuous(breaks = seq( 0, max(tabella$Abundance), 2) ) +
#   labs(x="Experiment day", y="Percent abundance",
#        color="",
#        title = "genera positely correlated with time\n according to Kendall correlation on CLR counts")
# ggsave(file="Results/Time_series_in_R1/Genera_positively_correlated_with_time_R1.png",width=7,height=5, dpi=300) 
# dev.off()
# 
# 
# # THE 12 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
# corr_top<-significative_abundances_R1_neg[(length(significative_abundances_R1_neg)-11):length(significative_abundances_R1_neg)]
# suppressWarnings(rm(top, others, tabella))
# 
# # resetting to proportions
# samples_selected<-sample_names(data.genus.temp)
# data.genus.temp <- subset_samples(data.genus.prop, Sample_name %in% samples_selected)
# { 
#   tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update) 
#   prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% corr_top)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   tabella<-psmelt(prune.dat_top)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
# }
# tabella<-tabella[order(tabella$Sample_name) , ] # days
# ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
#   theme_classic(base_size =12) + 
#   geom_point(aes(color=Genus), size =2) +
#   geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
#   geom_line(aes(group=Genus), linewidth= 0.4, alpha= 0.4) +
#   scale_color_manual(values= fill_color_19) +
#   theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
#         axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
#         panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
#         panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
#         legend.key.size = unit(0.4, "cm"),
#         legend.text = element_text ( size = 8.5 ),
#         legend.position="bottom",
#         legend.margin = margin(1,36,0,0)) +
#   guides(color=guide_legend(nrow=3)) + 
#   #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
#   scale_y_continuous(breaks = seq( 0, max(tabella$Abundance), 2) ) +
#   labs(x="Experiment day", y="Percent abundance",
#        color="",
#        title = "genera negatively correlated with time\n according to Kendall correlation on CLR counts")
# ggsave(file="Results/Time_series_in_R1/Genera_negatively_correlated_with_time_R1.png",width=7,height=5, dpi=300) 
# dev.off()
# 
# suppressWarnings(rm(fill_color_6, prune.dat_top, tax_selected, tabella, corr_6_top))
# 
# 
# 
# # doing this in R2 would require to remove the samples before the change in F/F ... but this does not mean that the resulting genus are increased in specifically R2! (Different, maybe due to the different point number...)
# # the difference between F/F ratio has to be explored through DESeq2





############## DA WITH DESEQ2 _ R1 versus R2 #############

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering steps yet?\n", fill = T)
}

# NB: THE OBJECT "DATA" ALREADY CONTAINS ONLY THE BIOMASSES AFTER THE START OF THE EXPERIMENT (--> SIMILAR PERIOD ESCLUDED)

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DESeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

# system(" echo 'In case of pair analysis, the log2foldchange displayes the change from R1 to R2, see https://support.bioconductor.org/p/105981/ \n' > Results/DA_DESeq2/NB.txt ")  
# system(" echo 'Every result with a baseMean value lower than 50 has been arbitrarly filtered.' >> Results/DA_DESeq2/NB.txt ")  

Table_tot<-NULL
Res_tot<-NULL

# correcting also the time (with more time, the difference increases)
sample_data(data_pruned)$time<-as.factor(gsub("^[0-9][0-9]_","",sample_data(data_pruned)$Sampling_date))

### --> error if included in the model: too much auto-correlations among the variables (the time itself is the paired factor already!)
# then, the three samples collected in april (low difference) will be discarded to avoid this bias
data_pruned2<- subset_samples(data_pruned, ! time %in% c("04_23") )


for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned2, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # if the paired samples are seen as numbers (it depends from on levels name)
  sample_data(d)[["Sampling_date"]]<-as.character(sample_data(d)[["Sampling_date"]])
  
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sampling_date +Sample_Type) # interaction: more time passed = even more differences!
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sample_Type", "R1", "R2"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_R1_vs_R2.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Sample_Type<-factor(Table_tot$Sample_Type, levels = unique(Table_tot$Sample_Type))
Table_tot<-Table_tot[order(match(Table_tot$Sample_Type, levels(Table_tot$Sample_Type))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Sample_Type)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Sample_Type)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

### to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Sample_Type)


# marking the PHA accumulators
# lista_PHA<-table_PHA[,1]     # they derive from papers


# plotting
unique(tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_ f ","\nf_",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_ o ","\no_",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_marine_group","\nmarine_group",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("aceae","ac.",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("lales","al.",tabella_g$Bacteria)
# to display only certain days
tabella_g$Experiment_day2 <- as.character(tabella_g$Experiment_day)
tabella_g$Experiment_day2[ !tabella_g$Experiment_day2 %in% "95" ] <- ""
tabella_g$Experiment_day2[ tabella_g$Experiment_day2 %in% "95" ] <- "last"

plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Sample_Type)) +
  theme_classic(base_size = 8.85) +
  scale_color_manual(values = c("R2"="royalblue","R1"="coral")) +
  facet_wrap2(nrow=5,~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Sample_Type), size=0.35, alpha=1 ) +
  geom_point(aes(color=Sample_Type), size=1.8, alpha=0.35) +
  theme(strip.text.x=element_text(size=6.1,colour="black"),
        strip.switch.pad.wrap = unit(3,"line"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.y = element_text(size=9.5)
  )+
  scale_x_discrete(labels=unique(levels(tabella_g$Sample_Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=6),
        panel.grid.minor.y= element_blank(),
        legend.margin=margin(0, 0, 0, 0), legend.position = "none")
plot_1 +
  geom_line(aes(group=Sampling_date), linewidth=0.1, color="darkgray") +
  labs(title= "Differently abundant genera",
       y="Percent abundance", fill=" ", x="") +
    geom_text(aes(label=Experiment_day2), size=1.95, color="darkgray", alpha= 1,
              hjust=ifelse(tabella_g$Sample_Type=="R1",1.25,-0.25)
    )
ggsave(filename = "Results/DA_DESeq2/Plot_genera_DeSeq2_Type_AFTER_EXP_STARTED.png",
       width = 4, height = 4.9, dpi=300)
# again but with name
plot_1 +
  geom_text(aes(label=Experiment_day), size=2.3, color="black", alpha= 0.8,
            hjust=ifelse(tabella_g$Sample_Type=="R1",1,0)
  ) +
  geom_line(aes(group=Sampling_date), size=0.2, col="darkgray") +
  labs(title= "Differently abundant genera",
       y="Percent abundance", fill="Sample_Type", x="")
ggsave(filename = "Results/DA_DESeq2/Plot_genera_DeSeq2_Type_with_name.png",
       width = 5, height = 6, dpi=300)
dev.off()


##### removing redundant results
{ 
  # 1) if redundant then same abundances
  auto_Max<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, max )
  auto_Min<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, min )
  auto_Median<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, median )
  # reordering the names (lower ranks first --> the redundants are those at higher levels)
  ordered_names <- c( unique(Table_tot[Table_tot$Taxa=="Genus", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Family", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Order", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Class", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Phylum", "Bacteria" ])
  )
  ordered_names[!is.na(ordered_names)]
  auto_Max<-auto_Max[ ordered_names  ]
  auto_Min<-auto_Min[ ordered_names  ]
  auto_Median<-auto_Median[ ordered_names ] # the table it self follows the correct taxonomic order because it is built in this way
  # same abundances
  auto_redund<-names(auto_Max)[ duplicated(auto_Min) & duplicated(auto_Max) & duplicated(auto_Median)  ]
  
  # 2) if redundant then same ASV
  auto_ASV<-Res_tot$ASV[duplicated(Res_tot$ASV)]
  auto_ASV_Names<-unique(Table_tot$Bacteria[Table_tot$OTU %in% auto_ASV])
  auto_redund<-auto_redund[auto_redund %in% auto_ASV_Names]
}

Redund<-auto_redund
#Redund<-c("") # to select redundants (same ASV, same results at different levels)
tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)

tabella_3<-subset(tabella_2, Taxa %in% c("Phylum","Class"))
tabella_3$Bacteria <- gsub("Gammaproteobacteria","Gammaproteobact.", tabella_3$Bacteria)

plot_3<-ggplot(tabella_3, aes(x= Xaxis, y=Abundance, fill=Sample_Type)) +
  theme_classic(base_size = 9) +
  scale_color_manual(values = c("R2"="royalblue","R1"="coral")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Sample_Type), size=0.7, alpha=1) +
  geom_point(aes(color=Sample_Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=8.4,colour="black"), 
        strip.switch.pad.wrap = unit(2,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_3$Sample_Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  labs(title= "Differently abundant families and orders",
       y="Percent Abundance", fill=" ", x="") +
  geom_line(aes(group=Sampling_date), size=0.2) 

tabella_4<-subset(tabella_2, ! Taxa %in% c("Phylum","Class"))
tabella_4$Bacteria<-gsub("-","\n",tabella_4$Bacteria, fixed = F)
# tabella_4$Bacteria<-gsub("Peptostreptococcales","Peptostreptococc.",tabella_4$Bacteria)
plot_4<-ggplot(tabella_4, aes(x= Xaxis, y=Abundance, fill=Sample_Type)) +
  theme_classic(base_size = 9) +
  scale_color_manual(values = c("R2"="royalblue","R1"="coral")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Sample_Type), size=0.7, alpha=1) +
  geom_point(aes(color=Sample_Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=7.8,colour="black"), 
        strip.switch.pad.wrap = unit(2,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_4$Sample_Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  labs(title= "Differently abundant families and orders",
       y="Percent Abundance", fill=" ", x="") +
  geom_line(aes(group=Sampling_date), size=0.2) 

# now an unique plot (no dates)
head_plot<-plot_4 + # refining first plot with the title
  labs(title= "Differently abundant families, orders, class and phyla", y="Percent abundance", x="") + 
  #theme(plot.title = element_text(size=12.5))
  png(filename = "Results/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants.png",
      width = 2750, height = 1500, res=300)
grid.arrange(head_plot,
             plot_3 + labs(title= "", y="Percent abundance", x=""))
dev.off()

# again, an unique plot (with dates)
head_plot<-plot_4 + # refining first plot with the title
  labs(title= "Differently abundant families, orders, class and phyla", y="Percent abundance", x="") + 
  theme(plot.title = element_text(size=12.5)) + 
  geom_text(aes(label=Experiment_day), size=2.3, color="black", alpha= 0.8,
            hjust=ifelse(plot_4$data$Sample_Type=="R1",1,0)
  )
plot3_bis<- plot_3 +
  geom_text(aes(label=Experiment_day), size=2.3, color="black", alpha= 0.8,
            hjust=ifelse(plot_3$data$Sample_Type=="R1",1,0)
  )
png(filename = "Results/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants_WITH_DATE.png",
    width = 2750, height = 1500, res=300)
grid.arrange(head_plot,
             plot3_bis + labs(title= "", y="Percent abundance", x=""))
dev.off()



suppressWarnings(rm(tabella_2, plot_2, plot_1, head_plot))

