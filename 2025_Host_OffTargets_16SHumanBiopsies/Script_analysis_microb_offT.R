#25/05/2024



##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggvenn")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("vegan")
  library("ggpubr")
  library("DESeq2")
  library("ALDEx2")
  # utilities
  library("reshape")
  library("stringr")
  library("xlsx")  
  library("qiime2R")
  library("compositions")
}

{ 
  dir.create("Data_check")
  dir.create("Data_check/PCoA_test")
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Beta_div")
}

options(scipen = 100) # disable scientific annotation


### choosing colors for the barplots  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_25<-c("wheat3","deepskyblue2",
                 "deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "cyan","brown",
                 "grey15","lightgreen","orange2",
                 "darkorchid", "blue1", "violet", "lightblue1",
                 "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                 "lightgrey", "red","darkslategray3")

status_colors<-c("Adenoma"="chartreuse2",
                 "CRC"="chartreuse4",
                 "Chron"="brown4",
                 "Ulcerative_colitis"="coral2",
                 "Healthy_tissue"="deepskyblue2",
                 "Healthy_subject"="blue",
                 "Churg_Strauss"="gold"
                 )
status_colors2<-status_colors 
names(status_colors2)<-gsub("_"," ",names(status_colors)) # to display the labels in the plots

# preset_levels<- c("Ulcerative_colitis", "Chron", "AP","CRC","Healthy_tissue","Healthy_subject")
# preset_levels2<-gsub("_"," ", preset_levels) # to display the labels in the plots

     


####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
for(x in dir("raw_FASTQ_collection_in_folders/") ) {  # each sub folder is the project name

  assign( x , qza_to_phyloseq(features=paste0("QIIME/table_",x,".qza"),       # assign(var_name, value)
                              taxonomy=paste0("QIIME/rep_seqs_decontam/taxonomy_",x,".qza"),
                              # tree=paste0("QIIME/rep_seqs_decontam/rooted-tree_",x,".qza")
                              ) 
  )
  
  assign( paste0(x,"_contam") , qza_to_phyloseq(features=paste0("QIIME/table_",x,".qza"),
                              taxonomy=paste0("QIIME/rep_seqs_with_offtargets/contam_taxonomy_",x,".qza"),
                              # tree=paste0("QIIME/rep_seqs_with_offtargets/rooted-tree_with_contam_",x,".qza")
                              ) 
  )
  
}

### creating unique objects ... 

{
decont_objs<- ls() [ grepl("PRJ",ls()) | grepl("GSE",ls()) ]
contam_objs <- decont_objs [ grepl("_contam", decont_objs) ]
decont_objs <- decont_objs [ ! grepl("_contam", decont_objs) ]
}

data<-get(decont_objs[ 1 ])  # NB: without "get", the loop will try to add a character, not the object named through such character
# this type of loop cannot be performed on a "NULL" object, in case of phyloseq merge function ... then using the first element
for(x in decont_objs [ -1 ]) {
data<-merge_phyloseq(data, get(x) )
}

data_contam<-get(contam_objs[ 1 ])
for(x in contam_objs [ -1 ]) {
  data_contam<-merge_phyloseq(data_contam, get(x) ) # NB: without "get", the loop will try to add a character, not the object named through such character
}


# to match names with the metadata
sample_names(data)<-gsub("_S.*","", sample_names(data))
sample_names(data_contam)<-gsub("_S.*","", sample_names(data_contam))

original_names<-sample_names(data)
sample<-sample_names(data)
#Metadata <- as.data.frame(read.delim(file = "Complete_Metadata.txt", sep="\t") )
Metadata <- as.data.frame(read.delim(file = "Complete_Metadata.txt", sep=",") )
row.names(Metadata)<-Metadata$Run # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(Metadata$Run[!is.na(Metadata$Run)])
Metadata$Run [!Metadata$Run %in% sample ]

# The download of the sample SRR8295839 (PRJNA507548) was not successfull, not even after a second try (?)
# It is a crypt CRC sample, which will be discarded afterwards regardless... then it is just eliminated from the metadata
Metadata<-Metadata[Metadata$Run!="SRR8295839", ]
Metadata<-Metadata[Metadata$Run!="SRR3667009", ] # moreover, there are problem in converting this sra in fastq...
original_names<- original_names[!original_names %in% "SRR8295839" ]
original_names<- original_names[!original_names %in% "SRR3667009" ]
                                  
Metadata<-Metadata[sample, ]
if( ! identical( as.numeric(length(Metadata$Run[!is.na(Metadata$Run)])) ,as.numeric(original_length)-2 ) ){
  stop("\nSomething went wrong with the sample names, or, at least, there are less samples than expected...")
  # Metadata[!Metadata$Run %in% original_names, "Run"]
}
sample_data(data)<-Metadata
sample_data(data_contam)<- Metadata

#head(sample_data(data), n=2)



### removing the colon crypt samples from the same patient where the mucosa is also sampled (the mucosa is the principal area sampled across the works)
Crypt_to_remove<- Metadata$Run[ grepl("Crypt", Metadata$Extra_infos) ]
data<-subset_samples(data, ! Run %in% Crypt_to_remove )
data_contam<-subset_samples(data_contam, ! Run %in% Crypt_to_remove )
Metadata <- Metadata[ ! Metadata$Run %in% Crypt_to_remove , ]


cleaning<-ls()
cleaning<-cleaning [ ! cleaning %in% c("data","data_contam", "Crypt_to_remove",
                                       "Metadata","original_names",
                                       "decont_objs","contam_objs", "status_colors","status_colors2",
                                       "fill_color_5","fill_color_25")
                     ]
rm(list=cleaning)
rm("cleaning")

# save.image("data_pre_noise_filters_gen2024.RData")




################### OFF-TARGETS PROPORTIONS IN EACH PROJECT ###################

suppressWarnings(rm(top, others, tabella))
only_off <- taxa_names(data_contam) [ ! taxa_names(data_contam) %in% taxa_names(data) ]  # only in the contam obj --> they are the off-targets removed
data_contam_modified<-  data_contam
tax_table(data_contam_modified)[ only_off ,"Phylum"]  <- "Host derived"
also_this<- tax_table(data_contam_modified)[ ,"Kingdom"] %in% c("d__Eukaryota", "Unassigned")
tax_table(data_contam_modified)[ also_this ,"Phylum"]  <- "Host derived"
tax_table(data_contam_modified)[ ! taxa_names(data_contam_modified) %in% only_off ,"Phylum"]  <- "Prokaryote derived"
data_contam_modified <- tax_glom(data_contam_modified, "Phylum", NArm = F) # also this is required to speed up the melting
# If the projects are the x axis then their percentage abundances (depending on the number of samples) will be summed... the percentages has to be computed in this step, for each bioproject!
tabella<- rbind("Prokaryote derived","Host derived")
for( x in unique(sample_data(data)[["BioProject"]]) ) {
  temp<-subset_samples(data_contam_modified, BioProject == x)  # work in progress...
  sum_project<-sum(otu_table(temp))
  host_prject<-sum(otu_table(subset_taxa(temp, Phylum=="Host derived")))
  prok_prject<-sum(otu_table(subset_taxa(temp, Phylum!="Host derived")))
  table_temp<-rbind.data.frame(prok_prject, host_prject)
  table_temp<- round( table_temp/sum_project * 100 , 3 )
  tabella<-cbind.data.frame(tabella,table_temp)
}
colnames(tabella) <- c("Source", unique(sample_data(data)[["BioProject"]]) )
tabella2<-melt(tabella)
ggplot(data=tabella2, aes(x=variable, y=value, fill=Source)) +
  #facet_grid2( BioProject ~ . , scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =10) +
  scale_fill_manual(values=c("Prokaryote derived"="deepskyblue3",
                             "Host derived"="grey35")) +
  theme(axis.text.x=element_text(size=9, angle=35, vjust=1, hjust=1),  
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.5 ),
        legend.position="bottom",
        legend.margin = margin(0,0,5,0),
  ) +
  scale_y_continuous(breaks = seq(0,100,10)) +
  guides(fill=guide_legend(nrow=1)) +
  labs(x="BioProject", y="Percent abundance") +
  coord_flip()
ggsave(file="Data_check/Off-targetPercentages.png",width=5,height=2.5, dpi=300)
dev.off()

write.csv2(tabella, file="Data_check/Off-targetPercentages.csv", row.names = F, quote = F)


# extra: checking each sample
data_contam_modified2 <- data_contam_modified
tax_table(data_contam_modified2)[,"Kingdom"] <-"OVERWRITTEN"
data_contam_modified2 <- tax_glom(data_contam_modified2, taxrank = "Phylum", NArm = F)
proportion_H_on_P_row <- as.numeric(otu_table(subset_taxa(data_contam_modified2, Phylum =="Host derived"))) / as.numeric(otu_table(subset_taxa(data_contam_modified2, Phylum !="Host derived")))
write.csv2( rbind(
                  cbind(
                        otu_table(data_contam_modified2),
                        tax_table(data_contam_modified2)[,"Phylum"]
                        ),
                  c( round(proportion_H_on_P_row,2) ,"proportion_H_on_P"),
                  c( sample_data(data_contam_modified2)$Work_symbol , "Project_ID")
                  ),
           file="Data_check/OffTargetCountsForEachSample.csv",
           row.names = F, quote = F
           )


rm(tabella, tabella2, also_this, only_off, data_contam_modified, data_contam_modified2, temp, table_temp, host_prject, prok_prject, sum_project , x )




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
unfiltered_data<-data
unfiltered_contam<-data_contam
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)
# updating the tax table to distinguish the uncultured and NAs
taxa_temp<-as.data.frame(tax_table(data.genus.temp))
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
  tax_table(data.genus.temp)<-as.matrix(taxa_temp)
}


###### cutting under 0.01% to remove noises/contaminants, too conservative but also safe cutoff, see  PMID: 31164452,  PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.01, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_001_max_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) max(x) <= 0.01, TRUE)))[["Genus"]]
filtered<-filtered[! is.na(filtered) ]
data<-subset_taxa(data, ! Genus %in% filtered)
data<-prune_taxa(taxa_sums(data) > 0, data) 


###### Applying also prevalence filter: at least 1 read of that genus at least in 1% of the samples
who<-as.data.frame(otu_table(data.genus.temp) )
who<-apply(who, MARGIN=2, function(x) ifelse(x>1, 1, 0)) # if at least one read in the sample then it gets "a point"
how_many<- round( length(sample_names(data)) * 1 / 100  ,0) # --> 12 samples
who<-who[!rowSums(who)>how_many, ] # more than x "points"
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
who<-who[!is.na(who)]
who # just a check --> ok!
write.csv(who, "Data_check/Filtered_genera_that_were_NOT_AT_LEAST_in_1percent_of_samples.csv", row.names = F)
data <-subset_taxa(data, ! Genus %in% who)


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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Unassigned_domain_checking_after_Bowtie2.csv", row.names = T, quote = F)
data<- subset_taxa(data, ! Kingdom %in% c("Unassigned","d__Eukaryota") ) # removing those unassigned reads
data<- subset_taxa(data, ! is.na(Phylum) ) # removing also unassigned at Phylum Level (too strange)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


# removing chloroplasts and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria\n\n")
  data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mitochondria"))
}


### Removing every contaminant also in the object with human ASVs, excluding "Unassigned" ASVs (human) and eukaryotes
data_contam <-subset_taxa(data_contam, ! Genus %in% who)
data_contam <-subset_taxa(data_contam, ! Genus %in% filtered)
# removing chloroplast and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data_contam))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria from the human contaminated obj \n\n")
  data_contam<-subset_taxa(data_contam, ! Genus %in% c("Chloroplast","Mitochondria"))
}

proof1<- "Marker of the filtering, it is required for the script"
suppressWarnings( rm(filtered, how_many, data.genus.temp, who, taxa_temp) )


# save.image("data_POST_noise_filters_26dic2024.RData")




############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

if(! "proof1" %in% ls() ){  # if proof 2 is present then same sample is already removed...
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

# Raw number of reads
Original_number_table<-NULL
for(x in decont_objs ) {
  temp<- read.table( paste0("QIIME/Original_number_of_reads_",x,".qza.tsv"), sep= "\t", header=T, row.names = 1)
  Original_number_table<-rbind.data.frame(Original_number_table, temp)
}

# DADA2 initial reads (after cutadapt) and processed reads
DADA2_read_number<-NULL
for(x in decont_objs ) {
   temp<- read.table( paste0("QIIME/Percentuals_Initial_and_processed_reads_",x,".tsv"), sep= "\t", header=T, row.names = 1)
   DADA2_read_number<-rbind.data.frame(DADA2_read_number, temp)
}

row.names(Original_number_table) <- gsub("_S.*","", row.names(Original_number_table) )
row.names(DADA2_read_number) <- gsub("_S.*","", row.names(DADA2_read_number) )
Original_number_table <- Original_number_table [ row.names(Original_number_table) %in% Metadata$Run , ]
DADA2_read_number <- DADA2_read_number[row.names(DADA2_read_number) %in% Metadata$Run , ]

Original_number<-sum( c(Original_number_table$forward.sequence.count, Original_number_table$reverse.sequence.count))

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
after_removing_offtargets<-as.data.frame(colSums(otu_table(unfiltered_data))) # unfiltered is the data post Bowtie but without decontaminants removed!
after_filter_number<-as.data.frame(colSums(otu_table(data)))
# same order of rows
Original_number_table<-Original_number_table[row.names(DADA2_read_number), ]
after_removing_offtargets<-after_removing_offtargets[row.names(DADA2_read_number), ]
after_filter_number<-after_filter_number[row.names(DADA2_read_number), ]
# creating the table
table<-cbind.data.frame(Original=Original_number_table$forward.sequence.count, # merged --> just one column, not both of them
                        After_quality_filter=DADA2_read_number$non.chimeric,
                        After_offTargets_filter=after_removing_offtargets,
                        After_contam_filter=after_filter_number)
# adding factors to re-group the table
table$factor<-as.data.frame(sample_data(data))[row.names(Original_number_table), ]$BioProject
table$Samples<-row.names(Original_number_table)
# plotting
ggplot(aes(x=Samples, y=Original), data=table) +
  facet_grid2( factor ~.,  scales = "free", independent=T ) +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Quality read filters"), alpha= 0.8,
            width = 0.6, stat = "identity") +
  geom_bar( aes(y= After_offTargets_filter, fill="Off Target filter"), alpha= 0.8,
            width = 0.45, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Contaminant filters"),
            width= 0.3, stat = "identity") +
  theme_classic( base_size = 10.8) +
  theme(axis.text.x = element_blank(),
        #axis.text.x = element_text(size=4, angle = 90, vjust=0.5),
        axis.text.y = element_text(size=3),
        strip.text = element_text(size=7.5),
        panel.grid.major.y = element_line(linewidth=0.1, color="grey"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size=9.2),
        legend.title = element_text(size=9.8),
        legend.key.height = unit(0.4,"cm")
  ) +
  scale_fill_manual(name='Number of reads:  ',
                    breaks=c('Original number', 'Off Target filter', 'Quality read filters', 'Contaminant filters'),
                    values=c('Original number'='green3', 'Quality read filters'='darkgreen', 'Off Target filter'='coral', 'Contaminant filters'='red4')) +
  scale_y_sqrt(breaks = c(10000, seq(0, max(table$Original), 75000))) +
  labs(y="Reads abundance", x="FASTQ")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 8, height = 9, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)

rm(table,  Original_number_table, DADA2_read_number, after_filter_number, after_removing_offtargets)




############################ RAREFACTION ANALYSIS ################################

if(! "proof1" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as generated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of the original samples (and in the same order!!!)
  sat<-0
  unsaturated_list<-NULL
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
      cat(i,names[i],slope,check,"\n")
      unsaturated_list<-rbind.data.frame(unsaturated_list, cbind(names[i], check))
      #			text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],labels=slope,col=check,cex=0.5)
      #			points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col=check,pch=16,cex=1)
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  unsaturated_list<<-unsaturated_list
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

png(file="Data_check/Rarefaction_curve.png",width=3500,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
evalslopes(r,sample_data(data)$Sample_name,lim=0.001,cex=0.4)
dev.off()


### excluding the unsaturated samples
unsaturated_list<-unsaturated_list[unsaturated_list$check=="red","V1"]
pairs_unsaturated<-Metadata[Metadata$Pair_ID!="no_pair" & Metadata$Sample_name %in% unsaturated_list , "Sample_name"]
# NB: the pairs are removed ONLY for the paired part of this script
{
  data<-subset_samples(data, ! Sample_name %in% unsaturated_list )
  unfiltered_data<-subset_samples(unfiltered_data, ! Sample_name %in% unsaturated_list )
  data_contam<-subset_samples(data_contam, ! Sample_name %in% unsaturated_list)
  Metadata_no_unsat<- Metadata[ ! Metadata$Sample_name %in% unsaturated_list , ]
  proof2<-paste0("A total of ",length(unsaturated_list)," unsaturated samples, of which ", length(pairs_unsaturated)," were in pairs with other samples in the dataset, have been removed from the microbiota analysis")
}

rm(r)



################# CHECKING THE COMPOSITION AFTER FILTERING (UNFILTERED_vs_FILTERED) ####################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\n Wait! Did you perform the filtering steps??? \n\n")
}

data.versus<-tax_glom(unfiltered_data, taxrank = "Genus", NArm=F)   # not yet filtered BUT decontaminated
data.versus<-transform_sample_counts(data.versus, function(x) (x/sum(x))*100)
data.genus<-tax_glom(data, taxrank = "Genus", NArm=F)   # not yet filtered BUT decontaminated
data.genus.prop<-transform_sample_counts(data.genus, function(x) (x/sum(x))*100)



################# BRAY HELLINGER

##### unfiltered BRAY
suppressWarnings(rm(data.versus.labels, data.sqrt_prop))
data.versus.labels<-data.versus
{data.sqrt_prop<-transform_sample_counts(data.versus.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  #ordBC$vectors[ ,"Axis.1"]<-   - ordBC$vectors[ ,"Axis.1"]  # to mirror the positions and compare them to the other plot
  #ordBC$vectors[ ,"Axis.2"]<-   - ordBC$vectors[ ,"Axis.2"]  # to mirror the positions and compare them to the other plot
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors) +
  guides(color="none") +
  geom_point(size=3, alpha= 0.5) + 
  theme_bw(base_size = 13) + 
  stat_ellipse(linewidth=0.16) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name),
            color="black", size=1, show.legend = FALSE) +
  labs(title="Unfiltered", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.versus.labels, data.sqrt_prop))
data.versus.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.versus.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors) +
  theme(legend.margin = margin(0,2,0,-2)) +
  geom_point(size=3, alpha= 0.5) +
  theme_bw(base_size = 13) +
  stat_ellipse(linewidth=0.16) + 
  geom_text(aes(label=Sample_name), 
            color="black", size=1, show.legend = FALSE) +
  labs(title="Filtered", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_Filtered_vs_Unfiltered.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1, widths = c(1,1.32))
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.versus.labels))



################# EUCLIDEAN HELLINGER

##### unfiltered EUCLIDEAN
suppressWarnings(rm(data.versus.labels, data.sqrt_prop, p1))
data.versus.labels<-data.versus
{data.sqrt_prop<-transform_sample_counts(data.versus.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#  ordBC$vectors[ ,"Axis.1"]<-   - ordBC$vectors[ ,"Axis.1"]  # to mirror the positions and compare them to the other plot
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors) +
  guides(color="none") +
  geom_point(size=3, alpha= 0.5) + theme_bw(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=Sample_name), 
            color="black", size=1, show.legend = FALSE) +
  labs(title="Unfiltered", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))


##### filtered EUCLIDEAN
suppressWarnings(rm(data.versus.labels, data.sqrt_prop))
data.versus.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.versus.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors) +
  geom_point(size=3, alpha= 0.5) +
  theme_bw(base_size = 13) +
  stat_ellipse(linewidth=0.16) + 
  theme(legend.margin = margin(0,2,0,-2)) +
  geom_text(aes(label=Sample_name), 
            color="black", size=1, show.legend = FALSE) +
  labs(title="Filtered", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_HELLINGER_Filtered_vs_Unfiltered.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1, widths = c(1,1.32))
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.versus.labels))




############# CHECKING THE COMPOSITION AFTER FILTERING (WITH_vs_WITHOUT HUMAN) ###############

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\n Wait! Did you perform the filtering steps??? \n\n")
}

if(! "data_contam_G_prop" %in% ls( )) {  # to avoid losing time due to how much big is this object!
  data_contam_G_prop<-tax_glom(data_contam, taxrank = "Genus", NArm=F)
  data_contam_G_prop<-transform_sample_counts(data_contam_G_prop, function(x) (x/sum(x))*100)
}


##### contam HELLINGER
suppressWarnings(rm(data.versus.labels, data.sqrt_prop, p1))
data.versus.labels<-data_contam_G_prop
{data.sqrt_prop<-transform_sample_counts(data.versus.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  #ordBC$vectors[ ,"Axis.2"]<-   - ordBC$vectors[ ,"Axis.2"]  # to mirror the positions and compare them to the other plot
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1_project<-plot_ordination(data.sqrt_prop, ordBC, color = "BioProject") +
  #scale_color_manual(values=c("AAAA"="coral","BBBB"="chartreuse")) +
  guides(color="none") +
  geom_point(size=3, alpha= 0.5) + theme_bw(base_size = 13) + stat_ellipse(size=0.2) + 
  # geom_text(aes(label=Sample_name), 
  #           color="black", size=1, show.legend = FALSE) +
  labs(
    #title="With human ASVs", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
p1_status<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors) +
  guides(color="none") +
  geom_point(size=3, alpha= 0.5) + 
  theme_bw(base_size = 13) + 
  stat_ellipse(linewidth=0.16) + 
  # geom_text(aes(label=Sample_name), 
  #           color="black", size=1, show.legend = FALSE) +
  labs(
    #title="With human ASVs", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")
       )


##### filtered HELLINGER
suppressWarnings(rm(data.versus.labels, data.sqrt_prop))
data.versus.labels<-data.genus.prop
sample_data(data.versus.labels)$Sample_Type<-gsub("_"," ", sample_data(data.versus.labels)$Sample_Type )
{data.sqrt_prop<-transform_sample_counts(data.versus.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  ordBC$vectors[ ,"Axis.2"]<-   - ordBC$vectors[ ,"Axis.2"]  ### ! to flip the Y axis, in order to mirror the positions and compare them to the other plot
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2_project<-plot_ordination(data.sqrt_prop, ordBC, color = "BioProject") +
  #scale_color_manual(values=c("AAAA"="coral","BBBB"="chartreuse")) +
  #guides(color="none") +
  geom_point(size=3, alpha= 0.5) + theme_bw(base_size = 13) + stat_ellipse(size=0.2) + 
  # geom_text(aes(label=Sample_name), 
  #           color="black", size=1, show.legend = FALSE) +
  labs(
    #title="Without human ASVs", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
p2_status<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors2) +
  #guides(color="none") +
  geom_point(size=3, alpha= 0.5) + 
  theme_bw(base_size = 13) +
  theme(legend.text = element_text(size=12.5),
        legend.margin = margin(0,1,0,0)
        ) +
  stat_ellipse(linewidth=0.16) + 
  # geom_text(aes(label=Sample_name), 
  #           color="black", size=1, show.legend = FALSE) +
  labs(
    #title="Without human ASVs", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))


png(filename = "Data_check/PCoA_test/PCoA_HELLINGER_ContamHuman_vs_Decontam_ABOUT_PROJECT.png", width = 3300, height = 1800, res=300)
ggarrange(p1_project,p2_project, nrow = 1, widths = c(1,1.32))
dev.off()

png(filename = "Data_check/PCoA_test/PCoA_HELLINGER_ContamHuman_vs_Decontam_ABOUT_STATUS.png", width = 3300, height = 1800, res=300)
ggarrange(p1_status,p2_status, nrow = 1, widths = c(1,1.32) )
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.versus.labels))



# ################# wUNIFRAC PROP
# 
# ##### unfiltered wUNIFRAC
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
# data.prop.labels<-data_contam_G_prop
# #{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
# {DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
#   ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p1<-plot_ordination(data.prop.labels, ordBC, color = "Sample_Type") +
#   scale_color_manual(values=c("AAAA"="coral","BBBB"="chartreuse")) +
#   guides(color="none", shape="none") +
#   geom_point(size=3, alpha= 0.5) + theme_bw(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=Sample_name), 
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="Contaminated", 
#        color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# ##### filtered wUNIFRAC
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
# data.prop.labels<-data.prop
# #{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
# {DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
#   ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p2<-plot_ordination(data.prop.labels, ordBC, color = "Sample_Type") +
#   scale_color_manual(values=c("AAAA"="coral","BBBB"="chartreuse")) +
#   theme(legend.margin = margin(0,2,0,-2)) +
#   geom_point(size=3, alpha= 0.5) + theme_bw(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label= Sample_name), 
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="Unfiltered", 
#        color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_Contaminated_vs_UnfilteredBUTdecontam.png", width = 3200, height = 1800, res=300)
# ggarrange(p1,p2, nrow = 1)
# dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))



############# PLOT TO VISUALIZE THE CONTINGENCY TABLE OF DISEASE AND WORKS ###########

### Using TOP 5 Phyla as background
data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
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
tabella$Sample_Type<-gsub("_"," ",tabella$Sample_Type)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) +
  facet_grid2( Sample_Type ~ BioProject, scales = "free",independent = T) + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =11.5) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 13.2 ),
        legend.position="bottom",
        legend.margin = margin(-2,0,0,0),
        axis.ticks = element_blank()
  ) +
  guides(fill=guide_legend(nrow=1)) +
  labs(x="Patients", y="Percent abundance", fill="")
ggsave(file="Data_check/Distribution_of_disease_and_works_AFTER_REMOVING_UNSATURATED.png",width=9,height=9, dpi=300)

# writing it down as a true contingency table
table(Metadata$BioProject, Metadata$Sample_Type)
write.csv ( file="Data_check/Distribution_of_disease_and_works_INCLUDING_UNSATURATED.csv", table(Metadata$BioProject, Metadata$Sample_Type) )
write.csv ( file="Data_check/Distribution_of_disease_and_works_AFTER_REMOVING_UNSATURATED.csv", table(Metadata_no_unsat$BioProject, Metadata_no_unsat$Sample_Type) )

suppressWarnings( rm (table, tabella, top, prune.dat_top, others, prune.data.others))



####################### PREPARATION OF THE DATA #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\n Wait! Did you perform the filtering steps??? \n\n")
}

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
#data.class = tax_glom(data, taxrank = "Class", NArm = F)
#data.order = tax_glom(data, taxrank = "Order", NArm = F)
#data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
# data.class.prop <- transform_sample_counts(data.class, function(x) x/sum(x)*100)
# data.order.prop <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
# data.fam.prop <- transform_sample_counts(data.fam, function(x) x/sum(x)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
#Taxa.fam<-as.data.frame(tax_table(data.fam))
Taxa.phy<-as.data.frame(tax_table(data.phy))
#Taxa.class<-as.data.frame(tax_table(data.class))
#Taxa.order<-as.data.frame(tax_table(data.order))
}


# adding informations to missing genera names
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


### unfilt and contam object
if(! "data_contam_G" %in% ls( )) {  # to avoid losing time due to how much big is this object!
  data_contam_G<-tax_glom(data_contam, taxrank = "Genus", NArm=F)
  data_contam_G_prop<-transform_sample_counts(data_contam_G, function(x) (x/sum(x))*100)
}


### Only Healthy and CRC (works with paired samples only)
#data_sub<-subset_samples(data, Sample_Type %in% c("Healthy_tissue","CRC") )
data_sub<-subset_samples(data, BioProject %in% c("PRJNA507548","PRJEB57580","PRJNA743150","PRJNA995580","PRJNA325650") )
data_sub<-subset_samples(data_sub, ! Pair_ID %in% "no_pair" ) # only the paired samples are here analysed
data_sub<-subset_samples(data_sub, ! Sample_name %in% pairs_unsaturated ) # the pairs of these sample were unsaturated and then removed
data_sub<-tax_glom(data_sub, taxrank = "Genus", NArm=F)
data_sub.prop<-transform_sample_counts(data_sub, function(x) (x/sum(x))*100 )


### Only CRC tissues
#data_CRC<-subset_samples(data, Sample_Type %in% c("CRC") )


# save.image(file="data_post_filters_every_object_prepared_gen2024.Rdata")



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
  write.csv2(cbind(as(otu_table(data_contam),"matrix"),as(tax_table(data_contam),"matrix")),file="Results/Abundances/Raw_counts/Contaminated_offTarget_obj_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
  }

options(scipen = 100)
dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data_contam_G_prop),"matrix"),as(tax_table(data_contam_G_prop),"matrix")),file="Results/Abundances/Relative_abundances/Contaminated_obj_Percentcount_genus.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
}



###################### ABUNDANCES BAR PLOT ##########################

### TOP 5 Phyla
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
tabella$Sample_Type<-gsub("_","\n",tabella$Sample_Type)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) +
  facet_grid2( Sample_Type ~ . , scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_blank(),  
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 13 ),
        legend.position="bottom",
        legend.margin = margin(0,0,5,0),
        ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Percent abundance",  
       title = "Five most abundant phyla",
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/TOP_5_phyla.png",width=8.2,height=7.5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)



### TOP Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:24]
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
tabella$Genus<-gsub("[","",tabella$Genus, fixed=T)
tabella$Genus<-gsub("]","",tabella$Genus, fixed=T)
tabella$Genus<-gsub("_group","",tabella$Genus, fixed=T)
# Unassigned_phylum and Others have to be the last ones and Ruminococcus_torques (the one that is not in TOP no more if off target are not excluded) should have a distinctive color...
tabella$Genus<-factor(tabella$Genus,
                      levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% c("Others","Ruminococcus_torques")],
                                 "Ruminococcus_torques","Others"
                      ))
}
tabella$Sample_Type<-gsub("_","\n",tabella$Sample_Type)
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
  facet_grid2( Sample_Type ~ . , scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =10.9) +
  scale_fill_manual(values=fill_color_25) +
  theme(axis.text.x=element_blank(),  
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=8), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.24, "cm"),
        legend.text = element_text ( size = 9.2 ),
        legend.position="bottom", 
        legend.margin = margin(-5,0,0,3),
        title=element_blank(),
        ) + 
  guides(fill=guide_legend(nrow=5)) + 
  labs(x="Patients", y="Percent abundance",
       # title = "Most abundant genera",
       fill="",
       caption = " 'Others' includes every genus below rank 24 ")
ggsave(file="Results/Abundances/TOP_genera.png",width=7.8,height=7,dpi=300)
dev.off()

# means with Contamintants
data_top_contam<- prune_taxa( taxa_names(prune.dat_top), data_contam_G_prop)# , Genus %in% tax_table(prune.dat_top)[,"Genus"] )
mwc<-round( as.numeric(apply( otu_table(data_top_contam) ,1,mean)),2)


# means of TOP genera
to_save<- cbind.data.frame( "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]] , 
                            "Genus"= as.data.frame(tax_table(data_top_contam))[["Genus"]] ,
                            "Average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2),
                            "Average with contaminants"= mwc 
                            )
to_save<-to_save[order(to_save$Average, decreasing=T), ]

write.xlsx(file = "Results/Abundances/TOP_genera_Average_abundances.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, to_save))



# again, but geometric means (required by the reviewer 1)
count_no_proportions <- as(otu_table(data.genus),"matrix")
count_no_proportions <- count_no_proportions+0.5   # to solve the zeros, otherwise every log will be infinite and the ratio 0!
#geom_means <- geometricmeanRow( count_no_proportions )
# each geometric mean is equal to the these of the contaminated obj, as the geom is computed on the row itself without intereference from other rows
count_no_proportions <- clr(t(count_no_proportions))  # NB: translating, as the functions expect the obs to be columns (see the example data of the compositional package)
clr_means <- colMeans(count_no_proportions)
clr_means_no_contam <- round( as.numeric(clr_means[top]),3 )  # "top" sorts and selects the most abundant genera abund
# now on the contaminated obj
count_no_proportions <- as(otu_table(data_contam_G),"matrix")
count_no_proportions <- count_no_proportions+0.5   # to solve the zeros, otherwise every log will be infinite and the ratio 0!
count_no_proportions <- clr(t(count_no_proportions))  # NB: translating, as the functions expect the obs to be columns (see the example data of the compositional package)
clr_means <- colMeans(count_no_proportions)
clr_means_with_contam <- round( as.numeric(clr_means[top]),3 )  # "top" sorts and selects the most abundant genera abund
to_save<- cbind.data.frame( "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "AverageCLR"= clr_means_no_contam,
                            "AverageCLR with contaminants"= clr_means_with_contam 
                            )
to_save<-to_save[order(to_save$AverageCLR, decreasing=T), ]

write.xlsx(file = "Results/Abundances/TOP_genera_Average_with_CLR.xlsx", row.names = F, to_save)




########################### PCoA  ########################

### on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Sample_Type<-gsub("_"," ",sample_data(data.prop.labels)$Sample_Type)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
### Project
plot_ordination(data.sqrt_prop, ordBC, color = "BioProject") +
  # scale_color_manual(values=c("BBBB"="chartreuse","AAAA"="coral")) +
  # geom_point(size=3, alpha= 0.3) + 
  theme_classic(base_size = 10.2) +
  stat_ellipse(linewidth=0.18) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name),
            color="black",
            size=1.7,
            show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Hellinger_on_Genera_PROJECT.png", width = 6.5, height = 5, dpi=300)
### Sample_Type
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors2) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10.2) +
  stat_ellipse(linewidth=0.18) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Work), color="black", size=2, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="The letter in each points indicates the project of provenience, see Methods")
ggsave(file="Results/Beta_div/PCoA_Hellinger_on_Genera_DISEASE.png", width = 6.5, height = 5, dpi=300)

rm(eigval, DistBC, ordBC,data.prop.labels)



######## AGAIN, BUT FOCUSING ON HEALTHY VS CRC
ASV.genus.prop<-as.data.frame(otu_table(data_sub.prop))
sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t

# NB: there are other samples beside then the paired ones
# distance
perm_g<- vegan::adonis(sample_OTU ~Sample_Type + Pair_ID , data=as(sample_data(data_sub),"data.frame"), permutations = 9999, method="euclidean")
perm_g$aov.tab

#dispersion
{dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  disper<-vegan::betadisper(dist,as(sample_data(data_sub),"data.frame")$Sample_Type)
  disp_genus<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  b<-disp_genus$tab
}
b

# plot
data.prop.labels<-data_sub.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
### Sample_Type and Work
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type", shape = "BioProject") +
  geom_point(size=3, alpha= 0.3) + theme_classic(base_size = 10.2) + 
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", 
       color="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Hellinger_Genera_Healthy_vs_CRC_v1.png", width = 5, height = 4, dpi=300)
### Sample_Type only
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=status_colors) +
  geom_point(size=3, alpha= 0.3) +
  theme_classic(base_size = 10) + 
  theme(title = element_text(size= 8.8),
        legend.margin = margin(0,0,0,-5)
        ) +
  stat_ellipse(linewidth=0.10) +
  geom_line(aes(group=Original_Sample_name),
            linewidth=0.14,
            color="gray",
            alpha=0.8
            ) +
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title=paste0("'~ Condition + SampleID ' Pr(>F): ", perm_g$aov.tab$`Pr(>F)`[1], " and ", perm_g$aov.tab$`Pr(>F)`[2], "\n", "Betadipersion Pr(>F):",b$`Pr(>F)`[1]),"\n\n", 
       color="", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Hellinger_Genera_Healthy_vs_CRC_v2.png", width = 5, height = 4, dpi=300)


rm(eigval, DistBC, ordBC,data.prop.labels)



#################### ALPHA DIV (CHECKING THE WORKS) #######################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="BioProject", color="BioProject")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Sample_name, obs$Sample_name) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$BioProject, New_data$Sample_name),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
# with points only
pAlpha + theme_classic() +
  labs(x="BioProjects", title="Genera (after decontamination and filters)") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=10)) +
  #  stat_compare_means(aes(group = BioProject), label="p.format", method = "wilcox.test", paired = T,
  #                    label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Data_check/Alfa_diversity_between_BioProjects.png", width = 6,height =5.6, dpi=300)

suppressWarnings(rm(pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor))



####################### ALPHA DIV (HEALTHY VS CRC) #######################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

# the dataset has to be even (x CRC vs x Healthy) to compute a pair wilcoxon ...
data.temp<-subset_samples(data_sub, ! grepl("left", as.data.frame(sample_data(data_sub))[["Extra_infos"]] ) )
c<-sample_data(subset_samples(data.temp, Sample_Type == "CRC" ) )[["Pair_ID"]] # the CRC list is somehow shorter than Healthy list (maybe few sequences are missing... check with the newer versions of the script) 28/12/2024
h<-sample_data(subset_samples(data.temp, Sample_Type != "CRC" ) )[["Pair_ID"]] # the CRC list is somehow shorter than Healthy list (maybe few sequences are missing... check with the newer versions of the script) 28/12/2024
without_pairs<- h[! h %in% c]
without_pairs<- c( c[! c %in% h] , without_pairs )
data.temp<-subset_samples(data.temp, ! Pair_ID %in% without_pairs )

pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Sample_Type", color="Sample_Type")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Original_Sample_name, obs$Original_Sample_name) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Sample_Type, New_data$Original_Sample_name),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
# with points only
pAlpha +
  theme_classic() +
  scale_color_manual(values = status_colors) +
  geom_line(aes(group = pAlpha$data$Pair_ID),
            col="grey",
            linewidth=0.08) +
  labs(x="",) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=22, vjust=1, hjust=1, size=9.8)) +
  stat_compare_means(aes(group = Sample_Type,
                         label = sprintf("p = %4.3g",
                                         as.numeric(..p.format..))),
                     method = "wilcox.test", paired = T,
                     label.x= 1.4, size=3.2, label.y.npc = "top",
                     vjust=-0.62, hjust=0.4) +
  theme(panel.grid.major.y = element_line(linewidth=0.4), panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Alfa_div_GENUS_paired_HealthyTiss_vs_CRC.png", width = 5.8,height =4.8, dpi=300)

suppressWarnings(rm(pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor))



#################### VEEN DIAGRAM (CRC _vs Healthy) ##########################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

data.venn<-data_sub.prop
# no_venn<-as.data.frame(tax_table(filter_taxa(data.venn, function(x) max(x) <= 0.01, TRUE)))[["Genus"]]
# data.venn<-subset_taxa(data.venn, ! Genus %in% no_venn)

taxa_temp<-as.data.frame(tax_table(data.venn))
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
  tax_table(data.venn)<-as.matrix(taxa_temp)
}


### abundance AND prevalence filter
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0, 1, 0)) # if more than 0 (just present) --> "a point"
#length(sample_names(subset_samples(data.venn, Sample_Type=="CRC")))
#length(sample_names(subset_samples(data.venn, Sample_Type!="CRC")))
# same length (150) in the two groups, being paired ...
how_many_venn<- length(sample_names(subset_samples(data.venn, Sample_Type=="CRC")))*10 / 100
# 10 % --> at least 30
who<-who[!rowSums(who)>=how_many_venn, ] # more than x "points" --> at least in x samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


CRC<-subset_samples(data.venn, Sample_Type=="CRC")
CRC<-as.character(tax_table(prune_taxa(taxa_sums(CRC)>0, CRC))[,"Genus"])

Healthy<-subset_samples(data.venn, Sample_Type=="Healthy_tissue")
Healthy<-as.character(tax_table(prune_taxa(taxa_sums(Healthy)>0, Healthy))[,"Genus"])

ONLY_IN_Healthy<- Healthy[! Healthy %in% CRC]
### Now collapsing in one vector for the script  
ONLY_IN_Healthy<- paste(ONLY_IN_Healthy, collapse = ", ")
head(ONLY_IN_Healthy)

ONLY_IN_CRC<- CRC[! CRC %in% Healthy]
ONLY_IN_CRC<- paste(ONLY_IN_CRC, collapse = ", ")
head(ONLY_IN_CRC)

con<-file("Results/Beta_div/HEALTHY_vs_CRC_exclusive_bacteria.txt")
sink(con, append=TRUE)
cat("\n\nONLY IN CRC", fill=TRUE)
cat(ONLY_IN_CRC)
cat("\n\nONLY IN Healthy", fill=TRUE)
cat(ONLY_IN_Healthy)
sink()
close(con)

# version 1
x<-list('CRC                 '=CRC,
        "          Healthy tissue"=Healthy)
ggvenn(x, stroke_size = 0.6,
       set_name_size = 5.2,
       show_percentage = F,
       stroke_color = "darkgray",
       text_size = 6,
       text_color = "black",
       fill_color = c("lightgray","coral")) +
  theme(plot.margin = margin(-5,-15,0,-15)
  )
ggsave(filename = "Results/Beta_div/Venn_Healthy_vs_CRC_onlyGeneraAbove10percentSamples.png", width = 4.5, height = 4.5, dpi=300, bg = "white")
dev.off()

rm(data.venn)



#################### VEEN DIAGRAM (Contaminated _vs Decontaminated) ##########################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


data.venn1<-unfiltered_contam   # contaminanted and unfiltered
taxa_temp<-as.data.frame(tax_table(data.venn1))
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
  tax_table(data.venn1)<-as.matrix(taxa_temp)
}

data.venn2<-unfiltered_data   # decontaminated but yet unfiltered
taxa_temp<-as.data.frame(tax_table(data.venn2))
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
  tax_table(data.venn2)<-as.matrix(taxa_temp)
}

Contaminated<-as.character(tax_table(prune_taxa(taxa_sums(data.venn1)>0, data.venn1))[,"Genus"])

Decontaminated<-as.character(tax_table(prune_taxa(taxa_sums(data.venn2)>0, data.venn2))[,"Genus"])

ONLY_IN_Decontaminated<- Decontaminated[! Decontaminated %in% Contaminated]
# removing every NA_oNA (unassigned, searching for bacteria!)
ONLY_IN_Decontaminated<-ONLY_IN_Decontaminated[!grepl("NA_ o NA",ONLY_IN_Decontaminated)]
### Now collapsing in one vector for the script  
ONLY_IN_Decontaminated<- paste(ONLY_IN_Decontaminated, collapse = ", ")
head(ONLY_IN_Decontaminated)

ONLY_IN_Contaminated<- Contaminated[! Contaminated %in% Decontaminated]
ONLY_IN_Contaminated<- ONLY_IN_Contaminated[!grepl("NA_ o NA",ONLY_IN_Contaminated)]
ONLY_IN_Contaminated<- unique(ONLY_IN_Contaminated)
ONLY_IN_Contaminated<- paste(ONLY_IN_Contaminated, collapse = ", ")
head(ONLY_IN_Contaminated)
# Contaminated[ !(grepl("NA_", Contaminated)) & !(Contaminated %in% Decontaminated) ]

con<-file("Results/Beta_div/Decontaminated_vs_Contaminated_exclusive_bacteria.txt")
sink(con, append=TRUE)
cat("\n\nONLY IN Contaminated", fill=TRUE)
cat(ONLY_IN_Contaminated)
cat("\n\nONLY IN Decontaminated", fill=TRUE)
cat(ONLY_IN_Decontaminated)
sink()
close(con)




######### DA WITH DESEQ2 (CRC vs Healthy Tissue, no human reads) #####################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
otu_table(data_pruned)<-otu_table(data_pruned)+1 # avoids the error "every gene contains at least one zero"
# NB: rows with small counts can led to warnings as "larger maxit argument with nbinomWaldTest" https://support.bioconductor.org/p/65091/


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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sample_Type + Pair_ID)
  DE<-DESeq(DEseq_data) 
  res<-results(DE, contrast= c("Sample_Type", "Healthy_tissue", "CRC"))
  res = res[order(res$padj, na.last=NA), ]
  fc_threshold<- abs(res$log2FoldChange) - abs(res$lfcSE) # minimum fc
  res<-res[(res$padj < 0.05) & fc_threshold>1,]
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Healthy_tissue_vs_CRC.csv"), row.names = F, quote=F, na = "")
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
# note from DESeq:
# "note: fitType='parametric', but the dispersion trend was not well captured by the function: y = a/x + b, and a local regression fit was automatically substituted."

# View(Res_tot)
write.csv(Res_tot, file="Results/Every_DA_CRC_vs_H_PAIRED_according_to_DESeq2.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Bacteria<-gsub("_group","",Table_tot$Bacteria)
Table_tot$Bacteria<-gsub("[","",Table_tot$Bacteria, fixed=T)
Table_tot$Bacteria<-gsub("]","",Table_tot$Bacteria, fixed=T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Sample_Type)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Healthy_tissue"="deepskyblue2","CRC"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() 
  ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa", y="Percent Abundance", 
       fill="Condition", x="")
#ggsave(filename = "Results/Every_DA_CRC_vs_H_according_to_DESeq2.png", width = 10, height = 7, dpi=300)
dev.off()

# Redund<- c("Fusobacteriia","Fusobacteriaceae","Campylobacteraceae",
#            "Campylobacterales","Campilobacterota","Fusobacteriota",
#            "Fusobacteriales","Campylobacteria", "Comamonadaceae")
Table_tot2<-subset(Table_tot, Taxa %in% "Genus") # to remove redundant results
#Table_tot2<-subset(Table_tot, !Bacteria %in% Redund) # to remove redundant results
set.seed(1994)
DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Sample_Type)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.9, linewidth= 0.4, alpha= 0.15, 
               outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.3),
             aes(color=Sample_Type), size= 0.4, alpha= 0.65) +  
  theme_classic(base_size = 10) +
  scale_fill_manual(values=c("Healthy_tissue"="deepskyblue2","CRC"="coral")) +
  scale_color_manual(values=c("Healthy_tissue"="deepskyblue2","CRC"="coral")) +
  theme(strip.text.x=element_text(size=11,colour="black")) + 
  guides( fill=guide_legend(nrow=1) , color="none") +
  theme(legend.margin=margin(-15, 0, 0, 0),
        legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=10.5),
        axis.text.x = element_text(angle = 25,
                                   vjust=1, hjust=1,
                                   size=9), 
        axis.text.y = element_text(size=8.5), 
        plot.title= element_text(size=12),
        panel.grid.major.y = element_line(size=0.1, color="gray"),
        plot.margin =  margin(t=10,r=5,b=5, l=10) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,6,2), seq(10,max(Table_tot$Abundance)+5,5))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(y="Percent Abundance", 
       fill="Condition", x="")

DE_plot
ggsave(filename = "Results/DA_CRC_vs_H_PAIRED_according_to_DESeq2.png", width = 5.2, height = 4.2, dpi=300)
dev.off()



######### DA WITH DESEQ2 (CRC vs Healthy Tissue, WITH OFFTARGETS) #####################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_contam_sub<-subset_samples(data_contam, BioProject %in% c("PRJNA507548", "PRJEB57580","PRJNA743150","PRJNA995580", "PRJNA325650") )
data_contam_sub<-subset_samples(data_contam_sub, ! Pair_ID %in% "no_pair" ) 
data_contam_sub<-subset_samples(data_contam_sub, ! Sample_name %in% pairs_unsaturated ) # the pairs of these sample were unsaturated and then removed
data_pruned<- prune_taxa(taxa_sums(data_contam_sub) > 10, data_contam_sub) 
otu_table(data_pruned)<-otu_table(data_pruned)+1 # avoids the error "every gene contains at least one zero"

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sample_Type + Pair_ID)
  DE<-DESeq(DEseq_data) 
  res<-results(DE, contrast= c("Sample_Type", "Healthy_tissue", "CRC"))
  res = res[order(res$padj, na.last=NA), ]
  fc_threshold<- abs(res$log2FoldChange) - abs(res$lfcSE) # minimum fc
  res<-res[(res$padj < 0.05) & fc_threshold>1,]
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Healthy_tissue_vs_CRC.csv"), row.names = F, quote=F, na = "")
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
# obtained this note from DESeq again:
# "note: fitType='parametric', but the dispersion trend was not well captured by the function: y = a/x + b, and a local regression fit was automatically substituted."

#View(Res_tot)
write.csv(Res_tot, file="Results/NO_DA_taxa_according_to_DESeq2_WITH_HumanASVs.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)




######### DA WITH ALDEx2 (CRC vs Healthy Tissue, no human reads) #####################

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
sample_data(data_pruned)[["Sample_Type"]]<-factor(sample_data(data_pruned)[["Sample_Type"]], levels = c("Healthy_tissue","CRC"))
# Healthy_tissue will be the base level --> Denominator of fold change

# removing eventual not duplicated samples
paired_samples <- sample_data(data_sub)[,"Original_Sample_name", drop=T]
paired_samples <- paired_samples$Original_Sample_name
paired_samples<-paired_samples[duplicated(paired_samples)]
data_pruned <- subset_samples(data_pruned, Original_Sample_name %in% paired_samples)

BiocParallel::register(BiocParallel::MulticoreParam(58))   # threads to use (requires BiocParallel package to be installed)
# NB: specify the package from which take these functions is required in our machine!

Table_tot<-NULL
Res_tot<-NULL
### Starting the analysis on each taxon level
for( t in c("Genus","Family","Order") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # using the "manual" pipeline of Aldex2 to better specify the design ...
  mm <- model.matrix(~ Sample_Type + Pair_ID , as(sample_data(d),"data.frame") ) # name of the sample data for last
  set.seed(1994)
  aldx <- aldex.clr(otu_table(d),  mm,
                    mc.samples=128,   #"DMC"
                    # according to ALDEx2 vignette, the number of samples in the smallest group multiplied by the number of DMC be equal at least to 1000
                    # length(which(sample_data(data)[["Sample_Type"]]=="CRC"))
                    denom="all", # every feature as denominator
                    verbose=T,
                    useMC = T
  )
  # aldx@reads # the original counts _ NB: added 0.5 to every number
  
  aldx2 <- aldex.glm(aldx, verbose = T)
  
  #### THEN, the summary of the results
  aldx3 <- aldex.glm.effect(aldx, # calculates the effect size for every binary variable in the mm
                            useMC = T,
                            CI = F # confidence interval of the effect size
  ) 
  aldx_final <- data.frame(aldx2,aldx3) 
  
  p_val<-aldx_final$Sample_TypeCRC.pval
  # p_adj<-p.adjust(p_val, method="holm")
  p_adj<-p.adjust(p_val, method="BH")
  
  #### Moreover, the vignette suggests an effect size cutoff of 1 or greater be used when analyzing HTS 
  # eff_size<-aldx_final$Sample_TypeCRC.effect 
  # Enough_eff_size <- eff_size>1
  res2<-aldx_final
  res2$p_adj_BH<-p_adj
  res2<-res2[res2$p_adj_BH<0.05 , ]
  # res2<-res2[res2$p_adj_BH<0.05 & Enough_eff_size, ]
  results<-row.names(res2)
  if(length(res2[,1])>0){ # if there are results...
    cat(paste( length(length(res2[,1])) ,"results for the",t,"level\n"))
    Sys.sleep(1)
    r<-as.data.frame(res2)
    Taxa.d<-as.data.frame(tax_table(d))
    Taxa.d$Tax_ID<-row.names(Taxa.d)
    r<- cbind.data.frame(Taxa.d[row.names(r), ] , r )
    r<-r[ r[,t] !="none", ] # removing the taxa labeled "none", they are artifact from glomming, composed of different taxa with no name for that level in NCBI
    if( length(r[,1])>0 ){ # if there are still results after removing the "none"
      # assign(paste(t,"results",sep="_"), r)
      r_level<-r
      r_level[, "Taxon"]<- rep(t)
      r_level <- r_level[r_level[,t] != "none" , ]
      Res_tot<-rbind.data.frame(Res_tot,r_level)
      ### single box plots
      target<-r[[t]]
      colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
      target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
      Table_DE<-psmelt(target)
      colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
      # assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
      ### appending to a unique box plot ...
      index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t) # from : to
      index<- index[-length(index)] # removing the last index, regarding the taxa of interest
      Table_DE[,index]<-NULL
      Table_DE$Taxa<-t
      colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
      Table_tot<-rbind.data.frame(Table_tot, Table_DE)
    } # closing the nested "if r > 0" (after the "none" removal)
  } else {  # closing the initial "if r > 0"
    cat("Any results for the",t,"level\n")
    Sys.sleep(1)
  }
}

columns_to_remove<- grepl("Intercept",colnames(Res_tot)) | grepl("Pair_ID",colnames(Res_tot)) | grepl("pval.holm",colnames(Res_tot)) 
colnames(Res_tot)<-gsub("Sample_TypeCRC","Sample_Type_CRC",colnames(Res_tot))
Res_tot<-Res_tot[ , !columns_to_remove ]
#View(Res_tot)
write.csv(Res_tot, file="Results/Every_DA_CRC_vs_H_PAIRED_according_to_Aldex2.csv", row.names = F)


# plot
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Sample_Type)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("CRC"="coral","Healthy_tissue"="deepskyblue")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1, seq(2,max(Table_tot$Abundance+3),1))) +
  labs(title= "Differently abundant taxa", y="Proportional Abundance", 
       fill="Sample_Type", x="")
#ggsave(filename = "Results/DA_Aldex2/DA_Sample_Type_every_result.png", width = 8, height = 7, dpi=300)
dev.off()

# removing redundant results
Table_tot2<-subset(Table_tot, ! Bacteria %in% c("uncultured","Fusobacteriaceae","Tannerellaceae","Leptotrichiaceae","Campylobacteraceae","Fusobacteriales","Erysipelotrichales") )
Table_tot2<-subset(Table_tot2, ! is.na(Table_tot2$Bacteria) )
Table_tot2$Sample_Type <- factor(Table_tot2$Sample_Type, levels = c("CRC","Healthy_tissue") )
Table_tot2$Bacteria<-gsub("[Eubacterium]_hallii_group","Eubacterium_hallii", Table_tot2$Bacteria, fixed = T)
DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Sample_Type)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
  geom_boxplot(width=0.9, linewidth= 0.4, alpha= 0.15, 
               outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Sample_Type), size= 0.02, alpha= 0.3) +
  scale_fill_manual(values=c("CRC"="coral","Healthy_tissue"="deepskyblue")) +
  scale_color_manual(values=c("CRC"="coral","Healthy_tissue"="deepskyblue2")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme_classic(base_size = 10) +
  theme(strip.text.x=element_text(size=11,colour="black"),
        legend.margin=margin(-15, 0, 0, 0),
        legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=10.5),
        axis.text.x = element_text(angle = 35,
                                   vjust=1, hjust=1,
                                   size=7.4), 
        axis.text.y = element_text(size=7), 
        plot.title= element_text(size=12),
        panel.grid.major.y = element_line(size=0.1, color="gray"),
        plot.margin =  margin(t=10,r=5,b=5, l=10)
  ) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  # scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,6,2), seq(10,max(Table_tot$Abundance)+5,5))) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,6,2), seq(10, 65, 5), seq(75, max(Table_tot$Abundance)+5,10) )) +
  labs(y="Percent abundance", 
       fill="Condition", color="Condition",
       x="")

DE_plot
ggsave(filename = "Results/DA_CRC_vs_H_PAIRED_according_to_Aldex2.png", width = 6.5, height = 4.2, dpi=300)
dev.off()

# to avoid re-performing the analysis just to plot it again (Aldex2 takes a whole day!)
write.csv(Table_tot2, file="Results/DA_CRC_vs_H_PAIRED_according_to_Aldex2_TABLE_TO_PLOT.csv", row.names = F)

# system(" echo 'According to many papers (e.g. PMID: 33986544  or  https://doi.org/10.1186/s12859-021-04212-6  or  https://doi.org/10.1186/s12864-022-08803-2 ) many low abundance taxa (e.g. about 0.1%) are potentially false positives (bad calls)! ' > Results/DA_Aldex2/WATCH_OUT_FOR_LOW_ABUNDANCES.txt ")




######### DA WITH ALDEx2 (CRC vs Healthy Tissue, WITH OFFTARGETS) #####################

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_contam_sub<-subset_samples(data_contam, BioProject %in% c("PRJNA507548", "PRJEB57580","PRJNA743150","PRJNA995580", "PRJNA325650") )
data_contam_sub<-subset_samples(data_contam_sub, ! Pair_ID %in% "no_pair" ) 
data_contam_sub<-subset_samples(data_contam_sub, ! Sample_name %in% pairs_unsaturated ) # the pairs of these sample were unsaturated and then removed
data_pruned<- prune_taxa(taxa_sums(data_contam_sub) > 10, data_contam_sub) 
sample_data(data_pruned)[["Sample_Type"]]<-factor(sample_data(data_pruned)[["Sample_Type"]], levels = c("Healthy_tissue","CRC"))
# Healthy_tissue will be the base level --> Denominator of fold change

# removing eventual not duplicated samples
paired_samples <- sample_data(data_contam_sub)[,"Original_Sample_name", drop=T]
paired_samples <- paired_samples$Original_Sample_name
paired_samples<-paired_samples[duplicated(paired_samples)]
data_pruned <- subset_samples(data_pruned, Original_Sample_name %in% paired_samples)

BiocParallel::register(BiocParallel::MulticoreParam(58))   # threads to use (requires BiocParallel package to be installed)
# NB: specify the package from which take these functions is required in our machine!

Table_tot<-NULL
Res_tot<-NULL
### Starting the analysis on each taxon level
for( t in c("Genus","Family","Order") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # using the "manual" pipeline of Aldex2 to better specify the design ...
  mm <- model.matrix(~ Sample_Type + Pair_ID , as(sample_data(d),"data.frame") ) # name of the sample data for last
  set.seed(1994)
  aldx <- aldex.clr(otu_table(d),  mm,
                    mc.samples=128,   #"DMC"
                    # according to ALDEx2 vignette, the number of samples in the smallest group multiplied by the number of DMC be equal at least to 1000
                    # length(which(sample_data(data)[["Sample_Type"]]=="CRC"))
                    denom="all", # every feature as denominator
                    verbose=T,
                    useMC = T
  )
  # aldx@reads # the original counts _ NB: added 0.5 to every number
  
  aldx2 <- aldex.glm(aldx, verbose = T)
  
  #### THEN, the summary of the results
  aldx3 <- aldex.glm.effect(aldx, # calculates the effect size for every binary variable in the mm
                            useMC = T,
                            CI = F # confidence interval of the effect size
  ) 
  aldx_final <- data.frame(aldx2,aldx3) 
  
  p_val<-aldx_final$Sample_TypeCRC.pval
  # p_adj<-p.adjust(p_val, method="holm")
  p_adj<-p.adjust(p_val, method="BH")
  
  #### Moreover, the vignette suggests an effect size cutoff of 1 or greater be used when analyzing HTS 
  # eff_size<-aldx_final$Sample_TypeCRC.effect 
  # Enough_eff_size <- eff_size>1
  res2<-aldx_final
  res2$p_adj_BH<-p_adj
  res2<-res2[res2$p_adj_BH<0.05 , ]
  # res2<-res2[res2$p_adj_BH<0.05 & Enough_eff_size, ]
  results<-row.names(res2)
  if(length(res2[,1])>0){ # if there are results...
    cat(paste( length(length(res2[,1])) ,"results for the",t,"level\n"))
    Sys.sleep(1)
    r<-as.data.frame(res2)
    Taxa.d<-as.data.frame(tax_table(d))
    Taxa.d$Tax_ID<-row.names(Taxa.d)
    r<- cbind.data.frame(Taxa.d[row.names(r), ] , r )
    r<-r[ r[,t] !="none", ] # removing the taxa labeled "none", they are artifact from glomming, composed of different taxa with no name for that level in NCBI
    if( length(r[,1])>0 ){ # if there are still results after removing the "none"
      # assign(paste(t,"results",sep="_"), r)
      r_level<-r
      r_level[, "Taxon"]<- rep(t)
      r_level <- r_level[r_level[,t] != "none" , ]
      Res_tot<-rbind.data.frame(Res_tot,r_level)
      ### single box plots
      target<-r[[t]]
      colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
      target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
      Table_DE<-psmelt(target)
      colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
      # assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
      ### appending to a unique box plot ...
      index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t) # from : to
      index<- index[-length(index)] # removing the last index, regarding the taxa of interest
      Table_DE[,index]<-NULL
      Table_DE$Taxa<-t
      colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
      Table_tot<-rbind.data.frame(Table_tot, Table_DE)
    } # closing the nested "if r > 0" (after the "none" removal)
  } else {  # closing the initial "if r > 0"
    cat("Any results for the",t,"level\n")
    Sys.sleep(1)
  }
}

columns_to_remove<- grepl("Intercept",colnames(Res_tot)) | grepl("Pair_ID",colnames(Res_tot)) | grepl("pval.holm",colnames(Res_tot)) 
colnames(Res_tot)<-gsub("Sample_TypeCRC","Sample_Type_CRC",colnames(Res_tot))
Res_tot<-Res_tot[ , !columns_to_remove ]
#View(Res_tot)
write.csv(Res_tot, file="Results/Every_DA_CRC_vs_H_PAIRED_according_to_Aldex2_WITH_OFFT.csv", row.names = F)


# plot
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Sample_Type)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("CRC"="coral","Healthy_tissue"="deepskyblue")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1, seq(2,max(Table_tot$Abundance+3),1))) +
  labs(title= "Differently abundant taxa", y="Proportional Abundance", 
       fill="Sample_Type", x="")
#ggsave(filename = "Results/DA_Aldex2/DA_Sample_Type_every_result.png", width = 8, height = 7, dpi=300)
dev.off()

# removing redundant results
Table_tot2<-subset(Table_tot, ! Bacteria %in% c("uncultured","Fusobacteriaceae","Tannerellaceae","Leptotrichiaceae","Campylobacteraceae","Fusobacteriales","Erysipelotrichales") )
Table_tot2<-subset(Table_tot2, ! is.na(Table_tot2$Bacteria) )
Table_tot2$Sample_Type <- factor(Table_tot2$Sample_Type, levels = c("CRC","Healthy_tissue") )
Table_tot2$Bacteria<-gsub("[Eubacterium]_hallii_group","Eubacterium_hallii", Table_tot2$Bacteria, fixed = T)
DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Sample_Type)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
  geom_boxplot(width=0.9, linewidth= 0.4, alpha= 0.15, 
               outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Sample_Type), size= 0.02, alpha= 0.3) +
  scale_fill_manual(values=c("CRC"="coral","Healthy_tissue"="deepskyblue")) +
  scale_color_manual(values=c("CRC"="coral","Healthy_tissue"="deepskyblue2")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme_classic(base_size = 10) +
  theme(strip.text.x=element_text(size=11,colour="black"),
        legend.margin=margin(-15, 0, 0, 0),
        legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=10.5),
        axis.text.x = element_text(angle = 35,
                                   vjust=1, hjust=1,
                                   size=7.4), 
        axis.text.y = element_text(size=7), 
        plot.title= element_text(size=12),
        panel.grid.major.y = element_line(size=0.1, color="gray"),
        plot.margin =  margin(t=10,r=5,b=5, l=10)
  ) + 
  scale_x_discrete(expand=c(-0.2, 1)) +
  # scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,6,2), seq(10,max(Table_tot$Abundance)+5,5))) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,6,2), seq(10, 65, 5), seq(75, max(Table_tot$Abundance)+5,10) )) +
  labs(y="Percent abundance", 
       fill="Condition", color="Condition",
       x="")

DE_plot
ggsave(filename = "Results/DA_CRC_vs_H_PAIRED_according_to_Aldex2_WITH_OFFT.png", width = 6.5, height = 4.2, dpi=300)
dev.off()

# to avoid re-performing the analysis just to plot it again (Aldex2 takes a whole day!)
write.csv(Table_tot2, file="Results/DA_CRC_vs_H_PAIRED_according_to_Aldex2_TABLE_WITH_OFFT.csv", row.names = F)

# system(" echo 'According to many papers (e.g. PMID: 33986544  or  https://doi.org/10.1186/s12859-021-04212-6  or  https://doi.org/10.1186/s12864-022-08803-2 ) many low abundance taxa (e.g. about 0.1%) are potentially false positives (bad calls)! ' > Results/DA_Aldex2/WATCH_OUT_FOR_LOW_ABUNDANCES.txt ")




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
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
