##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggvenn")
  library("egg")
  library("ggh4x")
  library("ggbreak") 
  # analysis packages
  library("vegan")
  # utilities
  library("stringr")
  library("xlsx")  
  library("qiime2R")
  #
  library("decontam") # check required by a reviewer
}

{ 
  dir.create("Data_check")
  dir.create("Data_check/PCoA_test")
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Abundances/Extra_subsets/")
  dir.create("Results/Beta_div")
}

options(scipen = 100) # disable scientific annotation


# used in bar plots
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","lightgray") # "others" will be setted as the last one




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_vsearch.qza" )#, tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample
sample<-gsub("00000_ID3003_[0-9][0-9]-","",sample)
sample<-gsub("00000_ID3003_[0-9]-","",sample)
sample<-gsub("-A1-[A-Z][0-9]","",sample)
sample_names(data)<-sample # update

Metadata <- as.data.frame(read.csv(file="Metadata.csv", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name

head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))
sample_data(data)<-Metadata


# settings orders and names
order_temp<-paste0("Blood", sample_data(data)$FASTQ_Order)
order_temp[sample_data(data)$Donor_Condition=="Coli_colture"]<-gsub("Blood","E_Coli",order_temp[sample_data(data)$Donor_Condition=="Coli_colture"])
order_temp[sample_data(data)$Donor_Condition=="Mice_cell"]<-gsub("Blood","Mice_DNA",order_temp[sample_data(data)$Donor_Condition=="Mice_cell"])
Names_ordered<-c( order_temp[! grepl("[0-9][0-9]",order_temp) & grepl("Blood",order_temp) ] ,
                  order_temp[ grepl("[0-9][0-9]",order_temp) & grepl("Blood",order_temp) ],
                  order_temp[ grepl("Coli",order_temp) | grepl("Mice",order_temp)] ) 
Names_ordered[grepl("d[0-9]$", Names_ordered)]<- sort(Names_ordered[grepl("d[0-9]$", Names_ordered)])
sample_data(data)$Names_ordered<-factor(order_temp, levels = Names_ordered)
sample_names(data)<-sample_data(data)$Names_ordered
sample_data(data)$Sample_Type<-factor(sample_data(data)$Sample_Type, levels=c("Contaminant_trap","Healthy")) # decides the order in plots

# View(sample_data(data))



########## IMPORTING ALSO THE CONTAMINATED ONE (FROM HOST DNA)

# the following method allows to avoid the already perform (slow) re-aligning to SILVA 16S with VSearch, just by adding the "unassigned" rows of the host ASVs removed during the processing

otu_table_contam<-qza_to_phyloseq(features="QIIME/table.qza")
# adding the removed off-targets ASVs to this tax table to maintein them in the related phyloseq object
off_targets_tax<-read.table(file="QIIME/eucar_contaminants_ASV/Every_OFF_TARGETS_in_ASVs.txt")
off_targets_tax<-as.character(off_targets_tax[grepl(">",off_targets_tax$V1), ]) # only identifiers
off_targets_tax<-gsub(">","",off_targets_tax)
fake_rows_off_targets<-c("Unassigned",NA,NA,NA,NA,NA,NA) # this exact output would be obtained if those ASVs are matched in SILVA 16S with Vsearch
temp<-matrix(nrow = length(off_targets_tax), ncol = length(fake_rows_off_targets))
row.names(temp)<-off_targets_tax
colnames(temp)<-colnames(tax_table(data))
for(i in row.names(temp)){
  temp[i, ]<-fake_rows_off_targets
}
tax_table_contam<-as.matrix(rbind(tax_table(data),temp))
if( length(row.names(otu_table_contam))==length(row.names(tax_table_contam)) ){
  cat("\nSame number of rows --> 'host contaminated' tax table completed! \n\n")
  data_contam<-phyloseq(otu_table(otu_table_contam, taxa_are_rows = T), tax_table(tax_table_contam))
  rm(temp, tax_table_contam, otu_table_contam)
} else {
  stop("\nAn error has occurred during the building of contaminanted tax table !!! \n\n")
}


if( identical(original_names,sample_names(data_contam))){
  sample_names(data_contam)<-sample # original modified names
  if( identical(sample_data(data)$FASTQ_ID,sample_names(data_contam))){
    sample_names(data_contam)<-sample_names(data)
    sample_data(data_contam)<-sample_data(data)
    cat("\nEverything looks fine!\n\n")
  }
}

write.csv2(cbind(otu_table(data_contam),tax_table(data_contam)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


suppressWarnings( rm(order_temp, original_length, original_names, sample, fake_rows_off_targets, i ) )

# save.image("data_backup_before_filters.RData")




#################### CHECKING THE HOST CONTAMINATION PROPORTION ###################

data_contam_k<- tax_glom(data_contam, taxrank = "Kingdom", NArm = F) # without Bowtie cleaning
data_contam_k<- transform_sample_counts(data_contam_k, function(x) (x/sum(x))*100)
data_contam_k<- subset_taxa(data_contam_k, is.na(tax_table(data_contam_k)[,"Kingdom"]) | as.character(tax_table(data_contam_k)[,"Kingdom"])=="Unassigned")
write.csv2(file="Data_check/Host_Contamin_proportion_in_raw_FASTQ.csv", cbind.data.frame(tax_table(data_contam_k)[,c("Kingdom","Phylum")],otu_table(data_contam_k) ))

if("Unassigned" %in% as.character (tax_table(data)[,"Kingdom"])){
  cat ("\nAlso in the decontaminated dataset there are still reads unassigned at Domain level! --> Removing them... \n\n")
}
Domain_clean_proof<-"Proof of performing this step, required for the script"

rm(data_contam_k)


#################### UPDATING THE TAXONOMY TABLE #################
# adding details to "uncultured" or "NA" taxa names to avoid confusing them during the filters

if(!"update_proof" %in% ls()){
  data.temp<-tax_glom(data, taxrank = "Genus")
  not_updated_g_taxonomy<-as.data.frame(tax_table(data.temp)) # required later for checking the number of NA in SILVA
  rm(data.temp)
}

taxa_temp<-as.data.frame(tax_table(data_contam))
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ c",taxa_temp[which(taxa_temp$Genus=="NA_ o NA")[1],"Class"])}
  tax_table(data_contam)<-as.matrix(taxa_temp)
}
taxa_temp<-as.data.frame(tax_table(data))
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ c",taxa_temp[which(taxa_temp$Genus=="NA_ o NA")[1],"Class"])}
  tax_table(data)<-as.matrix(taxa_temp)
}

update_proof<-"This objects allows the filter step"



########################## EXTRA: CHECK WITH DECONTAM ###############################

# NB: It is better to perform this check BEFORE using the actual filter (rel abundance)

# frequency method
sample_data(data)$ng.ul  # column with the concentration after the extraction
contam.test <- isContaminant(data, 
                             normalize=TRUE, 
                             method="frequency", 
                             threshold=0.3,
                             conc= "ng.ul",
                             #neg="Control"
                             )
Only.contaminants <- contam.test[which(contam.test$contaminant==TRUE), ]
decont_freq<-unique(tax_table(data)[row.names(Only.contaminants), "Genus"])

# prevalence method
sample_data(data)$Control <- sample_data(data)$Sample_Type=="Contaminant_trap"  # logical vector
contam.test <- isContaminant(data, 
                             normalize=TRUE, 
                             method="prevalence", 
                             threshold=0.3,
                             #conc= "ng.ul",
                             neg="Control"
)
Only.contaminants <- contam.test[which(contam.test$contaminant==TRUE), ]
decont_prev<-unique(tax_table(data)[row.names(Only.contaminants), "Genus"])
not_contam_according_to_dec_LIST<-as.character(unique(tax_table(subset_taxa(data, ! Genus %in%  as.character(decont_prev) ))[, "Genus"]))

con <- file("Data_check/EXTRA_CHECK_WITH_DECONTAM.txt")
sink(con, append = TRUE)
cat("Contaminants according to Decontam frequency method", fill=T)
cat(decont_freq)
cat("\n", "\n", fill=TRUE)
cat("Contaminants according to Decontam prevalence method (threshold 0.3)", fill=T)
cat(as.character(decont_prev))
cat("\n", "\n", fill=TRUE)
cat("Genera remaining if the above contaminants are excluded", fill=T)
cat(not_contam_according_to_dec_LIST)
cat("\n \nNB: The remaining genera are Burkholderia, Leifsonia, Flavobacterium... and other genera which are VERY abundant in controls... \nThere is also Pseudomonas, featured in controls and renowned as potential contaminant... and also Helicobacter (!)... this method does not appear to work properly on this dataset, as expected")
sink()
close(con)
suppressWarnings(rm(con))




#################### FILTERING NOISES FROM THE 'TRUE BLOOD' DATA SET ####################

if(!"update_proof" %in% ls()){
  stop("\n\n   Update the genus names before going further!\n\n")
}

# using the data *with* Mice DNA to better estimate the contaminants through proportions
filtering_blood_data<-subset_samples(data_contam, Sample_Type=="Healthy") 
filtering_blood_data_g<-tax_glom(filtering_blood_data, taxrank = "Genus", NArm = F) # NB: created WHILE the object is STILL contaminated!


###### cutting under 0.01% to remove noises/contaminants, too conservative but also safe cutoff, see  PMID: 31164452 and   DOI: 10.1128/mSystems.00290-19
# moreover, in PMC4882852 showed that the contaminants in blood DNA extractions kit may have an abundances far beyond 0.1% (not 0.01) !
suppressWarnings(rm(data.genus.temp))
data.genus.temp<-transform_sample_counts(filtering_blood_data_g, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.01, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_001_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.01, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured" & !is.na(filtered)] # to avoid the removal of other eventual uncultured genera (missed by the tax table update)
data<-subset_taxa(data, ! Genus %in% filtered)
rm(filtered, data.genus.temp)

 
###### Applying also prevalence filter: at least 1 read of that genus at least in 3 "true" blood samples (on 11)
data.temp<-filtering_blood_data_g
who<-as.data.frame(otu_table(data.temp) )
who<-apply(who, MARGIN=2, function(x) ifelse(x>1, 1, 0)) # if at least one read in the sample then it gets "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
who<-who[!is.na(who) & who!="uncultured"] # otherwise EVERY other NA and uncultured would be also eliminated
who # just a check --> ok!
write.csv(who, "Data_check/Filtered_genera_that_were_NOT_AT_LEAST_in_THREE_blood_sample.csv", row.names = F)
data <-subset_taxa(data, ! Genus %in% who)
rm(who, data.temp)


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
head(tax_table(data))
write.csv2(tax_table(data), file="Data_check/Every_filtered_ASV_and_taxonomic_assignment.csv", row.names = T)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


# removing chloroplast and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast...\n\n")
  data<-subset_taxa(data, Genus !="Chloroplast")
}
if("Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Mitochondria...\n\n")
  data<-subset_taxa(data, Genus !="Mitochondria")
}

proof2<- "Marker of the filtering, it is required for the script"



############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

if(! "proof2" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
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


### importing the percentuals of contaminants (for the plot below)
percent_fastq_human_count<-read.csv("QIIME/Host_off_target_counts/HUMAN_off_target_for_each_FASTQ.csv", header = F, row.names = 1)
percent_fastq_mice_count<-read.csv("QIIME/Host_off_target_counts/MICE_off_target_for_each_FASTQ.csv", header = F, row.names = 1)
# changing the mice rows in the output of GRCh39
percent_fastq_human_count[row.names(percent_fastq_mice_count), "V2"] <- percent_fastq_mice_count$V2

# reformatting
percent_fastq<-data.frame("FASTQ_Names"=row.names(percent_fastq_human_count),
                          "percentual"=percent_fastq_human_count$V2 )
rm(list=c("percent_fastq_mice_count", "percent_fastq_human_count"))
percent_fastq$percentual<-gsub("% overall alignment rate","",percent_fastq$percentual)
percent_fastq$percentual<-gsub(" ","",percent_fastq$percentual)
percent_fastq$percentual<-as.numeric(percent_fastq$percentual)
percent_fastq$FASTQ_Names<-gsub("_S1_L001_R[0-9]_001.fastq.gz ","",percent_fastq$FASTQ_Names)
percent_sample<-NULL
row<-NULL
for(i in unique(percent_fastq$FASTQ_Names) ){
  target_percent<-percent_fastq[percent_fastq$FASTQ_Names==i, "percentual"]
  row<-cbind(sample=i, percentual=as.numeric(mean(target_percent[1],target_percent[2])) )
  percent_sample<-rbind.data.frame(percent_sample, row)
}
percent_sample<-as.data.frame(percent_sample)
temp<-sample_data(data)
row.names(temp)<-temp$FASTQ_original_names
percent_sample$percentual<-round(as.numeric(percent_sample$percentual), 1)
percent_sample$Names_ordered<-temp[percent_sample$sample, ]$Names_ordered
percent_sample$Sample_code<-temp[percent_sample$sample, ]$Sample_code
percent_sample$Kit<-temp[percent_sample$sample, "Kit"]
percent_sample$Subject<-temp[percent_sample$sample, ][["Donor_Condition"]] 
percent_sample$Subject<-gsub("Healthy","Subject_", percent_sample$Subject)
percent_sample$Subject<-gsub("Coli_colture","ENC", percent_sample$Subject)
percent_sample$Subject<-gsub("Mice_cell","OTC", percent_sample$Subject)
rm(temp)


### barplot of reads
Original_read_number<-as.data.frame(read.table(file="QIIME/Original_number_of_reads_for_sample.tsv", header = T, sep="\t", row.names = 1))
DADA2_read_number<-as.data.frame(read.table(file="QIIME/DADA2_Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
after_filter_number<-as.data.frame(colSums(otu_table(data)))
# updating names from FASTQ codes to sample names
temp<-sample_data(data)
row.names(temp)<-temp$FASTQ_original_names
row.names(Original_read_number)<-as(temp[row.names(Original_read_number),"Names_ordered"],"vector")[[1]] # this [[1]] is because R still belives is a factor despite the coercition
row.names(DADA2_read_number)<-as(temp[row.names(DADA2_read_number),"Names_ordered"],"vector")[[1]] # this [[1]] is because R still belives is a factor despite the coercition
# same order of rows
Original_read_number<-Original_read_number[row.names(DADA2_read_number), ]
after_filter_number<-after_filter_number[row.names(DADA2_read_number), ]
row.names(percent_sample)<-percent_sample$Names_ordered
percent_sample<-percent_sample[row.names(DADA2_read_number), ]
# creating the table
table<-cbind.data.frame(Original=Original_read_number$forward.sequence.count, # same as reverse
                        percentual_host=percent_sample$percentual,
                        After_quality_filter=DADA2_read_number$non.chimeric,
                        After_contam_filter=after_filter_number,
                        Kit=percent_sample$Kit,
                        Subject=percent_sample$Sample_code)
table$Samples<-row.names(Original_read_number)
# computing the abundances of off-targets from the original number of reads
table$abundance_host<-0
for (i in 1:length(table$Original)) {
  table$abundance_host[i]<-round( as.numeric(table[i, "After_quality_filter" ])*as.numeric(table[i, "percentual_host" ])/100 , 0) 
}

# plotting
table$percentual_host<-paste0(table$percentual_host,"%") # now it can be numeric no more!
table$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",table$Kit)
table$Kit<-gsub("Microbiome_kit"," QIAamp® DNA Microbiome Kit",table$Kit)
table$Subject<-gsub("Coli_colture","ENC",table$Subject)
table$Subject<-gsub("Mice_cell","OTC",table$Subject)
ggplot(aes(x=Subject, y=Original), data=table) +
  facet_grid( ~ Kit, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.6) +
  geom_bar( aes(y= After_quality_filter, fill="After quality filters"), width = 0.75, alpha=0.75, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="After contaminant filters"), width= 0.5, alpha= 1, stat = "identity") +
  geom_text( aes( label= percentual_host ), vjust= 0, hjust=0.47, size= 2) +
  theme_classic( base_size = 9) +
  theme(axis.text.x = element_text(size=6.5, angle = 25, vjust=1, hjust=1),
        axis.text.y = element_text(size=5),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6.5),
        panel.grid.major.y = element_line(size=0.05, color="grey"),
        legend.position = "right",
        legend.margin = margin(0, -2, 0, -5),
        legend.key.height = unit(0.4,"cm")
        ) +
  scale_fill_manual(name='Number of reads',
                     breaks=c('Original number', 
                              'After quality filters',
                              'After contaminant filters'),
                     values=c('Original number'='green3',
                              'After quality filters'='coral',
                              "After contaminant filters"='red3')) +
  scale_y_continuous(breaks = c(10000, seq(0, max(table$Original), 25000))) +
  labs(y="Read abundance",
       caption="The percentual of host off-targets reads are reported on the top of each bar",
       x="Sample"
       )
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 5.8, height = 3.3, dpi=300)
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters_pdf.pdf", width = 5.8, height = 3.3, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)
  

# to compute percentages
table$Abundance<-table$After_quality_filter-table$abundance_host # residual before abundances filters
target<- !(table$Samples%in%c("E_Coli3","E_Coli4")) & !grepl("DNeasy® Blood & Tissue Kit",table$Kit)
target<- !(table$Samples%in%c("E_Coli3","E_Coli4"))
(sum( table$Abundance[ target ] ) / sum( table$Original[ target ] )) * 100



################# CHECKING THE COMPOSITION AFTER FILTERING ####################

#data.unf.prop<-transform_sample_counts(filtering_blood_data, function(x) (x/sum(x))*100) # if BOWTIE STEP IS NOT PERFORMED
data.unf.prop<-transform_sample_counts(data_contam, function(x) (x/sum(x))*100)


if(! "proof2" %in% ls()){
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
p1<-plot_ordination(data.sqrt_prop, ordBC,
                    color = "Sample_Type"
                    ) +
  scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) +
  geom_text(aes(label=Sample_Type),
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data",
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=Sample_Type), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

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
p1<-plot_ordination(data.prop.labels, ordBC, color = "Sample_Type") +
  scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=Sample_Type), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA sqrt Bray-Curtis (on proportional ASV)\n\n UNfiltered data", 
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.prop.labels))
data.prop.labels<-data.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "bray")
  DistBC = sqrt(DistBC)
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Sample_Type") +
  scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=Sample_Type), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=Sample_Type), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=Sample_Type), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


# ################# wUNIFRAC PROP
# 
# ##### unfiltered wUNIFRAC
# suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
# data.prop.labels<-data.unf.prop
# #{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
# {DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
#   ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# p1<-plot_ordination(data.prop.labels, ordBC, color = "Sample_Type") +
#   scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
#   guides(color="none", shape="none") +
#   geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) +
#   geom_text(aes(label=Sample_Type),
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data",
#        color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
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
#   scale_color_manual(values=c("Healthy"="coral","Contaminant_trap"="gray")) +
#   geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) +
#   geom_text(aes(label= Sample_Type),
#             color="black", size=2.5, show.legend = FALSE) +
#   labs(title="\n \n filtered data",
#        color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# 
# png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_test.png", width = 3200, height = 1800, res=300)
# ggarrange(p1,p2, nrow = 1)
# dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))




####################### PREPARATION OF THE DATA #######################

if(!"update_proof" %in% ls() | ! "proof2" %in% ls() ){
  stop("\n\n   Wait! Perform the filter steps before\n\n")
}

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

rm(taxa_temp)




#################### % ASSIGNED IN SILVA #########################

not_updated_g_taxonomy<-not_updated_g_taxonomy[taxa_names(data.genus), ] # this table has to be created before the filters --> filtering it now using ASVs codes

{a<-cbind(length(not_updated_g_taxonomy$Genus),length(which(!is.na(not_updated_g_taxonomy$Genus))),length(which(!is.na(not_updated_g_taxonomy$Genus)))/length(not_updated_g_taxonomy$Genus),"Genus")
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



########################### COUNTS EXPORT ##########################################

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
  }

options(scipen = 100)
dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}

# for the supplementary section
data.temp<-data.genus.prop
sample_data(data.temp)$Sample_code<-gsub("Mice_cell","OTC",sample_data(data.temp)$Sample_code)
sample_data(data.temp)$Sample_code<-gsub("Coli_colture","ENC",sample_data(data.temp)$Sample_code)
sample_names(data.temp)<-sample_data(data.temp)$Sample_code
write.xlsx(cbind(as(otu_table(data.temp),"matrix"),as(tax_table(data.temp),"matrix")),file="Results/Abundances/Relative_abundances/Table_supplementary_genus_abundances.xlsx",showNA = F, col.names = T, row.names = F)


############### TAXONOMIC PRESENCE/ABSENCE AMONG SAMPLES ####################

### two negative controls as E.coli e Mice
ENC_coli<-subset_samples(data.genus, Donor_Condition=="Coli_colture")
ENC_coli<-as.data.frame(tax_table( prune_taxa(taxa_sums(ENC_coli)>0, ENC_coli) ))
ENC_mice<-subset_samples(data.genus, Donor_Condition=="Mice_cell")
ENC_mice<-as.data.frame(tax_table( prune_taxa(taxa_sums(ENC_mice)>0, ENC_mice) ))
blood_micr_kit<-subset_samples(data.genus, Sample_Type=="Healthy" & Kit=="Microbiome_kit")
blood_micr_kit<-as.data.frame(tax_table( prune_taxa(taxa_sums(blood_micr_kit)>0, blood_micr_kit) ))
blood_blood_kit<-subset_samples(data.genus, Sample_Type=="Healthy" & Kit=="Blood_Tissue_kit")
blood_blood_kit<-as.data.frame(tax_table( prune_taxa(taxa_sums(blood_blood_kit)>0, blood_blood_kit) ))

complete_taxa<-as.data.frame(tax_table(data.genus))
complete_taxa<-complete_taxa[ , c("Phylum","Genus")]
complete_taxa$'ENC E.coli'<-ifelse(complete_taxa$Genus %in% ENC_coli$Genus , "*" ,"")
complete_taxa$'ENC mice DNA'<-ifelse(complete_taxa$Genus %in% ENC_mice$Genus , "*" ,"")
complete_taxa$'Blood_Microbiome kit'<-ifelse(complete_taxa$Genus %in% blood_micr_kit$Genus , "*" ,"")
complete_taxa$'Blood_Blood & Tissue Kit'<-ifelse(complete_taxa$Genus %in% blood_blood_kit$Genus , "*" ,"")

complete_taxa$Genus<-gsub("[","",complete_taxa$Genus, fixed = T)
complete_taxa$Genus<-gsub("]","",complete_taxa$Genus, fixed = T)
complete_taxa$Genus<-gsub("_group","",complete_taxa$Genus, fixed = T)

write.csv2(complete_taxa, file="Results/Report_of_presence_absence_among_samples_with_Ecoli_and_Mice.csv", row.names = F, quote=F)

rm(complete_taxa,blood_micr_kit, blood_blood_kit, ENC_coli, ENC_mice)


### two negative controls as neg Kit1 and neg Kit2

ENC_BloodT_K<-subset_samples(data.genus, Donor_Condition %in% c("Coli_colture","Mice_cell") & Kit=="Blood_Tissue_kit")
ENC_BloodT_K<-as.data.frame(tax_table( prune_taxa(taxa_sums(ENC_BloodT_K)>0, ENC_BloodT_K) ))
ENC_Microbiome_kit<-subset_samples(data.genus, Donor_Condition %in% c("Coli_colture","Mice_cell") & Kit=="Microbiome_kit")
ENC_Microbiome_kit<-as.data.frame(tax_table( prune_taxa(taxa_sums(ENC_Microbiome_kit)>0, ENC_Microbiome_kit) ))
blood_micr_kit<-subset_samples(data.genus, Sample_Type=="Healthy" & Kit=="Microbiome_kit")
blood_micr_kit<-as.data.frame(tax_table( prune_taxa(taxa_sums(blood_micr_kit)>0, blood_micr_kit) ))
blood_blood_kit<-subset_samples(data.genus, Sample_Type=="Healthy" & Kit=="Blood_Tissue_kit")
blood_blood_kit<-as.data.frame(tax_table( prune_taxa(taxa_sums(blood_blood_kit)>0, blood_blood_kit) ))

complete_taxa<-as.data.frame(tax_table(data.genus))
complete_taxa<-complete_taxa[ , c("Phylum","Genus")]
complete_taxa$'Neg Ctrls Blood_Tissue kit'<-ifelse(complete_taxa$Genus %in% ENC_BloodT_K$Genus , "*" ,"")
complete_taxa$'Neg Ctrls Microbiome kit'<-ifelse(complete_taxa$Genus %in% ENC_Microbiome_kit$Genus , "*" ,"")
complete_taxa$'Blood_Microbiome kit'<-ifelse(complete_taxa$Genus %in% blood_micr_kit$Genus , "*" ,"")
complete_taxa$'Blood_Blood & Tissue Kit'<-ifelse(complete_taxa$Genus %in% blood_blood_kit$Genus , "*" ,"")

complete_taxa$Genus<-gsub("[","",complete_taxa$Genus, fixed = T)
complete_taxa$Genus<-gsub("]","",complete_taxa$Genus, fixed = T)
complete_taxa$Genus<-gsub("_group","",complete_taxa$Genus, fixed = T)

write.csv2(complete_taxa, file="Results/Report_of_presence_absence_among_samples_with_Kit_Neg_Ctrls.csv", row.names = F, quote=F)

rm(complete_taxa,blood_micr_kit, blood_blood_kit, ENC_BloodT_K, ENC_Microbiome_kit)




###################### ABUNDANCES BAR PLOTS ##########################


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
tabella$Sample_Type<-gsub("Contaminant_trap","Controls",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
# tabella$Donor_Condition<-gsub("Healthy","Subject_", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Coli_colture","ENC2", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Mice_cell","ENC1", tabella$Donor_Condition)
tabella$Sample_code<-gsub("Mice_cell","OTC",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Phylum)) +
  facet_grid2( Kit ~ Sample_Type, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack", size=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,10,-16),
        legend.position="bottom",
        plot.margin = margin(3,1,3,-7)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant phyla", 
       caption = "'Others' includes every phylum below rank 5 ") #+
  #coord_flip() # usefull to flip everything with just a click!
#ggsave(file="Results/Abundances/TOP_5_phyla.png",width=5.85,height=3.7, dpi=300)
ggsave(file="Results/Abundances/TOP_5_phyla_v2.png",width=3.9,height=5.8, dpi=300)
#ggsave(file="Results/Abundances/TOP_5_phyla.pdf",width=5.85,height=3.7, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


### TOP 5 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_top$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.-Paraburk.", tabella_top$Genus)
  tabella_top$Genus<-gsub("Escherichia-Shigella","Escherichia", tabella_top$Genus)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_Type<-gsub("Contaminant_trap","Controls",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
# tabella$Donor_Condition<-gsub("Healthy","Subject_", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Coli_colture","ENC2", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Mice_cell","ENC1", tabella$Donor_Condition)
tabella$Sample_code<-gsub("Mice_cell","OTC",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Genus)) +
  facet_grid2( Kit ~ Sample_Type, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack", size=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,3,-9),
        legend.position="bottom",
        plot.margin = margin(3,2,3,-7)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant genera",
       caption = " 'Others' includes every genus below rank 5 ") #+
  # coord_flip() # usefull to flip everything with just a click!
#ggsave(file="Results/Abundances/TOP_5_genera.png",width=5.85,height=3.7, dpi=300)
ggsave(file="Results/Abundances/TOP_5_genera_v2.png",width=3.9,height=5.8, dpi=300)
# ggsave(file="Results/Abundances/TOP_5_genera_pdf.pdf",width=5.85,height=3.7, dpi=300)
dev.off()

rm(top, tabella, tabella_top)

### means of TOP5 genera
write.xlsx(file = "Results/Abundances/TOP_5_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




###################### ABUNDANCES BAR PLOTS (WITHOUT CONTROLS IN THE CALCULATION) ##########################


### TOP 5 Phyla
suppressWarnings(rm(top, others, tabella))
trueblood<-subset_samples(data.phy, ! Donor_Condition %in% c("Coli_colture","Mice_cell") )
trueblood<-transform_sample_counts(trueblood, function(x) (x/sum(x))*100 )
{top <- names(sort(taxa_sums(trueblood), decreasing=TRUE))[1:5]
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
tabella$Sample_Type<-gsub("Contaminant_trap","Controls",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
tabella$Sample_code<-gsub("Mice_cell","OTC",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Phylum)) +
  facet_grid2( Sample_Type ~ Kit, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack", size=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,3,-25),
        legend.position="bottom",
        plot.margin = margin(3,1,3,-7)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant genera in blood samples", 
       caption = "'Others' includes every phylum below rank 5 ") +
  coord_flip() # usefull to flip everything with just a click!
# ggsave(file="Results/Abundances/Extra_subsets/TOP_5_phyla_WITHOUT_Ctrls_in_the_calculation.png",width=5.85,height=2.9, dpi=300)
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_phyla_WITHOUT_Ctrls_in_the_calculation.png",width=5.85,height=3.7, dpi=300)
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_phyla_WITHOUT_Ctrls_in_the_calculation_pdf.pdf",width=5.85,height=3.7, dpi=300)
dev.off()

write.xlsx(file = "Results/Abundances/Extra_subsets/TOP_5_Phyla_WITHOUT_Ctrls_in_the_calculation.xlsx", row.names = F,
           cbind.data.frame("mean (no_ctrls_included)"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, trueblood)


### TOP 5 Genera
trueblood<-subset_samples(data.genus, ! Donor_Condition %in% c("Coli_colture","Mice_cell") )
trueblood<-transform_sample_counts(trueblood, function(x) (x/sum(x))*100 )
{top <- names(sort(taxa_sums(trueblood), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.-Paraburk.", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_Type<-gsub("Contaminant_trap","Controls",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
tabella$Sample_code<-gsub("Mice_cell","OTC",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Genus)) +
  facet_grid2( Sample_Type ~ Kit, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack") + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,3,-10),
        legend.position="bottom",
        plot.margin = margin(3,2,3,-7)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant genera in blood samples",
       caption = " 'Others' includes every genus below rank 5 ") +
  coord_flip() # usefull to flip everything with just a click!
# ggsave(file="Results/Abundances/Extra_subsets/TOP_5_genera_WITHOUT_Ctrls_in_calculation.png",width=5.85,height=2.9, dpi=300)
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_genera_WITHOUT_Ctrls_in_calculation.png",width=5.85,height=3.7, dpi=300)
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_genera_WITHOUT_Ctrls_in_calculation_pdf.pdf",width=5.85,height=3.7, dpi=300)
dev.off()

write.xlsx(file = "Results/Abundances/Extra_subsets/TOP_5_Genera_WITHOUT_Ctrls_in_calculation.xlsx", row.names = F,
           cbind.data.frame("mean (no_ctrls_included)"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, trueblood))



###################### ABUNDANCES BAR PLOTS WITH ONLY CONTROLS ##########################


### TOP 5 Phyla
suppressWarnings(rm(top, others, tabella))
noblood<-subset_samples(data.phy, Donor_Condition %in% c("Coli_colture","Mice_cell") )
noblood<-transform_sample_counts(noblood, function(x) (x/sum(x))*100 )
{top <- names(sort(taxa_sums(noblood), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,noblood)
  others<-taxa_names(noblood)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,noblood)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Sample_Type<-gsub("Contaminant_trap","ENC",tabella$Sample_Type)
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
tabella$Sample_code<-gsub("Mice_cell","ENC1",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC2",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Phylum)) +
  facet_grid2( Sample_Type ~ Kit, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack", size=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,3,-25),
        legend.position="bottom",
        plot.margin = margin(3,1,3,-7)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant genera in control samples", 
       caption = "'Others' includes every phylum below rank 5 ") +
  coord_flip() # usefull to flip everything with just a click!
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_phyla_ONLY_WITH_Ctrls.png",width=5.85,height=2.9, dpi=300)
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_phyla_ONLY_WITH_Ctrls_pdf.pdf",width=5.85,height=2.9, dpi=300)
dev.off()

write.xlsx(file = "Results/Abundances/Extra_subsets/TOP_5_Phyla_ONLY_WITH_Ctrls.xlsx", row.names = F,
           cbind.data.frame("mean (no_blood_included)"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, noblood)


### TOP 5 Genera
noblood<-subset_samples(data.genus, Donor_Condition %in% c("Coli_colture","Mice_cell") )
noblood<-transform_sample_counts(noblood, function(x) (x/sum(x))*100 )
{top <- names(sort(taxa_sums(noblood), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,noblood)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(noblood)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,noblood)
  tabella_top<-psmelt(prune.dat_top)
  tabella_top$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.-Paraburk.", tabella_top$Genus)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_Type<-gsub("Contaminant_trap","ENC",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
tabella$Sample_code<-gsub("Mice_cell","ENC1",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC2",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Genus)) +
  facet_grid2( Sample_Type ~ Kit, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack") + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,3,-10),
        legend.position="bottom",
        plot.margin = margin(3,2,3,-7)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant genera in control samples",
       caption = " 'Others' includes every genus below rank 5 ") +
  coord_flip() # usefull to flip everything with just a click!
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_genera_ONLY_WITH_Ctrls.png",width=5.85,height=2.9, dpi=300)
ggsave(file="Results/Abundances/Extra_subsets/TOP_5_genera_ONLY_WITH_Ctrls_pdf.pdf",width=5.85,height=2.9, dpi=300)
dev.off()

write.xlsx(file = "Results/Abundances/Extra_subsets/TOP_5_Genera_ONLY_WITH_Ctrls.xlsx", row.names = F,
           cbind.data.frame("mean (no_blood_included)"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, noblood))



########################### PCoA  #########################

### on Genera
data.prop.labels<-data.genus.prop
{
sample_data(data.prop.labels)$Donor_Condition2<-gsub("Healthy","",sample_data(data.prop.labels)$Donor_Condition)       
sample_data(data.prop.labels)$Donor_Condition2<-gsub("Coli_colture","E.coli",sample_data(data.prop.labels)$Donor_Condition2)       
sample_data(data.prop.labels)$Donor_Condition2<-gsub("_cell","",sample_data(data.prop.labels)$Donor_Condition2)
sample_data(data.prop.labels)$Donor_Condition3<-gsub("E.coli","ENC",sample_data(data.prop.labels)$Donor_Condition2)       
sample_data(data.prop.labels)$Donor_Condition3<-gsub("Mice","OTC",sample_data(data.prop.labels)$Donor_Condition3)   
sample_data(data.prop.labels)$Sample_Type<-gsub("Contaminant_trap","ENC",sample_data(data.prop.labels)$Sample_Type)
sample_data(data.prop.labels)$Sample_Type<-gsub("Healthy","Healthy blood",sample_data(data.prop.labels)$Sample_Type)
sample_data(data.prop.labels)$Kit<-gsub("Microbiome_kit","DNA Microbiome Kit",sample_data(data.prop.labels)$Kit) 
sample_data(data.prop.labels)$Kit<-gsub("Blood_Tissue_kit","Blood & Tissue Kit ",sample_data(data.prop.labels)$Kit)  
}
# Donor condition 2 = Mice and Ecoli    Donor Condition 3 = ENC e OTC     

{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type",
                shape = "Kit") +
  stat_ellipse(size=0.15) +
  scale_color_manual(values=c("Healthy blood"="coral","ENC"="gray")) +
  geom_point(size=5, alpha= 0.2) +
  geom_point(size=3.5, alpha= 0.6) +
  theme_classic(base_size = 10) +
  theme(legend.margin = margin(0,0,0,-2),
        legend.text = element_text(size=8.2) 
  ) +
  geom_text(aes(label=Donor_Condition2),
            color="black", size=2.3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
    color="Sample Type", x=paste("PC1: ",eigval[1],"% variation"),
    y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 5, height = 3.6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type",
                shape = "Kit") +
  scale_color_manual(values=c("Healthy blood"="coral","ENC"="gray")) +
  geom_point(size=5.4, alpha= 0.25) +
  geom_point(size=3.2, alpha= 0.88) +
  theme_classic(base_size = 10) +
  theme(legend.margin = margin(0,0,0,-2),
        legend.text = element_text(size=8.2) 
        ) +
  geom_text(aes(label=Donor_Condition3),
            color="black", size=2.3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       color="Sample Type", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_Genera_no_ellipse.png", width = 5, height = 3.6, dpi=300)
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_Genera_no_ellipse_pdf.pdf", width = 5, height = 3.6, dpi=300)
# without ellipses no names
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type", shape = "Kit") +
  scale_color_manual(values=c("Healthy blood"="coral","ENC"="gray")) +
  geom_point(size=5, alpha= 0.2) +
  geom_point(size=3.5, alpha= 0.6) +
  geom_point(size=0.5, alpha= 1) +
  theme_classic(base_size = 10) +
  theme(legend.margin = margin(0,0,0,-2)) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       color="Sample_Type", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_Genera_no_ellipse_no_names.png", width = 5, height = 3.6, dpi=300)


suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



#################### VEEN DIAGRAM : HEALTHY vs TRAPS #########################

data.venn<-data.genus.prop

Contaminant_trap<-subset_samples(data.venn, Sample_Type=="Contaminant_trap")
Contaminant_trap<-as.character(tax_table(prune_taxa(taxa_sums(Contaminant_trap)>0, Contaminant_trap))[,"Genus"])

Healthy<-subset_samples(data.venn, Sample_Type=="Healthy")
Healthy<-as.character(tax_table(prune_taxa(taxa_sums(Healthy)>0, Healthy))[,"Genus"])

ONLY_IN_Healthy<- Healthy[! Healthy %in% Contaminant_trap]
healthy_list<-ONLY_IN_Healthy
### further exploration for the plot
Blood_Kit<-Healthy[Healthy %in% as.character(tax_table(subset_samples(data.genus, Kit=="Blood_Tissue_kit"))[,"Genus"]) ]
Microbiome_Kit<-Healthy[Healthy %in% as.character(tax_table(subset_samples(data.genus, Kit!="Blood_Tissue_kit"))[,"Genus"]) ]
### Now collapsing in one vector for the script  
ONLY_IN_Healthy<- paste(ONLY_IN_Healthy, collapse = ", ")
head(ONLY_IN_Healthy)

ONLY_IN_Contaminant_trap<- Contaminant_trap[! Contaminant_trap %in% Healthy]
ONLY_IN_Contaminant_trap<- paste(ONLY_IN_Contaminant_trap, collapse = ", ")
head(ONLY_IN_Contaminant_trap)

con<-file("Results/Beta_div/HEALTHY_vs_TRAPS_exclusive_bacteria.txt")
sink(con, append=TRUE)
cat("\n\nONLY IN Contaminant_trap", fill=TRUE)
cat(ONLY_IN_Contaminant_trap)
cat("\n\nONLY IN Healthy", fill=TRUE)
cat(ONLY_IN_Healthy)

sink()
close(con)

# version 1
x<-list('Controls                 '=Contaminant_trap,
        "          Healthy blood"=Healthy)
ggvenn(x, stroke_size = 0.6,
       set_name_size = 5.2,
       show_percentage = F,
       stroke_color = "darkgray",
       text_size = 6,
       text_color = "black",
       fill_color = c("lightgray","coral")) +
  theme(plot.margin = margin(-5,-15,0,-15)
  )
ggsave(filename = "Results/Beta_div/Venn_Diagramm_Healthy_vs_TRAPS.png", width = 4.5, height = 4.5, dpi=300, bg = "white")
dev.off()


# version 2
x<-list('Controls         '=Contaminant_trap,
        "      DNA Microbiome Kit"=Microbiome_Kit,
        " Blood & Tissue Kit"=Blood_Kit)
ggvenn(x, stroke_size = 0.8,
       set_name_size = 5,
       show_percentage = F,
       stroke_color = "darkgray",
       text_size = 6,
       fill_alpha = 0.6,
       text_color = "black",
       fill_color = c("lightgray","deepskyblue3","coral")) +
  theme(plot.margin = margin(0,-15,0,-15)
  )
ggsave(filename = "Results/Beta_div/Venn_Diagramm_Kits_vs_TRAPS.png", width = 4.5, height = 4.5, dpi=300, bg = "white")
ggsave(filename = "Results/Beta_div/Venn_Diagramm_Kits_vs_TRAPS_pdf.pdf", width = 4.5, height = 4.5, dpi=300, bg = "white")
dev.off()


####### ABUNDANCES BAR PLOT OF THOSE "UNIQUE" BACTERIA
data.selected.healthy<-subset_samples(data.genus.prop, Sample_Type=="Healthy")
data.selected.healthy<-subset_taxa(data.selected.healthy, Genus %in% healthy_list) # the list comes from above, they are the 'exclusive' bacteria!
tabella<-psmelt(data.selected.healthy)
# to count in how many sample...
tabella_count<-tabella[tabella$Abundance>0, ]
number_of_subjects<-list()
number_of_samples<-list()
for(i in unique(tabella_count$Genus)){
  count_subj<-unique(tabella_count[tabella$Genus==i, "Donor_Condition"]) 
  number_of_subjects[i]<-length( count_subj[! is.na(count_subj)] )
  count_sampl<-unique(tabella_count[tabella$Genus==i, "Sample_code"]) 
  number_of_samples[i]<-length( count_sampl[! is.na(count_sampl)] )
  rm(count_subj, count_sampl)
}
# adding those counts to the table
tabella$number_of_subjects<-rep("temp")
tabella$number_of_samples<-rep("temp")
for( i in unique(tabella$Genus)) {
tabella[tabella$Genus==i, "number_of_subjects"] <- paste("In",as.numeric(number_of_subjects[i]),"subjects")
tabella[tabella$Genus==i, "number_of_samples"] <- paste("In",as.numeric(number_of_samples[i]),"samples")
}

# graphical settings
tabella$Genus<-factor(tabella$Genus, levels = unique(tabella$Genus)[order(unique(tabella$Genus))] )
max_y_height_count<-tapply(tabella$Abundance, tabella$Genus, max)
tabella$max_y_height_count<-0 # adding those y axis heights to the table
for( i in unique(tabella$Genus)) {
  tabella[tabella$Genus==i, "max_y_height_count"] <- as.numeric(max_y_height_count[i])
}

# plotting
ggplot(data=tabella, aes(x=Genus, y=Abundance, fill=Donor_Condition)) +
  geom_bar(stat="identity", position="stack") + 
  theme_classic(base_size =9.5) +
  theme(axis.text.x=element_text(angle=10, 
                                 vjust=1, hjust=1,
                                 size=9.5), 
        axis.title.x=element_text(size=0.2), # (more space!)
        axis.title.y=element_text(size=8.6), 
        plot.caption = element_text(size=2),
        panel.grid.major.y = element_line(size=0.18, color="grey"),
        panel.grid.minor.y = element_line(size=0.05, color="grey"),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.text = element_text ( size = 10 ),
        legend.spacing.x = unit(0.5,"cm"),
        legend.position="bottom",
        plot.margin = margin(-5,3,-10,3),
        legend.margin = margin(0,8,-8,0)
  ) + 
  # scale_y_sqrt(limits = c(0,100*(length(unique(tabella$FASTQ_ID))+1)), # each sample has a cap of 100 (max percent abundance) but here are considered the abundances for each sample
  #              breaks = c(0, 5, 25, 100, 200, 350,
  #                         seq(500, 100*(length(unique(tabella$FASTQ_ID))+1), 100) )
  #                    ) +
  #scale_y_break(c(250, 1050), scales=0.15, space = 0.3) +
  scale_y_break(c(50, 1200), scales=0.05, space = 0.5) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    breaks = seq(0,1200,5) ,
    limits = c(0,1200)
    ) +
  theme(axis.line.y.right = element_blank(), # to hide the dual axis
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank()
        ) +
  guides(fill=guide_legend(nrow=2)) +
  #geom_text(aes(label=number_of_subjects, y= max_y_height_count), vjust=-2.8) +
  #geom_text(aes(label=number_of_samples, y= max_y_height_count), vjust=-1.5) +
  geom_text(aes(label=number_of_subjects, y= max_y_height_count+7.4), vjust=-1) +
  geom_text(aes(label=number_of_samples, y= max_y_height_count+5), vjust=-1) +
  labs(x="",
       y="             Sum of percent abundances among blood samples", 
       #caption = "The DNA has been extracted from the blood of each subject through two different kits\n--> there are 6 subjects and 11 samples in total (the sequencing of 1 'true blood' sample has been unsuccessfull)",
       fill="")
ggsave(filename = "Results/Beta_div/Healthy_vs_TRAP_exclusive_reads.png", width = 5.8, height = 4.8, dpi=300)
ggsave(filename = "Results/Beta_div/Healthy_vs_TRAP_exclusive_reads_PDF.pdf", width = 5.8, height = 4.8, dpi=300)

write.csv(tabella, "Results/Beta_div/TABELLA_Healthy_vs_TRAP_exclusive_reads.csv", row.names = F, quote = F)

suppressWarnings(rm(ONLY_IN_Healthy, ONLY_IN_Contaminant_trap, x, con, Healthy, Contaminant_trap, data.venn, who))


# Lactobacillus is a common contaminants of those kits (PMC4882852)
# Only Phascolarctobacterium is NOT a KNOWN contaminant (yet) and isolated from human gut  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5585883/
# anyways, according to PMC4882852, many contaminants detected in those blood sequencings have a relative abundances of <1 % and included many species associated the GUT microbiota



################ WHAT IF : E_coli is removed from ctrl samples? ###################
# E_coli compose most of the negative controls except one... what if it is removed from them???

data_new1<-subset_samples(data, ! Donor_Condition %in% c("Coli_colture"))
data_new2<-subset_samples(data, Donor_Condition %in% c("Coli_colture"))
data_new2<-subset_taxa(data_new2, Genus!="Escherichia-Shigella") # removing E.coli only from neg ctrls
data_new<-merge_phyloseq(data_new1, data_new2)
data_new.genus<-tax_glom(data_new, taxrank = "Genus")
data_new.genus.prop<-transform_sample_counts(data_new.genus, function(x) (x/sum(x))*100 )
data_new.phyla<-tax_glom(data_new, taxrank = "Phylum")
data_new.phyla.prop<-transform_sample_counts(data_new.phyla, function(x) (x/sum(x))*100 )


######### Rarefaction curve
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
  # legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

small_labels<-sample_data(data_new)$Sample_code
small_labels<-gsub("Coli_colture","ENC",small_labels)
small_labels<-gsub("Mice_cell","OTC",small_labels)
# to un-hide their names
{ small_labels<-gsub("OTC_M","OTC_M   ",small_labels)
  small_labels<-gsub("OTC_B","  OTC_B\n",small_labels)
  small_labels<-gsub("ENC_M","  ENC_M\n",small_labels)
}

png(file="Data_check/Rarefaction_curve_after_Ecoli_removal_from_EcoliCtrls.png",width=2000,height=1800, res=300)
r<-rarecurve(t(as(otu_table(data_new.genus),"matrix")),
             label=F, ylab = "Number of different genera after filters", xlab="Sample depth\n(after removing E.coli from E.coli ctrls)")
evalslopes(r,small_labels,lim=0.0001,cex=0.75)
dev.off()

pdf(file="Data_check/Rarefaction_curve_after_Ecoli_removal_from_EcoliCtrls_pdf.pdf",width=10,height=8)
r<-rarecurve(t(as(otu_table(data_new.genus),"matrix")),
             label=F, ylab = "Number of different genera after filters", xlab="Sample depth\n(after removing E.coli from E.coli ctrls)")
evalslopes(r,small_labels,lim=0.0001,cex=0.75)
dev.off()

rm(r,small_labels)



######### TOP 5 Genera with no E.coli in neg ctrls
{top <- names(sort(taxa_sums(data_new.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data_new.genus.prop)
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
others<-taxa_names(data_new.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data_new.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_top$Genus<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.-Paraburk.", tabella_top$Genus)
tabella_top$Genus<-gsub("Escherichia-Shigella","Escherichia", tabella_top$Genus)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_Type<-gsub("Contaminant_trap","Controls",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
# tabella$Donor_Condition<-gsub("Healthy","Subject_", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Coli_colture","ENC2", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Mice_cell","ENC1", tabella$Donor_Condition)
tabella$Sample_code<-gsub("Mice_cell","OTC",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Genus)) +
  facet_grid2( Kit ~ Sample_Type, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack", size=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,3,-9),
        legend.position="bottom",
        plot.margin = margin(3,2,3,-7)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant genera (removing E.coli)",
       caption = " 'Others' includes every genus below rank 5 ") #+
  #coord_flip() # usefull to flip everything with just a click!
#ggsave(file="Results/Abundances/TOP_5_genera_BUT_REMOVING_Ecoli_FROM_NEG_CTRLS.png",width=5.85,height=3.7, dpi=300)
ggsave(file="Results/Abundances/TOP_5_genera_BUT_REMOVING_Ecoli_FROM_NEG_CTRLS_v2.png",width=3.9,height=5.8, dpi=300)
#ggsave(file="Results/Abundances/TOP_5_genera_BUT_REMOVING_Ecoli_FROM_NEG_CTRLS_pdf.pdf",width=5.85,height=3.7, dpi=300)
dev.off()

write.xlsx(file = "Results/Abundances/TOP_5_Genera_Averages_WITHOUT_Ecoli_in_Ctrls.xlsx", row.names = F,
           cbind.data.frame("mean (no_Ecoli_ctrls)"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

rm(top, tabella, tabella_top)




######### TOP 5 Phyla with no E.coli in neg ctrls

{top <- names(sort(taxa_sums(data_new.phyla.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data_new.phyla.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_new.phyla.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_new.phyla.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_top$Phylum<-gsub("Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Caball.-Paraburk.", tabella_top$Phylum)
  tabella_top$Phylum<-gsub("Escherichia-Shigella","Escherichia", tabella_top$Phylum)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Sample_Type<-gsub("Contaminant_trap","Controls",tabella$Sample_Type)
tabella$Sample_Type<-gsub("Healthy","Healthy blood",tabella$Sample_Type) 
tabella$Kit<-gsub("Microbiome_kit","QIAamp® DNA Microbiome Kit",tabella$Kit) 
tabella$Kit<-gsub("Blood_Tissue_kit","DNeasy® Blood & Tissue Kit ",tabella$Kit)  
# tabella$Donor_Condition<-gsub("Healthy","Subject_", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Coli_colture","ENC2", tabella$Donor_Condition)
# tabella$Donor_Condition<-gsub("Mice_cell","ENC1", tabella$Donor_Condition)
tabella$Sample_code<-gsub("Mice_cell","OTC",tabella$Sample_code)
tabella$Sample_code<-gsub("Coli_colture","ENC",tabella$Sample_code)
ggplot(data=tabella, aes(x=Sample_code, y=Abundance, fill=Phylum)) +
  facet_grid2( Kit ~ Sample_Type, scales = "free", independent = T) + 
  geom_bar(stat="identity", position="stack", size=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=6.5),
        axis.text.y=element_text(angle=0, size=7), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 9.8 ),
        legend.margin = margin(0,0,10,-15),
        legend.position="bottom",
        plot.margin = margin(3,2,3,-7)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="", y="Abundance %", 
       fill="",
       title = "Five most abundant phyla (removing E.coli)",
       caption = " 'Others' includes every phyla below rank 5 ") # +
  #coord_flip() # usefull to flip everything with just a click!
#ggsave(file="Results/Abundances/TOP_5_phyla_BUT_REMOVING_Ecoli_FROM_NEG_CTRLS.png",width=5.85,height=3.7, dpi=300)
ggsave(file="Results/Abundances/TOP_5_phyla_BUT_REMOVING_Ecoli_FROM_NEG_CTRLS_v2.png",width=3.9,height=5.8, dpi=300)
#ggsave(file="Results/Abundances/TOP_5_phyla_BUT_REMOVING_Ecoli_FROM_NEG_CTRLS_pdf.pdf",width=5.85,height=3.7, dpi=300)
dev.off()

write.xlsx(file = "Results/Abundances/TOP_5_Phyla_Averages_WITHOUT_Ecoli_in_Ctrls.xlsx", row.names = F,
           cbind.data.frame("mean (no_Ecoli_ctrls)"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, tabella, tabella_top)


########## PCoA WITHOUT E.coli
data.prop.labels<-data_new.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       
{
  sample_data(data.prop.labels)$Donor_Condition2<-gsub("Healthy","",sample_data(data.prop.labels)$Donor_Condition)       
  sample_data(data.prop.labels)$Donor_Condition2<-gsub("Coli_colture","E.coli",sample_data(data.prop.labels)$Donor_Condition2)       
  sample_data(data.prop.labels)$Donor_Condition2<-gsub("_cell","",sample_data(data.prop.labels)$Donor_Condition2)       
  sample_data(data.prop.labels)$Donor_Condition3<-gsub("E.coli","ENC",sample_data(data.prop.labels)$Donor_Condition2)       
  sample_data(data.prop.labels)$Donor_Condition3<-gsub("Mice","OTC",sample_data(data.prop.labels)$Donor_Condition3)       
  sample_data(data.prop.labels)$Sample_Type<-gsub("Contaminant_trap","ENC",sample_data(data.prop.labels)$Sample_Type)
  sample_data(data.prop.labels)$Sample_Type<-gsub("Healthy","Healthy blood",sample_data(data.prop.labels)$Sample_Type)
  sample_data(data.prop.labels)$Kit<-gsub("Microbiome_kit","DNA Microbiome Kit",sample_data(data.prop.labels)$Kit) 
  sample_data(data.prop.labels)$Kit<-gsub("Blood_Tissue_kit","Blood & Tissue Kit ",sample_data(data.prop.labels)$Kit)  
}
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type",
                shape = "Kit") +
  # stat_ellipse(size=0.15) +
  scale_color_manual(values=c("Healthy blood"="coral","ENC"="gray")) +
  geom_point(size=5.4, alpha= 0.25) +
  geom_point(size=3.2, alpha= 0.88) +
  theme_classic(base_size = 10) +
  theme(legend.margin = margin(0,0,0,-2),
        legend.text = element_text(size=8.2) 
  ) +
  geom_text(aes(label=Donor_Condition3),
            color="black", size=2.3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
    color="Sample Type", x=paste("PC1: ",eigval[1],"% variation"),
    y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_AFTER_REMOVING_Ecoli_FROM_CTRLS.png", width = 5, height = 3.6, dpi=300)
ggsave(file="Results/Beta_div/PCoA_Beta_AFTER_REMOVING_Ecoli_FROM_CTRLS_pdf.pdf", width = 5, height = 3.6, dpi=300)




##### PCoA IF E.Coli is removed ALSO in that strange Mice sample ...
data_new1<-subset_samples(data, ! Sample_code %in% c("Coli_colture_B","Coli_colture_M","Mice_cell_M"))
data_new2<-subset_samples(data, Sample_code %in% c("Coli_colture_B","Coli_colture_M","Mice_cell_M") )
data_new2<-subset_taxa(data_new2, Genus!="Escherichia-Shigella") # removing E.coli only from neg ctrls
data_new_again<-merge_phyloseq(data_new1, data_new2)
data_new.genus_again<-tax_glom(data_new_again, taxrank = "Genus")
data_new.genus.prop_again<-transform_sample_counts(data_new.genus_again, function(x) (x/sum(x))*100 )

data.prop.labels<-data_new.genus.prop_again
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       
{
  sample_data(data.prop.labels)$Donor_Condition2<-gsub("Healthy","",sample_data(data.prop.labels)$Donor_Condition)       
  sample_data(data.prop.labels)$Donor_Condition2<-gsub("Coli_colture","E.coli",sample_data(data.prop.labels)$Donor_Condition2)       
  sample_data(data.prop.labels)$Donor_Condition2<-gsub("_cell","",sample_data(data.prop.labels)$Donor_Condition2)       
  sample_data(data.prop.labels)$Donor_Condition3<-gsub("E.coli","ENC",sample_data(data.prop.labels)$Donor_Condition2)       
  sample_data(data.prop.labels)$Donor_Condition3<-gsub("Mice","OTC",sample_data(data.prop.labels)$Donor_Condition3)       
  sample_data(data.prop.labels)$Sample_Type<-gsub("Contaminant_trap","ENC",sample_data(data.prop.labels)$Sample_Type)
  sample_data(data.prop.labels)$Sample_Type<-gsub("Healthy","Healthy blood",sample_data(data.prop.labels)$Sample_Type)
  sample_data(data.prop.labels)$Kit<-gsub("Microbiome_kit","DNA Microbiome Kit",sample_data(data.prop.labels)$Kit) 
  sample_data(data.prop.labels)$Kit<-gsub("Blood_Tissue_kit","Blood & Tissue Kit ",sample_data(data.prop.labels)$Kit)  
}
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type",
                shape = "Kit") +
  # stat_ellipse(size=0.15) +
  scale_color_manual(values=c("Healthy blood"="coral","ENC"="gray")) +
  geom_point(size=5.4, alpha= 0.25) +
  geom_point(size=3.2, alpha= 0.88) +
  theme_classic(base_size = 10) +
  theme(legend.margin = margin(0,0,0,-2),
        legend.text = element_text(size=8.2) 
  ) +
  geom_text(aes(label=Donor_Condition3),
            color="black", size=2.3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
    color="Sample Type", x=paste("PC1: ",eigval[1],"% variation"),
    y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_AFTER_REMOVING_Ecoli_ALSO_FROM_MICE_B.png", width = 5, height = 3.6, dpi=300)
ggsave(file="Results/Beta_div/PCoA_Beta_AFTER_REMOVING_Ecoli_ALSO_FROM_MICE_B_pdf.pdf", width = 5, height = 3.6, dpi=300)



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
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
