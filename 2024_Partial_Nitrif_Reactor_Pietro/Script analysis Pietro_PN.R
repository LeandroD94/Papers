##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  # analysis packages
  library("vegan")
  # utilities
  library("reshape")
  library("xlsx")  
  library("qiime2R")
}

{dir.create("Data_check_PartialNitrif")
dir.create("Results_PartialNitrif")
dir.create("Results_PartialNitrif/Abundances")
dir.create("Results_PartialNitrif/Correlations_along_time/")
}
options(scipen = 100) # disable scientific annotation



### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue2","firebrick3","gold","gray","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group




################ DEFINING THE GROUPS OF ORGANISMS ################

### common PAO GAO list (PMID:38524765 ,  PMID: 22827168,  PMID: 15774634,  PMID: 33404187, MIDAS website)
PAO_list<-c("Candidatus_Phosphoribacter","Ca_Phosphoribacter","Phosphoribacter",
            "Tetrasphaera",
            "Candidatus_Accumulimonas","Accumulimonas","Ca_Accumulimonas",
            "Candidatus_Accumulibacter","Accumulibacter","Ca_Accumulibacter",
            "Candidatus_Microthrix","Microthrix","Ca_Microthrix",
            "Candidatus_Dechloromonas","Dechloromonas",
            "Malikia",
            "Quatrionicoccus",
            "Beggiatoa",
            "Gemmatimonas",
            "Friedmaniella",
            "Tessaracoccus",
            "Azonexus",
            # "Pseudomonas",   # too aspecific!
            "Candidatus_Accumulimonas","Accumulimonas","Ca_Accumulimonas",
            "Microlunatis","Microlunatus",
            "Candidatus_Lutibacillus","Lutibacillus","Ca_Lutibacillus")
GAO_list<-c("Candidatus_Competibacter","Candidatus_Contendobacter","Candidatus_Proximibacter",
            "Ca_Competibacter","Ca_Contendobacter","Ca_Proximibacter",
            "Competibacter","Contendobacter",
            # "Proximibacter", # may also be a PAO!
            "Defluviicoccus","Micropruina","Propionivibrio")

# AOBs derive from from MIDAS, also in other articles I can't find other names...
AOB_list <- c("Nitrosomonas", "Nitrosococcus",   # synonymous, see MIDAS
              "Nitrosospira", "Nitrosolobus", "Nitrosovibrio")   # synonymous, see MIDAS
AOB_list <- c(AOB_list, "Nitrososphaera", "Nitrosocaldus", "Nitrosotalea")  # adding also AOA, from PMID: 24559743
NOB_list <- c("Nitrospira",  # alcuni da MIDAS, ma tutti elencati in DOI: 10.1016/S0076-6879(11)86005-2
              "Nitrospina",
              "Nitrobacter",
              "Nitrococcus",
              "Nitrotoga", "Candidatus Nitrotoga")

# NO2_r_list<-c("Thauera","Rhodoferax","Dokdonella",  # these are from MIDAS
#               "Ca Competibacter","Candidatus Competibacter", "Competibacter",
#               # "Nitrosomonas",   # PMID: 11948173 and PMID: 17298375 ... it looks like it can also reduce NO2- ... but including it in this list would cause problems in the plot with AOB
#               # "Nitrosospira",   # PMID: 17298375 , same as above
#               # "Nitrospira", https://doi.org/10.1016/j.tim.2016.05.004 ... this is a NOB, yet same as above
#               "Ca Accumulibacter","Candidatus Accumulibacter","Accumulibacter",
#               "Haliangium",
#               "Ca Promineofilum","Candidatus Promineofilum","Promineofilum",
#               "Rhodoplanes", "Bradyrhizobium", "Iamia", "Thiothrix", "Sulfuritalea",
#               "Zoogloea", "Thermogutta", "Diaphorobacter" , "Simplicispira",
#               "Ca Phosphitivorax", "Candidatus Phosphitivorax", "Phosphitivorax",
#               "Steroidobacter" ,"Sterolibacterium", "Azospira", "Acidovorax",
#               "Pseudomonas",
#               "Uruburuella", "Anaerolinea", "Corynebacterium",
#               "Luteimonas", "Hyphomonas", "Halioglobus", "Azoarcus",
#               "Ca Solibacter", "Candidatus Solibacter", "Solibacter",
#               "Ca Brocadia", "Candidatus Brocadia", "Brocadia",
#               "Skermania" , "Pannonibacter", "Sulfurimonas" , "Desulfitibacter" ,
#               "Thiopseudomonas" , "Azonexus" ,
#               "Ca Proximibacter", "Candidatus Proximibacter", "Proximibacter",
#               #the following denitrifiers are from  PMID: 24559743
#               "Enterobacter", "Micrococcus", "Spirillus","Proteus", "Aerobacter","Flavobacterium",
#               # the following are from the FISH hand book  (Vedi anche tabella 3.1 del relativo capitolo)
#               "Curvibacter", "Paracoccus", "Rhizobium", "Delftia", "Simplicispira","Dechloromonas", "Halomonas", "Thermomonas", "Aminomonas", "Methylophilaceae",
#               "Azovibrio","Ralstonia", "Shinella","Ochrobactrum","Hyphomicrobium","Blastobacter", "Pseudoxanthomonas", "Stenotrophomonas", "Castellaniella",
#               "Herbaspirillum","Citrobacter","Methylophilaceae","Alicycliphilus","Ottowia","Diaphorobacter", "Rhodobacter", "Microvirgula","Aquaspirillum", "Vogesella"
# )




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza", tree = "QIIME/rooted-tree.qza")

# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

metadata <- as.data.frame(read.table(file="Metadata_Pietro.tsv", sep="\t", header = T))
row.names(metadata)<-metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(metadata)
original_length<-length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])
original_order<-metadata$Sample_name # to maintein this sample order in plots
metadata<-metadata[sample, ]
identical(as.numeric(length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)

sample_names(data)<-sample_data(data)$Sample_name
sample_data(data)$Sample_name<-factor(sample_data(data)$Sample_name,
                                      levels = original_order)
sample_data(data)$Exp_Day<-factor(sample_data(data)$Exp_Day,
                                         levels = unique(metadata$Exp_Day))





#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check_PartialNitrif/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.01% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
target_filter<- filter_taxa(data.genus.temp, function(x) max(x) <= 0.01, TRUE)
filtered<-taxa_names( target_filter )
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check_PartialNitrif/Filtered_genus_under_001_max_cutoff.csv")

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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check_PartialNitrif/Bacteria_Archaea_proportion_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)



# removing every Unassigned Phylum (probable contaminant)
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


con<-file("Data_check_PartialNitrif/DETAILS_ABOUT_ORIGINAL_AND_PROCESSED_READS_NUMBER.txt")
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
table <- cbind.data.frame(Original=Original_read_number$forward.sequence.count, # merged --> just one column, not both of them
                        After_quality_filter=DADA2_read_number$non.chimeric,
                        After_contam_filter=after_filter_number[,1])
table$FASTQ<-as.character(sample_data(data)$FASTQ_ID)
# re-naming and adding groups
if(identical( table$FASTQ , as.character(sample_data(data)$FASTQ_ID) ) ){ # same order
  table$Samples<-sample_data(data)$Sample_name
  #table$Experiment_day<-sample_data(data)$Experiment_day
  # table$factor <- sample_data(data)$Reactor_Type
  # table$factor <- gsub("N_limitation","N limiting",table$factor)
  # table$factor <- gsub("P_limitation","P limiting",table$factor)
  # table$factor <- as.factor(table$factor)
  cat("\n\nOK!\n\n")
}
# plotting
table$Samples<- factor( paste0(sample_data(data)$Exp_Day, "th day") , levels= paste0(sample_data(data)$Exp_Day, "th day") )
ggplot(aes(x=Samples, y=Original), data=table) +
  #facet_grid( ~ factor, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Read quality filters"), alpha= 0.8,
            width = 0.65, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
            width= 0.35, stat = "identity") +
  theme_classic( base_size = 10.2) +
  theme(axis.text.x = element_text(size=6.5, angle = 35, vjust=1, hjust=1),
        axis.text.y = element_text(size=5),
        panel.grid.major.y = element_line(linewidth =0.1, color="grey"),
        legend.position = "bottom",
        legend.margin = margin(0, 30, 0, 0),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9.5),
        legend.key.height = unit(0.4,"cm")
  ) +
  scale_fill_manual(name='Number of reads:  ',
                    breaks=c('Original number', 'Read quality filters', 'Relative abundance filters'),
                    values=c('Original number'='green3', 'Read quality filters'='coral', 'Relative abundance filters'='red3')) +
  scale_y_continuous(breaks = c(10000, seq(0, max(table$Original+10000), 25000))) +
  labs(y="Reads abundance", x="FASTQ")
ggsave(file="Data_check_PartialNitrif/Number_of_reads_pre_and_post_filters.png", width = 6, height = 3.8, dpi=300)
dev.off()
# saving also the table itself
write.csv2(table, file="Data_check_PartialNitrif/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)


suppressWarnings( rm(table, table2, Original_read_number, DADA2_read_number, after_filter_number) )



 
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
write.csv2(assigned,file="Data_check_PartialNitrif/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)




############################ RAREFACTION ANALYSIS ################################


png(file="Data_check_PartialNitrif/Rarefaction_curve.png",width=2500,height=1500, res=300)
data_temp<-data.genus
sample_names(data_temp)<-paste0( sample_data(data.genus)$Exp_Day, "th day")
rarecurve(t(as(otu_table(data_temp),"matrix")),
          step=100,label=T,
          cex=0.55,
          ylab = "Genera", xlab= "Reads amount")
dev.off()
rm(data_temp)


proof2<-"No unsaturated sample to remove"




#################### \\\\\\ STARTING THE ANALYSIS \\\\\ #######################


### everything ready, let's start with the analysis
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




########################### COUNTS EXPORT ##########################################

dir.create("Results_PartialNitrif/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_PartialNitrif/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results_PartialNitrif/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_order.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_PartialNitrif/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}




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
  tabella$Phylum<-gsub ("Candidatus","Ca.", tabella$Phylum)
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
tabella$Exp_Day<-factor(tabella$Exp_Day, levels = unique(metadata$Exp_Day) )
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  #facet_grid(~ Project_description, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=32,
                                 vjust=1,
                                 hjust=1,
                                 size= 9 
                                 ),
  axis.title =element_text(size=10),
  legend.key.width = unit(0.25, "cm"),
  legend.key.height = unit(0.38, "cm"),
  legend.text = element_text ( size = 10.5 ),
  legend.position="bottom",
  legend.margin = margin(-3,5,9,5),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percentual abundance of clades",
       title = "",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results_PartialNitrif/Abundances/TOP_microbial_phyla.png",width=6,height=5, dpi=300)
dev.off()

# means of TOP phyla
write.xlsx(file = "Results_PartialNitrif/Abundances/TOP_microb_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)




### TOP Genera (19)
suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:19]
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
  tabella$Genus<-gsub ("Candidatus","Ca.", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Exp_Day<-factor(tabella$Exp_Day, levels = unique(metadata$Exp_Day) )
#fill_color_modified<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","red","pink3", "blue","gray","gold","darkgreen","violet", "chartreuse3","deepskyblue2","wheat3","firebrick3","darkslategray3")
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_19) +
  theme(axis.text.x=element_text(angle=32,
                                 vjust=1,
                                 hjust=1,
                                 size= 9
  ),
  axis.title =element_text(size=10),
  #strip.text = element_text(size=6.7),
  legend.key.width = unit(0.22, "cm"),
  legend.key.height = unit(0.38, "cm"),
  legend.text = element_text ( size = 9 ),
  legend.position="bottom",
  legend.margin = margin(-10,42,5,1),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="", y="Percentual abundance of clades",
       title = "",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results_PartialNitrif/Abundances/TOP_microbial_genera_19.png",width=6,height=5, dpi=300)
dev.off()


# means of TOP genera
write.xlsx(file = "Results_PartialNitrif/Abundances/TOP_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



### TOP Genera (10)
suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:10]
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
  tabella$Genus<-gsub ("Candidatus","Ca.", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Exp_Day<-factor(tabella$Exp_Day, levels = unique(metadata$Exp_Day) )
fill_color_modified<-c("darkblue","brown4","springgreen2","wheat","darkgreen","coral","yellow3","darkmagenta","red","chartreuse2", "darkslategray3")
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_modified) +
  theme(axis.text.x=element_text(angle=32,
                                 vjust=1,
                                 hjust=1,
                                 size= 9
  ),
  axis.title =element_text(size=10),
  #strip.text = element_text(size=6.7),
  legend.key.width = unit(0.22, "cm"),
  legend.key.height = unit(0.38, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin(-3,42,9,1),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percentual abundance of clades",
       title = "",
       fill="",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results_PartialNitrif/Abundances/TOP_microbial_genera_10.png",width=6,height=5, dpi=300)




############################ ALPHA DIVERSITY ##############################

data.genus.temp<-data.genus
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Exp_Day
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
mix_alpha$Sample_name<-gsub(".","/",mix_alpha$Sample_name, fixed = T)
mix_alpha$Sample_name<-gsub("X","",mix_alpha$Sample_name, fixed = T)
mix_alpha$Sample_name<-factor(mix_alpha$Sample_name, levels=unique(metadata$Exp_Day))
mix_alpha$ID<-rep("unique_ID_for_a_line")
ggplot(data=mix_alpha, aes(y=value, x=Sample_name)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  geom_point(size =1.8, color="green") +
  geom_point(size =3.25, alpha=0.4, color="green") +
  #geom_line(aes(group=ID), size= 0.4, alpha= 0.4) + 
  geom_path(aes(group=as.character(ID)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic2(base_size = 9) +
  labs(title="Alpha diversity along the time", y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.035,0)) + # to have more space on the borders
  # geom_text(aes(label=Sample_name), color="black", size=1.8, vjust=-0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_text(angle=32,
                                   vjust=1,
                                   hjust=1,
                                   size= 7.35
  ),
  axis.text.y= element_text(angle=0, size=5.5),
  axis.title.x = element_blank(),
  panel.grid.major.y = element_line(linewidth=0.4),
  panel.grid.minor.y = element_line(linewidth=0.25),
  plot.margin = margin(5,1,5,1)
  )
ggsave(file="Results_PartialNitrif/Alfa_diversity_along_the_time_computed_on_genera.png", width = 5.5,height =4.5, dpi=300)


suppressWarnings(rm(mix_alpha) )




########################### PCoA BETA DIV ##########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Exp_Day<-as.numeric(sample_data(data.prop.labels)$Exp_Day)
sample_data(data.prop.labels)$group<-rep("Same_group_same_line")
sample_names(data.prop.labels)<- paste0(sample_data(data.prop.labels)$Exp_Day , "th day")
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color="Exp_Day") + # , color = "Sample_Type") +
  geom_path(aes(group=group), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  ) +
  scale_color_gradient2(low="darkgreen", mid = "chartreuse3", high="green", lim=c(NA,NA),
                        breaks=c(1,150,188),
                        midpoint = 5) +
  geom_point(size=3.8, alpha=0.5,  show.legend = F) +
  guides(color="none") +
  theme_classic(base_size = 12.5) + 
  theme(plot.title = element_blank()) + 
  #geom_text(aes(label=sample_data(data.sqrt_prop)$Exp_Day), color="gray20", size=2.5, vjust= -0.2, show.legend = FALSE) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="gray20", size=2.5, vjust= -0.2, show.legend = FALSE) +
  labs(
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results_PartialNitrif/Beta_divers_Hellinger_on_genera_WITH_DATES.png", width = 6.5, height = 5, dpi=300)




########## CORRELATIONS ABUNDANCES vs TIME _ N limiting (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
# data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
#data.genus.temp <- subset_samples(data.genus, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
data.genus.temp<- data.genus
data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sample_name
order_of_samples<-metadata[levels(sample_data(data)$Sample_name), "Sample_name"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]


Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "kendall") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")

write.csv2(Corr_results, "Results_PartialNitrif/Correlations_along_time/Correlations_CLR_bacteria_abundances_with_time_Kendall_algorithm.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )

con <- file("Results_PartialNitrif/Correlations_along_time/NNNNNOTHING_THERE.txt")
sink(con, append = T)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the samples (-->",min_prevalence,"samples) have been correlated with the passing of time.", fill=T)
cat("Among these, none was significantly increased or decreased along the passing of time after the BH adjustment of the p-value... maybe due to the low number of samples.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))




############################ FOCUS ON PAO AND GAO ########################

data.genus.temp<-data.genus.prop
PAO_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% PAO_list )
# with PAO65 off-targets (excluding Ca Competibacter)
# PAO_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% c(PAO_list,"Dechlorosoma", "Methyloglobulus", "Probionivibrio","Thauera") )
GAO_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% GAO_list)
tabella_PAO<-psmelt(PAO_identified_phylo)
tabella_GAO<-psmelt(GAO_identified_phylo)
tabella_PAO$GroupPG <- "PAO"
tabella_GAO$GroupPG <- "GAO"
tabella<- rbind.data.frame(tabella_PAO, tabella_GAO)
tabella$GroupPG<-factor(tabella$GroupPG, levels=c("PAO","GAO"))
tabella$Exp_Day<-factor(tabella$Exp_Day, levels=unique(metadata$Exp_Day))

ggplot(data=tabella, aes( x=Exp_Day, y=Abundance, fill=GroupPG)) + 
  #facet_grid2( ~Type, scales = "free", space="free",
  #             strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("PAO"="red2","GAO"="chartreuse2")) +
  theme(axis.text.x=element_text(angle=38, hjust=1,vjust=1, size=6.5),
        axis.title.y = element_text(size=9), 
        axis.text.y = element_text(size=7), 
        strip.text.x = element_text(size=8.5), 
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.8 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,1,2,1),
        legend.margin =  margin(-10,2,0,0)
  ) +
  scale_y_continuous(lim=c(0,50), breaks=seq(0,50,5)) +
  labs(x="", y="Percent abundance in the sample", 
       fill="",
       title = paste("PAO / GAO abundances in the samples"))
ggsave(file="Results_PartialNitrif/Abundances/PAO_GAO_comparison.png",width=6,height=5,dpi=300)
dev.off()


# Focus on Accumulibacter and Competibacter
tabella$GroupPG2<- "temp" # extra line to allow a quick reset of this part (the factors hate the g_subs )
tabella$GroupPG2[ tabella$Genus %in% c("Candidatus_Accumulibacter","Candidatus_Competibacter") ] <- tabella$Genus [tabella$Genus %in% c("Candidatus_Accumulibacter","Candidatus_Competibacter") ]
tabella$GroupPG2[ ! tabella$Genus %in% c("Candidatus_Accumulibacter") &  tabella$Genus %in% PAO_list ] <- "Other PAOs"
tabella$GroupPG2[ ! tabella$Genus %in% c("Candidatus_Competibacter") &  tabella$Genus %in% GAO_list ] <- "Other GAOs"
tabella$GroupPG2<- gsub("Candidatus_","Ca. ", tabella$GroupPG2)
tabella$GroupPG2<- factor( tabella$GroupPG2 , levels= c("Ca. Accumulibacter", "Other PAOs", "Ca. Competibacter", "Other GAOs") )
ggplot(data=tabella, aes( x=Exp_Day, y=Abundance, fill=GroupPG2)) + 
  # facet_grid2( ~Type, scales = "free", space="free",
  #              strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("Ca. Accumulibacter"="coral", "Ca. Competibacter"="chartreuse2",
                             "Other PAOs"="red3","Other GAOs"="chartreuse4")) +
  theme(axis.text.x=element_text(angle=38, hjust=1,vjust=1, size=6.5),
        axis.title.y = element_text(size=9), 
        axis.text.y = element_text(size=7), 
        strip.text.x = element_text(size=8.5), 
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,1,2,1),
        legend.margin =  margin(-10,2,0,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  scale_y_continuous(lim=c(0,50), breaks=seq(0,50,5)) +
  labs(x="", y="Percent abundance in the sample", 
       fill="",
       title = paste("PAO / GAO abundances in the samples"))
ggsave(file="Results_PartialNitrif/Abundances/PAO_GAO_comparison_version2.png",width=6,height=5,dpi=300)
dev.off()  



# exporting the abundances
PAO_identified_phylo<- phyloseq(otu_table(PAO_identified_phylo),tax_table(PAO_identified_phylo)) # removing the two trees
GAO_identified_phylo<- phyloseq(otu_table(GAO_identified_phylo),tax_table(GAO_identified_phylo))
tot <- merge_phyloseq(PAO_identified_phylo, GAO_identified_phylo)
# means of every PAO and GAO
tot<- cbind.data.frame("average_among_samples"=as.numeric(apply(otu_table(tot),1,mean)), "genus"= as.data.frame(tax_table(tot))[["Genus"]])
tot$Type <- ifelse( tot$genus %in% PAO_list, "PAO", "GAO")
write.xlsx(file = "Results_PartialNitrif/Abundances/PAO_GAO_names_and_averages_MIXEDLIQUORS.xlsx", row.names = F,
           tot
)


suppressWarnings(rm(tot, tabella, tabella_GAO, tabella_PAO) )



####################### FOCUS ON AOB, NOB and N reducer ########################

data.genus.temp<-data.genus.prop
AOB_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% AOB_list )
NOB_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% NOB_list)
# Nred_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% NO2_r_list)
tabella_AOB<-psmelt(AOB_identified_phylo)
#tabella_AOB<-aggregate( Abundance~Exp_Day, data=tabella_AOB, FUN=sum)
tabella_NOB<-psmelt(NOB_identified_phylo)
#tabella_NOB<-aggregate( Abundance~Exp_Day, data=tabella_NOB, FUN=sum)
# tabella_Nred<-psmelt(Nred_identified_phylo)
tabella_AOB$GroupPG <- "AOB (Nitrosomonas)  "
tabella_NOB$GroupPG <- "NOB (Nitrospira)"
# tabella_Nred$GroupPG <- "N_reducer"  # from MIDAS only Nitrite red, but added also other denitrifiers from articles
# tabella<- rbind.data.frame(tabella_AOB, tabella_NOB, tabella_Nred)
# tabella$GroupPG<-factor(tabella$GroupPG, levels=c("AOB","NOB", "N_reducer"))
tabella<- rbind.data.frame(tabella_AOB, tabella_NOB)
tabella$GroupPG<-factor(tabella$GroupPG, levels=c("AOB (Nitrosomonas)  ","NOB (Nitrospira)"))
tabella$Exp_Day<-factor(tabella$Exp_Day, levels=unique(metadata$Exp_Day))

ggplot(data=tabella, aes( x=Exp_Day, y=Abundance, fill=GroupPG)) + 
  # facet_grid2( ~Type, scales = "free", space="free",
  #              strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 1) +
  # geom_bar(stat="identity", position="stack", width = 0.9, alpha= 1) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("AOB (Nitrosomonas)  "="lightblue",
                             "NOB (Nitrospira)"="deepskyblue4"
                             #"N_reducer"="darkblue" )
                    ) ) +
  theme(axis.text.x=element_text(angle=38, hjust=1,vjust=1, size=9),
        axis.title.y = element_text(size=9), 
        axis.text.y = element_text(size=7), 
        strip.text.x = element_text(size=8.5), 
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.8 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,1,2,1),
        legend.margin =  margin(-5,2,0,0)
  ) +
  scale_y_continuous(lim=c(0,40), breaks=seq(0,50,5)) +
  labs(x="", y="Percent abundance in the sample", 
       fill="",
       title = paste("AOB / NOB relative abundance")
       )
ggsave(file="Results_PartialNitrif/Abundances/AOB_NOB_comparison.png",width=5,height=4.2,dpi=300)
dev.off()


# exporting the abundances
AOB_identified_phylo<- phyloseq(otu_table(AOB_identified_phylo),tax_table(AOB_identified_phylo)) # removing the two trees
NOB_identified_phylo<- phyloseq(otu_table(NOB_identified_phylo),tax_table(NOB_identified_phylo))
#Nred_identified_phylo<- phyloseq(otu_table(Nred_identified_phylo),tax_table(Nred_identified_phylo))
#tot <- merge_phyloseq(AOB_identified_phylo, NOB_identified_phylo, Nred_identified_phylo)
tot <- merge_phyloseq(AOB_identified_phylo, NOB_identified_phylo)
# means of every AOB, NOB and denitrif
tot<- cbind.data.frame("average among samples"=as.numeric(apply(otu_table(tot),1,mean)),
                       otu_table(tot), 
                       "genus"= as.data.frame(tax_table(tot))[["Genus"]]
)
tot$Type <- tot$genus
tot$Type[tot$Type %in% AOB_list] <-"AOB"
tot$Type[tot$Type %in% NOB_list] <-"NOB"
# tot$Type[tot$Type %in% NO2_r_list] <-"N_reducer"
write.xlsx(file = "Results_PartialNitrif/Abundances/AOB_NOB_names_and_averages_MIXEDLIQUORS.xlsx", row.names = F,
           tot 
)


suppressWarnings(rm(tot, tabella, tabella_NOB, tabella_AOB) )




##################### \\\\\\ R AND PACKAGES VERSION \\\\\\ #########################


### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results_PartialNitrif/R_version_and_packages.txt")
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

