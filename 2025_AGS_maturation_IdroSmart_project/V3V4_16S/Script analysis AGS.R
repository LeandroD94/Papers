# SRT e Daily temperature sono abbastanza inutili nelle correlazioni... sostiturire con altri due utili (o anche di più... poco male, meno risultati!)
# chiedi round migliore per meas


##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  library("ggvenn")
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
  dir.create("Results/Correlations/")
  dir.create("Results/Checking_Measurements")
  dir.create("Results/Extra_edited_figures_for_presentations")
}
options(scipen = 100) # disable scientific annotation



### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
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

NO2_r_list<-c("Thauera","Rhodoferax","Dokdonella",  # these are from MIDAS
             "Ca Competibacter","Candidatus Competibacter", "Candidatus_Competibacter", "Competibacter",
             # "Nitrosomonas",   # PMID: 11948173 and PMID: 17298375 ... it looks like it can also reduce NO2- ... but including it in this list would cause problems in the plot with AOB
             # "Nitrosospira",   # PMID: 17298375 , same as above
             # "Nitrospira", https://doi.org/10.1016/j.tim.2016.05.004 ... this is a NOB, yet same as above
             "Ca Accumulibacter","Candidatus Accumulibacter","Candidatus_Accumulibacter", "Accumulibacter",
             "Haliangium",
             "Ca Promineofilum","Candidatus Promineofilum","Candidatus_Promineofilum","Promineofilum",
             "Rhodoplanes", "Bradyrhizobium", "Iamia", "Thiothrix", "Sulfuritalea",
             "Zoogloea", "Thermogutta", "Diaphorobacter" , "Simplicispira",
             "Ca Phosphitivorax", "Candidatus Phosphitivorax",  "Candidatus_Phosphitivorax", "Phosphitivorax",
             "Steroidobacter" ,"Sterolibacterium", "Azospira", "Acidovorax",
             "Pseudomonas",
             "Uruburuella", "Anaerolinea", "Corynebacterium",
             "Luteimonas", "Hyphomonas", "Halioglobus", "Azoarcus",
             "Ca Solibacter", "Candidatus Solibacter", "Candidatus_Solibacter", "Solibacter",
             "Ca Brocadia", "Candidatus Brocadia","Candidatus_Brocadia", "Brocadia",
             "Skermania" , "Pannonibacter", "Sulfurimonas" , "Desulfitibacter" ,
             "Thiopseudomonas" , "Azonexus" ,
             "Ca Proximibacter", "Candidatus Proximibacter","Candidatus_Proximibacter", "Proximibacter",
             #the following denitrifiers are from  PMID: 24559743
             "Enterobacter", "Micrococcus", "Spirillus","Proteus", "Aerobacter","Flavobacterium",
             # the following are from the FISH hand book  (Vedi anche tabella 3.1 del relativo capitolo)
             "Curvibacter", "Paracoccus", "Rhizobium", "Delftia", "Simplicispira","Dechloromonas", "Halomonas", "Thermomonas", "Aminomonas", "Methylophilaceae",
             "Azovibrio","Ralstonia", "Shinella","Ochrobactrum","Hyphomicrobium","Blastobacter", "Pseudoxanthomonas", "Stenotrophomonas", "Castellaniella",
             "Herbaspirillum","Citrobacter","Methylophilaceae","Alicycliphilus","Ottowia","Diaphorobacter", "Rhodobacter", "Microvirgula","Aquaspirillum", "Vogesella"
             )




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza", tree = "QIIME/rooted-tree.qza")

# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

# NB: the order of the row in the Metadata is employed also to set the order of the samples in the plots
metadata <- as.data.frame(read.table("AGS_metadata.tsv", header = T, sep="\t"))
row.names(metadata)<-metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])
original_order<-metadata$Sample_name # to maintein this sample order in plots
metadata<-metadata[sample, ]
identical(as.numeric(length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])),as.numeric(original_length))

# metadata$Sample_name<- factor(metadata$Sample_name, levels = unique(metadata$Sample_name) )
metadata$Type<- factor(metadata$Type, levels = c("Mixed liquor","Direct extract","Influent","Out") )

sample_data(data)<-metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

# restoring the original row order of the metadata (it sets also the order in the plots!)
row.names(metadata)<-metadata$Sample_name
metadata<-metadata[original_order, ]

sample_names(data)<-sample_data(data)$Sample_name
sample_data(data)$Sample_name<-factor(sample_data(data)$Sample_name,
                                      levels = original_order)
sample_data(data)$Time_point<-factor(sample_data(data)$Time_point,
                                     levels = unique(metadata$Time_point))
sample_data(data)$Experiment_day<-factor(sample_data(data)$Experiment_day,
                                     levels = unique(metadata$Experiment_day))
rm(original_length)




##### Importing the measurements metadata
meas <- read.delim("Reactor_measurements.tsv")
row.names(meas)<- meas$Sample_name
meas$Sample_name<-NULL
meas<-round(meas, digits = 1)
# adjusting few names for better aesthetics...
colnames(meas)<-gsub("VSS.TSS", "VSS/TSS", colnames(meas))
colnames(meas)<-gsub("SVI30.5", "SVI30/SVI5", colnames(meas))
colnames(meas)<-gsub(".1", "", colnames(meas), fixed = T) # this is due a little importing issue...
colnames(meas)<-gsub(".", " ", colnames(meas), fixed = T)




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
target_filter<- filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE)
filtered<-taxa_names( target_filter )
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_max_cutoff.csv")


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


# save.image("AGS_16S_after_filters.RData")




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
  #table$Experiment_day<-sample_data(data)$Experiment_day
  table$factor <- sample_data(data)$Type
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
  theme(axis.text.x = element_text(size=6, angle = 40, vjust=1, hjust=1),
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
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 6.5, height = 4.5, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)


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
r<-rarecurve(t(as(otu_table(data.genus),"matrix")), step=100,label=F, ylab = "Genera", xlab= "Reads amount")
evalslopes(r,sample_names(data.genus),lim=0.001,cex=1)
dev.off()
rm(r)



proof2<-"No unsaturated sample to remove"




#################### \\\\\\ STARTING THE ANALYSIS \\\\\ #######################


if(! "proof1" %in% ls() | ! "proof2" %in% ls() ){
  stop("\nDid you perform the filtering steps yet?\n")
}


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

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
    #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
}




###################### ABUNDANCES BAR PLOT (EVERY SAMPLE) ##########################


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
{
  levels(tabella$Type) <- gsub ("Influent","Influ", levels(tabella$Type))
  levels(tabella$Type) <- gsub ("Direct extract","Direct", levels(tabella$Type))
}

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
                                ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=4),
  axis.title =element_text(size=10),
  strip.text = element_text(size=9.5),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.44, "cm"),
  legend.text = element_text ( size = 9.8 ),
  legend.position="bottom",
  legend.margin = margin(0 ,5, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results/Abundances/TOP_phyla_EverySample.png",width=7.1,height=4.6, dpi=300)
ggsave(file="Results/Extra_edited_figures_for_presentations/TOP_phyla_EverySample_v2.png",width=7.1,height=4.5, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_phyla_averages_EverySample.xlsx", row.names = F, to_save)

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
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
{
  levels(tabella$Type) <- gsub ("Influent","Influ", levels(tabella$Type))
  levels(tabella$Type) <- gsub ("Direct extract","Direct", levels(tabella$Type))
}
fill_color_modified<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","red","pink3", "blue","gray","gold","darkgreen","violet", "chartreuse3","deepskyblue2","wheat3","firebrick3","darkslategray3")
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_modified) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
                                ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=4),
  axis.title =element_text(size=10),
  strip.text = element_text(size=9.5),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.385, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(-5,28,0,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       #title = "Most abundant identified genera (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_EverySample.png",width=7.1,height=4.6, dpi=300)
ggsave(file="Results/Extra_edited_figures_for_presentations/TOP_genera_EverySample_v2.png",width=7.1,height=4.5, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame( "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2),
                            "Average in mature pilot"= round( as.numeric(apply( otu_table(subset_samples(prune.dat_top,Experiment_day%in%c("151","154","157","160",
                                                                                                                                           "164","169","176","184"))) ,1,mean)) ,2)
)
                           
                           
to_save<-to_save[order(to_save$Average, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_genera_Averages_EverySample.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




###################### ABUNDANCES BAR PLOT (Only Mixed Liquor) ##########################

data_mix<-subset_samples(data.phy.prop, Type=="Mixed liquor")

### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data_mix)
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

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9.9 ),
  legend.position="bottom",
  legend.margin = margin(0 ,5, 0 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified phyla",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results/Abundances/TOP_phyla_onlyMixedLiquor.png",width=7,height=5.4, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_phyla_averages_onlyMixedLiquor.xlsx", row.names = F, to_save )

# averages of mature reactor only
prune.dat_mature<- subset_samples(prune.dat_top, Time_point %in% as.character(seq(32,42,1)))  # last 10 samples
prune.dat_mature<- subset_samples(prune.dat_mature, ! Sample_name %in% c("other_AGS"))
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_phyla_averages_onlyMatureAGS.xlsx", row.names = F, to_save )

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, prune.dat_mature)




### TOP Genera
data_mix<-subset_samples(data.genus.prop, Type=="Mixed liquor")

suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data_mix), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data_mix)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_mix)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_mix)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
fill_color_modified<-c("darkblue","brown4","springgreen2","wheat","red","yellow3","coral","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","chartreuse3","lightcoral","darkslategray3")
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_modified) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=3),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(-5,28,0,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified genera",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_onlyMixedLiquor.png",width=7,height=5.4, dpi=300)
ggsave(file="Results/Extra_edited_figures_for_presentations/TOP_genera_onlyMixedLiquor_v2.png",width=6.8,height=4.5, dpi=300)
dev.off()
# poster version
tabella2<-tabella[!grepl("other",tabella$Experiment_day), ]
ggplot(data=tabella2, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =8.5) +
  scale_y_continuous(expand=c(0,1) ) +
  scale_fill_manual(values=fill_color_modified) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=3),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 9 ),
  legend.position="bottom",
  legend.margin = margin(-7,28,0,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percentual abundance of clades", fill="")
ggsave(file="Results/Extra_edited_figures_for_presentations/TOP_genera_no_other_samples.png",width=6.7,height=4.5, dpi=300)



# means of TOP genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_genera_Averages_onlyMixedLiquor.xlsx", row.names = F, to_save)


# averages of mature reactor only
prune.dat_mature<- subset_samples(prune.dat_top, Time_point %in% as.character(seq(32,42,1)))
prune.dat_mature<- subset_samples(prune.dat_mature, ! Sample_name %in% c("other_AGS"))
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_genera_averages_onlyMatureAGS.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, prune.dat_mature))





###################### ABUNDANCES BAR PLOT (Only Influents and outputs) ##########################


### TOP Genera in influents
data_mix<-subset_samples(data.genus.prop, Type=="Influent")

suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data_mix), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data_mix)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_mix)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_mix)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_19) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=-3),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(9,28,0,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified genera (only influents)",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_onlyInfluents.png",width=7,height=5.4, dpi=300)
dev.off()




### TOP Genera in influents
data_mix<-subset_samples(data.genus.prop, Type=="Out")

suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data_mix), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data_mix)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_mix)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_mix)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_19) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=-3),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 9.1 ),
  legend.position="bottom",
  legend.margin = margin(9,28,1,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified genera (only discarded biomass)",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_onlyDiscarded.png",width=7,height=5.4, dpi=300)
dev.off()



######################## ALPHA DIVERSITY (EVERY SAMPLE) ##############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Type", color="Type")
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
  geom_boxplot(data=pAlpha$data, aes(x=Type, y=value, color=NULL, fill=Type), alpha=0.15, linewidth=0.15 ) +
  theme_classic2(base_size = 9) +
  scale_color_manual(values = c("Mixed liquor"="deepskyblue",
                                "Influent"="lightblue3", "Out"="coral",
                                "Direct extract"="green2")
  ) +
  scale_fill_manual(values = c("Mixed liquor"="deepskyblue",
                                "Influent"="lightblue3", "Out"="coral",
                                "Direct extract"="green2")
  ) +
  geom_text(aes(label=Sample_name), color="black", size=2.05, show.legend = FALSE) +
  labs(title="Alpha diversity (every sample)") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=6),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25)
        )
ggsave(file="Results/Alfa_diversity_GENUS_BOXPLOT.png", width = 5.2,height =4, dpi=300)


suppressWarnings(rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor) )




######################## ALPHA DIVERSITY (Only Mixed Liquor) ##############################

data.genus.temp<-subset_samples(data.genus, Type=="Mixed liquor")
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Experiment_day
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
mix_alpha$Sample_name<- gsub("X","",mix_alpha$Sample_name)  # added by the melt function
mix_alpha$ID <- gsub("[1-9].*","mixed",mix_alpha$Sample_name) # to connect the reactor points (without other reactors)
mix_alpha$Sample_name2 <- gsub("_","\n",mix_alpha$Sample_name)  # to avoid the cut in the display of the "other_" labels
mix_alpha$Time<- gsub("other_CAS","0",mix_alpha$Sample_name)
mix_alpha$Time<- gsub("other_AGS","185",mix_alpha$Time) # the other AGS has been distinguished with one day more than the last experiment day to maintain an order in the plot
mix_alpha$Time_color<- as.numeric(mix_alpha$Time)    # invisible x axis
mix_alpha$Time_axisX<- factor(as.character(mix_alpha$Time), levels =c (0, unique(mix_alpha$Time)[!unique(mix_alpha$Time) %in% c(0,185)], 185) )  # continuous value to print a shade of color

ggplot(data=mix_alpha, aes(y=value, x=Time_axisX, color=Time_color)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(0,50,100,150,184),
                        midpoint = 50) +
  geom_point(size =2, aes(shape=ID)) +
  geom_point(size =3.8, alpha=0.4, aes(shape=ID)) +
  #geom_line(aes(group=ID), size= 0.4, alpha= 0.4) + 
  geom_path(aes(group=as.character(ID)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic2(base_size = 9) +
  labs(y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  geom_text(aes(label=Sample_name2), color="black", size=1.8, vjust=-0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_blank(),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  )
ggsave(file="Results/Alfa_diversity_GENUS_MixedLiquor_along_the_time.png", width = 5.5,height =4, dpi=300)


suppressWarnings(rm(mix_alpha) )




########################### PCoA BETA DIV ##########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Time_point<-as.numeric(as.character(sample_data(data.prop.labels)$Time_point))
sample_data(data.prop.labels)$ID <- "this_project"
sample_data(data.prop.labels)$ID[grepl("other_",sample_data(data.prop.labels)$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 


# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color="Type", shape="ID") + # , color = "Sample_Type") +
  scale_color_manual(values = c("Mixed liquor"="deepskyblue",
                                "Influent"="lightblue3", "Out"="coral",
                                "Direct extract"="green2")
                     ) +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none", shape="none") +
  theme_classic(base_size = 12.5) +
  theme(title=element_text(size=10)) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_divers_Hellinger_on_genera_WITH_NAMES.png", width = 6.8, height = 5, dpi=300)


#### Only Mixed liquors
data_mix<-subset_samples(data.genus.prop, Type=="Mixed liquor")
data.prop.labels<-data_mix
sample_data(data.prop.labels)$ID <- "this_project"
sample_data(data.prop.labels)$ID[grepl("other_",sample_data(data.prop.labels)$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 
#sample_data(data.prop.labels)$Time_point<-as.numeric(as.character(sample_data(data.prop.labels)$Time_point))
sample_data(data.prop.labels)$Time_point<-as.character(sample_data(data.prop.labels)$Experiment_day)
sample_data(data.prop.labels)$Time_point<-gsub("other_CAS","0",sample_data(data.prop.labels)$Time_point)
sample_data(data.prop.labels)$Time_point<-gsub("other_AGS","184",sample_data(data.prop.labels)$Time_point)
sample_data(data.prop.labels)$Time_point<-as.numeric(sample_data(data.prop.labels)$Time_point)
sample_data(data.prop.labels)$Type2<-as.character(sample_data(data.prop.labels)$Type)
sample_data(data.prop.labels)$Type2[! sample_data(data.prop.labels)$Sample_name %in% c("other_AGS","other_CAS")] <- "Uniti"
sample_data(data.prop.labels)$Type2[sample_data(data.prop.labels)$Sample_name %in% c("other_AGS")] <- "Staccato1"
sample_data(data.prop.labels)$Type2[sample_data(data.prop.labels)$Sample_name %in% c("other_CAS")] <- "Staccato2"
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# with shade of color AND also exp day
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point", shape="ID") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,40,80,120,150,184),
                        midpoint = 25) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=as.character(Type2)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=8.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9.2),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none") +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day), color="gray30", vjust= -0.28, size=2.2, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/Beta_divers_Hellinger_on_genera_shade_and_TimePointText.png", width = 5.4, height = 4.5, dpi=300)

# again but only with shades of colors
sample_data(data.sqrt_prop)$Sample_name2<-gsub("L.*","",sample_data(data.sqrt_prop)$Sample_name) # to maintain infos about only the "others" sample
new_plot<-plot_ordination(data.sqrt_prop, ordBC, color = "Time_point", shape="ID") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,40,80,120,150,184),
                        midpoint = 25) +
  geom_point(size=4.5, alpha=0.3) +
  geom_point(size=2.5, alpha=1) +
  geom_path(aes(group=as.character(Type2)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=8.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9.2),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none") +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
new_plot + labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances")
ggsave(file="Results/Beta_divers_Hellinger_on_genera_shade.png", width = 5.4, height = 4.5, dpi=300)
new_plot + theme(plot.margin = margin(1,-1,1,1))
ggsave(file="Results/Extra_edited_figures_for_presentations/Beta_divers_Hellinger_on_genera_shade.png", width = 5, height = 4.2, dpi=300)
rm(new_plot)


#### Only Mixed liquor again, but without "Others"
data_mix<-subset_samples(data.genus.prop, Type=="Mixed liquor" & !Sample_name %in% c("other_CAS","other_AGS"))
data.prop.labels<-data_mix
sample_data(data.prop.labels)$ID <- "this_project"
sample_data(data.prop.labels)$ID[grepl("other_",sample_data(data.prop.labels)$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 
sample_data(data.prop.labels)$Time_point<-as.character(sample_data(data.prop.labels)$Experiment_day)
sample_data(data.prop.labels)$Time_point<-gsub("other_CAS","0",sample_data(data.prop.labels)$Time_point)
sample_data(data.prop.labels)$Time_point<-gsub("other_AGS","184",sample_data(data.prop.labels)$Time_point)
sample_data(data.prop.labels)$Time_point<-as.numeric(sample_data(data.prop.labels)$Time_point)
sample_data(data.prop.labels)$Type2<-as.character(sample_data(data.prop.labels)$Type)
sample_data(data.prop.labels)$Type2[! sample_data(data.prop.labels)$Sample_name %in% c("other_AGS","other_CAS")] <- "Uniti"
sample_data(data.prop.labels)$Type2[sample_data(data.prop.labels)$Sample_name %in% c("other_AGS")] <- "Staccato1"
sample_data(data.prop.labels)$Type2[sample_data(data.prop.labels)$Sample_name %in% c("other_CAS")] <- "Staccato2"
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
sample_data(data.sqrt_prop)$Sample_name2<-gsub("L.*","",sample_data(data.sqrt_prop)$Sample_name) # to maintain infos about only the "others" sample
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point", shape="ID") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,40,80,120,150,184),
                        midpoint = 25) +
  geom_point(size=4.5, alpha=0.3) +
  geom_point(size=2.5, alpha=1) +
  geom_path(aes(group=as.character(Type2)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day), color="gray30", vjust= -0.5, hjust= -0.2, size=2.5, show.legend = FALSE) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=8.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9.2),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none", color="none") +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
  ) + theme(plot.margin = margin(1,-1,1,1))
ggsave(file="Results/Extra_edited_figures_for_presentations/Beta_divers_Hellinger_only_this_reactor_ML.png", width = 4, height = 3, dpi=300)
rm(new_plot)



#### Direct vs storage
subset_target <- subset_samples(data.genus, Type %in% c("Mixed liquor","Direct extract") )
sample_data(subset_target)$Type<-gsub("Mixed liquor", "After storage", sample_data(subset_target)$Type) 
sample_data(subset_target)$Type<-gsub("Direct extract", "After collection", sample_data(subset_target)$Type) 
sample_data(subset_target)$Time_point_char <- as.character(sample_data(subset_target)$Time_point)
#subset_target<- subset_samples(subset_target,  Time_point_char %in% c("5","6","9","39","40"))
subset_target<- subset_samples(subset_target,  ! grepl("other", Sample_name))
data.prop.labels <- subset_target
sample_data(data.prop.labels)$Time_point<-as.numeric(as.character(sample_data(data.prop.labels)$Time_point))
sample_data(data.prop.labels)$Type2<-as.character(sample_data(data.prop.labels)$Type)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# with shade of color AND also exp day
plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  scale_color_manual(values = c("After storage"="deepskyblue",
                                "After collection"="green2")
  ) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_line(aes(group=as.character(Time_point_char)), col= "darkgray", linewidth = 0.3) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=8.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9.2),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       color="DNA extraction ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/Beta_divers_Hellinger_on_genera_Direct_vs_Storage.png", width = 5.4, height = 4.5, dpi=300)





#################### CORRELATING ABUNDANCES vs TIME (CENTERED LOG RATIO) ########################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
data.genus.temp <- subset_samples(data.genus.temp, Type %in% c("Mixed liquor") )
data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("other_AGS","other_CAS") )
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
order_of_samples<- metadata[levels(sample_data(data)$Sample_name), "Sample_name"]   # the sample data has been ordered already during the preparation
order_of_samples <- order_of_samples[order_of_samples %in% sample_names(data.genus.temp) ]
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
# Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")

write.csv2(Corr_results, "Results/Correlations/Kendall_correlation_CLR_bacteria_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )

con <- file("Results/Correlations/README_Only_certain_genera_have_been_tested_in_correlation.txt")
sink(con, append = T)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))



########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE 12 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_top<-significative_abundances_pos[1:12] # the table is already ordered
suppressWarnings(rm(top, others, tabella))


data.genus.temp <- subset_samples(data.genus.prop, Type %in% c("Mixed liquor") )
data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("other_AGS","other_CAS") )
{ tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% corr_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella$Genus<-gsub("Subgroup_7","Holophagae_Subgroup7",tabella$Genus)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella<-tabella[order(tabella$Sample_name) , ] # days
fill_color_modified<-c("springgreen2","firebrick2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","brown4","gray30","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_modified) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.3 ),
        legend.position="bottom",
        legend.margin = margin(1,36,0,0)) +
  guides(color=guide_legend(nrow=3)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  scale_y_continuous(breaks = seq( 0, max(tabella$Abundance), 2) , limits = c(0, max(tabella$Abundance)) ) +
  labs(x="Experiment day", y="Percent abundance",
       color="",
       title = "genera positely correlated with time\n according to Kendall correlation on CLR counts")
ggsave(file="Results/Correlations/Genera_CLR_positively_correlated_with_time.png",width=7,height=4.6, dpi=300) 
dev.off()



# THE 12 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_top<-significative_abundances_neg[(length(significative_abundances_neg)-11):length(significative_abundances_neg)]
suppressWarnings(rm(top, others, tabella))
data.genus.temp <- subset_samples(data.genus.prop, Type %in% c("Mixed liquor") )
data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("other_AGS","other_CAS") )
{
  tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update) 
  prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% corr_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella<-tabella[order(tabella$Sample_name) , ] # days
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_19) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 8.5 ),
        legend.position="bottom", legend.margin = margin(1,36,0,0)) +
  guides(color=guide_legend(nrow=3)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  scale_y_continuous(breaks = seq( 0, max(tabella$Abundance), 2) , limits = c(0, max(tabella$Abundance)) ) +
  labs(x="Experiment day", y="Percent abundance",
       color="",
       title = "genera negatively correlated with time \n according to Kendall correlation on CLR counts")
ggsave(file="Results/Correlations/Genera_CLR_negatively_correlated_with_time.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(prune.dat_top, tax_selected, tabella, corr_top))





#################### CORRELATING BACTERIA (CLR) vs MIXED LIQUOR MEASUREMENTS ########################

# STILL TOO MANY!
# # selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
# data.genus.temp<-data.genus.prop
# data.genus.temp <- subset_samples(data.genus.temp, Type %in% c("Mixed liquor") )
# data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("other_AGS","other_CAS") )
# min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
# tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
# ### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
# who<-as.data.frame(otu_table(data.genus.temp))
# who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
# who<-who[!rowSums(who)>=min_prevalence,]
# who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# only the most abundant genera
top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:19]
# switching to log ratio transformation
#data.genus.temp <- subset_samples(data.genus, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
data.genus.temp<-data.genus # NB: not prop this time
data.genus.temp <- subset_samples(data.genus.temp, Type %in% c("Mixed liquor") )
data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("other_AGS","other_CAS") )
data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
# data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
data.genus.temp<- prune_taxa(top,data.genus.temp)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])
row.names(abundances)<- gsub ("Candidatus_","Ca. ", row.names(abundances) )
row.names(abundances)<- gsub ("uncultured","uncult.", row.names(abundances) )

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sample_name
order_of_samples<- metadata[levels(sample_data(data)$Sample_name), "Sample_name"]   # the sample data has been ordered already during the preparation
order_of_samples <- order_of_samples[order_of_samples %in% sample_names(data.genus.temp) ]
abundances<-t(abundances)
abundances<-abundances[ order_of_samples, ]
meas<-meas[order_of_samples, ] # subsetting and re-orderdering the samples
# identical(row.names(meas), colnames(abundances))



#### Computing the correlations...
corr<-NULL
for( x in 1:length( colnames(abundances)) ){
  for( y in 1:length( colnames(meas)) ) {
    loop<- suppressWarnings( cor.test( as.numeric(abundances[, x]) , meas[[y]] , method = "kendall" ) )
    # NB: warnings due to the ties... this is normal!
    new_row<-c(colnames(abundances)[x], colnames(meas)[y] , round(as.numeric(loop$estimate),2), round(as.numeric(loop$p.value),6) )
    corr<-rbind(corr, new_row)
  }
}

corr<-as.data.frame(corr)
row.names(corr)<-NULL
colnames(corr)<-c("Bacteria","Measurements","Corr","pvalue")
corr$padj_holm<-p.adjust(corr$pvalue, method = "holm")
corr$Corr<-as.numeric(corr$Corr)
corr$Sign<-corr$padj_holm
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = Bacteria, y = Measurements, fill = Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.y=element_text(angle = -35, hjust = 1, vjust=1 , size= 9),
        axis.text.x=element_text(angle = 40, hjust = 1, vjust=1 , size= 9.5),
        plot.title = element_text(size=11.2),
        legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(0.7, 'cm'), legend.title = element_text(size=12)
  ) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =10, vjust=0.8) +
  labs(title = "Kendall correlations (CLR abundances vs measurements)", y= "ML measurements", x= "ML most abundant genera", 
       caption= "\n adjusted p-value (holm) lower than 0.05 are displayed through * sign")
ggsave(file="Results/Correlations/Correlations_Bacteria_VS_measurements_Kendall_Holm_usingCLRabundances.png", dpi=300, width = 7, height = 6.5)
write.csv2(corr, file="Results/Correlations/Correlations_Bacteria_VS_measurements_Kendall_Holm_table_usingCLRabundances.csv", quote = F, na="", row.names = F)




### testing "visually" these correlations ...

if(! identical( sample_names(data.genus.temp),row.names(abundances)) ){
  stop("Wait! Now there is a different sample order here... Maybe it's better to check why... ")
}

signif <- corr[corr$Sign=="*" , ]

pdf(file = "Results/Correlations/EXTRA_to_visually_test_correl_Bacteria_vs_Measurem_using_True_Values.pdf")
par(mfrow = c(2,3))
for(x in 1:length(row.names(signif))){
  target_bacteria<-signif[x, "Bacteria"]
  target_bacteria2<- gsub("uncult.","uncultured", target_bacteria, fixed=T )  # it's necessary to reset the names to search in the phyloseq object
  target_bacteria2<- gsub("Ca. ","Candidatus_", target_bacteria2, fixed=T )
  target_measur<-signif[x, "Measurements"]
  # a<-as.numeric( abundances[ ,target_bacteria ] )  # this are post clr... it's better to show the proportions!
  a <- suppressWarnings( subset_taxa(data.genus.temp, Genus == target_bacteria2 ) ) # the warning is about the tree
  a <- round( as.numeric(otu_table(a)), 2 )
  b<-meas [ ,target_measur ]
  plot( a~b , ylab= paste(target_bacteria,"%") , xlab=target_measur, col="red")
  title("")
  abline( lm(a~b) )
}
dev.off()




########################## FOCUS ON PAO AND GAO ###########################

data.genus.temp<-subset_samples(data.genus.prop, Type=="Mixed liquor")
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
tabella$Sample_name<- factor(tabella$Sample_name, levels = c("other_CAS", paste0("L", as.character(seq(0,42,1))), "other_AGS"))

ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=GroupPG)) + 
  facet_grid2( ~Type, scales = "free", space="free",
                strip = strip_nested(size="constant"))+
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
ggsave(file="Results/Abundances/PAO_GAO_comparison_MIXEDLIQUORS_v1.png",width=5.6,height=4.2,dpi=300)
dev.off()


# Focus on Accumulibacter and Competibacter
tabella$GroupPG2<- "temp" # extra line to allow a quick reset of this part (the factors hate the g_subs )
tabella$GroupPG2[ tabella$Genus %in% c("Candidatus_Accumulibacter","Candidatus_Competibacter") ] <- tabella$Genus [tabella$Genus %in% c("Candidatus_Accumulibacter","Candidatus_Competibacter") ]
tabella$GroupPG2[ ! tabella$Genus %in% c("Candidatus_Accumulibacter") &  tabella$Genus %in% PAO_list ] <- "Other PAOs"
tabella$GroupPG2[ ! tabella$Genus %in% c("Candidatus_Competibacter") &  tabella$Genus %in% GAO_list ] <- "Other GAOs"
tabella$GroupPG2<- gsub("Candidatus_","Ca. ", tabella$GroupPG2)
tabella$GroupPG2<- factor( tabella$GroupPG2 , levels= c("Ca. Accumulibacter", "Other PAOs", "Ca. Competibacter", "Other GAOs") )
ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=GroupPG2)) + 
  facet_grid2( ~Type, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("Ca. Accumulibacter"="coral", "Ca. Competibacter"="chartreuse2",
                             "Other PAOs"="red3","Other GAOs"="chartreuse4")) +
  theme(axis.text.x=element_text(angle=42, hjust=1,vjust=1, size=6.5),
        axis.title.y = element_text(size=9), 
        axis.text.y = element_text(size=7), 
        strip.text.x = element_text(size=8.5), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.18, "cm"),
        legend.text = element_text ( size = 7.6 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        #legend.position="bottom",
        legend.position="right",
        plot.margin = margin(1,1,-8,1),
        legend.margin =  margin(-5,2,0,-5)
  ) +
  scale_y_continuous(lim=c(0,40), breaks=seq(0,50,5)) +
  labs(x="", y="Percent abundance in the sample", 
       fill=""
       #title = paste("PAO / GAO abundances in the samples")
       )
ggsave(file="Results/Abundances/PAO_GAO_comparison_MIXEDLIQUORS_v2.png",width=5.6,height=3.2,dpi=300)
dev.off()  
  


# exporting the abundances
PAO_identified_phylo<- phyloseq(otu_table(PAO_identified_phylo),tax_table(PAO_identified_phylo)) # removing the two trees
GAO_identified_phylo<- phyloseq(otu_table(GAO_identified_phylo),tax_table(GAO_identified_phylo))
tot <- merge_phyloseq(PAO_identified_phylo, GAO_identified_phylo)
# means of every PAO and GAO
tot<- cbind.data.frame("mean in Mix Liq"=as.numeric(apply(otu_table(tot),1,mean)), "genus"= as.data.frame(tax_table(tot))[["Genus"]])
tot$Type <- ifelse( tot$genus %in% PAO_list, "PAO", "GAO")
write.xlsx(file = "Results/Abundances/PAO_GAO_names_and_averages_MIXEDLIQUORS.xlsx", row.names = F,
           tot
           )


suppressWarnings(rm(tot, tabella, tabella_GAO, tabella_PAO) )




####################### FOCUS ON AOB, NOB and N reducer ########################

data.genus.temp<-subset_samples(data.genus.prop, Type=="Mixed liquor")
AOB_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% AOB_list )
NOB_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% NOB_list)
Nred_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% NO2_r_list)
tabella_AOB<-psmelt(AOB_identified_phylo)
tabella_NOB<-psmelt(NOB_identified_phylo)
tabella_Nred<-psmelt(Nred_identified_phylo)
tabella_AOB$GroupPG <- "AOB"
tabella_NOB$GroupPG <- "NOB"
tabella_Nred$GroupPG <- "N_reducer"  # from MIDAS only Nitrite red, but added also other denitrifiers from articles
tabella<- rbind.data.frame(tabella_AOB, tabella_NOB, tabella_Nred)
tabella$GroupPG<-factor(tabella$GroupPG, levels=c("AOB","NOB", "N_reducer"))
tabella$Sample_name<- factor(tabella$Sample_name, levels = c("other_CAS", paste0("L", as.character(seq(0,42,1))), "other_AGS"))

ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=GroupPG)) + 
  facet_grid2( ~Type, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("AOB"="lightblue",
                             "NOB"="deepskyblue4",
                             "N_reducer"="darkblue" )) +
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
  scale_y_continuous(lim=c(0,40), breaks=seq(0,50,5)) +
  labs(x="", y="Percent abundance in the sample", 
       fill="",
       title = paste("AOB / NOB / denitrifiers abundances in the samples"))
ggsave(file="Results/Abundances/AOB_NOB_denitrif_comparison_MIXEDLIQUORS.png",width=5.6,height=4.2,dpi=300)
dev.off()


# exporting the abundances
AOB_identified_phylo<- phyloseq(otu_table(AOB_identified_phylo),tax_table(AOB_identified_phylo)) # removing the two trees
NOB_identified_phylo<- phyloseq(otu_table(NOB_identified_phylo),tax_table(NOB_identified_phylo))
Nred_identified_phylo<- phyloseq(otu_table(Nred_identified_phylo),tax_table(Nred_identified_phylo))
tot <- merge_phyloseq(AOB_identified_phylo, NOB_identified_phylo, Nred_identified_phylo)
# means of every AOB, NOB and denitrif
tot<- cbind.data.frame("mean in Mix Liq"=as.numeric(apply(otu_table(tot),1,mean)),
                       # otu_table(tot), 
                       "genus"= as.data.frame(tax_table(tot))[["Genus"]]
)
tot$Type <- tot$genus
tot$Type[tot$Type %in% AOB_list] <-"AOB"
tot$Type[tot$Type %in% NOB_list] <-"NOB"
tot$Type[tot$Type %in% NO2_r_list] <-"N_reducer"
write.xlsx(file = "Results/Abundances/AOB_NOB_names_and_averages_MIXEDLIQUORS.xlsx", row.names = F,
           tot 
)


suppressWarnings(rm(tot, tabella, tabella_NOB, tabella_AOB) )




#################### VEEN DIAGRAM ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-subset_samples(data.genus.prop, Sample_name!="other_CAS")
data.venn<-subset_samples(data.genus.prop, Type!="Out")

sample_data(data.venn)$Type <- as.character(sample_data(data.venn)$Type) # removing the "factor" status
sample_data(data.venn)$Type[sample_data(data.venn)$Sample_name=="other_AGS"] <- "other_AGS"

### abundance AND prevalence filter (0.5% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.5, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)
data.venn <- prune_taxa(taxa_sums(data.venn)>0, data.venn)

Influent<-subset_samples(data.venn, Type=="Influent")
Influent<-as.character(tax_table(prune_taxa(taxa_sums(Influent)>0, Influent))[,"Genus"])

Mixed_liquor<-subset_samples(data.venn, Type=="Mixed liquor")
Mixed_liquor<-as.character(tax_table(prune_taxa(taxa_sums(Mixed_liquor)>0, Mixed_liquor))[,"Genus"])

Direct_extract<-subset_samples(data.venn, Type=="Direct extract")
Direct_extract<-as.character(tax_table(prune_taxa(taxa_sums(Direct_extract)>0, Direct_extract))[,"Genus"])

other_AGS<-subset_samples(data.venn, Type=="other_AGS")
other_AGS<-as.character(tax_table(prune_taxa(taxa_sums(other_AGS)>0, other_AGS))[,"Genus"])


ONLY_IN_Mixed_liquor<- Mixed_liquor[! Mixed_liquor %in% Direct_extract & ! Mixed_liquor %in% Influent & ! Mixed_liquor %in% other_AGS ]
ONLY_IN_Mixed_liquor<- paste(ONLY_IN_Mixed_liquor, collapse = ", ")
head(ONLY_IN_Mixed_liquor)

ONLY_IN_Influent<- Influent[! Influent %in% Direct_extract & ! Influent %in% Mixed_liquor & ! Influent %in% other_AGS ]
ONLY_IN_Influent<- paste(ONLY_IN_Influent, collapse = ", ")
head(ONLY_IN_Influent)

ONLY_IN_Direct_extract<- Direct_extract[! Direct_extract %in% Influent & ! Direct_extract %in% Mixed_liquor & ! Direct_extract %in% other_AGS]
ONLY_IN_Direct_extract<- paste(ONLY_IN_Direct_extract, collapse = ", ")
head(ONLY_IN_Direct_extract)

ONLY_IN_other_AGS<- other_AGS[! other_AGS %in% Influent & ! other_AGS %in% Mixed_liquor & ! other_AGS %in% Direct_extract]
ONLY_IN_other_AGS<- paste(ONLY_IN_other_AGS, collapse = ", ")
head(ONLY_IN_other_AGS)


IN_other_AGS_AND_Mixed_liquor_AND_Influent <- Mixed_liquor[Mixed_liquor %in% other_AGS & Mixed_liquor %in% Influent ] 
# head(IN_other_AGS_AND_Mixed_liquor_AND_Influent)


con<-file("Results/Venn_Diagr_List_of_exclusive_bacteria_and_intersection.txt")
sink(con, append=TRUE)
cat("ONLY IN Direct_extract:", fill=TRUE)
cat(ONLY_IN_Direct_extract)
cat("\n\nONLY IN Influent:", fill=TRUE)
cat(ONLY_IN_Influent)
cat("\n\nONLY IN Mixed liquor:", fill=TRUE)
cat(ONLY_IN_Mixed_liquor)
cat("\n\n\nIntersection (In other AGS, Mixed liquor and Influent):", fill=TRUE)
cat(paste(IN_other_AGS_AND_Mixed_liquor_AND_Influent, collapse = ",  "))
sink()
close(con)


x<-list(Influent=Influent,"Mixed liquor "=Mixed_liquor," Direct extract"=Direct_extract,"Other AGS"=other_AGS)
ggvenn(x, stroke_size = 0.4, set_name_size = 2.8, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue","gold")) +
  theme(plot.title = element_text(size=8), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram \n(only genera with minimal abundance > 0.5% at least in three samples)\n")
ggsave(filename = "Results/Venn_Diagramm_EveryGroup.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()

x<-list(Influent=Influent,"Mixed liquor"=Mixed_liquor,"Other AGS"=other_AGS)
ggvenn(x, stroke_size = 0.4, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue","gold")) +
  theme(plot.title = element_text(size=8.2), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram \n(only genera with minimal abundance > 0.5% at least in three samples)\n")
ggsave(filename = "Results/Venn_Diagramm_3groups.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()



suppressWarnings(rm(ONLY_IN_Mixed_liquor, ONLY_IN_Influent, ONLY_IN_Direct_extract, ONLY_IN_other_AGS, 
                    x, con, Mixed_liquor, other_AGS, Influent, Direct_extract, data.venn, who))




############# ABUNDANCES OF GENERA ALWAYS FEATURED (VENN INTERSECTION) ######################

data_temp<-subset_taxa(data.genus.prop, Genus %in% IN_other_AGS_AND_Mixed_liquor_AND_Influent) # from the Venn part of this script

suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data_temp), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data_temp)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_19) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=6.5),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(-10,28,5,1),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant genera always featured in both current and other AGS reactor and influents",
       fill="",
       caption = " 'Others' includes every common genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera_AlwaysFeatured_Venn_Intersection.png",width=7.1,height=5.4, dpi=300)
ggsave(file="Results/Extra_edited_figures_for_presentations/TOP_genera_AlwaysFeatured_Venn_Intersection_v2.png",width=7.1,height=4.5, dpi=300)
dev.off()


rm(top, tabella_top, data_temp, tabella)



######################## DA WITH DESEQ2 (Type Direct extract vs Mixed liquor) ###########################

if(! "unfiltered_data" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned, subset_target))

subset_target <- subset_samples(data.genus, Type %in% c("Mixed liquor","Direct extract") )
sample_data(subset_target)$Type <-gsub (" ","_", sample_data(subset_target)$Type)
sample_data(subset_target)$Time_point_char <- as.character(sample_data(subset_target)$Time_point)
subset_target <- subset_samples(subset_target,  Time_point_char %in% c("5","6","9","39","40"))
data_pruned<- prune_taxa(taxa_sums(subset_target) > 10, subset_target) 
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Time_point_char + Type)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Type", "Direct_extract", "Mixed_liquor"))
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Direct_extract_vs_Mixed_liquor.csv"), row.names = F, quote=F, na = "")
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

# write.csv(Res_tot, file="Results/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/DA_differ_abundances_Storage_vs_Direct_DESeq2.xlsx", showNA = F, col.names = T)



############ PLOTTING THESE RESULTS ...

Table_tot$Type<-factor(Table_tot$Type, levels = unique(Table_tot$Type))
Table_tot<-Table_tot[order(match(Table_tot$Type, levels(Table_tot$Type))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Type)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Type)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

## to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Type)


#unique(tabella_g$Bacteria)
tabella_g$Type <- gsub("_"," ",tabella_g$Type)
tabella_g$Bacteria<-gsub("_marine","\nmarine",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("Candidatus_","Candidatus\n",tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Type)) +
  theme_classic2(base_size = 12) +
  scale_fill_manual(values = c("Mixed_liquor"="coral","Direct_extract"="chartreuse")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = "Genus")~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Type), size=1.8, alpha=1) +
  geom_point(aes(color=Type), size=3, alpha=0.5) +
  geom_text(aes(label=Time_point_char), size=2.5, vjust= -0.3, color="darkgray", show.legend = F) +
  geom_line(aes(group=Time_point_char), linewidth=0.2) +
  scale_x_discrete(labels=unique(levels(tabella_g$Type)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=10,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=11.5),
        axis.text.y = element_text(size=8),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt"),
        panel.grid.minor.y= element_blank()
  ) +
  theme(legend.margin=margin(0, 0, 0, 0), 
        legend.position = "none") +
  guides(fill="none")
# plot_1 +
#   labs(title= "Differently abundant genera between directly extracted and stored samples",
#        y="Percentual abundance", fill="", color="", x="") +
#   theme(legend.margin=margin(-5, 0, 0, 0), # overwriting
#         legend.position = "bottom")
# ggsave(filename = "Results/DA_Diff_abund_genera_Direct_vs_Storage_DeSeq2.png", width = 7.2, height = 5, dpi=300)
# dev.off()



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
tabella_2$Type <- gsub("_"," ",tabella_2$Type)
tabella_2$Bacteria<-gsub("_marine","\nmarine",tabella_2$Bacteria)

# plot_2<-ggplot(tabella_2, aes(x= Xaxis, y=Abundance, fill=Type)) +
#   theme_classic2(base_size = 12) +
#   scale_fill_manual(values = c("Mixed_liquor"="coral","Direct_extract"="chartreuse")) +
#   facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
#   geom_point(aes(color=Type), size=1.8, alpha=1) +
#   geom_point(aes(color=Type), size=3, alpha=0.5) +
#   geom_text(aes(label=Time_point_char), size=2.5, vjust= -0.3, color="darkgray", show.legend = F) +
#   geom_line(aes(group=Time_point_char), linewidth=0.2) +
#   scale_x_discrete(labels=unique(levels(tabella_g$Type)), expand=c(0,0.5)) +
#   # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
#   theme(strip.text.x=element_text(size=10,colour="black"), 
#         strip.switch.pad.wrap = unit(10,"line"),
#         # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.text = element_text(size=11.5),
#         axis.text.y = element_text(size=8),
#         panel.grid.major.y = element_line(size=0.2, color="gray"),
#         plot.title= element_text(size=12),
#         panel.spacing.x = unit(1, "pt"),
#         panel.grid.minor.y= element_blank()
#   ) +
#   theme(legend.margin=margin(0, 0, 0, 0), 
#         legend.position = "none") +
#   guides(fill="none")
# plot_2 +
#   labs(title= "Differently abundant genera between directly extracted and stored samples",
#        y="Percentual abundance", fill="", color="", x="") +
#   theme(legend.margin=margin(-5, 0, 0, 0), # overwriting
#         legend.position = "bottom")
# #ggsave(filename = "Results/DA_Diff_abund_families_Direct_vs_Storage_DeSeq2.png", width = 7.2, height = 5, dpi=300)
# dev.off()


# # now a unique plot
# head_plot<-plot_1 + # refining first plot with the title
#   labs(title= "Differently abundant genera between directly extracted and stored samples",
#        y="Percentual abundance", fill="", color="", x="") +
#   theme(plot.title = element_text(size=12))
# png(filename = "Results/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants.png", width = 2900, height = 2600, res=300)
# grid.arrange(head_plot,
#              plot_2 + labs(title= "", y="Proportional abundance", x=""),
#              nrow=1)
# # + plot3 ecc
# dev.off()

tabella3 <-rbind.data.frame(tabella_g, tabella_2)   # few results, it's better to re-join them ...
tabella3$Type<-gsub("Mixed liquor", "Storage", tabella3$Type) 
ggplot(tabella3, aes(x= Xaxis, y=Abundance, fill=Type)) +
  theme_classic2(base_size = 12) +
  scale_color_manual(values = c("Storage"="deepskyblue",
                               "Direct extract"="green2") 
  ) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus"))~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Type), size=1.8, alpha=1) +
  geom_point(aes(color=Type), size=3.8, alpha=0.5) +
  geom_text(aes(label=paste0(Experiment_day,"th")), size=2.5, vjust= -0.3, color="gray38", show.legend = F) +
  geom_line(aes(group=Time_point_char), linewidth=0.2) +
  scale_x_discrete(labels=unique(levels(tabella_g$Type)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=8.2,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=11.8),
        axis.text.y = element_text(size=7.2),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        panel.grid.minor.y= element_blank(),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt")
  ) +
  guides(fill="none") +
  labs(title= "Differently abundant genera between directly extracted and stored samples",
       y="Percentual abundance", fill="", color="", x="") +
  theme(legend.margin=margin(-5, 0, 0, 0),
        legend.position = "bottom")
ggsave(filename = "Results/DA_Diff_abund_bacteria_Direct_vs_Storage_DeSeq2.png", width = 7, height = 5, dpi=300)




######################## DA WITH ALDEx2 (Type Direct extract vs Mixed liquor) ###########################


if(! "proof1" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
}


suppressWarnings(rm(data_pruned, data.genus_pruned, subset_target))

subset_target <- subset_samples(data.genus, Type %in% c("Mixed liquor","Direct extract") )
sample_data(subset_target)$Type <-gsub (" ","_", sample_data(subset_target)$Type)
sample_data(subset_target)$Time_point_char <- as.character(sample_data(subset_target)$Time_point)
subset_target <- subset_samples(subset_target,  Time_point_char %in% c("5","6","9","39","40"))
data_pruned<- prune_taxa(taxa_sums(subset_target) > 10, subset_target) 

sample_data(data_pruned)[["Type"]]<-factor(sample_data(data_pruned)[["Type"]], levels = c("Mixed_liquor","Direct_extract"))
# Mixed liquor will be the base level --> Denominator of fold change

Table_tot<-NULL
Res_tot<-NULL
### Starting the analysis on each taxon level
for( t in c("Genus","Family","Order","Class","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # using the "manual" pipeline of Aldex2 to better specify the design ...
  mm <- model.matrix(~ Type + Time_point_char, as(sample_data(d),"data.frame") ) # name of the sample data for last
  set.seed(1)
  aldx <- aldex.clr(otu_table(d),  mm,
                    mc.samples=128,   #"DMC"
                    # according to ALDEx2 vignette, the number of samples in the smallest group multiplied by the number of DMC be equal at least to 1000
                    # length(which(sample_data(data)[["Type"]]=="Direct_extract"))
                    denom="all", # every feature as denominator
                    verbose=T
  )
  # aldx@reads # the original counts _ NB: added 0.5 to every number
  
  aldx2 <- aldex.glm(aldx, verbose = T)
  #### THEN, the summary of the results
  aldx3 <- aldex.glm.effect(aldx, # calculates the effect size for every binary variable in the mm
                            CI = T # confidence interval of the effect size
  ) 
  aldx_final <- data.frame(aldx2,aldx3) 
  
  p_val<-aldx_final$TypeDirect_extract.pval
  # p_adj<-p.adjust(p_val, method="holm")
  p_adj<-p.adjust(p_val, method="BH")
  
  #### Moreover, the vignette suggests an effect size cutoff of 1 or greater be used when analyzing HTS 
  # eff_size<-aldx_final$TypeDirect extract.effect 
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
      index<- which(colnames(Table_DE)=="Domain") : which(colnames(Table_DE)==t) # from : to
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

# columns_to_remove<- grepl("Intercept",colnames(Res_tot)) | grepl("sequencing_depth",colnames(Res_tot)) | grepl("pval.holm",colnames(Res_tot)) 
# colnames(Res_tot)<-gsub("TypeDirect extract","Type_Direct extract",colnames(Res_tot))
# Res_tot<-Res_tot[ , !columns_to_remove ]
# #View(Res_tot)
# write.csv(Res_tot, file="Results/DA_Aldex2/Every_result_Aldex2.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/DA_DESeq2/Every_result_Aldex2.xlsx", showNA = F, col.names = T)
# 
# # plot
# ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Type)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   scale_fill_manual(values=c("Direct extract"="coral","Mixed liquor"="deepskyblue")) +
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
#         axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
#         axis.text.y = element_text(size=12), plot.title= element_text(size=18),
#         panel.grid.minor.y= element_blank() ) +   
#   scale_x_discrete(expand=c(-0.2, 1)) +
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1, seq(2,max(Table_tot$Abundance+3),1))) +
#   labs(title= "Differently abundant taxa", y="Proportional Abundance", 
#        fill="Tumor localisation", x="")
# ggsave(filename = "Results/DA_Aldex2/DA_Type_every_result.png", width = 8, height = 7, dpi=300)
# dev.off()

# # removing redundant results
# Redund<- c("")
# Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
# DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Type)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus","Species")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.9, size= 0.25, alpha= 0.2, outlier.alpha = 0) +
#   geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.6),
#              aes(color=Type), size= 0.3, alpha= 0.65) +
#   scale_fill_manual(values=c("Direct extract"="coral","Mixed liquor"="deepskyblue")) +
#   scale_color_manual(values=c("Direct extract"="coral","Mixed liquor"="deepskyblue2")) +
#   guides( fill=guide_legend(nrow=1) ) +
#   theme_classic2(base_size = 11) + 
#   theme(strip.text.x=element_text(size=12,colour="black"),
#         legend.margin=margin(-10, 0, 0, 0), legend.position="bottom" ,
#         axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=9.4), 
#         axis.text.y = element_text(size=9),
#         plot.margin = margin(5,2,5,22),
#         plot.title= element_text(size=15) ,
#         legend.key.size=unit(0.8,"cm"), 
#         legend.text=element_text(size=14),
#         panel.grid.major.y = element_line(size=0.12, color="gray")
#   ) + 
#   scale_x_discrete(expand=c(-0.2, 1)) +
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1, seq(2,max(Table_tot$Abundance+3),1))) +
#   labs(title= "Differently abundant Taxa", y="Percent abundance", 
#        fill="Tumor localisation", x="")
# 
# DE_plot
# ggsave(filename = "Results/DA_Aldex2/DA_Type_WITHOUT_redundants.png", width = 7, height = 6, dpi=300)
# dev.off()

system(" echo 'No results for the comparison of Direct extract vs Storage according to ALDEx2' > Results/DA_Nessuna_diff_abbondanza_Storage_vs_Direct_ALDEx2.txt ")




######################## EXTRA: CHECKING THE MEASUREMENTS #####################

# hist(meas$`sCOD_removal efficiency`, breaks=10)
# hist(meas$`PO4_removal efficiency`, breaks=10)
# hist(meas$tCOD_influent, breaks=10)
# hist(meas$NH4_influent, breaks=10)

samples_ordered<-metadata[metadata$Type %in% "Mixed liquor", "Sample_name"]
samples_ordered<-samples_ordered[ !grepl("other", samples_ordered) ]
meas_check<-meas[ samples_ordered, ]


# removal efficiency
target<-meas_check[ , grepl("removal", colnames(meas_check))]
target<-target[! apply(target, MARGIN = 1, function(x) length(which(is.na(x))>0) ), ] # this removes the row with at least a NA
to_plot<-prcomp(target)
plot_stats<-as.data.frame(summary(to_plot)$importance)
colors_shade<-as.numeric(metadata[metadata$Sample_name %in% row.names(target), "Experiment_day"])
TypeLine <-metadata[metadata$Sample_name %in% row.names(target), "Type" ]
ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade) 
       ) +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,40,80,120,150,184),
                        midpoint = 25) +
  # to change the legend labels (this can't be done beforehand, due to the slashes in the name that can cause problem in DESeq2! Moreover... R1 is always before then R2 --> labels wrote accordingly)
  theme_classic( )  +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  # geom_path(aes(group=as.character(TypeLine)), col= "darkgray", linewidth = 0.3,
  #           arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  # ) +
  # geom_text(mapping = aes(label=gsub("R[1-9]_","",row.names(as.data.frame(to_plot$x)))), color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
  geom_text(mapping = aes(label=gsub("R[1-9]_","",paste0(colors_shade,"th"))), color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
  guides(color="none") +
  labs(
        x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
        y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
        title="PCA of removal efficiency measurements"
  )
ggsave(paste0("Results/Checking_Measurements/PCA_removal_efficiency.png"), width= 5.4, height = 4.5, dpi=300)


# influent characteristics
target<-meas_check[ , grepl("influent", colnames(meas_check))]
target<-target[! apply(target, MARGIN = 1, function(x) length(which(is.na(x))>0) ), ] # this removes the row with at least a NA
to_plot<-prcomp(target)
plot_stats<-as.data.frame(summary(to_plot)$importance)
colors_shade<-as.numeric(metadata[metadata$Sample_name %in% row.names(target), "Experiment_day"])
TypeLine <-metadata[metadata$Sample_name %in% row.names(target), "Type" ]
ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade) 
) +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,40,80,120,150,184),
                        midpoint = 25) +
  # to change the legend labels (this can't be done beforehand, due to the slashes in the name that can cause problem in DESeq2! Moreover... R1 is always before then R2 --> labels wrote accordingly)
  theme_classic( )  +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  # geom_path(aes(group=as.character(TypeLine)), col= "darkgray", linewidth = 0.3,
  #           arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  # ) +
  # geom_text(mapping = aes(label=gsub("R[1-9]_","",row.names(as.data.frame(to_plot$x)))), color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
  geom_text(mapping = aes(label=gsub("R[1-9]_","",paste0(colors_shade,"th"))), color="grey35", vjust= -0.3, size= 3, lineheight = 0.7, alpha= 0.45) +
  guides(color="none") +
  labs(
    x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
    y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
    title="PCA of influent measurements"
  )
ggsave(paste0("Results/Checking_Measurements/PCA_influent_characteristics.png"), width= 5.4, height = 4.5, dpi=300)



# LINE PLOT OF REMOVAL EFFICIENCIES THROUGHOUT THE TIME
target<-meas_check[ , grepl("removal", colnames(meas_check))]
target<-target[! apply(target, MARGIN = 1, function(x) length(which(is.na(x))>0) ), ] # this removes the row with at least a NA
Name<-paste0(metadata[metadata$Sample_name %in% row.names(target), "Experiment_day"],"th")
target$Name<-Name
tabella<-melt(target)
tabella$Name<-factor(tabella$Name, levels= Name)  # the "Names" are the days --> they have to be ordered along X axis --> as factor
tabella$variable<-gsub("_removal efficiency","",tabella$variable)
ggplot(data=tabella, aes(y=value, x=Name, color=variable)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=variable), size =2) +
  geom_point(aes(color=variable), size =3.8, alpha=0.4) +
  geom_line(aes(group=variable), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_5) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 8.5 ),
        legend.position="bottom",
        legend.margin = margin(1,36,0,0)) +
  guides(color=guide_legend(nrow=1)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  # scale_y_continuous(breaks = seq( 0, max(tabella$Abundance), 2) , limits = c(0, max(tabella$Abundance)) ) +
  labs(x="Experiment day", y="Efficiency",
       color="",
       title = "C,N and P compounds removal efficiency throughout time")
ggsave(file="Results/Checking_Measurements/Removal_efficiency_throughout_time.png",width=7,height=5, dpi=300) 
dev.off()





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

