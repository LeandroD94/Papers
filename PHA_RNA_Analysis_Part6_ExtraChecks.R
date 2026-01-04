################ PREPARING THE ENVIRONMENT ###############

{
  library("ggplot2")
  library("ggvenn")
  library("reshape2")
  library("ggh4x")
  library("ecodist")
  library("Hmisc")
  suppressPackageStartupMessages( library("vegan") ) # it has a verbose loading!
}

options(scipen=100)

dir.create("RNA_PHA_Results")

fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
fill_color_30<-c("wheat3","deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "aquamarine2",
                 "bisque2","cyan","yellow3","brown",
                 "springgreen4", "firebrick3",
                 "grey15","lightblue1","lightgreen","orange2",
                 "darkorchid", "darkslategray3","blue1", "violet", 
                 "yellow", "pink3","yellow4", "chocolate4", "red","darkgreen",
                 "lightgrey","coral2","deepskyblue2")
fill_color_30_BIS<-c("chartreuse","deeppink", "firebrick3",
                     "orange","darkmagenta", 
                     "darkblue","grey60", "grey95",
                     "bisque","cyan","yellow3","brown4",
                     "springgreen2", "grey20", 
                     "deepskyblue2","darkgreen","orange3",
                     "darkorchid", "lightblue1" , "violet", "blue",
                     "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                     "lightgrey", "red","darkslategray3")


load(file="Tabelle_RNA_PHA_after_filters.RData")




######################## COMPUTING THE PROPORTIONS ##############################

if(! "proof_filter" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
}

### Proportions computed after the filters and with the unmapped
prop_with_unmap <- apply(table_abs, MARGIN = 2, function(x) x/sum(x))
prop_with_unmap <- prop_with_unmap * 100

### Proportions computed after the filters BUT without the unmapped (for statistical tests or alike)
prop_with_NO_unmap <- apply(table_abs[! row.names(table_abs) %in% "UNMAPPED", ], MARGIN = 2, function(x) x/sum(x))
prop_with_NO_unmap <- prop_with_NO_unmap * 100

### Gene proportions
prop_GENE_with_unmap <- apply(table_gene, MARGIN = 2, function(x) x/sum(x))
prop_GENE_with_unmap <- as.data.frame(prop_GENE_with_unmap) * 100




############ EXTRA: CANDIDATUS ACCUMULIB AND EUKARYOTES ABUNDANT mRNAs (!?) #################

### What are Rotaria, Candidatus Acc. and Diploscapter doing???
table <- table_abs
sub_dict <- dictionary [ grepl("Accumulibacter|Rotaria|Diploscapter",dictionary$Taxon) , ]
table <- table[ sub_dict$Entry , ]
table <- cbind.data.frame( table, sub_dict )
table$Genus <-  gsub("^Accumuli","Ca_Accumuli", table$Taxon)
table$Genus <- gsub(" .*","", table$Genus)
table <- table[ order(rowSums(table[ ,!colnames(table)%in%c("Taxon","Genus","Gene","Entry","Description","Pathway")]),decreasing=TRUE) , ]
AVERAGES_R1<- round(apply( table[, grepl("SR[1-9]",colnames(table))], MARGIN = 1, mean) ,0)
AVERAGES_R2<- round(apply( table[, grepl("FR[1-9]",colnames(table))], MARGIN = 1, mean) ,0)
to_save <- cbind.data.frame( table, AVERAGES_R1, AVERAGES_R2 )

write.csv2( to_save , 
            file = "RNA_PHA_Results/Diplo_Rotaria_Accumuli_transcripts.csv", row.names = F, quote = F)
rm(table)




############### EXTRA: BAR PLOTS OF THE PHAGE ABUNDANT RNAs ##################

# table<-table_abs
table <- as.data.frame(prop_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

# subsetting only phages
phages<-dictionary [ grepl( "Tequatro|phage|virales|virus|viridae|viriae|viricota|viricetes" , dictionary$Taxon ) , "Entry"]

top <- names(sort(rowSums(table[row.names(table) %in% phages , ]), decreasing=TRUE))[1:20]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("DNA invertase Pin-like site-specific DNA recombinase","DNA invertase Pin-like recomb",infos,fixed=T)
infos<-gsub("TtsA-like Glycoside hydrolase fam 108 domain-containing prot","TtsA-like prot with glycoside hydrol 108 domain",infos,fixed=T)
infos<-gsub("Candidatus","Ca.",infos,fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
table_to_plot$Protein <- gsub("Uncharact", "Uncharacterized", table_to_plot$Protein)
table_to_plot$Protein <- gsub(" prot", "", table_to_plot$Protein)
table_to_plot$Protein <- gsub(" 108", "", table_to_plot$Protein)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_19) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        plot.title= element_text(size=6),
        legend.key.height = unit(0.14, "cm"),
        legend.key.width = unit(0.28, "cm"),
        legend.spacing.x = unit(0.29, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,3,-34),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=9)) +
  labs(x="", y="Percent counts",
       fill="",
       title = "Phages RNA",
       caption = "NB: different depth among samples --> more 'others' RNA if more depth!")
ggsave("RNA_PHA_Results/About_phages.png", width= 6.4, height = 5, dpi=300)




################### EXTRA: CHECKING THE EUKARYOTIC "CONTAMINATION" ########################

# taxonomy <- read.delim("/media/matteo/SSD4/reference/Bwt2_SILVA138_2_RNA_both_SSU_LLU/taxonomy.txt", sep=";", header = F)
taxonomy <- read.delim("databases_dictionary/taxonomy.txt", sep=";", header = F)
# txt with taxon code (space) species name (;) Domain  (obtained from SILVA fasta headers)

taxonomy$Domain <- gsub(".* ","",taxonomy$V1)
taxonomy$V1 <- gsub(" .*","",taxonomy$V1)
colnames(taxonomy) <- c("Code","Species","Domain")

# NB: the expected targets are txt files composed of one column (one for each paired file), hence if they are all compressed in a unique gz it will be necessary to decompress it before!
targets <- dir("FASTQ_check/BBDuk/rRNA_alignments/")
x=targets[1]

list_tables <- list()
all_taxa <- NULL

for( x in 1:length(targets) ){
  cat(paste("\rImporting sample",x,"on",length(targets),"..."))
  file_name <- targets[x]
  sample_name <- gsub("_SILVA.txt","", file_name)
  import <- read.delim(file=paste0("FASTQ_check/BBDuk/rRNA_alignments/",file_name))
  import_count <- table(import)
  import_table <- cbind.data.frame( "Taxon"=names(import_count) , "Count" = as.numeric(import_count) )
  list_tables[[x]] <- import_table
  names(list_tables)[x] <- sample_name
  all_taxa <- c(all_taxa, import_table$Taxon )
  all_taxa <- unique(all_taxa)
}

complete_table <- cbind.data.frame( "Taxon"= all_taxa )
for( x in 1: length(list_tables) ){
  cat(paste("\rParsing sample",x,"on",length(targets),"..."))
  taxa_in_target <- as.character(list_tables[[x]][["Taxon"]])
  missing_in_target <- all_taxa[!all_taxa%in%taxa_in_target]
  table_missing <- cbind.data.frame("Taxon"=missing_in_target , "Count"=rep(0))
  table_target <- rbind.data.frame(list_tables[[x]] , table_missing)   # same taxa vector
  row.names(table_target)<-table_target$Taxon
  table_target <- table_target[ all_taxa , ]    # same order for every sample
  complete_table <- cbind.data.frame( complete_table , "new"=table_target$Count )
  colnames(complete_table)[colnames(complete_table)=="new"] <- names(list_tables)[x]
}

complete_table$Average <- round(rowMeans( complete_table[!colnames(complete_table)%in%c("Taxon","Average")] ), 0)

selected_taxa <- taxonomy
selected_taxa <- selected_taxa[!duplicated(selected_taxa$Code) , ]  # few codes are replicated in the database (but related to the same taxon)
row.names(selected_taxa) <- selected_taxa$Code
selected_taxa <- selected_taxa[ complete_table$Taxon , ]
complete_table <- cbind.data.frame( complete_table , "Species"=selected_taxa$Species , "Domain"=selected_taxa$Domain )

complete_table <- complete_table[ order(complete_table$Average, decreasing = T), ]
# head(complete_table, n=20)

proper_names <- metadata$Sample_name
names(proper_names) <- metadata$FASTQ_ID
proper_names <- proper_names[colnames(complete_table) [colnames(complete_table)%in%metadata$FASTQ_ID] ]
colnames(complete_table) [colnames(complete_table)%in%metadata$FASTQ_ID] <- as.character(proper_names)

write.csv2(complete_table, file="FASTQ_check/Count_of_rRNA_assignments.csv", row.names = F )




################### EXTRA: CHECKING THE ncRNA besides the ribosomial #########################
# NB: this part of the script works only if the related optional step has been performed during the files processing

# taxonomy <- read.delim("/media/matteo/SSD4/reference/prokaryotes_ncRNA/Dictionary_ncRNA_microb_scfbio_24_01_24.csv", sep=";", header = F)
taxonomy <- read.delim("databases_dictionary/Dictionary_ncRNA_microb_scfbio_24_01_24.csv.gz", sep=",", header = T)
# csv with code, description and reference taxon as columns, manually build from the downloaded ScfBio ncRNA fasta

taxonomy$taxon <- gsub("uncultured ","uncultured_",taxonomy$taxon)
taxonomy$taxon <- gsub("Candidatus ","Candidatus_",taxonomy$taxon)
taxonomy$taxon <- gsub(" .*","",taxonomy$taxon) # going to Genus level, it is more accurate!

# NB: the expected targets are txt files composed of two columns (one for each paired file), with the second column being the alignment name in the db
targets <- dir("FASTQ_check/BBDuk/ncRNA_alignments/")
x=targets[1]

list_tables <- list()
all_taxa <- NULL

for( x in 1:length(targets) ){
  cat(paste("\rImporting sample",x,"on",length(targets),"..."))
  file_name <- targets[x]
  sample_name <- gsub("_ncRNA.tsv","", file_name)
  import <- read.delim(file=paste0("FASTQ_check/BBDuk/ncRNA_alignments/",file_name), 
                       sep="\t", header = F,
                       quote = ""   # because few times the ScfBio dicionary has quotes in the description!
  )
  import <- import[ , "V2", drop=F]
  import$V2 <- gsub("ncRNA_","", import$V2)
  import$V2 <- gsub("_.*","", import$V2)
  import_count <- table(import)
  import_table <- cbind.data.frame( "Code"=names(import_count) , "Count" = as.numeric(import_count) )
  list_tables[[x]] <- import_table
  names(list_tables)[x] <- sample_name
  all_taxa <- c(all_taxa, import_table$Code )
  all_taxa <- unique(all_taxa)
}

complete_table <- cbind.data.frame( "Code"= all_taxa )
for( x in 1: length(list_tables) ){
  cat(paste("\rParsing sample",x,"on",length(targets),"..."))
  taxa_in_target <- as.character(list_tables[[x]][["Code"]])
  missing_in_target <- all_taxa[!all_taxa%in%taxa_in_target]
  table_missing <- cbind.data.frame("Code"=missing_in_target , "Count"=rep(0))
  table_target <- rbind.data.frame(list_tables[[x]] , table_missing)   # same taxa vector
  row.names(table_target)<-table_target$Code
  table_target <- table_target[ all_taxa , ]    # same order for every sample
  complete_table <- cbind.data.frame( complete_table , "new"=table_target$Count )
  colnames(complete_table)[colnames(complete_table)=="new"] <- names(list_tables)[x]
}

complete_table$Average <- round(rowMeans( complete_table[!colnames(complete_table)%in%c("Code","Average")] ), 0)

selected_taxa <- taxonomy
selected_taxa <- selected_taxa[!duplicated(selected_taxa$TXNAME) , ]  
row.names(selected_taxa) <- selected_taxa$TXNAME
selected_taxa <- selected_taxa[ complete_table$Code , ]
complete_table <- cbind.data.frame( complete_table , "GENEID"=selected_taxa$GENEID , # "Type"=selected_taxa$UniProt_Swiss_Prot_ID ,
                                    "Description"=selected_taxa$protein_fuction , "Taxon"= selected_taxa$taxon )

complete_table <- complete_table[ order(complete_table$Average, decreasing = T), ]
# head(complete_table, n=20)

proper_names <- metadata$Sample_name
names(proper_names) <- metadata$FASTQ_ID
proper_names <- proper_names[colnames(complete_table) [colnames(complete_table)%in%metadata$FASTQ_ID] ]
colnames(complete_table) [colnames(complete_table)%in%metadata$FASTQ_ID] <- as.character(proper_names)

write.csv2(complete_table, file="FASTQ_check/Count_of_ncRNA_assignments.csv", row.names = F )



##### LET'S PLOT THEM
to_plot <- complete_table[ , grepl("SR|FR|Description|Taxon",colnames(complete_table)) ]
to_plot <- to_plot[to_plot$Description!="",]
#### solving few almost similar names before the glomming
# unique(complete_table$Description)
{
  to_plot$Description<-gsub("_"," ",to_plot$Description)
  to_plot$Description<-gsub("; predicted by Infernal","6S RNA",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("6S RNA regulator of RNA polymerase holoenzyme","",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("^6[S-s]","6S",to_plot$Description)
  to_plot$Description<-gsub("^6S$","6S RNA",to_plot$Description)
  to_plot$Description<-gsub("RNaseP bact a","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("Bacterial RNase P class A","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("bacterial RNase P, class A","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("[R-r]ibonuclease P","RNA component of RNaseP",to_plot$Description)
  to_plot$Description<-gsub("ribonuclease P RN.*","RNA component of RNaseP",to_plot$Description)
  to_plot$Description<-gsub("RNase P RNA component precursor RnpB","RNaseP's RNA precursor RnpB",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNase P, RNA component precursor RnpB","RNaseP's RNA precursor RnpB",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNase P" ,"RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("rnasep","RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("^RNaseP RNA$","RNA component of RNaseP",to_plot$Description)
  to_plot$Description<-gsub("RNaseP, RNA component; M1 RNA","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNaseP, M1 RNA component","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNA component of RNaseP RNA component","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNA component of RNaseP RNA","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNaseP RNA component","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("RNaseP RNA subunit","RNA component of RNaseP",to_plot$Description, fixed = T)
  to_plot$Description<-gsub("CrcZ","CrcZ sRNA - Carbon Catabolism Repression",to_plot$Description, fixed = T) # see doi: 10.4161/rna.23019
}
# unique(to_plot$Description)
to_plot$DescriptionTaxon <- paste0(to_plot$Description," (",to_plot$Taxon,")")
to_plot$Description<-NULL
to_plot$Taxon<-NULL
to_plot <- aggregate(.~DescriptionTaxon, to_plot, FUN=sum)
# top 20 most abundant
to_plot <- to_plot[ order(rowMeans(to_plot[,!colnames(to_plot)%in%"DescriptionTaxon"]) , decreasing = T) , ]
to_plot <- to_plot[1:20, ]
# computing percentages among themself
to_plot[,!colnames(to_plot)%in%"DescriptionTaxon"] <- apply(to_plot[,!colnames(to_plot)%in%"DescriptionTaxon"], MARGIN = 2, function(x) x/sum(x))
to_plot[,!colnames(to_plot)%in%"DescriptionTaxon"] <- to_plot[,!colnames(to_plot)%in%"DescriptionTaxon"] * 100
table_to_plot<-melt(to_plot, id="DescriptionTaxon")
table_to_plot$DescriptionTaxon <- gsub("tmRNA", "tmRNA (transfer-messenger RNA)", table_to_plot$DescriptionTaxon)
table_to_plot$Reactor <- gsub("^SR.*", "R1", table_to_plot$variable)
table_to_plot$Reactor <- gsub("^FR.*", "R2", table_to_plot$Reactor)
table_to_plot$variable <- gsub(".*R2.*", "63th day", table_to_plot$variable)
table_to_plot$variable <- gsub(".*R5.*", "81th day", table_to_plot$variable)
table_to_plot$variable <- gsub(".*R7.*", "95th day", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=DescriptionTaxon )) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =8.4) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values= c("navyblue","darkslategray3",fill_color_19[5:18],"orange","royalblue2","green","black") ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=9.5),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.23, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.45 ),
        legend.margin = margin(-15,0,2,-38),
        legend.position="bottom",
        plot.margin = margin(2,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       title = "Most abundant ncRNAs (excluding rRNA)"
  )
ggsave("RNA_PHA_Results/Most_abundant_NonCodingRNA.png", width= 4.92, height = 4, dpi=300)




############## EXTRA: DEFINING THE QUORUM SENSING RELATED ENTRIES #################

# PHA related genes, references: 1) doi 10.1016/j.watres.2023.119814 
#                                2) doi 10.3390/microorganisms9020239 
#                                3) doi 10.3390/microorganisms10061239 
#                                4) doi 10.1016/j.arabjc.2024.106050  
#                                5) doi 10.1016/j.jece.2025.116522
#                                6) doi 10.1128/JB.01304-08
#                                7) doi 10.1096/fasebj.2020.34.s1.05447
#                                8) doi.org/10.1007/s11274-025-04512-6
#                                9) doi.org/10.3389/fmars.2022.808132
#                               10) doi.org/10.1016/j.biotechadv.2025.108733
#                               11) doi.org/10.1016/j.enzmictec.2013.12.010 (penycillin acyclase)
only_QuorumGenes_IDs<-dictionary [ grepl("[Q-q]uorum", dictionary$Description ) |
                                     grepl("^QS", dictionary$Description ) |
                                     grepl(" QS ", dictionary$Description ) |
                                     grepl(" QS-", dictionary$Description, fixed=T ) |
                                     grepl("communicat", dictionary$Description ) |
                                     grepl("acyl-HSL", dictionary$Description , fixed=T ) |
                                     grepl("HSL", dictionary$Description , fixed=T ) |
                                     grepl("AHL", dictionary$Description , fixed=T ) |
                                     grepl("OHL", dictionary$Description , fixed=T ) |
                                     grepl("^Acyl carrier prot", dictionary$Description ) |   # linked to SAM, then AHL
                                     grepl(" acyl carrier prot", dictionary$Description , fixed=T ) |   # linked to SAM, then AHL
                                     grepl("homoserine-lacton", dictionary$Description , fixed=T ) |
                                     grepl("homoserine lacton", dictionary$Description ) |
                                     grepl("[A-a]utoinducer bind", dictionary$Description ) |
                                     grepl("[A-a]utoinducer-bind", dictionary$Description ) |
                                     grepl("[A-a]utoinducer 2", dictionary$Description ) |
                                     grepl("[A-a]utoinducer-2", dictionary$Description ) |
                                     grepl("CAI-1", dictionary$Description , fixed=T ) |   # also AI-2E are correct patterns!
                                     grepl("^AI-2", dictionary$Description ) |   # also AI-2E are correct patterns!
                                     grepl(" AI-2", dictionary$Description , fixed=T ) |   # also AI-2E are correct patterns!
                                     grepl("autoinducer 3", dictionary$Description, fixed=T ) |
                                     grepl("^AI-3", dictionary$Description) |
                                     grepl(" AI-3", dictionary$Description, fixed=T ) |
                                     grepl("AiiA", dictionary$Description, fixed=T ) |    # autoinducer lactonase
                                     grepl("AIP ", dictionary$Description, fixed=T ) |    # NB: non AIPR, which is related to fage resistance!!!
                                     grepl("hydroxyketon", dictionary$Description, fixed=T ) |
                                     grepl("signal factor", dictionary$Description, fixed=T ) |
                                     grepl("DSF", dictionary$Description, fixed=T ) |
                                     grepl("11-methyl-2-dodecenoic", dictionary$Description, fixed=T ) |   # a DSF
                                     grepl("BDSF", dictionary$Description, fixed=T ) |  
                                     grepl("cis-2-dodecenoic", dictionary$Description, fixed=T ) |   # alias BDSF
                                     grepl("BDE", dictionary$Description, fixed=T ) |
                                     grepl("S-adenosylmethionine", dictionary$Description, fixed=T ) |
                                     grepl("lactone acylase", dictionary$Description, fixed=T ) |
                                     grepl("[A-a]gracinopine", dictionary$Description, fixed=T ) |
                                     grepl("OH-PAME", dictionary$Description, fixed=T ) |
                                     grepl("OH-MAME", dictionary$Description, fixed=T ) |
                                     grepl("tetrahydroxytetrahydrofuran borate", dictionary$Description, fixed=T ) |    # (AI-2) , doi: 10.1128/mBio.02086-17
                                     grepl("hydroxytridecan-4-one", dictionary$Description, fixed=T ) |    # (AI-2) , doi: 10.1128/mBio.02086-17
                                     grepl("hydroxyketone", dictionary$Description, fixed=T ) |    
                                     grepl("hydroxy-ketone", dictionary$Description, fixed=T ) |    
                                     grepl("dimethylpyrazin-2-ol", dictionary$Description, fixed=T ) |    # (AI-2) , doi: 10.1128/mBio.02086-17
                                     grepl("1-acyl-sn-glycerol-3-phosphate acyltransf", dictionary$Description, fixed=T ) |    # AHLs , doi.org/10.1016/j.watres.2023.119814
                                     grepl("1-acyl-sn-glycerol-3-P acyltr", dictionary$Description, fixed=T ) |  # doi.org/10.1016/j.watres.2023.119814
                                     grepl("4,5-dihydroxy-2,3-pentanedione", dictionary$Description, fixed=T ) |    # see doi.org/10.1016/j.biotechadv.2025.108733
                                     grepl("3-oxo-octanoyl", dictionary$Description, fixed=T ) |  # see doi.org/10.1016/j.biotechadv.2025.108733
                                     grepl("3-hydroxy-4-quinolinone", dictionary$Description, fixed=T ) |  # see doi.org/10.1016/j.biotechadv.2025.108733
                                     grepl("3-hydroxyacyl-ACP", dictionary$Description, fixed=T ) |  # linked to QS through DSFs, see doi.org/10.1016/j.biotechadv.2025.108733
                                     grepl("2-alkylquinolines", dictionary$Description, fixed=T ) |  
                                     grepl("N-acyl homoserine", dictionary$Description, fixed=T ) |  
                                     grepl("N-acyl homoserine", dictionary$Description, fixed=T ) |  
                                     grepl("N-acyl-homoserine", dictionary$Description, fixed=T ) |  
                                     grepl("S-adenosylhomocysteine ", dictionary$Description , fixed=T ) |  # it is LuxS                                
                                     grepl("S-ribosylhomocysteine ", dictionary$Description , fixed=T ) |  # it is LuxS                                
                                     grepl("homoserine deaminase", dictionary$Description ) |                                  
                                     grepl("homoserine decarb", dictionary$Description ) |                                  
                                     grepl("homoserine lactonas", dictionary$Description ) |
                                     grepl("homoserine-lactonas", dictionary$Description , fixed=T) |                                  
                                     grepl("penicillin acylase", dictionary$Description ) |                                  
                                     grepl("dCACHE", dictionary$Description ) |                                  
                                     # grepl("CDP-diacylglycerol--glycerol", dictionary$Description ) |     # this was a gene downregulated by Lux in a paper, yet it is not linked with the quorum sensing, not directly at least                             
                                     grepl("Rpf", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("Lux ", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("Lux ", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("Lux$", dictionary$Description ) |    # NB: yes, in description
                                     grepl("HdtS", dictionary$Description ) |    # NB: yes, in description
                                     grepl("AttM", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("PQS", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("pqsR", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("CQS", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("Cqs-like", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     # grepl("Clp", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("Lsr", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     # grepl("Phn", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("SdiA", dictionary$Description, fixed=T ) |    # NB: yes, in description
                                     grepl("Blc", dictionary$Description, fixed=T ) |    # NB: yes, in description (outer membrane)
                                     # NB: do not search for "DGC", it would lead to other proteins not involved in QS
                                     
                                     grepl("^ahl", dictionary$Gene ) |  
                                     grepl("^aiiA", dictionary$Gene ) |  # yes, only A  
                                     grepl("^att", dictionary$Gene ) |   
                                     grepl("^blcC", dictionary$Gene ) |   # ref only about C (an outer membrane peptidase)
                                     grepl("^blp", dictionary$Gene ) |
                                     grepl("^cciR", dictionary$Gene ) |  # ref only for R
                                     # grepl("^clp", dictionary$Gene ) | 
                                     grepl("^cqs", dictionary$Gene ) | 
                                     grepl("^cepI", dictionary$Gene ) |  # ref only for I and L
                                     grepl("^cepL", dictionary$Gene ) |  # ref only for I and L
                                     grepl("^crp", dictionary$Gene ) | 
                                     grepl("^dgc", dictionary$Gene ) |   # NB: do not search for "DGC", it would lead to other proteins not involved in QS
                                     grepl("^esa", dictionary$Gene ) | 
                                     grepl("^eam", dictionary$Gene ) |   
                                     grepl("^hapR", dictionary$Gene ) |     # only R
                                     grepl("^hdt", dictionary$Gene ) |     # only R
                                     grepl("^lsr", dictionary$Gene ) |    
                                     grepl("^lux", dictionary$Gene ) |    
                                     grepl("nprX", dictionary$Gene ) |    
                                     # grepl("^phn", dictionary$Gene ) |   
                                     grepl("^phz", dictionary$Gene ) |   
                                     grepl("^papR", dictionary$Gene ) |     # ref only for R
                                     grepl("^psa", dictionary$Gene ) |   
                                     grepl("^pvd", dictionary$Gene ) |   
                                     grepl("^pqs", dictionary$Gene ) |   
                                     grepl("^tof", dictionary$Gene ) |   
                                     grepl("^tra", dictionary$Gene ) |   
                                     grepl("^trb", dictionary$Gene ) |   
                                     grepl("^vfm", dictionary$Gene ) |   
                                     grepl("^rcsA", dictionary$Gene ) |    # ref only for A
                                     grepl("^rhi", dictionary$Gene ) |   
                                     grepl("^rpf", dictionary$Gene ) |   
                                     grepl("^sdi", dictionary$Gene ) |   
                                     grepl("^xag", dictionary$Gene ) |
                                     grepl("^qseC", dictionary$Gene ) ,  
]

### Searching for other quorum genes using UniprotKB ...
## downloading results as table from https://www.uniprot.org/uniprotkb?query=quorum+sensing
# system("wget -r https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28quorum+sensing%29  -O  quorumS_results_UniProtKB.gz")
# system(gzip -d "quorumS_results_UniProtKB.gz")
# file.rename("quorumS_results_UniProtKB","quorumS_results_UniProtKB.tsv")
UniProt_KB_IDs<-read.delim("quorumS_results_UniProtKB.tsv", header = T, sep = "\t")
UniProt_KB_IDs$Gene.Names <- gsub(" .*", "", UniProt_KB_IDs$Gene.Names)
UniProt_KB_IDs<- UniProt_KB_IDs[ !duplicated(UniProt_KB_IDs$Gene.Names) , ]
# UniProt_KB_IDs[grepl("quorum", UniProt_KB_IDs)]  # just to check

### adding the "others" quorum related genes to the list
UniProt_Plus <- UniProt_KB_IDs [ UniProt_KB_IDs$Entry %in% dictionary$Entry , "Entry" ]
UniProt_Plus <- UniProt_Plus [! UniProt_Plus %in% only_QuorumGenes_IDs$Entry ]
length(UniProt_Plus)   # other 26 entries retrieved using UniProt ...
rm(UniProt_KB_IDs)


# removing certain genes ...
only_QuorumGenes_IDs <- only_QuorumGenes_IDs[only_QuorumGenes_IDs$Gene!="sppA" , ]  # it is a "signal peptidase" but this signal is about the proteins translocation!
only_QuorumGenes_IDs <- only_QuorumGenes_IDs[only_QuorumGenes_IDs$Gene!="queA" , ]  # it is related to tRNAs, not QS ...
only_QuorumGenes_IDs <- only_QuorumGenes_IDs[! grepl("QueA", only_QuorumGenes_IDs$Description) , ]  # it is related to tRNAs, not QS ...
only_QuorumGenes_IDs <- only_QuorumGenes_IDs[! grepl("di-iron", only_QuorumGenes_IDs$Description) , ]  
only_QuorumGenes_IDs <- only_QuorumGenes_IDs[! grepl("methylthioadenosine", only_QuorumGenes_IDs$Description, fixed = T) , ]   # Involved in SAM/SAH cycle BUT too much indirect for the quorum ...
only_QuorumGenes_IDs <- only_QuorumGenes_IDs[!grepl("Adineta|Rotaria|Diploscapter|Carchesium|Stentor",only_QuorumGenes_IDs$Taxon) , ] 


# REF for Clp Proteins regulating the quorum...
# doi.org/10.1371/journal.pbio.1002449




################ EXTRA: QS GENES *** GLOMMING THROUGH DESCRIPTION + GENUS *** ###########################

table_descr<-as.data.frame( table_abs[ c( only_QuorumGenes_IDs$Entry , UniProt_Plus ) ,  ] ) 
selected_infos<-dictionary[dictionary$Entry %in% row.names(table_descr), ]
row.names(selected_infos) <- selected_infos$Entry
selected_infos<-selected_infos[ row.names(table_descr), ]

# Only few (quick) minor modification to the descriptions (this is "just" an extra check)
{selected_infos$Description<-gsub("S-adenosylmethionine synthetase.*","S-adenosylmethionine synth",selected_infos$Description)
  selected_infos$Description<-gsub("LuxR fam transcript regulator.*","LuxR fam transcript regulator",selected_infos$Description)
  selected_infos$Description<-gsub("coenzyme A","CoA",selected_infos$Description,fixed=T)
  selected_infos$Description<-gsub("-coA ","-CoA",selected_infos$Description,fixed=T)
  selected_infos$Description<-gsub("Autoinducer 2","AI-2",selected_infos$Description,fixed=T)
  selected_infos$Description<-gsub("Autoinducer-2","AI-2",selected_infos$Description,fixed=T)
  selected_infos$Description<-gsub("autoinducer-2 (AI-2)","AI-2",selected_infos$Description,fixed=T)
  selected_infos$Description<-gsub("AI-2E fam transp.*","AI-2E fam transp",selected_infos$Description)
  selected_infos$Description<-gsub("3-hydroxy-5-phosphonooxypentane-2,4-dione thiolase LsrF","AI-2 thiolase LsrF",selected_infos$Description, fixed = T) # Thiolase/aldolase involved in AI2-E final modifications, 10.2210/pdb4p2v/pdb
  selected_infos$Description<-gsub("AI-2 aldolase","AI-2 thiolase LsrF",selected_infos$Description, fixed = T) # Thiolase/aldolase involved in AI2-E final modifications, 10.2210/pdb4p2v/pdb
  selected_infos$Description<-gsub("Acyl-homoserine lactone acylase QuiP","Acyl-homoserine lactone acylase",selected_infos$Description,fixed=T)
  selected_infos$Description<-gsub("Acyl carrier.*","Acyl carrier",selected_infos$Description)
  selected_infos$Description[grepl("RpfG", selected_infos$Description)] <- "Cyclic di-GMP phosphodiesterase response regulator RpfG"
  selected_infos$Description[grepl("TrbG", selected_infos$Description)] <- "P-type conjugative transfer TrbG"
}
# Extra modifications always needed before glomming
{
  selected_infos$Description<- gsub(" prot$", "", selected_infos$Description)
  selected_infos$Description<- gsub(" prot ", " ", selected_infos$Description)
  selected_infos$Description<- gsub("subunit","",selected_infos$Description)
  selected_infos$Description<- gsub("^[P-p]robable ", "", selected_infos$Description)  
  selected_infos$Description<- gsub("[P-p]utative ", "", selected_infos$Description)  
  selected_infos$Description<- gsub(", putative", "", selected_infos$Description, fixed=T)  
  selected_infos$Description<- gsub("'", "", selected_infos$Description, fixed=T)  
  selected_infos$Description<- gsub(" $", "", selected_infos$Description)  
  selected_infos$Description<- Hmisc::capitalize(selected_infos$Description)
}
#unique(selected_infos$Description[grepl("x",selected_infos$Description)] )

selected_infos$Taxon[selected_infos$Entry=="A0A918V768"] <- "Algibacter mikhailovii"

# from strain to Genus
{ table_descr$Taxa <- selected_infos$Taxon 
  table_descr$Taxa <- gsub("'", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub("[", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub("]", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub(" bact ", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub(" bact$", "", table_descr$Taxa)
  table_descr$Taxa <- gsub("uncult ", "uncultured ", table_descr$Taxa)
  # to genus level ...
  table_descr$Taxa <- gsub("uncultured ", "", table_descr$Taxa)
  table_descr$Taxa <- gsub("unclassified ", "", table_descr$Taxa)
  table_descr$Taxa <- gsub("unidentified  ", "", table_descr$Taxa)
  table_descr$Taxa <- gsub("Candidatus ", "Candidatus_", table_descr$Taxa)
  table_descr$Taxa <- gsub("candidate ", "candidate_", table_descr$Taxa)
  table_descr$Taxa <- gsub("alpha ", "alpha_", table_descr$Taxa)
  table_descr$Taxa <- gsub("beta ", "beta_", table_descr$Taxa)
  table_descr$Taxa <- gsub("Alphaproteobacteria.*", "Alphaproteobacteria", table_descr$Taxa)
  table_descr$Taxa <- gsub("Betaproteobacteria.*", "Betaproteobacteria", table_descr$Taxa)
  table_descr$Taxa <- gsub("Bacteroidetes.*", "Bacteroidetes", table_descr$Taxa)
  table_descr$Taxa <- gsub("Paracoccaceae.*", "Paracoccus", table_descr$Taxa)
  table_descr$Taxa <- gsub( " .*", "" , table_descr$Taxa )
}

table_descr$descr <-paste0( selected_infos$Description," (", table_descr$Taxa,")" )
# Now "descr" is Description + Genus
table_descr$Taxa <- NULL
table_descr<-aggregate(.~descr, table_descr, FUN=sum)
# head( table_descr[order(table_descr$descr), ] )
row.names(table_descr)<-table_descr$descr
table_descr$descr<-NULL

# The unmapped are not included in the subsetted dictionary
prop_descr_withOUT_unmap <- apply(table_descr, MARGIN = 2, function(x) x/sum(x))  
prop_descr_withOUT_unmap <- as.data.frame(prop_descr_withOUT_unmap) * 100

table <- as.data.frame( prop_descr_withOUT_unmap ) 
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")
top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# Further modification to the descriptions here, but this time they are mostly for aesthetics
table_to_plot$Protein <- gsub("S-adenosylmethionine synth", "S-adenosylmethionine synthase", table_to_plot$Protein, fixed=T)
table_to_plot$Protein <- gsub("carrier (", "carrier protein (", table_to_plot$Protein, fixed=T)
table_to_plot$Protein <- gsub("fam transp (", "fam transporter (", table_to_plot$Protein, fixed=T)
table_to_plot$Protein <- gsub("Acyl-homoserine-lactone synth", "Acyl-homoserine-lactone synthase", table_to_plot$Protein, fixed=T)
table_to_plot$Protein <- gsub("1-acyl-sn-glycerol-3-P acyltransf", "1-acyl-sn-glycerol-3-P acyltransf [HdtS ?]", table_to_plot$Protein, fixed=T) # doi.org/10.1016/j.watres.2023.119814

colors_type <- c("deeppink1","deeppink3",
                 "wheat3",
                 "greenyellow","chartreuse1","darkgreen",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "cyan","lightskyblue","deepskyblue1","cadetblue2","royalblue3","deepskyblue3","skyblue2","skyblue1",
                 "tomato2","indianred1","red1","firebrick1","darkred","orangered1","tomato1","firebrick3",
                 "darkmagenta","magenta3",
                 "gray35","black",
                 "gold3","yellow","gold", # ,"yellow2","gold",
                 "pink2"
)
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =8.5) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values= colors_type ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
        axis.title.y=element_text(size=8.65),
        axis.ticks.y=element_blank(),
        strip.text.x=element_text(size=9,colour="black", 
                                  lineheight = 1,
                                  margin = margin(2,3,2,3, "pt")
        ),
        strip.background = element_rect(color = "black", linewidth = 0.45),
        # plot.title = element_text(size=6.5),
        legend.key.height = unit(0.05, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.45 ),
        legend.margin = margin(-16,0,2,-37),
        legend.position="bottom",
        plot.margin = margin(2,1,0.5,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       # title = "Most abundant mRNA of starvation related genes (transcript level)"
  )
ggsave("RNA_PHA_Results/GlommingThroughDescriptions_QsTAXON_Most_abundant.png", width= 5.45, height = 4.8, dpi=300)



 
################ EXTRA: CHECKING THE rRNA REGULATORY GENES ###########################
# 
# rRNA_amount_genes<-dictionary [ grepl("^rrn", dictionary$Gene ) |    # general ribosomial related genes, e.g. rrnB P1 is a ribosomal RNA promoter in E coli, see DOI: 10.1002/j.1460-2075.1990.tb07586.x 
#                                   # grepl("^rps", dictionary$Gene ) |  # ribosomal proteins (Small)
#                                   # grepl("^rpl", dictionary$Gene ) |  # ribosomal proteins (Large)
#                                   grepl("^fis", dictionary$Gene ) |  # Fis activates rrnB transcription, hence rRNA transcr DOI: 10.1002/j.1460-2075.1990.tb07586.x
#                                   grepl(" Fis", dictionary$Description ) |  # Fis activates rrnB transcription, hence rRNA transcr DOI: 10.1002/j.1460-2075.1990.tb07586.x
#                                   grepl("Fis$", dictionary$Description ) |  # Fis activates rrnB transcription, hence rRNA transcr DOI: 10.1002/j.1460-2075.1990.tb07586.x
#                                   grepl("[F-f]actor for inversion stimul", dictionary$Description ) |  
#                                   grepl("[F-f]actor inversion stimul", dictionary$Description ) |
#                                   grepl("HupA", dictionary$Description ) |
#                                   # *** level of Fis is highest during early-exponential growth phase but extremely low in stationary phase  https://doi.org/10.1016/j.ygeno.2019.07.013
#                                   # grepl("^bipA", dictionary$Gene ) |  # "BipA is a ribosome binding GTP required for Fis synthesis" (...or, at least, it is strongly correlated!) , see https://doi.org/10.1038/sj.emboj.7600343
#                                   # grepl("BipA", dictionary$Description ) |  # "BipA is a ribosome binding GTP required for Fis synthesis" (...or, at least, it is strongly correlated!) , see https://doi.org/10.1038/sj.emboj.7600343
#                                   # *** However, BipA is involved in too many other functions --> discarding it from the list
#                                   grepl("^hns", dictionary$Gene ) |  # gene of H-NS, against Fis activity, see https://www.uniprot.org/uniprotkb/P0ACF8/entry
#                                   grepl("H-NS", dictionary$Description , fixed=T) | # see above
#                                   grepl("histone-like nucleoid-structuring prot", dictionary$Description , fixed=T) | # see above
#                                   grepl("ppGpp", dictionary$Description ) | # it should be "(p)ppGpp" but the first p is omitted to avoid pattern misinterpretations 
#                                   grepl("guanosine 5'-diphosphate 3'-diphosph", dictionary$Description ) | # it should be "(p)ppGpp" but the first p is omitted to avoid pattern misinterpretations 
#                                   #  grepl("lrp", dictionary$Gene) | # collaborates with H-NS, see DOI: 10.1111/j.1365-2958.2005.04873.x
#                                   grepl("Lrp ", dictionary$Description , fixed=T) | # collaborates with H-NS, see DOI: 10.1111/j.1365-2958.2005.04873.x
#                                   grepl("Lrp$", dictionary$Description) | # collaborates with H-NS, see DOI: 10.1111/j.1365-2958.2005.04873.x
#                                   grepl("^relX", dictionary$Gene ) | # gene involved in ppGpp amount regulation , see https://pubmed.ncbi.nlm.nih.gov/342913/
#                                   grepl("^relA", dictionary$Gene ) | # gene involved in ppGpp amount regulation , see https://pubmed.ncbi.nlm.nih.gov/342913/
#                                   grepl("^relP", dictionary$Gene ) | # domain homology with (p)ppGpp synthetase, see https://doi.org/10.1128/iai.01439-09
#                                   grepl("^relQ", dictionary$Gene ) | # domain homology with (p)ppGpp synthetase, see https://doi.org/10.1128/iai.01439-09
#                                   grepl("RelX", dictionary$Description ) | 
#                                   grepl("RelA", dictionary$Description ) |
#                                   grepl("GTP pyrophosphokinase", dictionary$Description ) |   # relA coded protein, see the gene ref above
#                                   grepl("RelP", dictionary$Description ) | 
#                                   grepl("RelQ", dictionary$Description ) | 
#                                   grepl("^rsh", dictionary$Gene ) | # (p)ppGpp synthase, homolog of relA, see https://doi.org/10.1128/iai.01439-09
#                                   grepl(" RSH", dictionary$Description ) |
#                                   grepl("RSH$", dictionary$Description ) |
#                                   # *** RSH is responsible for (p)ppGpp synthesis in response to amino acid deprivation (lack of Leu/Val) or mupirocin treatment (see the ref above)
#                                   grepl("^spoT", dictionary$Gene ) | # (p)ppGpp synthase, homolog of relA and rsh, see https://doi.org/10.1128/iai.01439-09
#                                   grepl("SpoT", dictionary$Description, fixed=T ) | # (p)ppGpp synthase, homolog of relA and rsh, see https://doi.org/10.1128/iai.01439-09
#                                   grepl("Guanosine-3',5'-bis(DiP) 3'-pyrophosphohydrolas", dictionary$Gene , fixed=T ) | # description of spoT
#                                   # *** From "PMC7356644" : In E. coli, synthesis of (p)ppGpp is carried out by RelA, during amino acid starvation, 
#                                   # *** [...] and the bifunctional SpoT protein during other stress treatments [...]  SpoT is also responsible for degradation of (p)ppGpp .
#                                   grepl("DksA", dictionary$Gene ) |   # transcr fact, binds to RNApol augmentig ppGpp induced transcr initiation, see  www.ncbi.nlm.nih.gov/gene/944850 and https://journals.asm.org/doi/10.1128/iai.01439-09
#                                   grepl("DksA", dictionary$Description ) ,    # transcr fact, binds to RNApol augmentig ppGpp induced transcr initiation, see  www.ncbi.nlm.nih.gov/gene/944850 and https://journals.asm.org/doi/10.1128/iai.01439-09
#                                   ### *** The following are involved in biogenesis, not transcription regulation
#                                   # grepl("Ribosome biogenesis ", dictionary$Description ) |
#                                   # grepl("Ribosome degrad ", dictionary$Description ) |
#                                   # grepl("^ylqF", dictionary$Gene ) | # involved in 50S subunit biogenesis, see DOI: 10.1111/j.1365-2958.2005.04948.x
#                                   # grepl("^rbgA", dictionary$Gene ) , # same as ylqF, see DOI: 10.1111/j.1365-2958.2005.04948.x
#                                 # grepl("^rpoS", dictionary$Gene ) |  # it is a sigma factor related to stress and starvation, increased with high ppGpp and Hfq , see PMC7356644
#                                 # grepl("RpoS", dictionary$Description ) ,  # it is a sigma factor related to stress and starvation, increased with high ppGpp and Hfq , see PMC7356644
#                                 
# ]
# # The RNases may act also in rRNA maturation and quality-control related degradations, hence they have NOT been included in this list
# 
# 
# rRNA_amount_genes<-rRNA_amount_genes [ !grepl("DNA-bind Fis domain|Fis-type|54|Fis-like|Fis [F-f]am|NifA", rRNA_amount_genes$Description ) , ]  # Fis activates rrnB transcription, hence rRNA transcr DOI: 10.1002/j.1460-2075.1990.tb07586.x
# # *** There are many Fis-like proteins which have the same structure but are not actually Fis... For example, those with sigma 54 should not be Fis!
# # according to   DOI: 10.1046/j.1365-2958.1996.388920.x  sigma 54 is involved in metabolism related activation, e.g. nitrogen assimilation and fixation,
# # Infact, Fis is more related to sigmaS (38) and sigma70, see https://doi.org/10.1111/j.1365-2958.2006.05560.x
# 
# rRNA_amount_genes <- rRNA_amount_genes[!grepl("Adineta|Rotaria|Diploscapter|Carchesium|Stentor",rRNA_amount_genes$Taxon) , ]
# 
# 
# 
# ### TRANSCRIPT LEVEL
# table <- as.data.frame( prop_with_unmap[ rRNA_amount_genes$Entry ,  ] ) 
# colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")
# 
# backup_numbers <- table  # to count them afterwards
# 
# top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
# table_top<- table[top, ]
# table_top$Protein<-row.names(table_top)
# ### Converting IDs
# selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
# {infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
#   infos<- gsub("()", "(Unknown source)", infos, fixed=T)
#   infos<- gsub("unclassified", "unclass", infos, fixed=T)
#   infos<- gsub("Betaproteobacteria bact HGW-Betaproteobacteria-21", "HGW-Betaproteobacteria-21", infos, fixed=T)
#   infos<- gsub("Bifunctional (P)ppGpp synthetase/guanosine-3',5'-bis(DiP) 3'-pyrophosphohydrolase", "Bifunctional enzyme (p)ppGpp related SpoT", infos, fixed=T)
#   infos<- gsub("transcript", "transcr", infos)
#   infos<- gsub("prot H-NS", "H-NS", infos, fixed=T)
#   infos<- gsub(", Lrp", " Lrp", infos, fixed=T)
#   infos<- gsub("Ribosome biogenesis GTPase A", "Ribosome biogenesis GTPase A 'relA'", infos, fixed=T)
#   infos<- gsub("Small ribosomal subunit prot", "RPS (S30) protein", infos, fixed=T)
#   infos<- gsub("Large ribosomal subunit prot", "RPL (S50) protein", infos, fixed=T)
# }
# table_top$Protein<-infos
# only_them<-table_top
# only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
#                                                        MARGIN = 2, function(x) x/sum(x) )
# only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# # plotting ...
# table_to_plot<-only_them
# table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
# table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
# table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
#   geom_bar(stat="identity", position="stack", linewidth=0.5) +
#   theme_classic(base_size =9) +
#   facet_grid2( . ~ variable+Reactor,
#                scales = "free",
#                space="free",
#                switch = "x", # This places the grid at the bottom (as x axis)
#                strip = strip_nested(size="constant"))+
#   scale_fill_manual(values= fill_color_30_BIS ) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
#         axis.ticks.y=element_blank(),
#         strip.text = element_text(size=10),
#         plot.title = element_text(size=6.5),
#         legend.key.height = unit(0.11, "cm"),
#         legend.key.width = unit(0.23, "cm"),
#         legend.spacing.x = unit(0.1, "cm"),
#         legend.title = element_text ( size = 9 ),
#         legend.text = element_text ( size = 6.4 ),
#         legend.margin = margin(-14,0,2,-40),
#         legend.position="bottom",
#         plot.margin = margin(3,1,3,1)
#   ) +
#   scale_x_discrete (expand = c(0.01,0) ) +
#   scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
#   guides(fill=guide_legend(nrow=15)) +
#   labs(x="", y="Percentages in samples\n(excluding other RNAs)",
#        fill="",
#        title = "Most abundant mRNA of rRNA synthesis related genes (transcript level)"
#   )
# ggsave("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_rRNAsynthesisRelated_Most_abundant_transcr_AbsCounts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)
# 
# colSums(backup_numbers[ , !colnames(backup_numbers)%in%"Protein" ])  # proportions on whole community
# # 0.05737485  0.03417964  0.03709485  0.03751821  0.03820338  0.04128640   
# # ten folded if the rps and rpl protein genes are included!
# 
# 
# 
# ### GENE LEVEL
# these_genes <- dictionary[ dictionary$Entry %in% rRNA_amount_genes$Entry , "Gene" ]
# these_genes <- these_genes[these_genes!=""]  # Related *BUT* not linked to any gene --> Impossible to merge to other transcripts
# table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% these_genes , ] ) # used later when *only* PAO and GAO are required
# colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")
# 
# top <- names(sort(rowSums(table[row.names(table)!="" , ]), decreasing=TRUE))[1:30]   # "" = unmapped
# table_top<- table[top, ]
# table_top$Protein<-row.names(table_top)
# # Converting IDs
# selected_infos<-dictionary[dictionary$Gene %in% table_top$Protein, ]
# selected_infos<-selected_infos[! duplicated(selected_infos$Gene), ]
# row.names(selected_infos)<-selected_infos$Gene
# selected_infos<-selected_infos[table_top$Protein, ]
# selected_infos$Description[grepl("phage",selected_infos$Taxon)]<-paste(selected_infos$Description[grepl("phage",selected_infos$Taxon)], "(virus)")
# infos<-paste0( selected_infos$Description," (", selected_infos$Gene,")" )
# { infos<- gsub("()", "(no gene name)", infos, fixed=T)
#   infos<- gsub("Bifunctional (P)ppGpp synthetase/guanosine-3',5'-bis(DiP) 3'-pyrophosphohydrolase", "Bifunctional enzyme (p)ppGpp related SpoT", infos, fixed=T)
#   infos<- gsub("transcript", "transcr", infos)
#   infos<- gsub("prot H-NS", "H-NS", infos, fixed=T)
#   infos<- gsub(", Lrp", " Lrp", infos, fixed=T)
#   infos<- gsub("(P)", "(p)", infos, fixed=T)
#   infos<- gsub("Small ribosomal subunit prot", "RPS (S30) protein", infos, fixed=T)
#   infos<- gsub("Large ribosomal subunit prot", "RPL (S50) protein", infos, fixed=T)
#   infos<- gsub("Putative Fis-like DNA-bind prot", "Fis-like DNA-bind protein (putative)", infos, fixed=T)
#   # infos<- gsub("Sigma-54-dependent Fis fam transcr regulator", "σ-54-dependent Fis family prot", infos, fixed=T)
#   # infos<- gsub("Two component sigma-54 specific Fis fam transcr regulator", "Two component σ-54 specific Fis family prot", infos, fixed=T)
# }
# table_top$Protein<-infos
# only_them<-table_top
# only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
#                                                        MARGIN = 2, function(x) x/sum(x) )
# only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# # plotting ...
# table_to_plot<-only_them
# table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
# table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
# table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# base_plot <- ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
#   geom_bar(stat="identity", position="stack", linewidth=0.5) +
#   theme_classic(base_size =9) +
#   facet_grid2( . ~ variable+Reactor,
#                scales = "free",
#                space="free",
#                switch = "x", # This places the grid at the bottom (as x axis)
#                strip = strip_nested(size="constant"))+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
#         axis.ticks.y=element_blank(),
#         strip.text = element_text(size=10),
#         plot.title = element_text(size=6.5),
#         legend.key.height = unit(0.11, "cm"),
#         legend.key.width = unit(0.32, "cm"),
#         legend.spacing.x = unit(0.3, "cm"),
#         legend.title = element_text ( size = 9 ),
#         legend.text = element_text ( size = 7.5 ),
#         legend.margin = margin(-14,0,2,-39),
#         legend.position="bottom",
#         plot.margin = margin(1,2,1,1)
#   ) +
#   scale_x_discrete (expand = c(0.01,0) ) +
#   scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
#   guides(fill=guide_legend(nrow=15)) +
#   labs(x="", y="Percentages in samples\n(excluding other genes)",
#        fill="",
#        title = "Most abundant genes of rRNA synthesis related genes (gene level)"
#   )
# # COLORS FOR EACH TYPE OF GENE ...
# colors_type <- c("orangered","firebrick2","maroon",
#                  "red","coral","firebrick4","tomato2","coral3",
#                  "red2",
#                  "magenta",
#                  "black","gray30",
#                  # "green","chartreuse3","yellowgreen","darkgreen","greenyellow","forestgreen","green3","chartreuse","olivedrab","springgreen3",
#                  "cyan2","deepskyblue2","blue1","cyan4","skyblue","deepskyblue","royalblue3","dodgerblue","lightskyblue","cyan3","royalblue1","blue3","dodgerblue3",
#                  
#                  "darkmagenta","violet","magenta3",
#                  # "gray35","gray18","slategray4",
#                  # "royalblue1","deepskyblue","blue","cyan4","royalblue2","skyblue","deepskyblue2","dodgerblue","lightskyblue","cyan3","darkblue","cyan","dodgerblue3"
#                  "chartreuse3","green"
# )
# base_plot + scale_fill_manual(values= colors_type )
# ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_GENES_rRNAsynthesisRelated_Most_abundant_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)
