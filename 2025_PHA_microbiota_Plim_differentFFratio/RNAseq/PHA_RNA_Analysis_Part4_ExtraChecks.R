# WORK IN PROGRESS!!!


################ PREPARING THE ENVIRONMENT ###############

{
  library("ggplot2")
  library("ggvenn")
  library("reshape2")
  library("ggh4x")
  library("ecodist")
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
taxonomy <- read.delim("databases_dictionary/Dictionary_ncRNA_microb_scfbio_24_01_24.csv", sep=",", header = T)
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

