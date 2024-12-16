####################### PREPARING THE ENVIRONMENT #####################

{
  library("ggplot2")
  library("ggvenn")
  library("reshape2")
  library("ggh4x")
  library("ecodist")
  suppressPackageStartupMessages( library("vegan") ) # it has a verbose loading!
  library("DESeq2")
}

options(scipen=100)

dir.create("RNA_AGS_Results")
dir.create("RNA_AGS_Results/Most_abundant_RNAs")
dir.create("RNA_AGS_Results/Most_abundant_GENES")
dir.create("RNA_AGS_Results/Focus_on_PAO_and_GAO")
dir.create("RNA_AGS_Results/Correlations")


table_abs<-read.table("Diamond_output/Abs_counts_from_Diamond.tsv", row.names = 1, header = T)
# table_TPM<-read.table("AGS_RNA_Diamond_output/TPM_counts_from_Diamond.tsv", row.names=1, header = T)
dictionary<-read.delim("Diamond_output/IDs_infos.tsv", header = T, sep="\t")
metadata<-read.delim("AGS_RNA_metadata.tsv", header = T, sep="\t")


colnames(table_abs)<-gsub("^X","",colnames(table_abs))
# colnames(table_TPM)<-gsub("^X","",colnames(table_TPM))



#### Adding unmapped reads
mapping<-read.table("Diamond_output/Diamond_matching_percents.tsv", row.names = 1, header = T)
mapping<-mapping[colnames(table_abs),  ]
unmapped_count <- mapping$INPUT_SEQS_FROM_BOTH - mapping$CONCORD_MATCHES   # these are unmapped but still (theoretically) microbial
# NB: "inputs from both" does not mean summing up R1 and R2, moreover the "Concord_matches" includes also the unique R1 and R2 matches (therefore is the count of every successful not discarded mapping)
table_abs<-rbind(table_abs, UNMAPPED=unmapped_count)

# just a further sanity check ...
# not_micr<-as.data.frame(read.delim("Diamond_output/abs_NOT_microb_counts_from_Diamond.tsv", row.names = 1, header = T, sep="\t"))
# colnames(not_micr)<-gsub("^X","",colnames(not_micr))
# not_micr_maps<-colSums(not_micr[ ,  colnames(table_abs)] ) #  this step both orders and subsets
# mapping [1, "CONCORD_MATCHES"] - not_micr_maps[[1]] == colSums(table_abs[row.names(table_abs)!="UNMAPPED",])[[1]]    # the "not microbial" were removed from the total successful mappings in a last "manual" filter during the alignments parsing



#### Custom names
table_abs <- table_abs[ , metadata$FASTQ_ID] # this re-orders
colnames(table_abs)<-metadata$Sample_name
# table_TPM <- table_TPM[ , metadata$FASTQ_ID]
# colnames(table_TPM)<-metadata$Sample_name



fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
fill_color_20<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","firebrick3","springgreen4","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","darkgreen","lightgrey","coral","deepskyblue2") 
fill_color_30<-c("wheat3","deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "aquamarine2",
                 "bisque2","cyan","yellow3","brown",
                 "springgreen4", "firebrick3",
                 "grey15","lightblue1","lightgreen","orange2",
                 "darkorchid", "darkslategray3","blue1", "violet", 
                 "yellow", "pink3","yellow4", "chocolate4", "red","darkgreen",
                 "lightgrey","coral2","deepskyblue2")


# this overlaps with ssrA tm-RNA
# sequence from UniProt blasted with tblastn --> 67% of identity, 78% cover, 3 e-37 Evalue, overlaps on exact same area (seen on Genome Viewer)
dictionary[dictionary$Entry=="A0A011MM28", "Description" ] <- "Uncharact (~ssrA)"





##################### FILTERING ACCORDING TO THE ABS COUNT TABLE #####################

# this is an important step according to PMID 27441086 and PMID 27508061

if ( ! "table_abs_CPM" %in% ls() ) {
  table_abs_original <- table_abs
}
table_abs_CPM <- table_abs
table_abs_CPM <- apply( table_abs_CPM, MARGIN = 2, function(x) as.numeric(x) )
table_abs_CPM <- apply( table_abs_CPM, MARGIN = 2, function(x) x/sum(x) )
table_abs_CPM <- table_abs_CPM * 1000000 # NB: 1 million
row.names(table_abs_CPM) <- row.names(table_abs_original) # used afterward... 
CPM_min<-2
how_many_min<-1
# keep  <- apply(table_abs_CPM> 0.5) >= 2   # reference : PMID 27508061... but here there is a different sample (the other AGS) and also many of its abundant unique rows will be discarded if >2 (at least two samples) is maintenined!
keep  <- rowSums(table_abs_CPM> CPM_min) >= how_many_min   # ... therefore, just this threshold has to be reached in at least X number of sample (= minimum value)
# discarded <- table_abs[! keep, ]
# discarded
table_abs<-table_abs[keep, ]
remaining_message<- paste( length(which(keep)), "maintained on", length(table_abs_original[,1]) ,"filtered reads -->", round( (length(which(keep))/length(table_abs_CPM[,1]))*100 , 3) ,"% remaining")
remaining_message

# just a check
colSums(table_abs)/1000000
tail(table_abs_original)
tail(table_abs)

con <- file("RNA_AGS_Results/Further_discarderd_IDs_according_to_low_abund.txt")
sink(con, append = T)
cat("Every RNA that did not reach a CPM count of", CPM_min, "in at least", how_many_min, "samples have been discarded", fill=T)
cat("-->", remaining_message, fill=T)
cat("\nNB: the high threshold (CPM 2) has been chosen to slightly reduce the depth difference among samples, otherwise the sample LR39 (the one with the higher remaining depth) would have TOO MANY unique reads with low abundance", fill=T)
cat("\n\n\n")
cat("Total remaining_reads ABS divided 1million (= ANALYSED DEPTH DIFFERENCE)", fill=T)
print( colSums(table_abs)/1000000  ) 
sink()
close(con)


# table_abs<-table_abs_original # to reset the filter


### also the dictionary has to be re-filtered (to avoid bugs and NA during the re-subsettings)
dictionary <- dictionary[ dictionary$Entry %in% row.names(table_abs) , ]
length(unique(dictionary$Entry))
length(unique(dictionary$Taxon))

proof_filter<-"This is the proof of having performed the filtering steps"

# NB: the table "table_abs_CPM" still includes also the NOT filtered reads and it's already in CPM ... it will be used for further filters!


rm(remaining_message, unmapped_count, keep)

# rm(table_abs_CPM, table_abs_original, mapping)
# save.image("Backup_tabelle_RNA_AGS_after_abund_filters.RData")
# load("Backup_tabelle_RNA_AGS_after_abund_filters.RData")




###################### PREPARING THE GENE LEVEL OBJECT ##########################
# this is a slow step, hence it is performed only once, here... and after the filtering!

# glomming to gene level
table_gene<-cbind.data.frame(table_abs[dictionary$Entry, ], Gene=dictionary$Gene)
table_gene<-aggregate(.~Gene, table_gene, FUN=sum)
row.names(table_gene)<-table_gene$Gene
table_gene$Gene<-NULL

# length dim(table_gene)



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
prop_GENE_with_unmap <- prop_GENE_with_unmap * 100




############## SUBSETTING AND CHECKING PAOs AND GAOs SOURCED RNA #################

# common PAO GAO list (PMID:38524765 ,  PMID: 22827168,  PMID: 15774634,  PMID: 33404187, MIDAS website
only_PAO_IDs<-dictionary [ grepl("Accumulibacter", dictionary$Taxon ) | 
                             grepl("Phosphoribacter", dictionary$Taxon ) |
                             grepl("Tetrasphaera", dictionary$Taxon ) |
                             grepl("Accumulimonas", dictionary$Taxon ) |
                             grepl("Microthrix", dictionary$Taxon ) |
                             grepl("Dechloromonas", dictionary$Taxon ) |
                             grepl("Malikia", dictionary$Taxon ) |
                             grepl("Quatrionicoccus", dictionary$Taxon ) |
                             grepl("Beggiatoa", dictionary$Taxon ) |
                             grepl("Gemmatimonas", dictionary$Taxon ) |
                             grepl("Friedmaniella", dictionary$Taxon ) |
                             grepl("Tessaracoccus", dictionary$Taxon ) |
                             grepl("Azonexus", dictionary$Taxon ), 
]
only_GAO_IDs<-dictionary [ grepl("Contendobacter", dictionary$Taxon ) |
                             grepl("Competibacter", dictionary$Taxon ) |
                             grepl("Defluviicoccus", dictionary$Taxon ) |
                             grepl("Proximibacter", dictionary$Taxon ) |
                             grepl("Micropruina", dictionary$Taxon ) |
                             grepl("Propionivibrio", dictionary$Taxon ) ,
]
# unique(only_PAO_GAO_IDs$Taxon)



# Just a quick sanity check...
if ( length(which( grepl("virus", only_PAO_IDs$Taxon) | grepl("phage", only_PAO_IDs$Taxon) )) == 0 ){
  cat ("Ok, no RNA related to phage or virus among the selected GAO entries...", fill=T)
}
if ( length(which( grepl("virus", only_PAO_IDs$Taxon) | grepl("phage", only_PAO_IDs$Taxon) )) == 0 ){
  cat ("Ok, no RNA related to phage or virus among the selected PAO entries...", fill=T)
}


### Checking and plotting the ratio of functions with not PAO or GAO
suppressWarnings(rm(table))
table<-table_abs[!row.names(table_abs) %in% "UNMAPPED", ]
table$PAO_GAO<-"not PAO or GAO"
table[only_PAO_IDs$Entry, "PAO_GAO"] <-"other PAOs"
table[ only_PAO_IDs[ grepl("Accumulibacter", only_PAO_IDs$Taxon) , "Entry" ] , "PAO_GAO"] <- "Ca Accumulibacter"
table[only_GAO_IDs$Entry, "PAO_GAO"] <-"other GAOs"
table[ only_GAO_IDs[ grepl("Competibacter", only_GAO_IDs$Taxon) , "Entry" ] , "PAO_GAO"] <- "Ca Competibacter"
table_to_plot<-table
# to plot relative percentages
table_to_plot[ , ! colnames(table_to_plot) %in% c("Protein", "PAO_GAO") ] <- apply( table_to_plot[ ,! colnames(table_to_plot) %in% c("Protein", "PAO_GAO") ],
                                                                                    MARGIN = 2, function(x) as.numeric(x) )
table_to_plot[ , ! colnames(table_to_plot) %in% c("Protein", "PAO_GAO") ] <- apply( table_to_plot[ , ! colnames(table_to_plot) %in% c("Protein", "PAO_GAO") ],
                                                                                    MARGIN = 2, function(x) x/sum(x) )
table_to_plot[ , ! colnames(table_to_plot) %in% c("Protein", "PAO_GAO") ] <- table_to_plot[ , ! colnames(table_to_plot) %in% c("Protein", "PAO_GAO") ] * 100
# plotting ...
colnames(table_to_plot)<- c( paste0(metadata$Experiment_day, "th day") ,"PAO_GAO" ) # already set the same order beforehand ,
colnames(table_to_plot)<-gsub("AGSth day", "AGS", colnames(table_to_plot)) 
table_to_plot<-melt(table_to_plot, id.vars = "PAO_GAO")
table_to_plot<-aggregate( .~PAO_GAO+variable, data=table_to_plot, FUN=sum)
table_to_plot$PAO_GAO<-factor(table_to_plot$PAO_GAO, levels= c("Ca Accumulibacter","Ca Competibacter", "other PAOs", "other GAOs", "not PAO or GAO") )
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=PAO_GAO)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values= c("Ca Accumulibacter"="coral",
                              "Ca Competibacter"="chartreuse2",
                              "other PAOs"="red3",
                              "other GAOs"="chartreuse4",
                              "not PAO or GAO"="deepskyblue2"
                              )
                    )+
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7.5),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 10 ),
        legend.margin = margin(-10,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment day", y="absolute count % (only PAO and GAO RNAs)", 
       fill="",
       title = "Ratios of PAO/GAO related RNAs in the dataset")
ggsave("RNA_AGS_Results/Focus_on_PAO_and_GAO/PAO_GAO_RNA_ratios.png", width= 5.5, height = 4.8, dpi=300)



###### EXPORTING INFOs ABOUT EVERY PAO GAO ENTRIES FOUND

ordered_infos<-dictionary
table2<-table[table$PAO_GAO != "not PAO or GAO", ]
table2[grepl("Accumulibacter",table2$PAO_GAO), "PAO_GAO"]<-"PAO"
table2[table2$PAO_GAO=="other PAOs", "PAO_GAO"]<-"PAO (not Ca Accumuli)"
table2[grepl("Competibacter",table2$PAO_GAO), "PAO_GAO"]<-"GAO"
table2[table2$PAO_GAO=="other GAOs", "PAO_GAO"]<-"GAO (not Ca Competib)"
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table2), ]
# ordered_infos$Pathway<-NULL
write.csv2(file=paste0("RNA_AGS_Results/Focus_on_PAO_and_GAO/Every_PAO_GAO_inferred_RNAs_Info.csv"), cbind(table2, ordered_infos), row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_to_plot))




########### BAR PLOTS OF THE MOST ABUNDANT RNA (ABS_count) ###############

# table<-table_abs
table <- as.data.frame(prop_with_unmap) # prop abs
#colnames(table)<- paste0(metadata$Experiment_day, "th day")
#colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 
colnames(table)<- metadata$Experiment_day

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("Candidatus","Ca.",infos,fixed=T)
infos<-gsub("bacterium","",infos,fixed=T)
infos<-gsub(" )",") ",infos,fixed=T) # the space here allows to distinguish two transcripts with the same name
infos<-gsub("denitrificans Run_A_D11","denit RUN_A_D11",infos,fixed=T)
infos<-gsub("Trichosanthes kirilowii picorna-like virus","Trichosanthes kirilowii",infos,fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.42, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-16,0,3,-36),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="Experiment day", y="absolute count %", 
       fill="",
       title = "Most abundant mRNA in AGS mixed liquor",
       caption = "NB: different depth among samples --> if more depth then more 'others' RNA!")
ggsave("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)

# write.csv2(file="RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_WHOLE_DATASET_infos.csv", table_top, row.names=T, quote=F)




########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                               MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.45, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-5,0,-2,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="Experiment day", y="absolute count % (excluding other RNAs)", 
       fill="",
       title = "Most abundant mRNA in AGS mixed liquor",
       caption = "NB: these are ONLY the most abundant ones!")
ggsave("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)
# write.csv2(file="RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_RNA_computed_with_abs__PROPORTIONS_ON_THEMSELVES_infos.csv", table_top, row.names=T, quote=F)


# AGAIN, BUT DESIGNED FOR A POSTER
table_to_plot2<-table_to_plot[ table_to_plot$variable!="other_AGS" , ]
table_to_plot2$Protein<-gsub(" Run_A_D11","", table_to_plot2$Protein)
table_to_plot2$Protein<-gsub(" Run_B_J11","", table_to_plot2$Protein)
table_to_plot2$Protein<-gsub("_[0-9]","", table_to_plot2$Protein)
ggplot(data=table_to_plot2, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=6.8),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.2), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 8.45 ),
        legend.margin = margin(-12,0,3,-45),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % \n(excluding other RNAs)", 
       fill="")
ggsave("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_RNA_poster_version.png", width= 6.7, height = 4.5, dpi=300)



########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( table_top, MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( table_top[colnames(table_top)!="other_AGS"], MARGIN = 1, mean)
to_save <- cbind.data.frame( table_top[ , !grepl("^LR", colnames(table_top) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_RNAs/Table_Most_Abundant_RNA_ABSpropAverages_WHOLE_DATASET_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them[ colnames(only_them)!="Protein"], MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[ !colnames(only_them) %in% c("other_AGS","Protein")], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS,  ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_RNAs/Table_Most_Abundant_RNA_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




# ########### BAR PLOTS OF THE MOST ABUNDANT RNA (TPM_count) ###############
# 
# table<-table_TPM
# # table<-table[ as.numeric(apply(table, 1, FUN = max)) >0, , drop=F] # more than 1 read to be sums to the 'others'
# 
# top <- names(sort(rowSums(table), decreasing=TRUE))[1:20]
# table_top<- table[top, ]
# table_top$Protein<-row.names(table_top)
# ### Converting IDs
# selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
# infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
# infos<- gsub("()", "(Unknown source)", infos, fixed=T)
# infos<-gsub("domain-containing protein","domain-containing prot",infos,fixed=T)
# infos<-gsub("polymerase","pol",infos,fixed=T)
# infos<-gsub("Candidatus","Ca.",infos,fixed=T)
# table_top$Protein<-infos
# table_to_plot<-table_top
# # to plot relative percentages
# table_to_plot[ , colnames(table_to_plot)!="Protein"] <- apply( table_to_plot[ , colnames(table_to_plot)!="Protein"],
#                                                                MARGIN = 2, function(x) as.numeric(x) )
# table_to_plot[ , colnames(table_to_plot)!="Protein"] <- apply( table_to_plot[ , colnames(table_to_plot)!="Protein"],
#                                                                MARGIN = 2, function(x) x/sum(x) )
# table_to_plot[ , colnames(table_to_plot)!="Protein"] <- table_to_plot[ ,colnames(table_to_plot)!="Protein"] * 100
# # plotting ...
# colnames(table_to_plot)<- c( paste0(metadata$Experiment_day, "th day") ,"Protein" ) # already set the same order beforehand , 
# colnames(table_to_plot)<-gsub("AGSth day", "AGS", colnames(table_to_plot))
# table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# # bar plot
# ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
#   geom_bar(stat="identity", position="stack", linewidth=0.5) + 
#   theme_classic(base_size =10) +
#   scale_fill_manual(values=fill_color_20) +
#   theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
#         axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
#         axis.ticks.y=element_blank(), 
#         legend.key.height = unit(0.2, "cm"),
#         legend.key.width = unit(0.5, "cm"),
#         legend.spacing.x = unit(0.3, "cm"),
#         legend.title = element_text ( size = 9 ),
#         legend.text = element_text ( size = 7.5 ),
#         legend.margin = margin(-10,0,3,-35),
#         legend.position="bottom",
#         plot.margin = margin(3,1,3,1)
#   ) +
#   guides(fill=guide_legend(nrow=10)) +
#   labs(x="", y="absolute count %", 
#        fill="",
#        title = "Most abundant mRNA in AGS mixed liquor (TPM counts)",
#        caption = "NB: those are ONLY the most abundant ones")
# ggsave("RNA_AGS_Results/Most_abundant_mRNA_in_AGS_computed_with_TPM_counts.png", width= 7, height = 5, dpi=300)
# 
# write.csv2(file="RNA_AGS_Results/Most_abundant_mRNA_in_AGS_computed_with_TPM_counts_infos.csv", table_top, row.names=T, quote=F)


# ########## EXPORTING ALSO *EVERY* TPM COUNT WITH INFO
# 
# ordered_infos<-dictionary
# row.names(ordered_infos)<-ordered_infos$Entry
# ordered_infos$Entry<-NULL
# ordered_infos<-ordered_infos[row.names(table), ]
# ordered_infos$Pathway<-NULL
# write.csv2(file=paste0("RNA_AGS_Results/Every_TPM_count_above_threshold_with_infos.csv"), cbind(table,ordered_infos), row.names=T, quote=F)




########### BAR PLOTS OF THE MOST ABUNDANT RNA (ONLY PAO AND GAO, ABS COUNTS) ###############

table <- as.data.frame( prop_with_unmap[ c( only_PAO_IDs$Entry, only_GAO_IDs$Entry), ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("domain-containing protein","domain-containing prot",infos,fixed=T)
infos<-gsub("CheX-like domain-containing prot","with CheX-like domain",infos,fixed=T)
infos<-gsub("OmpA-like domain-containing prot","Prot with OmpA-like domain",infos,fixed=T)
infos<-gsub("Competibacteraceae bacterium","Competibacteraceae",infos,fixed=T)
infos<-gsub("denitrificans Run_A_D11","denitr.",infos,fixed=T)
infos<-gsub("polymerase","pol",infos,fixed=T)
infos<-gsub("Candidatus","Ca.",infos,fixed=T)
infos<-gsub("sp.","",infos,fixed=T)
infos<-gsub("system","",infos,fixed=T)
infos<-gsub("Carboxypeptidase regulatory-like domain-containing prot","Carboxypeptidase regulatory-like domain prot",infos,fixed=T)
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.29, "cm"),
        legend.spacing.x = unit(0.25, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % (excluding other RNAs)", 
       fill="",
       title = "Most abundant PAO and GAO related mRNA",
       caption = "NB: those are ONLY the most abundant ones related to PAO and GAO")
ggsave("RNA_AGS_Results/Focus_on_PAO_and_GAO/Most_abundant_PAO_GAO_RNA_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)

write.csv2(file="RNA_AGS_Results/Focus_on_PAO_and_GAO/Most_abundant_PAO_GAO_RNA_abs_counts_PROPORTIONS_ON_THEMSELVES_infos.csv", table_top, row.names=T, quote=F)



########## EXPORTING ALSO *EVERY* ABS COUNT WITH INFO

ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(only_them), ]
ordered_infos$Pathway<-NULL
only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them, MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[colnames(only_them)!="other_AGS"], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Focus_on_PAO_and_GAO/Table_most_abundant_RNA_of_PAO_and_GAO_ABScounts_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)


suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))



########### BAR PLOTS OF THE MOST ABUNDANT GENES (ABS_count) ###############

table <- as.data.frame(prop_GENE_with_unmap) # prop abs
#colnames(table)<- paste0(metadata$Experiment_day, "th day")
#colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 
colnames(table)<- metadata$Experiment_day

top <- names(sort(rowSums(table[row.names(table)!="" , ]), decreasing=TRUE))[1:30]   # "" = unmapped
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Gene %in% table_top$Protein, ]
selected_infos<-selected_infos[! duplicated(selected_infos$Gene), ]
row.names(selected_infos)<-selected_infos$Gene
selected_infos<-selected_infos[table_top$Protein, ]
selected_infos$Description[grepl("phage",selected_infos$Taxon)]<-paste(selected_infos$Description[grepl("phage",selected_infos$Taxon)], "(virus)")
infos<-paste0( selected_infos$Description," (", selected_infos$Gene,")" )
infos<- gsub("()", "(no gene name)", infos, fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,3,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %", 
       fill="",
       title = "Most abundant genes in AGS mixed liquor",
       caption = "NB: different depth among samples --> if more depth then more 'others' genes!")
ggsave("RNA_AGS_Results/Most_abundant_GENES/Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)

# write.csv2(file="RNA_AGS_Results/Most_abundant_GENES/Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_WHOLE_DATASET_infos.csv", table_top, row.names=T, quote=F)




########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.42, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.95 ),
        legend.margin = margin(-5,0,-2,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="Experiment day", y="absolute count % (excluding other genes)", 
       fill="",
       title = "Most abundant genes in AGS mixed liquor",
       caption = "NB: these are ONLY the most abundant ones!")
ggsave("RNA_AGS_Results/Most_abundant_GENES/Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)
# write.csv2(file="RNA_AGS_Results/Most_abundant_GENES/Most_abundant_RNA_computed_with_abs__PROPORTIONS_ON_THEMSELVES_infos.csv", table_top, row.names=T, quote=F)



########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

AVERAGES_EVERY_SAMPLE<-apply( table_top[ colnames(table_top)!="Protein"], MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( table_top[ !colnames(table_top) %in% c("other_AGS","Protein")], MARGIN = 1, mean)
to_save <- cbind.data.frame( table_top[ , !grepl("^LR", colnames(table_top) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_GENES/Table_Most_Abundant_GENES_ABSpropAverages_WHOLE_DATASET_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them[ colnames(only_them)!="Protein"], MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[ !colnames(only_them) %in% c("other_AGS","Protein")], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_GENES/Table_Most_Abundant_GENES_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




#################### VEEN DIAGRAM (SEARCHING FOR UNIQUE RNA) ##########################

# avoiding the RNAs featured only in traces
min_CPM<-50 # this corresponds to 0.00005% min
min_samples<- 1
keep  <- rowSums(table_abs_CPM [!row.names(table_abs_CPM) %in% "UNMAPPED" , ] > min_CPM) >= min_samples   # NB: this table still includes also the NOT filtered reads and it's already in CPM 
table_abs_mains<-table_abs[keep, ]
cat( length(which(keep)), "maintained on", length(table_abs[,1]) ,"-->", (length(which(keep))/length(table_abs_CPM[,1]))*100 ,"% remaining")

Inoculum<-table_abs_mains[ , colnames(table_abs_mains) %in% "LR1" , drop=F]
Inoculum<-Inoculum[ rowSums(Inoculum)!=0 , , drop=F]
Inoculum<-row.names(Inoculum)

AGS<-table_abs_mains[ , colnames(table_abs_mains) %in% c("LR13","LR19","LR25","LR28","LR39") , drop=F]
AGS<-AGS[ rowSums(AGS)!=0 , , drop=F]
AGS<-row.names(AGS)

other_AGS<-table_abs_mains[ , colnames(table_abs_mains) %in% "other_AGS" , drop=F]
other_AGS<-other_AGS[ rowSums(other_AGS)!=0 , , drop=F]
other_AGS<-row.names(other_AGS)


# searching for unique RNA
ONLY_IN_AGS<- AGS[! AGS %in% other_AGS & ! AGS %in% Inoculum ]
ONLY_IN_AGS_list <- ONLY_IN_AGS # used for the bar plot (NB: it's not really a "list")
ONLY_IN_AGS<- paste(ONLY_IN_AGS, collapse = ", ")
head(ONLY_IN_AGS)

ONLY_IN_Inoculum<- Inoculum[! Inoculum %in% other_AGS & ! Inoculum %in% AGS ]
ONLY_IN_Inoculum<- paste(ONLY_IN_Inoculum, collapse = ", ")
head(ONLY_IN_Inoculum)

ONLY_IN_other_AGS<- other_AGS[! other_AGS %in% AGS & ! other_AGS %in% Inoculum ]
ONLY_IN_other_AGS_list <- ONLY_IN_other_AGS # used for the bar plot (NB: it's not really a "list")
ONLY_IN_other_AGS<- paste(ONLY_IN_other_AGS, collapse = ", ")
head(ONLY_IN_other_AGS)


IN_other_AGS_AND_this_AGS_AND_Inoculum <- AGS[AGS %in% other_AGS & AGS %in% Inoculum ] 
#head(IN_other_AGS_AND_this_AGS_AND_Inoculum)

IN_other_AGS_AND_this_AGS <- AGS[AGS %in% other_AGS ] 
#head(IN_other_AGS_AND_this_AGS)


con<-file("RNA_AGS_Results/Venn_Diagr_List_of_exclusive_RNA_and_intersection.txt")
sink(con, append=TRUE)
cat("ONLY IN the other AGS (16/06/2023):", fill=TRUE)
cat(ONLY_IN_other_AGS)
cat("\n\nONLY IN Inoculum:", fill=TRUE)
cat(ONLY_IN_Inoculum)
cat("\n\nONLY IN the current AGS:", fill=TRUE)
cat(ONLY_IN_AGS)
cat("\n\n\nIntersection --> (In other AGS, this AGS and Inoculum):", fill=TRUE)
cat(paste(IN_other_AGS_AND_this_AGS_AND_Inoculum, collapse = ",  "))
cat("\n\n\nIn both the other AGS and in this AGS):", fill=TRUE)
cat(paste(IN_other_AGS_AND_this_AGS, collapse = ",  "))
cat("\n\n\n")
cat("NB: only RNA with more than", min_CPM, "in at least", min_samples, "samples have included as features of the Venn Diagramm, to avoid focusing on unuseful traces", fill=T)
cat("Therefore only", length(which(keep)), "RNAs on", length(table_abs[,1]), "(filtered) have been used")
sink()
close(con)


x<-list("Floc inoculum  "=Inoculum, "  Current AGS"=AGS, "other AGS"=other_AGS)
ggvenn(x, stroke_size = 0.4, set_name_size = 2.8, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue","gold")) +
  theme(plot.title = element_text(size=8), plot.caption = element_text(size=7) ) +
  labs(title = "RNA identified among different samples")
ggsave(filename = "RNA_AGS_Results/Venn_Diagramm_EverySample.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()

x<-list("Current AGS"=AGS, "other AGS"=other_AGS )
ggvenn(x, stroke_size = 0.4, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue","gold")) +
  theme(plot.title = element_text(size=8.2), plot.caption = element_text(size=7) ) +
  labs(title = "RNA identified among different samples")
ggsave(filename = "RNA_AGS_Results/Venn_Diagramm_CurrentAGS_vs_OtherAGS.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()


suppressWarnings(rm(ONLY_IN_Inoculum, x, con, AGS, other_AGS, Inoculum))





########## BAR PLOT OF THE MOST ABUNDANT UNIQUE RNA OF ONE SAMPLE (AGS S46) USING Abs count in prop   ##########

table<- as.data.frame( prop_with_unmap[ row.names(prop_with_unmap) %in% ONLY_IN_other_AGS_list, colnames(table_abs)%in% "other_AGS", drop=F] )
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", , drop=F]), decreasing=TRUE))[1:30]
table_top<- table[top, , drop=F]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("domain-containing protein","domain-containing prot",infos,fixed=T)
infos<-gsub("polymerase","pol",infos,fixed=T)
infos<-gsub("Candidatus","Ca.",infos,fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-melt(table_top, id.vars = "Protein")
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + # coord polar turns the bar plot values into relative values by itself
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0), breaks = c(0,0.2,0.4,0.6)) +
  theme( #axis.text.x=element_blank(),
    axis.text.y=element_blank(), 
    # axis.title.y=element_blank(),
    axis.ticks.y=element_blank(), 
    axis.line.y = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.spacing.x = unit(0.28, "cm"),
    legend.title = element_text ( size = 9 ),
    legend.text = element_text ( size = 8 ),
    legend.margin = margin(-14,-25,-10,0),
    legend.position="bottom",
    plot.margin = margin(2,0,3,-7)
  ) +
  # scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
  guides(fill=guide_legend(nrow=30)) +
  labs(y="", x="Percentual abundance \nin the whole sample", 
       fill="",
       title = paste0("Most abundant unique RNA of the other reactor\n (16/06/23)" )
  )
ggsave(paste0("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_UNIQUE_RNA_in_the_otherAGS_fromVennDiagram.png"), width= 5.5, height = 5.8, dpi=300)



# ########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO
# 
# ordered_infos<-dictionary
# row.names(ordered_infos)<-ordered_infos$Entry
# ordered_infos$Entry<-NULL
# ordered_infos<-ordered_infos[row.names(table_top), ]
# ordered_infos$Pathway<-NULL
# 
# table_top$Protein<-NULL
# table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), , drop=F]
# AVERAGES_IN_OTHER_AGS<-apply( table_top, MARGIN = 1, mean)
# AVERAGES_IN_OTHER_SAMPLES<- rep(0, length(AVERAGES_IN_OTHER_AGS))
# to_save <- cbind.data.frame( table_top[ , !grepl("^LR", colnames(table_top) ), drop=F] , AVERAGES_IN_OTHER_AGS, AVERAGES_IN_OTHER_SAMPLES, ordered_infos )
# write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_RNAs/Table_Most_Abundant_UNIQUE_RNA_in_the_otherAGS_fromVennDiagram_with_infos.csv"), to_save, row.names=T, quote=F)
# 
# suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_IN_OTHER_AGS, AVERAGES_IN_OTHER_SAMPLES))
# 



########## BAR PLOT OF THE MOST ABUNDANT UNIQUE RNA (NOT FEATURED IN THE OTHER AGS) USING Abs prop   #############

# table<-table_abs
table <- as.data.frame(prop_with_unmap[ row.names(prop_with_unmap) %in% ONLY_IN_AGS_list , ! colnames(prop_with_unmap) %in% "other_AGS"  ] ) # prop abs
colnames(table)<- paste0(metadata$Experiment_day[! metadata$Experiment_day %in% "other_AGS"] , "th day")

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", ]), decreasing=TRUE))[1:26] # NB: only 26 unique RNAs!
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("domain-containing protein","domain-containing prot",infos,fixed=T)
infos<-gsub("Candidatus","Ca. ",infos,fixed=T)
infos<-gsub("Prepilin-type N-term cleavage/methylation domain-containing prot","Prepilin-type prot with cleav/methylat domain",infos,fixed=T)
infos<-gsub(" bacterium","",infos,fixed=T)
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated

table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30[1:26]) + # only 26 RNAs
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=13)) +
  labs(x="", y="absolute count % (excluding other RNAs)",
       fill="",
       title = "Every unique RNA of this AGS reactor (not in the other AGS reactor)",
       caption = "NB: The inoculum (flocs) RNAs are excluded from this list"
       )
ggsave("RNA_AGS_Results/Most_abundant_RNAs/Table_Most_Abundant_UNIQUE_RNA_in_this_AGS_fromVennDiagram.png", width= 7, height = 5, dpi=300)
# NB: the warning is expected because these most abundant ones are not featured in the inoculum!




########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top[, !colnames(table_top)%in% "1th day"], MARGIN = 1, mean), decreasing = T), ]
AVERAGES_FROM_2nd_to_last_SAMPLE<-apply( table_top, MARGIN = 1, mean)
AVERAGES_IN_OTHER_AGS<-rep(0, length(AVERAGES_FROM_2nd_to_last_SAMPLE))

only_them$Protein<-NULL
# only_them$`1th day`<-rep(0) # to overwrite the NaN assigned during its creation (not important for the plot)
only_them <- only_them[order(apply( only_them [, !colnames(table_top)%in% "1th day"], MARGIN = 1, mean), decreasing = T), ]
AVERAGES_RECOMPUTED_USING_ONLY_THEM<-apply( only_them, MARGIN = 1, mean)
to_save <- cbind.data.frame( table_top[ , !grepl("^LR", colnames(table_top) ), drop=F] , AVERAGES_FROM_2nd_to_last_SAMPLE, AVERAGES_RECOMPUTED_USING_ONLY_THEM, AVERAGES_IN_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_RNAs/Table_Most_Abundant_UNIQUE_RNA_in_this_AGS_fromVennDiagram_with_infos.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_FROM_2nd_to_last_SAMPLE, AVERAGES_RECOMPUTED_USING_ONLY_THEM))




############ ORDINATION PLOTS ON BOTH TABLES (INCLUDING THE OTHER AGS) ##############

# Hereby a loop will be used because no custom setting are required for each element


list_tables<-list( table_abs [!row.names(table_abs) %in% "UNMAPPED", ]  ) # , table_TPM )
#                       1                                                               2

for (number in 1:length(list_tables) ) {
  
  table<-list_tables[[number]]
  
  name_type<- ifelse( number == 1 ,   # n1 = table abs
                      "ABS_counts",
                      "TPM_counts")
  cat("\n Analysing",name_type,"data... \n\n")
  
  # NB: ggplot will create the subfolders by itself...
  
  
  # Custom names
  colnames(table)<-paste0(metadata$Experiment_day, "th")  # already set the same order beforehand
  colnames(table)<-gsub("AGSth", "AGS", colnames(table))
  colnames(table)<-gsub("_", "\n", colnames(table))
  
  
  # Where the same sum among samples is required ...
  if( name_type == "ABS_counts"){
    prop<- prop_with_NO_unmap
  } else {
    prop <- table  # the TPM counts already have with the same total in each sample
  }
  
  
  # Preparing also the DESeq2 normalized objects ...
  if( name_type == "ABS_counts"){
    cts<-as.matrix(table[rowSums(table)>10, ])
    coldata<- cbind( Fake_groups=rep("none", length(colnames(table))) )
    row.names(coldata) <- colnames(table)
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ 1)               
    dds2 <- estimateSizeFactors(dds)
    DS2_counts <- counts(dds2, normalized=TRUE)
  }     
  
  
  # To color the plots
  colors_shade <- metadata$Experiment_day
  colors_shade[colors_shade=="other_AGS"] <- 165 # to set it as a strong color (mature reactor)
  colors_shade<-as.numeric(colors_shade)
  
  
  
  ### HELLINGER PCOA
  
  dist <- vegdist(t(sqrt(prop)),  method = "euclidean")   # sqrt of prop --> Hellinger
  PCOA<-pco(dist)
  ggplot(data=PCOA$vectors, aes(x=PCOA$vectors[,1],y=PCOA$vectors[,2] , color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs( x=paste("PC1:", round((PCOA$values[1]/sum(PCOA$values))*100 ,2) ,"%"),
          y=paste("PC2:", round((PCOA$values[2]/sum(PCOA$values))*100 ,2) ,"%"),
          title=paste0("PCoA on Hellinger distance (",name_type,")")
    )
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"/PCoA_Hellinger_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ### BRAY PCOA
  
  dist <- vegdist(t(prop),  method = "bray")
  PCOA<-pco(dist)
  ggplot(data=PCOA$vectors, aes(x=PCOA$vectors[,1],y=PCOA$vectors[,2] , color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( ) +
    labs( x=paste("PC1:", round((PCOA$values[1]/sum(PCOA$values))*100 ,2) ,"%"),
          y=paste("PC2:", round((PCOA$values[2]/sum(PCOA$values))*100 ,2) ,"%"),
          title=paste0("PCoA on Bray-Curtis distance (",name_type,")")
    ) +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none")
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"/PCoA_Bray_proportions_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ##### Bray-Curtis NMDS

  set.seed(1)
  x <-metaMDS(t(prop), dist="bray", trace=FALSE)  # trace=F --> quiet
  data.scores <- as.data.frame(scores(x)$sites)
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs(
      title=paste0("NMDS on Bray-Curtis dissimilarity index (",name_type,")")
    )
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"/NMDS_Bray_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ##### Hellinger NMDS
  set.seed(1)
  x <-metaMDS(t(sqrt(prop)), dist="euclidean", trace=FALSE)  # sqrt of prop --> Hellinger
  data.scores <- as.data.frame(scores(x)$sites)
  new_names <- gsub ("\n", " ", colnames(table))   # for aesthetics reasons
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=new_names),color="grey35", vjust= -0.3, hjust= 0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs( title=paste0("NMDS on Hellinger distance index (",name_type,")")
    )
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"/NMDS_Hellinger_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  
  ##### PCA with Deseq2 median of ratio normalization 
  if( name_type == "ABS_counts"){
    to_plot<-prcomp(t(DS2_counts))
    plot_stats<-as.data.frame(summary(to_plot)$importance)
    new_names <- gsub ("\n", " ", colnames(table))
    ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade)) +
      scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                            breaks=c(1,40,80,120,150,184),
                            midpoint = 25) +
      theme_classic( )  +
      geom_point(size=4.3, alpha=0.3) +
      geom_point(size=2.3, alpha=1) +
      geom_text(mapping = aes(label=new_names),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
      guides(color="none") +
      labs( title="PCA on DESeq2 median of ratios counts", 
            x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
            y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
      )
    ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"/PCA_DESeq2Normalization_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  }
  
  
  
  ##### PCA with Deseq2 VST normalization 
  if( name_type == "ABS_counts"){
    vst_data <- vst(dds, blind=TRUE )    # it accounts for different library sizes by itself, see https://support.bioconductor.org/p/9138134/
    # blind=T -> no within-group normalizat (unsupervised)
    to_plot<-prcomp(t(assay(vst_data)))
    plot_stats<-as.data.frame(summary(to_plot)$importance)
    new_names <- gsub ("\n", " ", colnames(table))
    ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade)) +
      scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                            breaks=c(1,40,80,120,150,184),
                            midpoint = 25) +
      theme_classic( )  +
      geom_point(size=4.3, alpha=0.3) +
      geom_point(size=2.3, alpha=1) +
      geom_text(mapping = aes(label=new_names),color="grey35", vjust= -0.3, hjust= 0.61, size= 3, lineheight = 0.7) +
      guides(color="none") +
      labs( title="PCA on vst normalized counts",
            x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
            y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
      )
    ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"/PCA_vst_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  }
  
}




############ ORDINATION PLOTS ON BOTH TABLES (EXCLUDING THE OTHER AGS) ##############


list_tables<-list( table_abs [!row.names(table_abs) %in% "UNMAPPED", ]  ) # , table_TPM )
#                       1                                                               2


for (number in 1:length(list_tables) ) {
  
  table<-list_tables[[number]]
  
  name_type<- ifelse( number == 1 ,   # n1 = table abs
                      "ABS_counts",
                      "TPM_counts")
  cat("\n Analysing",name_type,"data without the other AGS... \n\n")
  
  # NB: ggplot will create the subfolders by itself...
  
  
  
  # Removing here the other AGS sample
  table <- table[ , colnames(table)!= "other_AGS"]
  
  
  
  # Where the same sum among samples is required ...
  if( name_type == "ABS_counts"){
    prop<- prop_with_NO_unmap[ , ! colnames(prop_with_NO_unmap) %in% "other_AGS" ]
  } else {
    prop <- table  # the TPM counts already have with the same total in each sample
  }
  
  
  
  colnames(table)<-paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "th")  # already set the same order beforehand
  
  
  
  
  # Preparing also the DESeq2 normalized objects ...
  if( name_type == "ABS_counts"){
    cts<-as.matrix(table[rowSums(table)>10, ])
    coldata<- cbind( Fake_groups=rep("none", length(colnames(table))) )
    row.names(coldata) <- colnames(table)
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ 1)               
    dds2 <- estimateSizeFactors(dds)
    DS2_counts <- counts(dds2, normalized=TRUE)
  }     
  
  
  # To color the plots
  colors_shade <- metadata$Experiment_day[metadata$Experiment_day!="other_AGS"]
  colors_shade<-as.numeric(colors_shade)
  
  
  
  ### HELLINGER PCOA
  
  dist <- vegdist(t(sqrt(prop)),  method = "euclidean")   # sqrt of prop --> Hellinger
  PCOA<-pco(dist)
  ggplot(data=PCOA$vectors, aes(x=PCOA$vectors[,1],y=PCOA$vectors[,2] , color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs( x=paste("PC1:", round((PCOA$values[1]/sum(PCOA$values))*100 ,2) ,"%"),
          y=paste("PC2:", round((PCOA$values[2]/sum(PCOA$values))*100 ,2) ,"%"),
          title=paste0("PCoA on Hellinger distance (",name_type,")")
    )
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"_excluding_otherAGS/PCoA_Hellinger_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ### BRAY PCOA
  
  dist <- vegdist(t(prop),  method = "bray")
  PCOA<-pco(dist)
  ggplot(data=PCOA$vectors, aes(x=PCOA$vectors[,1],y=PCOA$vectors[,2] , color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( ) +
    labs( x=paste("PC1:", round((PCOA$values[1]/sum(PCOA$values))*100 ,2) ,"%"),
          y=paste("PC2:", round((PCOA$values[2]/sum(PCOA$values))*100 ,2) ,"%"),
          title=paste0("PCoA on Bray-Curtis distance (",name_type,")")
    ) +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none")
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"_excluding_otherAGS/PCoA_Bray_proportions_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ##### Bray-Curtis NMDS
  
  set.seed(1)
  x <-metaMDS(t(prop), dist="bray", trace=FALSE)  # trace=F --> quiet
  data.scores <- as.data.frame(scores(x)$sites)
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs(
      title=paste0("NMDS on Bray-Curtis dissimilarity index (",name_type,")")
    )
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"_excluding_otherAGS/NMDS_Bray_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ##### Hellinger NMDS
  set.seed(1)
  x <-metaMDS(t(sqrt(prop)), dist="euclidean", trace=FALSE)  # sqrt of prop --> Hellinger
  data.scores <- as.data.frame(scores(x)$sites)
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs( title=paste0("NMDS on Hellinger distance index (",name_type,")")
    )
  ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"_excluding_otherAGS/NMDS_Hellinger_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  
  
  ##### PCA with Deseq2 VST normalization 
  if( name_type == "ABS_counts"){
    vst_data <- vst(dds, blind=TRUE )  # it accounts for different library sizes by itself, see https://support.bioconductor.org/p/9138134/
    # blind=T -> no within-group normalizat (unsupervised)
    to_plot<-prcomp(t(assay(vst_data)))
    plot_stats<-as.data.frame(summary(to_plot)$importance)
    ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade)) +
      scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                            breaks=c(1,40,80,120,150,184),
                            midpoint = 25) +
      theme_classic( )  +
      geom_point(size=4.3, alpha=0.3) +
      geom_point(size=2.3, alpha=1) +
      geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
      guides(color="none") +
      labs( title="PCA on vst normalized counts",
            x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
            y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
      )
    ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"_excluding_otherAGS/PCA_vst_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  }
  
  
  
  ##### PCA with Deseq2 median of ratio normalization 
  if( name_type == "ABS_counts"){
    to_plot<-prcomp(t(DS2_counts))
    plot_stats<-as.data.frame(summary(to_plot)$importance)
    ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade)) +
      scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                            breaks=c(1,40,80,120,150,184),
                            midpoint = 25) +
      theme_classic( )  +
      geom_point(size=4.3, alpha=0.3) +
      geom_point(size=2.3, alpha=1) +
      geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
      guides(color="none") +
      labs( title="PCA on DESeq2 median of ratios counts",
            x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
            y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
      )
    ggsave(paste0("RNA_AGS_Results/ordinations_",name_type,"_excluding_otherAGS/PCA_DESeq2Normalization_on_",name_type,".png"), width= 5.2, height = 4.5, dpi=300)
  }
  
  
}




############ ORDINATION PLOTS ON BOTH TABLES (ONLY PAO AND GAO SOURCED RNAs) ##############

target <- table_abs[ c(only_GAO_IDs$Entry, only_PAO_IDs$Entry), ]
list_tables<-list( target ) #, table_TPM )   # only ABS counts
#                       1          2

for (number in 1:length(list_tables) ) {
  
  table<-list_tables[[number]]
  
  name_type<- ifelse( number == 1 ,   # n1 = table abs
                      "ABS_counts",
                      "TPM_counts")
  cat("\n Analysing",name_type,"data without the non PAO and GAO sourced RNAs... \n\n")
  
  # NB: ggplot will create the subfolders by itself...
  
  
  # Custom names
  colnames(table)<-paste0(metadata$Experiment_day, "th")  # already set the same order beforehand
  colnames(table)<-gsub("AGSth", "AGS", colnames(table))
  colnames(table)<-gsub("_", "\n", colnames(table))
  
  
  
  # Where the same sum among samples is required ...
  if( name_type == "ABS_counts"){
    prop<- prop_with_NO_unmap
  } else {
    prop <- table  # the TPM counts already have with the same total in each sample
  }
  
  
  
  # Preparing also the DESeq2 normalized objects ...
  if( name_type == "ABS_counts"){
    cts<-as.matrix(table[rowSums(table)>10, ])
    coldata<- cbind( Fake_groups=rep("none", length(colnames(table))) )
    row.names(coldata) <- colnames(table)
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ 1)               
    dds2 <- estimateSizeFactors(dds)
    DS2_counts <- counts(dds2, normalized=TRUE)
  }     
  
  
  # To color the plots
  colors_shade <- metadata$Experiment_day
  colors_shade[colors_shade=="other_AGS"] <- 165 # to set it as a strong color (mature reactor)
  colors_shade<-as.numeric(colors_shade)
  
  
  
  ### HELLINGER PCOA
  
  dist <- vegdist(t(sqrt(prop)),  method = "euclidean")   # sqrt of prop --> Hellinger
  PCOA<-pco(dist)
  ggplot(data=PCOA$vectors, aes(x=PCOA$vectors[,1],y=PCOA$vectors[,2] , color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( )  +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none") +
    labs( x=paste("PC1:", round((PCOA$values[1]/sum(PCOA$values))*100 ,2) ,"%"),
          y=paste("PC2:", round((PCOA$values[2]/sum(PCOA$values))*100 ,2) ,"%"),
          title=paste0("PCoA on Hellinger distance (only PAOs and GAOs)")
    )
  ggsave(paste0("RNA_AGS_Results/Focus_on_PAO_and_GAO/PCoA_Hellinger_on_",name_type,"_only_PAO_GAO.png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ### BRAY PCOA
  
  dist <- vegdist(t(prop),  method = "bray")
  PCOA<-pco(dist)
  ggplot(data=PCOA$vectors, aes(x=PCOA$vectors[,1],y=PCOA$vectors[,2] , color= colors_shade)) +
    scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                          breaks=c(1,40,80,120,150,184),
                          midpoint = 25) +
    theme_classic( ) +
    labs( x=paste("PC1:", round((PCOA$values[1]/sum(PCOA$values))*100 ,2) ,"%"),
          y=paste("PC2:", round((PCOA$values[2]/sum(PCOA$values))*100 ,2) ,"%"),
          title=paste0("PCoA on Bray-Curtis distance (only PAOs and GAOs)")
    ) +
    geom_point(size=4.3, alpha=0.3) +
    geom_point(size=2.3, alpha=1) +
    geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
    guides(color="none")
  ggsave(paste0("RNA_AGS_Results/Focus_on_PAO_and_GAO/PCoA_Bray_proportions_on_",name_type,"_only_PAO_GAO.png"), width= 5.2, height = 4.5, dpi=300)
  
  
  
  ##### PCA with Deseq2 median of ratio normalization 
  if( name_type == "ABS_counts"){
    to_plot<-prcomp(t(DS2_counts))
    plot_stats<-as.data.frame(summary(to_plot)$importance)
    ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade)) +
      scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                            breaks=c(1,40,80,120,150,184),
                            midpoint = 25) +
      theme_classic( )  +
      theme(plot.title = element_text(size=12.5)) +
      geom_point(size=4.3, alpha=0.3) +
      geom_point(size=2.3, alpha=1) +
      geom_text(mapping = aes(label=colnames(table)),color="grey35", vjust= -0.3, size= 3, lineheight = 0.7) +
      guides(color="none") +
      labs( title="PCA on DESeq2 median of ratios (only PAOs and GAOs)", 
            x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
            y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
      )
    ggsave(paste0("RNA_AGS_Results/Focus_on_PAO_and_GAO/PCA_DESeq2Normalization_on_",name_type,"_only_PAO_GAO.png"), width= 5.2, height = 4.5, dpi=300)
  }
  
  
  
  ##### PCA with Deseq2 VST normalization 
  if( name_type == "ABS_counts"){
    vst_data <- vst(dds, blind=TRUE )  # it accounts for different library sizes by itself, see https://support.bioconductor.org/p/9138134/
    # blind=T -> no within-group normalizat (unsupervised)
    to_plot<-prcomp(t(assay(vst_data)))
    plot_stats<-as.data.frame(summary(to_plot)$importance)
    new_names <- gsub ("\n", " ", colnames(table))
    ggplot(as.data.frame(to_plot$x), aes(x = PC1, y = PC2, color= colors_shade)) +
      scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                            breaks=c(1,40,80,120,150,184),
                            midpoint = 25) +
      theme_classic( )  +
      geom_point(size=4.3, alpha=0.3) +
      geom_point(size=2.3, alpha=1) +
      geom_text(mapping = aes(label=new_names),color="grey35", vjust= -0.3, hjust= 0.61, size= 3, lineheight = 0.7) +
      guides(color="none") +
      labs( title="PCA on vst normalized counts of PAOs and GAOs",
            x=paste("PC1:", round( (plot_stats[2,1])*100, 2) ,"%"),
            y=paste("PC2:", round( (plot_stats[2,2])*100, 2) ,"%"),
      )
    ggsave(paste0("RNA_AGS_Results/Focus_on_PAO_and_GAO/PCA_vst_on_",name_type,"_only_PAO_GAO.png"), width= 5.2, height = 4.5, dpi=300)
  }
  
  
}




#################### CORRELATING RNA counts vs TIME ########################

# selecting only the most abundant and present RNAs (transcripts)
min_CPM<-10 # this corresponds to 0.00001% min in CPM
min_samples<- 3
keep  <- rowSums(table_abs_CPM [ !row.names(table_abs_CPM) %in% "UNMAPPED" ,! colnames(table_abs_CPM) %in% "other_AGS"] > min_CPM) >= min_samples  # NB: this table still includes also the NOT filtered reads and it's already in CPM
# NB: keep is a logical from an object with different length... but its names (named vector) are always valid!
table_abs_mains<-as.data.frame(prop_with_NO_unmap)[row.names(prop_with_NO_unmap) %in% names(which(keep)), ] # NB: filter on CPM but classic proportions will be used
cat( length(which(keep)), "maintained on", length(table_abs[,1]) ,"-->", (length(which(keep))/length(table_abs_CPM[,1]))*100 ,"% remaining")

abundances<- table_abs_mains[ , ! colnames(table_abs_mains) %in% "other_AGS"]
colnames(abundances)  # checking that they are properly ordered...
abundances<-round(abundances, digits = 3) # NB: these are not CPM, they are percentuals... because it's no use keeping this level of "details" in the decimals (noise)

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  cat("\r ...testing the", i, "th row on",  length(row.names(abundances)), "...") 
  save<-cor.test( as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "kendall") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
# Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
# Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")
Corr_results$sign<-ifelse(Corr_results$pvalue<0.05,"*","")                # NB: classic p value!
#Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")



temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
# temp_for_save

temp_for_save
ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(temp_for_save), ]
# ordered_infos$Pathway<-NULL
temp_for_save<-cbind.data.frame(temp_for_save, ordered_infos)
write.csv2(temp_for_save, "RNA_AGS_Results/Correlations/Signif_Kendall_correlation_RNA_ABSproportialCounts_vs_time.csv", row.names = T, quote = F)

con <- file("RNA_AGS_Results/Correlations/README_Only_certain_RNA_have_been_tested_in_correlation.txt")
sink(con, append = T)
cat("Only RNA with more than", min_CPM, "CPM in at least", min_samples, "time points have been correlated with the passing of time", fill=T)
cat("Therefore only", length(which(keep)), "RNAs on", length(table_abs[,1]), "(filtered) have been tested")
sink()
close(con)

# used below...
significative_abundances_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )


suppressWarnings(rm(con, Corr_results, min_CPM, min_samples, save, new_row))




########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE RNA  MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_top<-significative_abundances_pos[1:12] # the table is already ordered

table_to_plot <- as.data.frame(abundances[corr_top, ])
infos <- temp_for_save[ corr_top, ]
table_to_plot$Entry <- paste0( infos$Description," (", infos$Taxon,")" )
table_to_plot$Entry<- gsub("()", "(Unknown source)", table_to_plot$Entry, fixed=T)
# table_to_plot$Entry<-gsub("Candidatus","Ca.",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Candidatus ","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("denitrificans","denitr.",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub(" bacterium","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub(" Run_A_D11","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("(RNA-pol)-subunit","(RNA-pol)",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("secretion prot IcmF","secretion IcmF",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Alpha","α",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("nitrate reductase","NO3- reduct",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Bifunctional ","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("phosphodiesterase","phosphodiest",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("cyclase","cycl",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Type VI","TypeVI",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("^digua","Digua",table_to_plot$Entry,fixed=F)
table_to_plot$Entry<-gsub("^prot","Protein",table_to_plot$Entry,fixed=F)
table_to_plot$Entry<-gsub("^h","H",table_to_plot$Entry,fixed=F)


#colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "th day") ,"Entry" ) # already set the same order beforehand
colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "") ,"Entry" ) # already set the same order beforehand
table_to_plot<-melt(table_to_plot, id.vars = "Entry")

ggplot(data=table_to_plot, aes(y=value, x=variable, color=Entry)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Entry), size =2) +
  geom_point(aes(color=Entry), size =3.8, alpha=0.4) +
  geom_line(aes(group=Entry), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_19[3:18]) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=8),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 8.5 ),
        legend.position="bottom",
        legend.margin = margin(1,42,0,0)
        ) +
  guides(color=guide_legend(nrow=6)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  scale_y_continuous(breaks = seq( 0, max(table_to_plot$value)+0.01, 0.01) ) +
  labs(x="Experiment day", y="ABS counts %",
       color="",
       title = "RNA most positely correlated with time\n according to Kendall correlation on ABS proportional counts")
ggsave(file="RNA_AGS_Results/Correlations/RNA_most_positively_correlated_with_time.png",width=6.2,height=5, dpi=300) 
dev.off()



# THE RNA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
# corr_top<-significative_abundances_neg[(length(significative_abundances_neg)-9):length(significative_abundances_neg)]  # 10 
corr_top<-significative_abundances_neg[1:12]

table_to_plot <- as.data.frame(abundances[corr_top, ])
infos <- temp_for_save[ corr_top, ]
table_to_plot$Entry <- paste0( infos$Description," (", infos$Taxon,")" )
table_to_plot$Entry<- gsub("()", "(Unknown source)", table_to_plot$Entry, fixed=T)
table_to_plot$Entry<-gsub("(Fragment)","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Candidatus","Ca.",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("denitrificans","denitr.",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("tube prot ","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub(" bacterium","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Replicase large subunit","Replicase subunit",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub(" rugose","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry[duplicated(table_to_plot$Entry)] <- paste( table_to_plot$Entry[duplicated(table_to_plot$Entry)] , row.names(table_to_plot)[duplicated(table_to_plot$Entry)] )

#colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "th day") ,"Entry" ) # already set the same order beforehand
colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "") ,"Entry" ) # already set the same order beforehand
table_to_plot<-melt(table_to_plot, id.vars = "Entry")

ggplot(data=table_to_plot, aes(y=value, x=variable, color=Entry)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Entry), size =2) +
  geom_point(aes(color=Entry), size =3.8, alpha=0.4) +
  geom_line(aes(group=Entry), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_19) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=8),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 8.6 ),
        legend.position="bottom",
        legend.margin = margin(1,42,0,0)) +
  guides(color=guide_legend(nrow=6)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  # scale_y_continuous(breaks = seq( 0, max(table_to_plot$value)+0.01, 0.01) , limits = c(0, 0.03) ) +
  labs(x="Experiment day", y="ABS counts %",
       color="",
       title = "RNA most negatively correlated with time\n according to Kendall correlation on ABS proportional counts")
ggsave(file="RNA_AGS_Results/Correlations/RNA_most_negatively_correlated_with_time.png",width=6.2,height=5, dpi=300) 
dev.off()


suppressWarnings(rm(prune.dat_top, tax_selected, tabella, corr_top))




#################### CORRELATING RNA counts vs MEASUREMENTS ########################

# selecting only the most abundant and present RNAs
min_CPM<-10 # this corresponds to 0.00001% min in CPM
min_samples<- 3
keep  <- rowSums(table_abs_CPM [ !row.names(table_abs_CPM) %in% "UNMAPPED" ,! colnames(table_abs_CPM) %in% "other_AGS"] > min_CPM) >= min_samples  # NB: this table still includes also the NOT filtered reads and it's already in CPM
table_abs_mains<-as.data.frame(prop_with_NO_unmap)[row.names(prop_with_NO_unmap) %in% names(which(keep)), ] # NB: filter on CPM but classic proportions will be used
cat( length(which(keep)), "maintained on", length(table_abs[,1]) ,"-->", (length(which(keep))/length(table_abs_CPM[,1]))*100 ,"% remaining")

abundances<- table_abs_mains[ , ! colnames(table_abs_mains) %in% "other_AGS"]
colnames(abundances)  # checking that they are properly ordered...
abundances<-round(abundances, digits = 3) # NB: these are not CPM, they are percertuals... no use keeping this level of "details" in the decimals (noise)
abundances<-t(abundances)

##### Importing the measurements metadata
meas <- read.delim("Reactor_measurements.tsv")
meas$Sample_name<-gsub("L","LR", meas$Sample_name)
row.names(meas)<- meas$Sample_name
meas<-meas[row.names(abundances), ] # subsetting and re-orderdering the samples
# identical(row.names(meas), colnames(abundances))
meas$Sample_name<-NULL
meas<-round(meas, digits = 1)
# adjusting few names for better aesthetics...
colnames(meas)<-gsub("VSS.TSS", "VSS/TSS", colnames(meas))
colnames(meas)<-gsub("SVI30.5", "SVI30/SVI5", colnames(meas))
colnames(meas)<-gsub(".1", "", colnames(meas), fixed = T) # this is due a little importing issue...
colnames(meas)<-gsub(".", " ", colnames(meas), fixed = T)



#### Computing the correlations...
corr<-NULL
for( x in 1:length( colnames(abundances)) ){
  
  for( y in 1:length( colnames(meas)) ) {
    if( length(meas[[y]] [ ! is.na( meas[[y]] ) ]) >=3 ){   # to avoid the error "non ci sono abbastanza osservazioni finite"
      loop<- suppressWarnings( cor.test( as.numeric(abundances[, x]) , meas[[y]] , method = "kendall" ) )
      # NB: warnings due to the ties... this is normal!
      new_row<-c(colnames(abundances)[x], colnames(meas)[y] , round(as.numeric(loop$estimate),2), round(as.numeric(loop$p.value),6) )
      corr<-rbind(corr, new_row)
    }
  }
  
}

corr<-as.data.frame(corr)
row.names(corr)<-NULL
colnames(corr)<-c("RNA","Measurements","Corr","pvalue")
corr <- corr[!is.na(corr$Corr), ] # NA obtained if standard deviation of a variable is 0
corr$padj_bh<-p.adjust(corr$pvalue, method = "BH")
corr$Corr<-as.numeric(corr$Corr)
corr$Sign<-corr$pvalue
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""

# # too many RNA to plot in an heatplot!
# ggplot(corr, aes(x = RNA, y = Measurements, fill = Corr)) +
#   geom_tile(color = "white", lwd = 0.5, linetype = 1) +
#   scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
#   theme_bw(base_size=12) +
#   theme(axis.text.y=element_text(angle = -35, hjust = 1, vjust=1 , size= 9),
#         axis.text.x=element_text(angle = 40, hjust = 1, vjust=1 , size= 9.5),
#         plot.title = element_text(size=11.2),
#         legend.text = element_text(size=10), legend.key.height= unit(2, 'cm'),
#         legend.key.width= unit(0.7, 'cm'), legend.title = element_text(size=12)
#   ) +
#   guides(fill= guide_colourbar(title ="rho")) +
#   geom_text(aes(label= Sign), color= "white", size =10, vjust=0.8) +
#   labs(title = "Kendall correlations (CLR abundances vs measurements)", y= "ML measurements", x= "ML most abundant genera", 
#        caption= "\n adjusted p-value (holm) lower than 0.05 are displayed through * sign")
# ggsave(file="Results/Correlations/Correlations_Bacteria_VS_measurements_Kendall_Holm_usingCLRabundances.png", dpi=300, width = 7, height = 6.5)


# exporting
ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[corr$RNA, ]
corr<-cbind.data.frame(corr, ordered_infos)
write.csv2(corr, file="RNA_AGS_Results/Correlations/Correlations_RNA_VS_measurements_Kendall_table_usingABSpropCounts.csv", quote = F, na="", row.names = F)



### testing "visually" these correlations ...

signif <- as.data.frame( corr[corr$Sign=="*" , ] )
signif$Description <- gsub("domain-containing protein","domain-containing prot", signif$Description ,fixed=T)
signif$Description <- gsub("polymerase","pol", signif$Description, fixed=T)
signif$Taxon <-gsub("Candidatus","Ca.", signif$Taxon, fixed=T)
signif$Taxon <-gsub("denitrificans","denitr.", signif$Taxon, fixed=T)
signif$Taxon <-gsub(" bacterium","", signif$Taxon, fixed=T)


pdf(file = "RNA_AGS_Results/Correlations/EXTRA_to_visually_test_signif_correlations_RNA_vs_Measurem.pdf")
par(mfrow = c(2,3))
for(x in 1:length(row.names(signif))){
  target_RNA<-signif[x, "RNA"]
  target_measur<-signif[x, "Measurements"]
  # a<-as.numeric( abundances[ ,target_bacteria ] )  # this are post clr... it's better to show the proportions!
  a <- abundances[, target_RNA]
  a <- round( as.numeric(a), 2 )
  b<-meas [ ,target_measur ]
  plot( a~b , ylab= paste0( signif[x, "Description"], " (",signif[x, "Taxon"],")" ," %") , xlab=target_measur, col="red", cex=0.65)
  title("")
  abline( lm(a~b) )
}
dev.off()




#################### CORRELATING GENE counts vs TIME ########################

# selecting only the most abundant and present GENEs
min_CPM<-0.005 # NB: much higher then the threshold used for the transcript!
min_samples<- 3
keep  <- rowSums(prop_GENE_with_unmap [  ,! colnames(prop_GENE_with_unmap) %in% "other_AGS"] > min_CPM) >= min_samples  # NB: this table still includes also the NOT filtered reads and it's already in %
# NB: keep is a logical from an object with different length... but its names (named vector) are always valid!
table_abs_mains<-as.data.frame(prop_GENE_with_unmap)[row.names(prop_GENE_with_unmap) %in% names(which(keep)), ] # NB: filter on CPM but classic proportions will be used
cat( length(which(keep)), "maintained on", length(prop_GENE_with_unmap[,1]) ,"-->", (length(which(keep))/length(prop_GENE_with_unmap[,1]))*100 ,"% remaining")

abundances<- table_abs_mains[ , ! colnames(table_abs_mains) %in% "other_AGS"]
colnames(abundances)  # checking that they are properly ordered...
abundances<-round(abundances, digits = 3) # NB: these are not CPM, they are percentuals... because it's no use keeping this level of "details" in the decimals (noise)

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  cat("\r ...testing the", i, "th row on",  length(row.names(abundances)), "...") 
  save<-cor.test( as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "kendall") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
# Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
# Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")
Corr_results$sign<-ifelse(Corr_results$pvalue<0.05,"*","")                # NB: classic p value!
#Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")



temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
# temp_for_save

temp_for_save
ordered_infos<-dictionary
ordered_infos<-ordered_infos[ordered_infos$Gene %in% row.names(temp_for_save), ]
ordered_infos<-ordered_infos[!duplicated(ordered_infos$Gene), ]
row.names(ordered_infos)<-ordered_infos$Gene
# ordered_infos$Pathway<-NULL
ordered_infos<-ordered_infos[row.names(temp_for_save), ]
temp_for_save<-cbind.data.frame(temp_for_save, ordered_infos)
write.csv2(temp_for_save, "RNA_AGS_Results/Correlations/Signif_Kendall_correlation_GENES_ABSproportialCounts_vs_time.csv", row.names = T, quote = F)

con <- file("RNA_AGS_Results/Correlations/README_Only_certain_GENES_have_been_tested_in_correlation.txt")
sink(con, append = T)
cat("Only Genes with more than", min_CPM, "% count in at least", min_samples, "time points have been correlated with the passing of time", fill=T)
cat("Therefore only", length(which(keep)), "Genes on", length(table_abs[,1]), "(filtered) have been tested")
sink()
close(con)

# used below...
significative_abundances_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )


suppressWarnings(rm(con, Corr_results, min_CPM, min_samples, save, new_row))




########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE RNA  MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_top<-significative_abundances_pos[1:12] # the table is already ordered

table_to_plot <- as.data.frame(abundances[corr_top, ])
infos <- temp_for_save[ corr_top, ]
table_to_plot$Entry <- paste0( infos$Description," (", infos$Gene,")" )
# table_to_plot$Entry<- gsub("()", "(Unknown source)", table_to_plot$Entry, fixed=T)
# table_to_plot$Entry<-gsub("Candidatus ","",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub("denitrificans","den.",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub(" bacterium","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("serine/threonine prot","ser/thr",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("(RNA-pol)-subunit","RNApol-sub",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub("associated domain","domain",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub("Putative ","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Type VI secretion prot","TypeVI secretion",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("Alpha","α",table_to_plot$Entry,fixed=T)

# colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "th day") ,"Entry" ) # already set the same order beforehand
colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "") ,"Entry" ) # already set the same order beforehand
table_to_plot<-melt(table_to_plot, id.vars = "Entry")

ggplot(data=table_to_plot, aes(y=value, x=variable, color=Entry)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Entry), size =2) +
  geom_point(aes(color=Entry), size =3.8, alpha=0.4) +
  geom_line(aes(group=Entry), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_19[3:17]) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=8),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 8.55 ),
        legend.position="bottom",
        legend.margin = margin(1,44,0,0)
  ) +
  guides(color=guide_legend(nrow=6)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  scale_y_continuous(breaks = seq( 0, max(table_to_plot$value)+0.01, 0.01) ) +
  labs(x="Experiment day", y="ABS counts %",
       color="",
       title = "Genes most positely correlated with time\n according to Kendall correlation on ABS proportional counts")
ggsave(file="RNA_AGS_Results/Correlations/Genes_most_positively_correlated_with_time.png",width=6.2,height=5, dpi=300) 
dev.off()



# THE RNA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
# corr_top<-significative_abundances_neg[(length(significative_abundances_neg)-9):length(significative_abundances_neg)]  # 10 
corr_top<-significative_abundances_neg[1:11] # there is no 12th

table_to_plot <- as.data.frame(abundances[corr_top, ])
infos <- temp_for_save[ corr_top, ]
table_to_plot$Entry <- paste0( infos$Description," (", infos$Gene,")" )
# table_to_plot$Entry<- gsub("()", "(Unknown source)", table_to_plot$Entry, fixed=T)
# table_to_plot$Entry<-gsub("(Fragment)","",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub("Candidatus ","",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub("denitrificans","den.",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub(" bacterium","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("serine/threonine prot","ser/thr",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("(RNA-pol)-subunit","RNApol-sub",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("associated domain","domain",table_to_plot$Entry,fixed=T)
# table_to_plot$Entry<-gsub("Putative ","",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("membrane prot","membrane",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("synthase","synth",table_to_plot$Entry,fixed=T)
table_to_plot$Entry<-gsub("gamma-","γ-",table_to_plot$Entry,fixed=T)
table_to_plot$Entry[duplicated(table_to_plot$Entry)] <- paste( table_to_plot$Entry[duplicated(table_to_plot$Entry)] , row.names(table_to_plot)[duplicated(table_to_plot$Entry)] )

# colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "th day") ,"Entry" ) # already set the same order beforehand
colnames(table_to_plot)<- c( paste0(metadata$Experiment_day[metadata$Experiment_day!="other_AGS"], "") ,"Entry" ) # already set the same order beforehand
table_to_plot<-melt(table_to_plot, id.vars = "Entry")

ggplot(data=table_to_plot, aes(y=value, x=variable, color=Entry)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Entry), size =2) +
  geom_point(aes(color=Entry), size =3.8, alpha=0.4) +
  geom_line(aes(group=Entry), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_19[5:19]) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=8),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 8.6 ),
        legend.position="bottom",
        legend.margin = margin(1,42,0,0)) +
  guides(color=guide_legend(nrow=6)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  # scale_y_continuous(breaks = seq( 0, max(table_to_plot$value)+0.01, 0.01) , limits = c(0, 0.03) ) +
  labs(x="Experiment day", y="ABS counts %",
       color="",
       title = "Genes most negatively correlated with time\n according to Kendall correlation on ABS proportional counts")
ggsave(file="RNA_AGS_Results/Correlations/Genes_most_negatively_correlated_with_time.png",width=6.2,height=5, dpi=300) 
dev.off()


suppressWarnings(rm(prune.dat_top, tax_selected, tabella, corr_top))




######################## SEARCHING FOR THE KNOWN AGS RELATED GENES ################

# AGS related genes from https://doi.org/10.1016/j.scitotenv.2023.169103 and https://doi.org/10.1016/j.watres.2023.120700
only_AGSgenes_IDs<-dictionary [ grepl("^ppk", dictionary$Gene ) |    # poly-phosphate kinase
                                  grepl("^actP", dictionary$Gene ) |    # proton acetate symporter
                                  grepl("^pit", dictionary$Gene ) |     # inorganic phosphate transporter
                                  grepl("^pta", dictionary$Gene ) |    # acetate use related gene
                                  grepl("^ackA", dictionary$Gene ) |   # acetate use related gene
                                  grepl("^ppa", dictionary$Gene ) |    # pyrophosphatase
                                  grepl("^glg", dictionary$Gene ) |    # glycogenin glucosyltransferase; EC 2.4.1.186, used by GAO
                                  grepl("^bglX", dictionary$Gene ) |   # beta-glucosidase like enzymes, used by GAO
# also AGS attachments genes from https://doi.org/10.1016/j.jes.2023.07.011
                                  grepl("^rbmA", dictionary$Gene ) |
                                  grepl("^bap1", dictionary$Gene ) |
                                  grepl("^rbmC", dictionary$Gene ) |
                                  grepl("^galU", dictionary$Gene ) |
                                  grepl("^pGK", dictionary$Gene ) |
                                  grepl("^pGM", dictionary$Gene ) |
                                  grepl("^rpfF", dictionary$Gene ) |
                                  grepl("^rmlA", dictionary$Gene ) ,
                                ]
only_AGSgenes_IDs<-only_AGSgenes_IDs[ ! grepl("Copper-transporting", only_AGSgenes_IDs$Description) ,]  # homonym



### TRANSCRIPT LEVEL
table <- as.data.frame( prop_with_unmap[ only_AGSgenes_IDs$Entry, ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("alpha","α",infos,fixed=T)
infos<-gsub("Alpha","α",infos,fixed=T)
infos<-gsub("α-1,4","1,4-α",infos,fixed=T)
infos<-gsub("GlgB","",infos,fixed=T)
infos<-gsub("enzyme","enz",infos,fixed=T)
infos<-gsub("phosphate","phosph",infos,fixed=T)
infos<-gsub("Glycogen","Glycog",infos,fixed=T)
infos<-gsub("Candidatus","Ca.",infos,fixed=T)
infos<-gsub("Competibacteraceae bacterium","Competibacteraceae",infos,fixed=T)
infos<-gsub("Cation acetate","Cation/acetate",infos,fixed=T)
infos<-gsub("Cation/acetate symporter actP acetate and glycolate permease","Cation/acetate symporter and permease",infos,fixed=T)
infos<-gsub("denitrificans Run_A_D11","denitr.",infos,fixed=T)
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-16,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % (excluding other RNAs)", 
       fill="",
       title = "Most abundant mRNA of AGS related genes (transcript level)",
       caption = "NB: those are ONLY the most abundant ones associated to AGS related genes ")
ggsave("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_AGS_related_RNA_genes_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)

# EXPORTING ALSO *EVERY* ABS COUNT WITH INFO
ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Pathway<-NULL
ordered_infos2<-ordered_infos # used below
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(only_them), ]
only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them, MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[colnames(only_them)!="other_AGS"], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_AGS_related_RNA_genes_abs_counts_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)

# extra: exporting every detected function with ref names, avoiding duplications
ordered_infos2<-ordered_infos2[ordered_infos2$Entry %in% only_AGSgenes_IDs$Entry, ]
ordered_infos2<-ordered_infos2[ !duplicated(paste(ordered_infos2$Taxon , ordered_infos2$Gene) ) , ]
write.csv(file="RNA_AGS_Results/Every_AGS_trascript_found.csv", ordered_infos2[order(ordered_infos2$Taxon), ] , quote = F, row.names = F)

suppressWarnings(rm(table, table2, ordered_infos, ordered_infos2, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))



### GENE LEVEL
table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% only_AGSgenes_IDs$Gene, ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

top <- names(sort(rowSums(table[row.names(table)!="" , ]), decreasing=TRUE))[1:30]   # "" = unmapped
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
# Converting IDs
selected_infos<-dictionary[dictionary$Gene %in% table_top$Protein, ]
selected_infos<-selected_infos[! duplicated(selected_infos$Gene), ]
row.names(selected_infos)<-selected_infos$Gene
selected_infos<-selected_infos[table_top$Protein, ]
selected_infos$Description[grepl("phage",selected_infos$Taxon)]<-paste(selected_infos$Description[grepl("phage",selected_infos$Taxon)], "(virus)")
infos<-paste0( selected_infos$Description," (", selected_infos$Gene,")" )
infos<- gsub("()", "(no gene name)", infos, fixed=T)
infos<-gsub("alpha","α",infos,fixed=T)
infos<-gsub("Alpha","α",infos,fixed=T)
infos<-gsub("α-1,4","1,4-α",infos,fixed=T)
infos<-gsub("enzyme","enz",infos,fixed=T)
infos<-gsub("phosphate","phosph",infos,fixed=T)
infos<-gsub("Glycogen","Glycog",infos,fixed=T)
infos<-gsub("Cation acetate","Cation/acetate",infos,fixed=T)

table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-15,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % (only displayed genes)", 
       fill="",
       title = "Most abundant genes of AGS related genes (gene level)",
       caption = "NB: those are ONLY the most abundant ones associated to AGS related genes ")
ggsave("RNA_AGS_Results/Most_abundant_GENES/Most_abundant_GENES_already_known_in_AGS_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)

#exporting gene infos...
ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them, MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[colnames(only_them)!="other_AGS"], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_GENES/Most_abundant_GENES_already_known_in_AGS_abs_counts_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




######################## SEARCHING FOR THE N OXIDATION RELATED GENES ################

# from PMID: 19830422 ; DOI:10.1016/j.scitotenv.2022.157513
only_Ngenes_IDs<-dictionary [ grepl("^amo", dictionary$Gene ) |   # ammonia monooxygenase
                           grepl("^hao", dictionary$Gene ) | # hydroxylamine oxidoreductase
                           grepl("^hzo", dictionary$Gene ) | # hydrazine oxidoreductase
                           grepl("^nar", dictionary$Gene ) |    # nitrate ox reductase
                           grepl("^n2or", dictionary$Gene ) |    # last step
                           grepl("^nos", dictionary$Gene ) |    # nitrous oxide reductase (synonimous of n2or)
                           grepl("^nir", dictionary$Gene ) |
                           grepl("^nor", dictionary$Gene ) | # nitric oxide reductase
                           grepl("nitrogen", dictionary$Description ) |
                           grepl("nitrate", dictionary$Description ) |
                           grepl("nitrite", dictionary$Description ) |
                           grepl("N2O", dictionary$Description ) |
                           grepl("ammonium", dictionary$Description ) ,
                         ]
only_Ngenes_IDs<-only_Ngenes_IDs[ ! grepl("drug", only_Ngenes_IDs$Description) ,]  # norM, which is not involved with N!



### TRANSCRIPT LEVEL
table <- as.data.frame( prop_with_unmap[ only_Ngenes_IDs$Entry, ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("alpha","α",infos,fixed=T)
infos<-gsub("Candidatus","Ca.",infos,fixed=T)
{infos<-gsub(" bacterium","",infos,fixed=T)
  infos<-gsub("uncultured)","uncultured bact)",infos,fixed=T)
}
infos<-gsub("reductase","reduct",infos,fixed=T)
infos<-gsub("sp.","",infos,fixed=T)
infos<-gsub("prot","",infos,fixed=T)
infos<-gsub("subunit","sub",infos,fixed=T)
infos<-gsub("biosynthesis","biosynth",infos,fixed=T)
infos<-gsub("Partial","",infos,fixed=T)
infos<-gsub("hydroxylamine","NH2OH",infos,fixed=T)
infos<-gsub("nitric oxide","NO",infos,fixed=T)
infos<-gsub("Nitric oxide","NO",infos,fixed=T)
infos<-gsub("^ ","",infos)
infos<-gsub("Nitrate","NO3-",infos,fixed=T)
infos<-gsub("nitrate","NO3-",infos,fixed=T)
infos<-gsub("nitrite","NO2-",infos,fixed=T)
infos<-gsub("Nitrite","NO2-",infos,fixed=T)
infos<-gsub("NO3-:","NO3-",infos,fixed=T)
infos<-gsub("NO2-:","NO2-",infos,fixed=T)
infos<-gsub("Ammonia","NH3",infos,fixed=T)
infos<-gsub("cytochrome","",infos,fixed=T)
infos<-gsub("assembly","",infos,fixed=T)
infos<-gsub("Putative ","",infos,fixed=T)
infos<-gsub("denitrificans Run_A_D11","denitr. Run_A_D11",infos,fixed=T)

table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.26, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-16,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % (excluding other RNAs)", 
       fill="",
       title = "Most abundant mRNA of N related genes (transcript level)",
       caption = "NB: those are ONLY the most abundant ones associated to N related genes ")
ggsave("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_N_related_RNA_genes_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)

# EXPORTING ALSO *EVERY* ABS COUNT WITH INFO
ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Pathway<-NULL
ordered_infos2<-ordered_infos # used below
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(only_them), ]
only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them, MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[colnames(only_them)!="other_AGS"], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_RNAs/Most_abundant_N_related_RNA_genes_abs_counts_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)

# extra: exporting every detected function with ref names, avoiding duplications
ordered_infos2<-ordered_infos2[ordered_infos2$Entry %in% only_Ngenes_IDs$Entry, ]
ordered_infos2<-ordered_infos2[ !duplicated(paste(ordered_infos2$Taxon , ordered_infos2$Gene) ) , ]
write.csv2(file="RNA_AGS_Results/Every_N_trascript_found.csv", ordered_infos2[order(ordered_infos2$Taxon), ] , quote = F, row.names = F)

suppressWarnings(rm(table, table2, ordered_infos, ordered_infos2, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))



### GENE LEVEL
table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% only_Ngenes_IDs$Gene[only_Ngenes_IDs$Gene!=""], ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

# length(unique(only_Ngenes_IDs$Gene))
top <- names(sort(rowSums(table[row.names(table)!="" , ]), decreasing=TRUE))[1:30]   # "" = unmapped
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
# Converting IDs
selected_infos<-dictionary[dictionary$Gene %in% table_top$Protein, ]
selected_infos<-selected_infos[! duplicated(selected_infos$Gene), ]
row.names(selected_infos)<-selected_infos$Gene
selected_infos<-selected_infos[table_top$Protein, ]
selected_infos$Description[grepl("phage",selected_infos$Taxon)]<-paste(selected_infos$Description[grepl("phage",selected_infos$Taxon)], "(virus)")
infos<-paste0( selected_infos$Description," (", selected_infos$Gene,")" )
infos<- gsub("()", "(no gene name)", infos, fixed=T)
infos<-gsub("nitric oxide","NO",infos,fixed=T)
infos<-gsub("Nitric oxide","NO",infos,fixed=T)
infos<-gsub("^ ","",infos)
infos<-gsub("subunit c552","c552",infos,fixed=T)
infos<-gsub("assembly","",infos,fixed=T)
infos<-gsub("synthetase","synthet",infos,fixed=T)
infos<-gsub("reductase","reduct",infos,fixed=T)
infos<-gsub("cytochrome","cytochr",infos,fixed=T)
infos<-gsub("nitrogen","N",infos,fixed=T)
infos<-gsub("Nitrate","NO3-",infos,fixed=T)
infos<-gsub("nitrate","NO3-",infos,fixed=T)
infos<-gsub("nitrite","NO2-",infos,fixed=T)
infos<-gsub("Nitrite","NO2-",infos,fixed=T)
infos<-gsub("NO3-:","NO3-",infos,fixed=T)
infos<-gsub("NO2-:","NO2-",infos,fixed=T)
infos<-gsub("Ammonia","NH3",infos,fixed=T)
infos<-gsub("Copper-containing NO3- reduct","Copper-containing NO2- reduct",infos,fixed=T)  # it was wrong, because it is related to nirK
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.14, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 8.2 ),
        legend.margin = margin(-10,0,8,-32),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % (only displayed genes)", 
       fill="",
       title = "Most abundant genes of N related genes (gene level)",
       caption = "NB: those are ONLY the most abundant ones associated to AGS related genes ")
ggsave("RNA_AGS_Results/Most_abundant_GENES/Most_abundant_GENES_related_to_N_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)

#exporting gene infos...
ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_EVERY_SAMPLE<-apply( only_them, MARGIN = 1, mean)
AVERAGES_WITHOUT_OTHER_AGS<-apply( only_them[colnames(only_them)!="other_AGS"], MARGIN = 1, mean)
to_save <- cbind.data.frame( only_them[ , !grepl("^LR", colnames(only_them) ), drop=F] , AVERAGES_EVERY_SAMPLE, AVERAGES_WITHOUT_OTHER_AGS, ordered_infos )
write.csv2(file=paste0("RNA_AGS_Results/Most_abundant_GENES/Most_abundant_GENES_related_to_N_abs_counts_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




############### EXTRA: BAR PLOTS OF THE PHAGE ABUNDANT RNAs ##################

# table<-table_abs
table <- as.data.frame(prop_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Experiment_day,
                         "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 
# subsetting only phages
phages<-dictionary [ grepl( "phage" , dictionary$Taxon ) |  grepl( "virus" , dictionary$Taxon ) , "Entry"]

top <- names(sort(rowSums(table[row.names(table) %in% phages , ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<-gsub("Picornavirus capsid domain-containing prot","Protein with Picornavirus capsid domain",infos,fixed=T)
infos<-gsub("Trichosanthes kirilowii picorna-like virus","Trichosanthes kirilowii",infos,fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-15,0,3,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %",
       fill="",
       title = "Phages RNA in AGS mixed liquor",
       caption = "NB: different depth among samples --> more 'others' RNA if more depth!")
ggsave("RNA_AGS_Results/About_phages.png", width= 7, height = 5, dpi=300)




################### EXTRA: CHECKING THE EUKARYOTIC "CONTAMINATION" ########################

table_euka<-read.table("FASTQ_check/Every_eukar_read_found_by_kraken2_or_Diamond.tsv", row.names = 1, header = T)
colnames(table_euka)<-gsub("^X","",colnames(table_euka))
table_euka <- table_euka[ , metadata$FASTQ_ID] # this re-orders
colnames(table_euka)<-metadata$Sample_name


euka_prop <- apply(table_euka, MARGIN = 2, function(x) x/sum(x))
euka_prop <- euka_prop * 100
table <- as.data.frame(euka_prop) # prop abs
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table))

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", ]), decreasing=TRUE))[1:19]
table_top<- table[top, ]
table_others<- t(as.data.frame(colSums(table[!row.names(table) %in% top, ])))
row.names(table_others)<-"Others"
# plotting ...
table_to_plot<-rbind.data.frame(table_top, table_others)
table_to_plot$Taxon<-row.names(table_to_plot)
table_to_plot<-melt(table_to_plot, id.vars = "Taxon")
table_to_plot$Taxon<-factor(table_to_plot$Taxon, levels = c(unique(table_to_plot$Taxon)[! unique(table_to_plot$Taxon) %in% "Others"],"Others"))
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Taxon)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_20) +
  scale_y_continuous(expand=c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 8.5 ),
        legend.margin = margin(-10,0,3,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="absolute count %",
       fill="",
       title = "Most abundant eukaryotic RNA (filtered)")
ggsave("FASTQ_check/Most_abundant_Eukar_contams_by_Kraken2_or_Diamond.png", width= 7, height = 5, dpi=300)




################### EXTRA: CHECKING THE ncRNA besides the ribosomial #########################
# NB: this part of the script works only if the optional Salmon step has been performed during the processing

# importing ncRNA
files_to_import<-dir("FASTQ_check/Salmon_ncRNAmapped_no_rRNA/")
files_to_import<- files_to_import[grep("ncRNA_found", files_to_import)]
if(length(files_to_import)==0){stop("\n\n Salmon files not found in 'FASTQ_check' folder! \nHave you performed the (optional) Salmon step? \n\n")}

files_list<-list()
for( i in 1:length(files_to_import) ) {
  files_list[[i]]<- as.data.frame(read.table( paste0("FASTQ_check/Salmon_ncRNAmapped_no_rRNA/",files_to_import[i]) , sep="\t", header = T, dec = ".", quote="") )
} 
names(files_list)<-gsub("ncRNA_found_in_", "", files_to_import)
names(files_list)<-gsub(".tsv", "", names(files_list), fixed=T)
updating_names<-metadata
row.names(updating_names)<-metadata$FASTQ_ID # required below to ensure the same order between these file and the metadata
names(files_list)<-paste0(updating_names$Experiment_day, "th day")
names(files_list)<-gsub("AGSth day", "AGS", names(files_list)) 

# joining every report to a temporary unique table (to search for unique ncRNA in every assignment and then creating the related row in the final object)
complete_table<-NULL
for(i in 1:length(files_list)){
  new_sample<-cbind.data.frame(sample=names(files_list)[i], files_list[[i]]) # sample name (repeated along the column) with the related results
  complete_table<-rbind.data.frame(complete_table, new_sample)
}
colnames(complete_table)[colnames(complete_table)=="NumReads"]<-"Abundance"  # selecting the estimation by bracken as abundance

# from list (each sample is repeated along a colum) to feature table (each sample is a column)
feature_table<-cbind("Temp_column"=unique(complete_table$Name))
for( s in 1:length(unique(complete_table$sample) )) {  # NB: the numeric index is required to deal with the column name afterward!
  sample<- unique(complete_table$sample)[s] # 1 sample at time
  table_4_that_sample<- complete_table[complete_table$sample==sample , ]
  if( length(which(duplicated(table_4_that_sample$Name)))>0 ) {  # sometime a different RNA (e.g. a different piece) may be mapped on the same transcript ID
    table_4_that_sample$Length<-NULL
    table_4_that_sample$sample<-NULL
    table_4_that_sample<-aggregate( .~Name, data=table_4_that_sample, FUN=sum)
  }
  abundances<-NULL # resets the vector, ready for the next code
  for( b in unique(complete_table$Name)) {
    # scan each code in order: if found in that sample then appends the correspective abundances, if absent paste a 0
    if( b %in% table_4_that_sample$Name ){
      abundances<-c(abundances, table_4_that_sample[ table_4_that_sample$Name==b, "Abundance"])
    } else { 
      abundances<-c(abundances, 0 ) # that bacterium is absent --> zero
    }
  }
  feature_table<-cbind.data.frame(feature_table, abundances)
  colnames(feature_table)[s+1] <- sample   # the +1 is for first column, which is temporarily the taxonomy column
}

Type<-sub(".*[0-9]_","",feature_table$Temp_column)
Type<-sub('"','',Type, fixed=T)
Type[Type==""]<-"no_info"
Code<-gsub("_.*","",feature_table$Temp_column)
feature_table$Temp_column<-NULL

feature_table<-feature_table[order(apply( feature_table, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_IN_THIS_REACTOR<-round(apply( feature_table[, !grepl("other",colnames(feature_table))], MARGIN = 1, mean),0)
to_save <- cbind.data.frame( feature_table , AVERAGES_IN_THIS_REACTOR, Type, Code )
write.csv2(file=paste0("FASTQ_check/Salmon_ncRNAmapped_no_rRNA/Feature_table_ncRNA.csv"), to_save, row.names=T, quote=F)




################ EXTRA: CHECKING THE rRNA REGULATORY GENES ###########################

rRNA_amount_genes<-dictionary [ grepl("^fis", dictionary$Gene ) |
                                  grepl("^hns", dictionary$Gene ) |  # gene of H-NS, see https://www.uniprot.org/uniprotkb/P0ACF8/entry
                                  grepl("^relX", dictionary$Gene ) | # gene involved in ppGpp amount regulation , see https://pubmed.ncbi.nlm.nih.gov/342913/
                                  grepl("^relA", dictionary$Gene ) | # gene involved in ppGpp amount regulation , see https://pubmed.ncbi.nlm.nih.gov/342913/
                                  grepl("ppGpp", dictionary$Description ) ,
                              ]

### GENE LEVEL
table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% rRNA_amount_genes$Gene, ] )
colnames(table)<- paste0(metadata$Experiment_day, "th day")
colnames(table)<-gsub("AGSth day", "AGS", colnames(table)) 

top <- names(sort(rowSums(table[row.names(table)!="" , ]), decreasing=TRUE))[1:30]   # "" = unmapped
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
# Converting IDs
selected_infos<-dictionary[dictionary$Gene %in% table_top$Protein, ]
selected_infos<-selected_infos[! duplicated(selected_infos$Gene), ]
row.names(selected_infos)<-selected_infos$Gene
selected_infos<-selected_infos[table_top$Protein, ]
selected_infos$Description[grepl("phage",selected_infos$Taxon)]<-paste(selected_infos$Description[grepl("phage",selected_infos$Taxon)], "(virus)")
infos<-paste0( selected_infos$Description," (", selected_infos$Gene,")" )
infos<- gsub("()", "(no gene name)", infos, fixed=T)
infos<- gsub("Candidatus", "Ca.", infos, fixed=T)
infos<- gsub("guanosine-3',5'-bis(Diphosphate) 3'-pyrophosphohydrolase", "diphosphate", infos, fixed=T)

table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) + 
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=7),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5), 
        axis.ticks.y=element_blank(), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.35, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-15,0,3,-35),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count % (excluding other RNAs)", 
       fill="",
       title = "Most abundant genes of Fis, H-NS and ppGpp (gene level)",
       caption = "NB: those are ONLY the most abundant ones associated rRNA transcr regulation")
ggsave("RNA_AGS_Results/About_Fis_HNS_ppGpp_PROPORTIONS_ON_THEMSELVES.png", width= 7, height = 5, dpi=300)




##################### \\\\\\ R AND PACKAGES VERSION , ANALYSIS NOTES \\\\\\ #########################


### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("RNA_AGS_Results/R_version_and_packages.txt")
sink(con, append = TRUE)

cat(package$R.version$version.string)
cat("   running on", package$running)
cat("\n", "\n", fill=TRUE)
package$otherPkgs$ggplot2[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$vegan[c(1,3)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$DESeq2[1:2]
cat("\n", "\n", fill=TRUE)
cat("\n", "\n", fill=TRUE)
print("vegan")
packageVersion("vegan")
cat("\n", "\n", fill=TRUE)
print("ecodist")
packageVersion("ecodist")
cat("\n", "\n", fill=TRUE)
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))


# Personal notes
con <- file("RNA_AGS_Results/Important_notes_about_the_analyses.txt")
sink(con, append = TRUE)
cat("* Every abundances filter has been computed including the unmapped microbial reads (filtered during the processing but not identified) in the dataset, using the raw abs table in CPM, to better identified which observation was a just a noise\n", fill=T)
cat("* The ordinations and the correlations have been computed on proportions not including the unmapped microbial reads, to avoid including into the analyses a 'missing information' represented by cumulative values (the sum of the unmapped read)' \n", fill=T)
cat("* Many barplots report the relative abundances of the displayed RNA only! Indeed, including the 'others' RNA or the unmapped count, would make impossible to see any color in the barplots!\n", fill=T)
cat("* The p values of the displayed correlations are not adjusted (no significant correlation after the adjustement due to the copious tests!) \n", fill=T)
sink()
close(con)
suppressWarnings(rm(con))

# Bye!