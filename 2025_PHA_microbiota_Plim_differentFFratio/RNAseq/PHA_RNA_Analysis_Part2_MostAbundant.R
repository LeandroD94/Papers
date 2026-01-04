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

suppressWarnings( dir.create("RNA_PHA_Results") )
dir.create("RNA_PHA_Results/Most_abundant_RNAs")
dir.create("RNA_PHA_Results/Most_abundant_GENES")


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
                     "darkblue","grey60", "darkolivegreen",
                     "slateblue3","cyan","yellow3","brown4",
                     "springgreen2", "grey15", 
                     "deepskyblue2","darkgreen","orange3",
                     "darkorchid", "lightblue1" , "violet", "blue",
                     "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                     "greenyellow", "red","darkslategray3")


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




########### BAR PLOTS OF THE MOST ABUNDANT RNA (ABS_count) ###############
# COMPUTATION ON THE WHOLE DATASET (BOTH R1 AND R2 SAMPLES)

# table<-table_abs
table <- as.data.frame(prop_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
selected_infos$Description[selected_infos$Description=="Uncharact"]<-paste( selected_infos$Description[selected_infos$Description=="Uncharact"], selected_infos$Entry[selected_infos$Description=="Uncharact"] )
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<- gsub("methanol/ethanol fam", "meth/eth fam", infos, fixed=T)
infos<- gsub("prot-sorting", "sorting", infos, fixed=T)
infos<- gsub("Uncharact ", "Uncharacterized ", infos, fixed=T)
infos<- gsub("Candidatus","Ca.",infos,fixed=T)
infos<- gsub("Ca_","Ca.",infos,fixed=T)
infos<- gsub("Phasin fam prot","Phasin fam protein",infos,fixed=T)
infos<- gsub("^phasin","Phasin",infos)
infos<- gsub(" fam "," family ",infos)
infos<- gsub("sp. 27","sp.27",infos, fixed=T)
infos<- gsub(" prot "," ",infos, fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =10) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,2,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %",
       fill="",
       title = "Most abundant mRNA in the entire dataset",
       caption = "NB: different depth among samples --> if more depth then more 'others' RNA!")
ggsave("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)

# write.csv2(file="RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_WHOLE_DATASET_infos.csv", table_top, row.names=T, quote=F)



########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30_BIS) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.26, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-40),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       title = "Most abundant mRNA in the entire dataset")
ggsave("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

# write.csv2(file="RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs__PROPORTIONS_ON_THEMSELVES_infos.csv", table_top, row.names=T, quote=F)




########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<-round( apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean) ,2)
AVERAGES_R2<-round( apply( table_top[, !grepl("R1",colnames(table_top))], MARGIN = 1, mean) ,2)
to_save <- cbind.data.frame( round(table_top,2) , AVERAGES_R1, AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Table_Most_Abundant_RNA_ABSpropAverages_WHOLE_DATASET_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<-round( apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean) , 2)
AVERAGES_R2<-round( apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean) ,2 )
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R1, AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Table_Most_Abundant_RNA_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




########### BAR PLOTS OF THE MOST ABUNDANT RNA ***IN R1*** (ABS_count) ###############

# table<-table_abs
table <- as.data.frame(prop_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", grepl("R1",colnames(table)) ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top <- table_top[ , grepl("R1",colnames(table_top))  ]  # Selecting R1 samples only

table_top$Protein<-row.names(table_top)

### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<- gsub("subunit", "", infos, fixed=T)
infos<- gsub("Uncharact", "Uncharacterized", infos, fixed=T)
infos<- gsub("PEP-CTERM prot-sorting domain", "PEP-CTERM sorting domain", infos, fixed=T)
infos<- gsub("methanol/ethanol", "meth/eth", infos, fixed=T)
infos<- gsub("Candidatus","Ca.",infos,fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =10) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.25, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,3,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="absolute count %",
       fill="",
       title = "Most abundant mRNA in R1 samples",
       caption = "NB: different depth among samples --> if more depth then more 'others' RNA!")
# ggsave("RNA_PHA_Results/Most_abundant_RNAs/R1_Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)
# about 4% on the whole dataset!

# write.csv2(file="RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_WHOLE_DATASET_infos.csv", table_top, row.names=T, quote=F)



########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.315, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,2,-44),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       title = "Most abundant mRNA in R1 samples")
ggsave("RNA_PHA_Results/Most_abundant_RNAs/R1_Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

# write.csv2(file="RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs__PROPORTIONS_ON_THEMSELVES_infos.csv", table_top, row.names=T, quote=F)




########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<-round(apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean),2)
# AVERAGES_R2<-apply( table_top[, !grepl("R1",colnames(table_top))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(table_top,2) , AVERAGES_R1, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/R1_Table_Most_Abundant_RNA_ABSpropAverages_WHOLE_DATASET_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<-round(apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean),2)
# AVERAGES_R2<-apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R1, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/R1_WHOLE_DATASET_Table_Most_Abundant_RNA_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

to_save_R1_transcr<-to_save   # to be use after

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




########### BAR PLOTS OF THE MOST ABUNDANT RNA ***IN R2*** (ABS_count) ###############

# table<-table_abs
table <- as.data.frame(prop_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", grepl("R2",colnames(table)) ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top <- table_top[ , grepl("R2",colnames(table_top))  ]  # Selecting R2 samples only

table_top$Protein<-row.names(table_top)

### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<- gsub("()", "(Unknown source)", infos, fixed=T)
infos<- gsub("sigma-", "σ", infos, fixed=T)
infos<- gsub("fam outer membrane", "outer membr", infos, fixed=T)
infos<- gsub("prot-sorting domain", "sorting domain", infos, fixed=T)
infos<- gsub("sorting domain", "sort domain", infos, fixed=T)
infos<- gsub("transcriptional regulator", "transcr regul", infos, fixed=T)
infos<- gsub("Insertion element IS402-like domain-containing prot", "Insertion IS402 containing prot", infos, fixed=T)
infos<- gsub("Candidatus","Ca.",infos,fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =10) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,2,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %",
       fill="",
       title = "Most abundant mRNA in R2 samples",
       caption = "NB: different depth among samples --> if more depth then more 'others' RNA!")
# ggsave("RNA_PHA_Results/Most_abundant_RNAs/R2_Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)
# about 2% on the whole dataset

# write.csv2(file="RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_WHOLE_DATASET_infos.csv", table_top, row.names=T, quote=F)



########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text = element_text(size=10),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-40),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       title = "Most abundant mRNA in R2 samples")
ggsave("RNA_PHA_Results/Most_abundant_RNAs/R2_Most_abundant_RNA_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

# write.csv2(file="RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_Most_abundant_RNA_computed_with_abs__PROPORTIONS_ON_THEMSELVES_infos.csv", table_top, row.names=T, quote=F)




########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R2<-round(apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean),2)
# AVERAGES_R2<-apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(table_top,2), AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/R2_Table_Most_Abundant_RNA_ABSpropAverages_WHOLE_DATASET_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R2<-round(apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean),2)
# AVERAGES_R2<-apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/R2_WHOLE_DATASET_Table_Most_Abundant_RNA_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

to_save_R2_transcr<- to_save # to be used after

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




########### BAR PLOTS OF THE MOST ABUNDANT GENES (ABS_count) ###############
# COMPUTATION ON THE WHOLE DATASET (BOTH R1 AND R2 SAMPLES)

table <- as.data.frame(prop_GENE_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

top <- names(sort(rowSums(table[!row.names(table)%in%c("","UNMAPPED"), ]), decreasing=TRUE))[1:30]   # "" = unmapped
# NB: "" removed, because they are not linked to any gene --> Impossible to merge with other transcripts
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
infos<- gsub("methanol/ethanol", "meth/eth", infos, fixed=T)
infos<- gsub("fam outer membrane", "outer membr", infos, fixed=T)
infos<- gsub("PEP-CTERM prot-sorting", "PEP-CTERM sorting", infos, fixed=T)
infos<- gsub("HPF/RaiA fam ", "HPF/RaiA ", infos, fixed=T)
infos<- gsub("Uncharact (", "Uncharacterized (", infos, fixed=T)
infos<- gsub(" fam "," family ",infos)
infos<- gsub(" prot "," ",infos, fixed=T)
infos<- gsub("Transmembrane (","Transmembrane protein (",infos, fixed=T)
infos<- gsub(" membr (","membr prot (",infos, fixed=T)
infos<- gsub(" dehydrog ("," dehydrogenase (",infos, fixed=T)
infos<- gsub("containing (","containing protein (",infos, fixed=T)
infos<- gsub("sorting domain-containing protein","sorting domain-containing prot",infos, fixed=T)
infos<- gsub(" attachment ("," attachment protein (",infos, fixed=T)
infos<- gsub(" lipoprot "," lipoprotein ",infos, fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =10) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,2,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %", 
       fill="",
       title = "Most abundant genes in the entire dataset",
       caption = "NB: different depth among samples --> if more depth then more 'others' genes!")
ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)

# write.csv2(file="RNA_PHA_Results/Most_abundant_GENES/Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_WHOLE_DATASET_infos.csv", table_top, row.names=T, quote=F)




########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30_BIS) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.4, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.29, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-40),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other genes)",
       fill="",
       title = "Most abundant genes in the entire dataset")
ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)



########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<-round(apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean) ,2)
AVERAGES_R2<-round(apply( table_top[, !grepl("R1",colnames(table_top))], MARGIN = 1, mean) ,2)
to_save <- cbind.data.frame( round(table_top,2), AVERAGES_R1, AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_Table_Most_Abundant_GENES_ABSpropAverages_WHOLE_DATASET_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<- round( apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean) ,2)
AVERAGES_R2<- round( apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean) ,2)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R1, AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_Table_Most_Abundant_GENES_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




########### BAR PLOTS OF THE MOST ABUNDANT GENES *** IN R1*** (ABS_count) ###############

table <- as.data.frame(prop_GENE_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")
table <- table[ , grepl("R1",colnames(table))]   # selecting R1

top <- names(sort(rowSums(table[!row.names(table)%in%c("","UNMAPPED"), ]), decreasing=TRUE))[1:30]   # "" = unmapped
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
infos<- gsub("iron-sulfur", "Fe-Sulph", infos, fixed=T)
infos<- gsub("prot ErpA", "prot", infos, fixed=T)
infos<- gsub("methanol/ethanol", "meth/eth", infos, fixed=T)
infos<- gsub("()", "(no gene name)", infos, fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =10) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,2,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %", 
       fill="",
       title = "Most abundant genes in R1 samples",
       caption = "NB: different depth among samples --> if more depth then more 'others' genes!")
# ggsave("RNA_PHA_Results/Most_abundant_GENES/R1_Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)

# write.csv2(file="RNA_PHA_Results/Most_abundant_GENES/Most_abundant_RNA_computed_with_abs_PROPORTIONS_ON_R1_infos.csv", table_top, row.names=T, quote=F)



########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-16),
        legend.position="bottom",
        plot.margin = margin(3,1,2,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other genes)",
       fill="",
       title = "Most abundant genes in R1 samples")
ggsave("RNA_PHA_Results/Most_abundant_GENES/R1_Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)



########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<- round( apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean) ,2)
#AVERAGES_R2<-apply( table_top[, !grepl("R1",colnames(table_top))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(table_top,2), AVERAGES_R1, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/R1_Table_Most_Abundant_GENES_ABSpropAverages_R1_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<-round( apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean) ,2)
# AVERAGES_R2<-apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R1, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/R1_Table_Most_Abundant_GENES_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

to_save_R1_gene<-to_save # to be used after

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




########### BAR PLOTS OF THE MOST ABUNDANT GENES *** IN R2*** (ABS_count) ###############

table <- as.data.frame(prop_GENE_with_unmap) # prop abs
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")
table <- table[ , grepl("R2",colnames(table))]   # selecting R2

top <- names(sort(rowSums(table[!row.names(table)%in%c("","UNMAPPED"), ]), decreasing=TRUE))[1:30]   # "" = unmapped
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
infos<- gsub("sigma-", "σ", infos, fixed=T)
infos<- gsub("fam outer membrane prot", "outer membr", infos, fixed=T)
infos<- gsub("transcriptional regulator", "transcr regul", infos, fixed=T)
infos<- gsub("transcript regulator", "transcr regul", infos, fixed=T)
infos<- gsub("σ54 specific Fis", "σ54-Fis", infos, fixed=T)
infos<- gsub("prot-sorting domain", "domain", infos, fixed=T)
infos<- gsub("sorting domain", "domain", infos, fixed=T)
infos<- gsub("Insertion element IS402-like", "Insertion IS402-like", infos, fixed=T)
table_top$Protein<-infos
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ Reactor + variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-12,0,2,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="absolute count %", 
       fill="",
       title = "Most abundant genes in R2 samples",
       caption = "NB: different depth among samples --> if more depth then more 'others' genes!")
# ggsave("RNA_PHA_Results/Most_abundant_GENES/R2_Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_WHOLE_DATASET.png", width= 7, height = 5, dpi=300)



########## AGAIN BUT PROPORTIONS RE-COMPUTED ON THEMSELVES ONLY
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ Reactor+variable,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.38, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-30),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other genes)",
       fill="",
       title = "Most abundant genes in R2 samples")
ggsave("RNA_PHA_Results/Most_abundant_GENES/R2_Most_abundant_GENES_computed_with_abs_counts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)



########## EXPORTING ALSO *MEANS* ABS COUNT WITH INFO

ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

table_top$Protein<-NULL
table_top <- table_top[order(apply( table_top, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R2<-round(apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean),2)
#AVERAGES_R2<-apply( table_top[, !grepl("R2",colnames(table_top))], MARGIN = 1, mean)
to_save <- cbind.data.frame( round(table_top,2), AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/R2_Table_Most_Abundant_GENES_ABSpropAverages_R2_with_infos.csv"), to_save, row.names=T, quote=F)

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R2<-round(apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean),2)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/R2_Table_Most_Abundant_GENES_ABSpropAverages_PROPORTIONS_ON_THEMSELVES_with_infos.csv"), to_save, row.names=T, quote=F)

to_save_R2_gene<-to_save # to be used after

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




############# INTERSECTION BETWEEN MOST ABUNDANT R1 AND R2 ##############

inters_transc<-to_save_R1_transcr[row.names(to_save_R1_transcr) %in% row.names(to_save_R2_transcr), c("Taxon","Description","AVERAGES_R1")]
temp <- to_save_R2_transcr[row.names(inters_transc), "AVERAGES_R2" ]
inters_transc<-cbind.data.frame(inters_transc, "AVERAGES_R2"=temp)

inters_gene<-to_save_R1_gene[row.names(to_save_R1_gene) %in% row.names(to_save_R2_gene), c("Gene","Description","AVERAGES_R1")]
temp <- to_save_R2_gene[row.names(inters_gene), "AVERAGES_R2" ]
inters_gene<-cbind.data.frame(inters_gene, "AVERAGES_R2"=temp)

# exporting
con <- file("RNA_PHA_Results/Most_abund_intersections_R1_R2.txt")
sink(con)
cat("The following are the observations which figure among the most abundant ones in both R1 and R2 dedicated computations",fill=T)
cat("\nTRANSCRIPTS WHICH ARE THE MOST ABUNDANT IN BOTH R1 AND R2 SAMPLE", fill=T)
print(inters_transc)
cat("\n\n\nGENES WHICH ARE THE MOST ABUNDANT IN BOTH R1 AND R2 SAMPLE", fill=T)
inters_gene
sink()
rm(con)




################# CHECKING TAXONOMIES (MAIN ACTORS) ##################

table <- table_abs
table$Taxa <- row.names(table)  # these are transcript IDs, but will be converted in taxa...
Taxa_ordered <- dictionary
row.names(Taxa_ordered) <- Taxa_ordered$Entry
table$Taxa <- Taxa_ordered[table$Taxa, "Taxon"]
table <- table[!grepl("RepID=", table$Taxa) , ]
table$Taxa <- gsub("'", "", table$Taxa, fixed=T)
table$Taxa <- gsub("[", "", table$Taxa, fixed=T)
table$Taxa <- gsub("]", "", table$Taxa, fixed=T)
table$Taxa <- gsub(" bact ", "", table$Taxa, fixed=T)
table$Taxa <- gsub(" bact$", "", table$Taxa)
table$Taxa <- gsub("uncult ", "uncultured ", table$Taxa)
# to genus level ...
table$Taxa <- gsub("uncultured ", "", table$Taxa)
table$Taxa <- gsub("unclassified ", "", table$Taxa)
table$Taxa <- gsub("unidentified  ", "", table$Taxa)
table$Taxa <- gsub("Candidatus ", "Candidatus_", table$Taxa)
table$Taxa <- gsub("candidate ", "candidate_", table$Taxa)
table$Taxa <- gsub("alpha ", "alpha_", table$Taxa)
table$Taxa <- gsub("beta ", "beta_", table$Taxa)
table$Taxa <- gsub("Alphaproteobacteria.*", "Alphaproteobacteria", table$Taxa)
table$Taxa <- gsub("Betaproteobacteria.*", "Betaproteobacteria", table$Taxa)
table$Taxa <- gsub("Bacteroidetes.*", "Bacteroidetes", table$Taxa)
table$Taxa <- gsub("Paracoccaceae.*", "Paracoccus", table$Taxa)
table$Taxa <- gsub( " .*", "" , table$Taxa )

table <- melt(table, id = "Taxa")
table <- aggregate( value ~ Taxa+variable, table, FUN=sum)
table <- reshape::cast(table, variable="Taxa" )
# row.names(table) <- table$Taxa

colnames(table)<- c("Taxa" , paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day") )

# to proportions (I can't use the prop table directly, don't know why ... )
table[ , !colnames(table)=="Taxa"] <- apply(table[ , !colnames(table)=="Taxa"], MARGIN = 2, function(x) x/sum(x))
table[ , !colnames(table)=="Taxa"] <- as.data.frame(table[ , !colnames(table)=="Taxa"]) * 100

top <-  order( rowSums( table[ , !colnames(table)%in%"Taxa"] ) , decreasing=TRUE)
table_top <- table[ top , ]
table_top$Taxa[ 20:length(table_top$Taxa) ] <- "Others"
# plotting ...
table_to_plot<-table_top
table_to_plot<-melt(table_to_plot, id.vars = "Taxa")
table_to_plot <- aggregate( value ~ variable + Taxa, table_to_plot, FUN=sum )
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
table_to_plot$Taxa <- factor( table_to_plot$Taxa , levels = 
                                c( rev(unique(table_top$Taxa[table_top$Taxa!="Others"])) , "Others" )
                              )
levels(table_to_plot$Taxa) <- gsub("Ca_", "Ca. ", levels(table_to_plot$Taxa))
levels(table_to_plot$Taxa) <- gsub("Accumulibacter", "Accumulibacter (*)", levels(table_to_plot$Taxa))

fill_color_taxa<-c("darkmagenta", "yellow3","wheat","indianred4","slateblue", "black", "chocolate","brown4","springgreen2",
                 "blue2","violet","wheat3","darkgreen",
                 "yellow", "red", "orange",
                  "navyblue", "green2","deepskyblue3","darkslategray3")
# bar plot
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Taxa)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =8.5) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values= fill_color_taxa ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.48, hjust=0.5, size=5.5),
        axis.ticks.y=element_blank(),
        strip.text.x=element_text(size=7.5,colour="black", 
                                  lineheight = 1,
                                  margin = margin(1.5,3,1.5,3, "pt")
        ),
        strip.background = element_rect(color = "black", linewidth = 0.45),        
        legend.key.height = unit(0.32, "cm"),
        legend.key.width = unit(0.24, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 6.2 ,face="italic" ),
        legend.margin = margin(-16,0,2,-28),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percent abundance",
       fill="")
ggsave("RNA_PHA_Results/Top_actors_of_the_theatre.png", width= 4, height = 3, dpi=300)

# NB: having overwritten few transcripts (see part 1) did not changed the plot beside showing Ca. accumulibacter instead of Azoarcus

