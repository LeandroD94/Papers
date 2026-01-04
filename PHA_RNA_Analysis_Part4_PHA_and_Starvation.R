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




################ MOST ABUNDANT PHA RELATED RNA and GENES ################

### TRANSCRIPT LEVEL
table <- as.data.frame( prop_with_unmap[ PHA_IDs , ] ) 
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

backup_numbers <- table  # to count them afterwards

top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
# Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
{infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
  infos<- gsub("()", "(Unknown source)", infos, fixed=T)
  infos<- gsub("Beta-", "β-", infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("poly(R)-hydroxyalkanoic", "poly(R)-hydroxyal.", infos, fixed=T)
  infos<- gsub("SWub3 = DSM 12120)", "SWub3)", infos, fixed=T)
  infos<- gsub("multifunctional regulatory prot", "regulatory prot", infos)
  infos<- gsub("unclassified", "unclassif", infos)
  infos<-gsub("Thauera linaloolentis.*","Thauera linaloolentis)",infos)
  infos<-gsub("Acetoacetyl-CoA reduct","Acetoacetyl-CoA reductase",infos)
  infos<-gsub("Acetoacetyl-CoA acetyltransf ","Acetoacetyl-CoA acetyltransferase ",infos)
  infos<-gsub("Acetyl-CoA acetyltransf ","Acetyl-CoA acetyltransferase ",infos)
  infos<-gsub("hydroxyal. acid synth ","hydroxyal. synth " ,infos , fixed=T)
  infos<-gsub("PHA depol ","PHA depolimerase ",infos)
  infos<- gsub("^phasin","Phasin",infos)
  infos<- gsub("Phasin fam prot","Phasin family",infos)
  infos<- gsub("Phasin prot","Phasin",infos)
  infos<- gsub("prot (PHASIN)","prot, phasin",infos, fixed=T)
  infos<- gsub("uncult bact","uncultured bact",infos)
  infos<- gsub("phbB PhbB2","phbB",infos)
  infos<- gsub("acetoacetyl-CoA reduct (Paracoccus ","Acetoacetyl-CoA reductase (Paracoccus ",infos , fixed=T)
  infos<- gsub("3-hydroxyacyl-CoA dehydrog NAD-bind domain-containing prot ","3-hydroxyacyl-CoA dehydrog ",infos , fixed=T)
  infos<- gsub("3-hydroxyacyl-CoA dehydrog","3-hydroxyacyl-CoA dehydrogenase",infos , fixed=T)
  infos<- gsub("synth fam prot ","synthase fam prot ",infos , fixed=T)
  infos<- gsub("dehydrog (","dehydrogenase (",infos , fixed=T)
  infos<- gsub(", class I"," class I",infos , fixed=T)
  infos<- gsub("Class I poly(R)-hydroxyal. synth","Poly(R)-hydroxyalkanoic acid synthetase class I",infos , fixed=T)
  infos<- gsub("Poly(R)-hydroxyalkanoic acid","PHA",infos , fixed=T)
  infos<- gsub("synth ","synthetase ",infos , fixed=T)
  infos<- gsub("granule regulatory prot","phasin",infos , fixed=T)
  infos<- gsub("(Bacteria)","(unclassified bacteria)",infos , fixed=T)
  infos<- gsub("coenzyme A","coA",infos , fixed=T)
  infos<- gsub("Acyl-coA dehydrogenase","Dehydrogenase of acyl-coenzyme A",infos , fixed=T)
}
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot$Protein <- gsub("(Thauera aminoaromatica) A0A418ZR87","(T. aminoaromatica) A0A418ZR87", table_to_plot$Protein , fixed = T )
table_to_plot$Protein <- gsub("PhaM fam PHA phasin","PhaM (PHB synthase activator)", table_to_plot$Protein , fixed = T )  #see doi: 10.1128/AEM.02935-13
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
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
        legend.key.width = unit(0.295, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.7 ),
        legend.margin = margin(-14,0,2,-39),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       # title = "Most abundant mRNA of PHA related genes (transcript level)"
  )
ggsave("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_PHA_related_Most_abundant_transcr_AbsCounts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

# AGAIN, BUT WITH CUSTOM COLORS FOR EACH FUNCTION(e.g., see table 1 of DOI 10.1080/21655979.2025.2458363 )
colors_type <- c("darkgreen" , "chartreuse4","chartreuse3", "darkseagreen4","green",
                 # "green2","darkolivegreen3", "chartreuse",
                 "orangered3","firebrick1","red3","coral2",
                 "yellow","darkmagenta",
                 #"blue3",
                 "gold","goldenrod2","orange","gold2","yellow2",
                 "cyan2","dodgerblue3",
                 "deepskyblue","aquamarine",
                 "blue3","dodgerblue4",
                 "royalblue1","cyan4", "turquoise","deepskyblue3","dodgerblue","lightskyblue",
                 "goldenrod1",
                 "darkviolet")
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
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
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.2),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.295, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.7 ),
        legend.margin = margin(-14,0,2,-39),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       title = "Most abundant mRNA of PHA related genes (transcript level)"
  )
ggsave("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_PHA_related_Most_abundant_transcr_AbsCounts_PROPORTIONS_ON_THEMSELVES2.png", width= 6, height = 5, dpi=300)

colSums(backup_numbers[ , !colnames(backup_numbers)%in%"Protein" ])  # proportions on whole dataset
# R1_63th day R2_63th day R1_81th day R2_81th day R1_95th day R2_95th day 
#  0.8127034   0.3068876   0.3405863   0.2233728   0.4273481   0.2517308 


# EXPORTING ALSO *EVERY* ABS COUNT WITH INFO
ordered_infos<-dictionary
row.names(ordered_infos)<-ordered_infos$Entry
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(only_them), ]
ordered_infos$Pathway<-NULL
only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<- round(apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean),2)
AVERAGES_R2<- round(apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean),2)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R1, AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_RNAs/Most_abundant_PHA_related_RNA_genes_abs_counts_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))



### GENE LEVEL
these_genes <- dictionary[ dictionary$Entry %in% PHA_IDs , "Gene" ]
these_genes <- these_genes[these_genes!=""]  # PHA related *BUT* not linked to any gene --> Impossible to merge to other transcripts
these_genes <- these_genes[!these_genes%in%c("thiG","ychF","cat","ilvD","ccrB")] # not related
table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% these_genes , ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

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
{ infos<- gsub("()", "(no gene name)", infos, fixed=T)
  infos<-gsub("alpha","α",infos,fixed=T)
  infos<-gsub("Phasin fam prot","Phasin family",infos,fixed=T)
  infos<-gsub("poly(R)-hydroxyalkanoic acid","poly(R)-hydroxyalk. acid",infos, fixed=T)
  infos<-gsub("Thauera linaloolentis.*","Thauera linaloolentis)",infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("^phasin","Phasin",infos)
  infos<- gsub("Phasin prot","Phasin",infos)
  infos<- gsub("phasin prot","phasin",infos)
  infos<- gsub("PHA synth (","PHA synthesis (",infos, fixed=T)
  infos<- gsub("depol","depolimerase",infos, fixed=T)
  infos<- gsub("depolimerase PhaZ","depolimerase",infos, fixed=T)
  infos<- gsub("ATPase YchF","ATPase",infos, fixed=T)
  infos<- gsub("reduct","reductase",infos, fixed=T)
  infos<- gsub("acetyltransf","acetyltransferase",infos, fixed=T)
  infos<- gsub("associated prot","associated protein",infos, fixed=T)
  infos<- gsub("dehydrog","dehydrogenase",infos)
  infos<- gsub("PaaC (paaC)","(paaC)",infos, fixed=T)
  infos<- gsub("synth (B","synthase (B",infos, fixed=T)
  infos<- gsub("synth, class I","synthase I",infos, fixed=T)
}
# infos<-gsub("enzyme","enz",infos,fixed=T)
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
base_plot <- ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.42, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-4),
        legend.position="bottom",
        plot.margin = margin(1,3,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other genes)",
       fill="",
       title = "Most abundant genes of PHA related genes (gene level)"
  )
base_plot + scale_fill_manual(values= rev(fill_color_30_BIS) )
# ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_GENES_PHA_associated_Most_abundant_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

# again, with colors for each type of function... (e.g., see table 1 of DOI 10.1080/21655979.2025.2458363 )
colors_type <- c("darkgreen" , "chartreuse4","darkslategray", "darkseagreen4",
                 # "green2","darkolivegreen3", "chartreuse",
                 "orangered2","darkorange","red","coral2",
                 "yellow","gold",
                 #"blue3",
                 "tomato2", "firebrick3",
                 "cyan2","dodgerblue3",
                 "deepskyblue","aquamarine","royalblue2","cyan4", "turquoise","deepskyblue3","dodgerblue","lightskyblue","aquamarine3",
                 "gray35","black","gray20",
                 "red3","orange2",
                 "chocolate4","brown3")
base_plot + scale_fill_manual(values= colors_type )
ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_GENES_PHA_associated_Most_abundant_PROPORTIONS_ON_THEMSELVES_2.png", width= 6, height = 5, dpi=300)


#exporting gene infos...
ordered_infos<-selected_infos # these are already ordered
row.names(ordered_infos)<-ordered_infos$Gene
ordered_infos$Entry<-NULL
ordered_infos<-ordered_infos[row.names(table_top), ]
ordered_infos$Pathway<-NULL
colnames(ordered_infos)[colnames(ordered_infos)=="Taxon"] <- "Example_of_ref_Taxon"

only_them$Protein<-NULL
only_them <- only_them[order(apply( only_them, MARGIN = 1, mean), decreasing = T), ]
AVERAGES_R1<- round(apply( only_them[, !grepl("R2",colnames(only_them))], MARGIN = 1, mean) ,2)
AVERAGES_R2<- round(apply( only_them[, !grepl("R1",colnames(only_them))], MARGIN = 1, mean) ,2)
to_save <- cbind.data.frame( round(only_them,2), AVERAGES_R1, AVERAGES_R2, ordered_infos )
write.csv2(file=paste0("RNA_PHA_Results/Most_abundant_GENES/Most_abundant__PHA_GENES_already_known_PROPORTIONS_ON_THEMSELVES.csv"), to_save, row.names=T, quote=F)

suppressWarnings(rm(table, table2, ordered_infos, table_top, table_to_plot, to_save, AVERAGES_EVERY_SAMPLE))




################ MOST ABUNDANT PHA RELATED GENES (Excluding phasins) ################

these_genes <- dictionary[ dictionary$Entry %in% PHA_IDs , c("Description","Gene") ]
these_genes <- these_genes[ !grepl("[P-p]hasin",these_genes$Description) , "Gene" ]
these_genes <- these_genes[these_genes!=""]
these_genes <- these_genes[!these_genes%in%c("thiG","ychF","cat","ilvD","ccrB")] # not related
these_genes <- these_genes[] # not related
table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% these_genes , ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

# backup_numbers <- table  # to count them afterwards

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
{ infos<- gsub("()", "(no gene name)", infos, fixed=T)
  infos<-gsub("alpha","α",infos,fixed=T)
  infos<-gsub("Phasin fam prot","Phasin family",infos,fixed=T)
  infos<-gsub("poly(R)-hydroxyalkanoic acid","poly(R)-hydroxyalk. acid",infos, fixed=T)
  infos<-gsub("Thauera linaloolentis.*","Thauera linaloolentis)",infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("^phasin","Phasin",infos)
  infos<- gsub("Phasin prot","Phasin",infos)
  infos<- gsub("phasin prot","phasin",infos)
  infos<- gsub("PHA synth (","PHA synthesis (",infos, fixed=T)
  infos<- gsub("depol","depolimerase",infos, fixed=T)
  infos<- gsub("depolimerase PhaZ","depolimerase",infos, fixed=T)
  infos<- gsub("ATPase YchF","ATPase",infos, fixed=T)
  infos<- gsub("reduct","reductase",infos, fixed=T)
  infos<- gsub("acetyltransf","acetyltransferase",infos, fixed=T)
  infos<- gsub("associated prot","associated protein",infos, fixed=T)
  infos<- gsub("dehydrog.*(fadJ)","dehydrogenase (fadJ",infos)
  infos<- gsub("dehydrog ","dehydrogenase ",infos)
  infos<- gsub("synth (B","synthase (B",infos, fixed=T)
  infos<- gsub("synth, class I","synthase I",infos, fixed=T)
  infos<- gsub("PHA","Poly-hydroxy-alkanoate",infos, fixed=T)
  infos<- gsub("depolimerase (SAMN","depolim (SAMN",infos, fixed=T)
  infos<- gsub("depolimerase, intracellular","depolim, intracellular",infos, fixed=T)
  infos<- gsub("fam prot","family",infos, fixed=T)
  infos<- gsub("genaseenase","genase",infos, fixed=T)
  infos<- gsub("coenzyme A","CoA",infos, fixed=T)
  infos<- gsub("Class I poly(R)-hydroxyalk. acid synth (phaC)","Poly(R)-hydroxyalk. acid synthase class I (phaC)",infos, fixed=T)
}
# infos<-gsub("enzyme","enz",infos,fixed=T)
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# bar plot
base_plot <- ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.45, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.5 ),
        legend.margin = margin(-14,0,2,-5),
        legend.position="bottom",
        plot.margin = margin(1,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other genes)",
       fill="",
       title = "Most abundant genes of PHA related genes (gene level), AFTER EXCLUDING PHASIN GENES"
  )
base_plot + scale_fill_manual(values= rev(fill_color_30_BIS) )
# ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_GENES_PHA_associated_Most_abundant_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

# again, with colors for each type of function...
colors_type <- c("darkgreen" ,"chartreuse4","darkslategray","darkolivegreen",
                 "darkolivegreen3","forestgreen", # ,"greenyellow", "limegreen",
                 "red","coral2","orangered","firebrick","orange","orangered3","red2","firebrick2",
                 "yellow", "gold",
                 "darkblue", 
                 "coral", 
                 "deepskyblue",
                 "violetred", # mch --> see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                 "black","royalblue4","gray30","navy","darkslateblue",
                 "goldenrod1","yellow2","gold2",
                 "chocolate4","brown3")
base_plot + scale_fill_manual(values= colors_type )
ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_GENES_PHA_NoPHASINS_Most_abundant_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)



################ PHA RELATED GENES * GLOMMING THROUGH DESCRIPTION * ###########################

# MODIFICATION TO DESCRIPTIONS HERE PERFORMED *BEFORE* THE COUNT GLOMMING
modification_PHA <- function(infos){
  # Enzyme types
  infos<- gsub("depol ","depolimerase ",infos, fixed=T)
  infos<- gsub("depol$","depolimerase",infos)
  infos<- gsub("depolimerase PhaZ","depolimerase",infos, fixed=T)
  infos<- gsub("ATPase YchF","ATPase",infos, fixed=T)
  infos<- gsub("reduct ","reductase ",infos, fixed=T)
  infos<- gsub("reduct$","reductase",infos)
  infos<- gsub("acetyltransf ","acetyltransferase ",infos, fixed=T)
  infos<- gsub("acetyltransf$","acetyltransferase",infos)
  infos<- gsub("dehydrog ","dehydrogenase ",infos)
  infos<- gsub("dehydrog, ","dehydrogenase ",infos)
  infos<- gsub("dehydrog$","dehydrogenase",infos)
  infos<-gsub("coenzyme A","CoA",infos,fixed=T)
  infos<-gsub("-coA ","-CoA",infos,fixed=T)
  infos<-gsub("family","fam",infos,fixed=T)
  # General PHA names
  infos<-gsub("Poly(R)-hydroxyalkanoic acid","PHA",infos, fixed=T)
  infos<-gsub("poly(R)-hydroxyalkanoic acid","PHA",infos, fixed=T)
  infos<-gsub("PolyR-hydroxyalkanoic acid","PHA",infos, fixed=T)
  infos<- gsub("[P-p]oly-hydroxyalkanoic acid", "PHA", infos)
  infos<-gsub("[P-p]olyhydroxyalkanoic acid","PHA",infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("Poly-beta-hydroxyalkanoate","PHA",infos, fixed=T)
  infos<- gsub("Poly(3-hydroxyalkanoate)", "PHA", infos, fixed = T)
  infos<- gsub("poly(3-hydroxyalkanoate)", "PHA", infos, fixed = T)
  # PHB
  infos<- gsub("[P-p]oly-beta-hydroxybutyrate", "PHB", infos)
  infos<- gsub("Poly-beta-hydroxybutyrate", "PHB", infos)
  infos<- gsub("Poly(3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("poly(3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("Poly (3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("poly (3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("Polyhydroxybutyrate", "PHB", infos, fixed = T)
  # Phasin related
  infos<-gsub("Phasin fam prot","Phasin family",infos,fixed=T)
  infos<-gsub("[P-]phasin fam$","Phasin family",infos)
  infos<- gsub("Phasin prot","Phasin",infos)
  infos<- gsub("phasin prot","phasin",infos)
  infos<- gsub("TIGR01841 fam phasin","Phasin family",infos)
  # Synthases
  infos<- gsub("PHA/PHB synth fam","PHA/PHB synthase",infos, fixed=T)
  infos<- gsub("3-hydroxyalkanoate synthetase","PHA/PHB synthase",infos, fixed=T) # see https://www.uniprot.org/uniprotkb/A0A133XN55/entry , search for the "InterPro" link
  infos<- gsub("PHA synth$","PHA synthase",infos)
  infos<- gsub("PHB synth$","PHB synthase",infos)
  infos<- gsub("PHA pol domain.*","PHA synthase",infos)
  infos<- gsub("PHA pol subunit.*","PHA synthase",infos)
  infos<- gsub("PHB pol family.*","PHB synthase",infos)
  infos<- gsub("PHB pol fam$","PHB synthase",infos)
  infos<- gsub("PHB pol domain.*","PHB synthase",infos)
  infos<- gsub("PHB pol subunit.*","PHB synthase",infos)
  infos<- gsub("PHB pol$","PHB synthase",infos)
  infos<- gsub("PHB pol N-term domain-containing","PHB synthase",infos)
  infos<- gsub("PhaM fam PHA granule multifunctional regulatory","PHB synthase activator",infos) #see doi: 10.1128/AEM.02935-13
  # Other enzymes
  infos<- gsub("Acetyl-CoA acetyltransf-1","Acetyl-CoA acetyltransferase",infos, fixed=T)
  infos<- gsub("Acetyl-CoA acetyltransfs","Acetyl-CoA acetyltransferase",infos, fixed=T)
  infos<- gsub("C-acetyltransferase","acetyltransferase",infos, fixed=T)
  infos<- gsub("AcetoAcetoacetyl","Acetoacetoacetyl",infos, fixed=T)
  infos<- gsub("(R)-specific enoyl-CoA hydratase","3-hydroxybutyryl-CoA dehydratase",infos, fixed=T)
  infos<- gsub("^[C-c]lass I ","",infos) # phaC
  infos<- gsub("synth, class I","synthase",infos, fixed=T) 
  infos<- gsub("synth class I","synthase",infos, fixed=T) 
  infos<- gsub("3-ketoacyl-CoA thiolase [fadN-fadA-fadE operon]","3-ketoacyl-CoA thiolase",infos, fixed=T) 
  infos<- gsub("[S-s]teroid 3-ketoacyl-CoA thiolase","3-ketoacyl-CoA thiolase",infos) 
  infos<- gsub("Beta-ketothiolase BktB","Beta-ketothiolase",infos) 
  infos<- gsub("Acyl-CoA dehydrogenaseFadE9","Acyl-CoA dehydrogenase",infos) 
  infos<- gsub("Short chain enoyl-CoA hydratase /3-hydroxyacyl-CoA dehydrogenase","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase/ enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T)  # fad J
  infos<- gsub("NAD-bind 3-hydroxyacyl-CoA dehydrogenase","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase NAD-bind.*","3-hydroxyacyl-CoA dehydrogenase",infos) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase NAD bind.*","3-hydroxyacyl-CoA dehydrogenase",infos) 
  infos<- gsub("(3R)-3-hydroxyacyl-CoA","3-hydroxyacyl-CoA",infos, fixed = T) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase type-2","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T) 
  infos<- gsub("(2S)-methylsuccinyl-CoA dehydro","Methylsuccinyl-CoA dehydro",infos, fixed = T) 
  # focus on depolimerases
  infos<- gsub("^[E-e]sterase/P","P",infos) 
  infos<- gsub("^[E-e]sterase P","P",infos) 
  infos<- gsub("^Esterase, P","P",infos) 
  infos<- gsub("[I-i]ntracellular ","",infos) 
  infos<- gsub("depolimerase fam esterase$","depolimerase",infos) 
  infos<- gsub("depolimerase fam esterase ","depolimerase",infos) 
  infos<- gsub("depolimerase fam$","depolimerase",infos) 
  infos<- gsub("PHB depolimerase family","depolimerase",infos) 
  infos<- gsub("depolimerase, intracellular","depolimerase",infos) 
  infos<- gsub("depol, intracellular","depolimerase",infos) 
  # Removing sporadic extras
  infos<-gsub("alpha","α",infos,fixed=T)
  infos<- gsub(" prot$", "", infos)  
  infos<- gsub("^[P-p]robable ", "", infos)  
  infos<- gsub("[P-p]utative ", "", infos)  
  infos<- gsub(", putative", "", infos, fixed=T)  
  infos<- gsub("'", "", infos, fixed=T)  
  infos<- gsub(" $", "", infos)  
  infos<- gsub("fam$", "family", infos)  
  infos<- Hmisc::capitalize(infos)
  infos_superassign <<- infos
}
# the descriptions are repeated along the dictionary, they have to be "glommed"
these_genes <- dictionary[ dictionary$Entry %in% PHA_IDs , ]
these_genes <- these_genes[!these_genes%in%c("thiG","ychF","cat","ilvD","ccrB") , ] # not related
PHA_IDs_unique <- these_genes
# NB: the starvation_genes obj was already created before, see the paragraph above
modification_PHA(infos = PHA_IDs_unique$Description)
PHA_IDs_unique$Description <- infos_superassign

### few descriptions are still not unique for the same gene, thus I'm solving them manually according to their Description and/or gene name..
{
  # Pha/phB
  PHA_IDs_unique[ grepl("[P-p]haC",PHA_IDs_unique$Description) , "Description"] <- "PHA synthase"
  PHA_IDs_unique[ grepl("[P-p]haC",PHA_IDs_unique$Description) , "Description"] <- "PHA synthase"
  PHA_IDs_unique[ grepl("^[P-p]haC",PHA_IDs_unique$Gene) , "Description"] <- "PHA synthase"
  PHA_IDs_unique[ grepl("[P-p]hbB",PHA_IDs_unique$Description) , "Description"] <- "Acetoacetyl-CoA reductase"
  PHA_IDs_unique[ grepl("^phbB",PHA_IDs_unique$Gene) , "Description"] <- "Acetoacetyl-CoA reductase"
  PHA_IDs_unique[ grepl("^phaJ",PHA_IDs_unique$Gene) , "Description"] <- "3-hydroxybutyryl-CoA dehydratase"
  PHA_IDs_unique[ grepl("^phaR",PHA_IDs_unique$Gene) , "Description"] <- "PHA synthesis repressor"
  PHA_IDs_unique[ grepl("PhaR",PHA_IDs_unique$Description) , "Description"] <- "PHA synthesis repressor"
  PHA_IDs_unique[ grepl("PhbR",PHA_IDs_unique$Description) , "Description"] <- "PHA synthesis repressor"
  # Phasins
  PHA_IDs_unique[ grepl("^phaP",PHA_IDs_unique$Gene) , "Description"] <- "Phasin"
  PHA_IDs_unique[ grepl("^phbP",PHA_IDs_unique$Gene) , "Description"] <- "Phasin"
  PHA_IDs_unique[ grepl("[P-p]haP",PHA_IDs_unique$Description) , "Description"] <- "Phasin"
  PHA_IDs_unique[ grepl("PHA granule-associated",PHA_IDs_unique$Description) , "Description"] <- "Phasin"
  PHA_IDs_unique[ grepl("PHA-granule associated",PHA_IDs_unique$Description) , "Description"] <- "Phasin"
  PHA_IDs_unique[ grepl("Granule-associated prot",PHA_IDs_unique$Description) , "Description"] <- "Phasin"
  # Metabolism
  PHA_IDs_unique[ grepl("^fadD",PHA_IDs_unique$Gene) , "Description"] <- "Long-chain acyl-CoA synthetase"
  PHA_IDs_unique[ grepl("^acsA",PHA_IDs_unique$Gene) , "Description"] <- "Acetyl-CoA synthetase"
  PHA_IDs_unique[ grepl("AcsA",PHA_IDs_unique$Description) , "Description"] <- "Acetyl-CoA synthetase"
  PHA_IDs_unique[ grepl("^paaC",PHA_IDs_unique$Gene) , "Description"] <- "3-hydroxyacyl-CoA dehydrogenase"
  PHA_IDs_unique[ grepl("PaaC",PHA_IDs_unique$Description) , "Description"] <- "3-hydroxyacyl-CoA dehydrogenase"
  PHA_IDs_unique[ grepl("PaaH",PHA_IDs_unique$Description) , "Description"] <- "3-hydroxyacyl-CoA dehydrogenase" # yes, the same as PaaC, see Uniprot descriptions
}
# View(PHA_IDs_unique[!duplicated(PHA_IDs_unique$Description) , ])

table_descr<-cbind.data.frame(table_abs[PHA_IDs_unique$Entry, ], descr=PHA_IDs_unique$Description)
# The abs counts are used in these as the abundance is not crucial "selection" rows, then the proportion will be re-computed afterwards
table_descr <- table_descr[ ! table_descr$descr=="" , ]  # fragments and proteins without name
table_descr <- table_descr[ ! is.na(table_descr$descr), ]

table_descr<-aggregate(.~descr, table_descr, FUN=sum)
# head( table_descr[order(table_descr$descr), ] )
row.names(table_descr)<-table_descr$descr
table_descr$descr<-NULL

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
table_to_plot$Protein <- gsub("PHA system","PHA system related protein", table_to_plot$Protein) # e.g. see N6XUS2 
table_to_plot$Protein <- gsub("cyclohydrolase$","cyclohydr.", table_to_plot$Protein) # e.g. see N6XUS2 
table_to_plot$Protein <- gsub("3-hydroxyacyl-CoA dehydrogenase / enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase", 
                              "3-hydroxyacyl-CoA dehydr/enoyl hydr/3-hydroxybutyryl epim.", table_to_plot$Protein, fixed = T)
# unique(table_to_plot$Protein)
colors_type <- c("chartreuse","green","green3","olivedrab3", "greenyellow","limegreen","chartreuse2", "darkgreen","green4",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "violet", 
                 "royalblue2","blue","deepskyblue",
                 "magenta","darkviolet",
                 "brown4",
                 "gray85","slategray4","slategray2",
                 "black","firebrick1","darkred","orange","red",
                 "yellow2","gold2",
                 "gray25",
                 "tomato", "orangered2",
                 "pink2"                 
)
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
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
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.245, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.62 ),
        legend.margin = margin(-14.5,0,2,-3),
        legend.position="bottom",
        plot.margin = margin(2,3,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       # title = "Most abundant mRNA of starvation related genes (transcript level)"
  )
ggsave("RNA_PHA_Results/GlommingThroughDescriptions_PHARelated_Most_abundant_PROPORTIONS_ON_THEMSELVES.png", width= 5.8, height = 5, dpi=300)




################ PHA GENES *** GLOMMING THROUGH DESCRIPTION + GENUS *** ###########################

# NB: this is the same dictionary built in the previous paragraph
table_descr<-cbind.data.frame(table_abs[PHA_IDs_unique$Entry, ], descr=PHA_IDs_unique$Description)
# The abs counts are used in these as the abundance is not crucial "selection" rows, then the proportion will be re-computed afterwards
table_descr <- table_descr[ ! table_descr$descr=="" , ]  # fragments and proteins without name
table_descr <- table_descr[ ! is.na(table_descr$descr), ]

# Genus
{ table_descr$Taxa <- PHA_IDs_unique$Taxon # NB: the table is already ordered according to this dictionary
  table_descr$Taxa <- gsub("'", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub("[", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub("]", "", table_descr$Taxa, fixed=T)
  table_descr$Taxa <- gsub("uncult bact", "Bacteria", table_descr$Taxa, fixed=T)  # same as another "Bacteria" taxon (both of which are very a-specific)
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

table_descr$descr <-paste0( table_descr$descr," (", table_descr$Taxa,")" )
# Now "descr" is Description + Genus
table_descr$Taxa <- NULL
table_descr<-aggregate(.~descr, table_descr, FUN=sum)
# head( table_descr[order(table_descr$descr), ] )
row.names(table_descr)<-table_descr$descr
table_descr$descr<-NULL

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
table_to_plot$Protein <- gsub("PHB synthase activator","PHB synthase activators", table_to_plot$Protein) # e.g. see N6XUS2 
table_to_plot$Protein <- gsub("cyclohydrolase$","cyclohydr.", table_to_plot$Protein) # e.g. see N6XUS2 
table_to_plot$Protein <- gsub("3-hydroxyacyl-CoA dehydrogenase / enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase", 
                              "3-hydroxyacyl-CoA dehydr/enoyl hydr/3-hydroxybutyryl epim.", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("domain$", "domain-containing", table_to_plot$Protein)
table_to_plot$Protein <- gsub("N-term domain-containing$", "N-term domain-contain", table_to_plot$Protein)

colors_type <- c("firebrick3", "orangered1","red3","indianred1",
                 "green3","green1", "darkgreen", "forestgreen","palegreen3","mediumspringgreen",
                 "brown3","firebrick2",
                 "green4",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "navyblue",
                 "darkgoldenrod2","gold1","yellow1",
                 "lightgoldenrod1","orange","goldenrod1",
                 "royalblue2","cadetblue3","deepskyblue1","royalblue2","cadetblue3",
                 "black","gray38",
                 "yellow2",
                 "darkmagenta","magenta"
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
        legend.key.height = unit(0.045, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.3 ),
        legend.margin = margin(-16,0,2,0),
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
ggsave("RNA_PHA_Results/GlommingThroughDescriptions_PHA_TAXON_Most_abundant.png", width= 5.45, height = 4.8, dpi=300)




################ CHECKING THE STARVATION REGULATORY RNA and GENES ###########################

### TRANSCRIPT LEVEL
table <- as.data.frame( prop_with_unmap[ starvation_genes$Entry ,  ] ) 
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

backup_numbers <- table  # to count them afterwards

top <- names(sort(rowSums(table), decreasing=TRUE))[1:30]
table_top<- table[top, ]
table_top$Protein<-row.names(table_top)
### Converting IDs
selected_infos<-dictionary[dictionary$Entry %in% table_top$Protein, ]
{infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
  infos<- gsub("()", "(Unknown source)", infos, fixed=T)
  # infos<- gsub("unclassified", "unclass", infos, fixed=T)
  infos<- gsub("Betaproteobacteria bact HGW-Betaproteobacteria-21", "HGW-Betaproteobacteria-21", infos, fixed=T)
  infos<- gsub("Bacteroidetes bact HGW Bacteroidetes 2", "HGW-Bacteroidetes-2", infos, fixed=T)
  infos<- gsub("communis SWub3 = DSM 12120", "SWub3=DSM12120", infos, fixed=T)
  infos<- gsub("ATP-dependent","ATP-dep",infos)
  infos<- gsub("ATP-bind ","",infos)
  infos<- gsub("proteolytic ","",infos)
  infos<- gsub("^Clp ","ATP-dep Clp ",infos)
  infos<- gsub("subunit","",infos)
  infos<- gsub(" prot Cl"," Cl",infos)
  infos<- gsub("PAS domain-containing prot","PAS domain",infos)
  infos<- gsub("PAC/Chase","PAC",infos , fixed=T)
  infos<- gsub("PAC and GAF","PAC-GAF",infos , fixed=T)
  infos<- gsub("^bi","Bi",infos)
  infos<- gsub("^diguanylate cyclase","diG-cyclase (c-di-GMP cyclase)",infos)
  infos<- gsub("sensor-containing diguanylate cyclase","sensor + diguanylate cyclase",infos)
  infos<- gsub("receiver modulated","",infos)
  infos<- gsub("phosphodiesterase","PDE",infos)
  infos<- gsub("uncult bact","uncultured bact",infos)
  infos<- gsub("diguanylate cyclase","diG-cyclase",infos)
  infos<- gsub("regulator  ","regulator ",infos)
  infos<- gsub("Clp protease ","protease ",infos)
  infos<- gsub("Universal stress (", "Universal stress protein (", infos, fixed=T)
  infos<- gsub("Nucleotide-bind universal stress, UspA fam", "Universal stress protein UspA fam", infos, fixed=T)
  infos<- gsub("  ", " ", infos, fixed=T)
  infos<- gsub("Heat-", "Heat ", infos, fixed=T)
  infos<- gsub(" like prot", "-like", infos, fixed=T)
  infos<- gsub("HSP", "Hsp", infos, fixed=T)
  infos<- gsub("^Hsp20 ", "Heat shock Hsp20", infos, fixed=T)
  infos<- gsub("fam", "family", infos, fixed=T)
  infos<- gsub(" prot ", " ", infos, fixed=T)
  # infos<- gsub("Heat shock Hsp", "Hsp", infos, fixed=T)
  infos<- gsub("Molecular chaper ", "", infos, fixed=T)
  infos<- gsub("^u", "U", infos)
  infos<- gsub("PDE (unclassified Thauera)", "PDE (Thauera)", infos, fixed=T)
  infos<- gsub("uncult ", "uncultured ", infos, fixed=T)
  infos<- gsub("Chaperonin GroEL ", "Hsp60-system chaperonin GroEL ", infos, fixed=T)
  infos<- gsub("stress ", "stress protein Usp ", infos, fixed=T)
  infos<- gsub("(Accumulibacter sp.)", "(Ca. Accumulibacter)*", infos, fixed=T)
  infos<- gsub("Nitrosomonas)", "Nitrosomonas)*", infos, fixed=T)
  infos<- gsub("Nitrosovibrio tenuis)", "Nitrosovibrio tenuis)*", infos, fixed=T)
  infos<- gsub("bact)", "bacterium)", infos, fixed=T)
}
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
table_to_plot$Protein <- gsub("Hsp20/alpha family (uncult Planktosalinus sp.)", "Hsp20/alpha family (uncultured Planktosalinus sp.)", table_to_plot$Protein, fixed=T)
colors_type <- c("darkgreen","chartreuse","forestgreen","greenyellow",
                 "darkolivegreen3","green2","springgreen2","springgreen3","chartreuse4",
                 "cyan2","lightskyblue", "blue", "dodgerblue3","deepskyblue2","deepskyblue3",
                 "deeppink2","yellow",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "indianred3", "firebrick2","tomato",
                 "red2","coral","firebrick","orangered", # "darkred",
                 "yellow2","gold2",
                  "gray32","black",
                 "violet","darkviolet"
                 # "gray35","gray18","slategray4",
                 # "royalblue1","deepskyblue","blue","cyan4","royalblue2","skyblue","deepskyblue2","dodgerblue","lightskyblue","cyan3","darkblue","cyan","dodgerblue3"
)
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
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
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.65 ),
        legend.margin = margin(-14,0,2,-40),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       title = "Most abundant mRNA of starvation related genes (transcript level)"
  )
ggsave("RNA_PHA_Results/Most_abundant_RNAs/WHOLE_DATASET_StarvationRelated_Most_abundant_transcr_AbsCounts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)

colSums(backup_numbers[ , !colnames(backup_numbers)%in%"Protein" ])  # proportions on whole community
# 1.394900    1.566607    1.175422    1.009808    1.296853    1.082531  



########### GENE LEVEL

these_genes <- dictionary[ dictionary$Entry %in% starvation_genes$Entry , "Gene" ]
these_genes <- these_genes[these_genes!=""]  # Quorum related *BUT* not linked to any gene --> Impossible to merge to other transcripts
table <- as.data.frame( prop_GENE_with_unmap[ row.names(prop_GENE_with_unmap) %in% these_genes , ] ) # used later when *only* PAO and GAO are required
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

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
{ infos<- gsub("()", "(no gene name)", infos, fixed=T)
  infos<- gsub("ATP-dependent","ATP-dep",infos)
  infos<- gsub("ATP-bind ","",infos)
  infos<- gsub("proteolytic ","",infos)
  infos<- gsub("^Clp ","ATP-dep Clp ",infos)
  infos<- gsub("subunit","",infos)
  infos<- gsub(" prot Cl"," Cl",infos)
  infos<- gsub("PAS domain-containing prot","PAS domain",infos)
  infos<- gsub("PAC/Chase","PAC",infos , fixed=T)
  infos<- gsub("PAC and GAF","PAC-GAF",infos , fixed=T)
  infos<- gsub("^bi","Bi",infos)
  infos<- gsub("^diguanylate cyclase","diG-cyclase (c-di-GMP cyclase)",infos)
  infos<- gsub("sensor-containing diguanylate cyclase","sensor + diguanylate cyclase",infos)
  infos<- gsub("receiver modulated","",infos)
  infos<- gsub("phosphodiesterase","PDE",infos)
  infos<- gsub("diguanylate cyclase","diG-cyclase",infos)
  infos<- gsub("regulator  ","regulator ",infos)
  infos<- gsub("Clp protease ","protease ",infos)
  infos<- gsub("Universal stress (", "Universal stress protein Usp (", infos, fixed=T)
  infos<- gsub("  ", " ", infos, fixed=T)
  infos<- gsub("Heat-", "Heat ", infos, fixed=T)
  infos<- gsub(" like prot", "-like", infos, fixed=T)
  infos<- gsub("HSP", "Hsp", infos, fixed=T)
  infos<- gsub("fam", "family", infos, fixed=T)
  infos<- gsub(" prot ", " ", infos, fixed=T)
  # infos<- gsub("Heat shock Hsp", "Hsp", infos, fixed=T)
  infos<- gsub("Molecular chaper ", "", infos, fixed=T)
  infos<- gsub("Chaperone DnaJ (dnaJ)", "Hsp chaperone DnaJ (dnaJ)", infos, fixed=T)
  infos<- gsub("Chaperone DnaJ (dnaJ)", "Hsp chaperone DnaJ (dnaJ)", infos, fixed=T)
  infos<- gsub("^Hsp20 ", "Heat shock Hsp20 ", infos)
  infos<- gsub("^u", "U", infos)
  infos<- gsub("Chaperonin GroEL ", "Hsp60-system chaperonin GroEL ", infos, fixed=T)
  infos<- gsub("Co-chaperonin GroES ", "Hsp60-system cochaperonin GroES ", infos, fixed=T)
  infos<- gsub("Chaperonin (groEL)", "Hsp60-system chaperonin GroEL (groEL)", infos, fixed=T)
  infos<- gsub("Chaperone DnaK ", "Hsp70-system chaperone DnaK ", infos, fixed=T)
  # infos<- gsub("Heat inducible transcription repressor HrcA", "Hsp GroEL transcriptional repressor", infos, fixed=T)
  infos<- gsub("Heat inducible transcription repressor HrcA", "GroEL's transcriptional repressor", infos, fixed=T)
  infos<- gsub("Protein GrpE", "Hsp70-system cochaperone GrpE", infos, fixed=T)
  infos<- gsub("Protease HtpX homolog (htpX)", "HSP protease HtpX homolog (htpX)", infos, fixed=T)  # yes, it is a HSP itself
  infos<- gsub("Universal stress ", "Universal stress protein Usp ", infos, fixed=T)  # yes, it is a HSP itself
}
table_top$Protein<-infos
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
# plotting ...
table_to_plot<-only_them
table_to_plot [ duplicated(table_to_plot$Protein) , "Protein" ] <- paste( table_to_plot[duplicated(table_to_plot$Protein),"Protein" ] , row.names(table_to_plot[duplicated(table_to_plot$Protein),]) ) # otherwise duplicated entries will be aggregated
table_to_plot<-melt(table_to_plot, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
colors_type <- c("chartreuse","forestgreen","greenyellow",
                 "darkolivegreen3","green2","springgreen3","chartreuse4","darkgreen",
                 "cyan2","blue","deepskyblue2","lightskyblue","royalblue1",
                 "yellow",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "darkviolet", 
                 "firebrick2","tomato","red3","indianred",
                 "red","coral","firebrick","orangered","coral2" , "darkred",
                 "yellow2","gold2",
                 # "black",
                 "grey50",
                 "black", "violet"
                 # "gray35","gray18","slategray4",
                 # "royalblue1","deepskyblue","blue","cyan4","royalblue2","skyblue","deepskyblue2","dodgerblue","lightskyblue","cyan3","darkblue","cyan","dodgerblue3"
)
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9) +
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
        axis.ticks.y=element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.11, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 7.65 ),
        legend.margin = margin(-14,0,2,-40),
        legend.position="bottom",
        plot.margin = margin(3,1,3,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other genes)",
       fill="",
       title = "Most abundant mRNA of starvation related genes (gene level)"
  )
ggsave("RNA_PHA_Results/Most_abundant_GENES/WHOLE_DATASET_StarvationRelated_Most_abundant_transcr_AbsCounts_PROPORTIONS_ON_THEMSELVES.png", width= 6, height = 5, dpi=300)




################ STARVATION GENES * GLOMMING THROUGH DESCRIPTION * ###########################

# MODIFICATION TO DESCRIPTIONS HERE PERFORMED *BEFORE* THE COUNT GLOMMING
modification_Stress <- function(infos){ 
  infos<- gsub("ATP-dependent","ATP-dep",infos)
  infos<- gsub("ATP-bind ","",infos)
  infos<- gsub("proteolytic ","",infos)
  infos<- gsub("^Clp ","ATP-dep Clp ",infos)
  infos<- gsub("subunit","",infos)
  infos<- gsub(" prot Cl"," Cl",infos)
  infos<- gsub("PAS domain-containing prot","PAS domain",infos)
  infos<- gsub("PAC/Chase","PAC",infos , fixed=T)
  infos<- gsub("PAC and GAF","PAC-GAF",infos , fixed=T)
  infos<- gsub("^bi","Bi",infos)
  # infos<- gsub("^diguanylate cyclase","diG-cyclase (c-di-GMP cyclase)",infos)
  infos<- gsub("sensor-containing diguanylate cyclase","sensor + diguanylate cyclase",infos)
  infos<- gsub("receiver modulated","",infos)
  infos<- gsub("phosphodiesterase","PDE",infos)
  infos<- gsub("diguanylate cyclase","diG-cyclase",infos)
  infos<- gsub("regulator  ","regulator ",infos)
  infos<- gsub("Clp protease N-term domain-containing prot","ATP-dep protease Clp N-term domain-containing",infos)
  infos<- gsub("Clp protease ","protease ",infos)
  infos<- gsub("Universal stress protein Usp ", "Universal stress protein ", infos, fixed=T)
  infos<- gsub("  ", " ", infos, fixed=T)
  infos<- gsub("'", "", infos, fixed=T)
  infos<- gsub("Heat-", "Heat ", infos, fixed=T)
  infos<- gsub(" like prot", "-like", infos, fixed=T)
  infos<- gsub("HSP", "Hsp", infos, fixed=T)
  infos<- gsub("fam", "family", infos, fixed=T)
  infos<- gsub(" prot ", " ", infos, fixed=T)
  # infos<- gsub("Heat shock Hsp", "Hsp", infos, fixed=T)
  infos<- gsub("Molecular chaper ", "", infos, fixed=T)
  infos<- gsub("Chaperone DnaJ", "Hsp chaperone DnaJ", infos, fixed=T)
  infos<- gsub("^Hsp20 ", "Heat shock Hsp20 ", infos)
  infos<- gsub("^u", "U", infos)
  infos<- gsub("Chaperonin GroEL ", "Hsp60-system chaperonin GroEL ", infos, fixed=T)
  infos<- gsub("Co-chaperonin GroES ", "Hsp60-system cochaperonin GroES ", infos, fixed=T)
  infos<- gsub("Co-chaper GroES", "Hsp60-system cochaperonin GroES", infos, fixed=T)
  infos<- gsub("Co-chaper GroES", "Hsp60-system cochaperonin GroES", infos)
  infos<- gsub("^GroES$", "Hsp60-system cochaperonin GroES ", infos)
  # infos<- gsub("GroES-like", "Hsp60-system cochaperonin GroES ", infos)
  infos<- gsub(".*GroES L", "Hsp60-system cochaperonin GroES ", infos)
  infos<- gsub("^Chaperonin GroEL", "Hsp60-system chaperonin GroEL", infos)
  infos<- gsub("Chaperonin (groEL)", "Hsp60-system chaperonin GroEL", infos, fixed=T)
  infos<- gsub(".*chaperonin.*GroES", "Hsp60-system cochaperonin GroES ", infos)  
  infos<- gsub("^GroEL$", "Hsp60-system chaperonin GroEL", infos)
  infos<- gsub("Chaperone DnaK ", "Hsp70-system chaperone DnaK ", infos, fixed=T)
  infos<- gsub("Heat shock 70.*", "Hsp70-system chaperone DnaK", infos)  
  infos<- gsub(".*Heat shock 70.*", "Hsp70-system chaperone DnaK", infos)  
  infos<- gsub(".*70 kDa heat*", "Hsp70-system chaperone DnaK", infos)  
  infos<- gsub("18 kDa heat shock", "Heat shock Hsp18", infos)  
  infos<- gsub("Heat shock Hsp20.*", "Heat shock Hsp20", infos)  
  infos<- gsub("Ribosome-associated heat shock Hsp15", "Heat shock Hsp15", infos)  
  # infos<- gsub("Heat inducible transcription repressor HrcA", "Hsp GroEL transcriptional repressor", infos, fixed=T)
  infos<- gsub("Heat inducible transcription repressor HrcA", "GroEL's transcriptional repressor", infos, fixed=T)
  infos<- gsub("Protein GrpE", "Hsp70-system cochaperone GrpE", infos, fixed=T)
  infos<- gsub("Protease HtpX homolog (htpX)", "HSP protease HtpX homolog", infos, fixed=T)  # yes, it is a HSP itself
  infos<- gsub("ATP-dep protease adapter ClpS", "ATP-dep protease ClpS", infos, fixed=T)  
  infos<- gsub("Universal stress protein Usp protein ", "Universal stress protein ", infos, fixed=T)  
  infos<- gsub("Endopeptidase Clp regulatory (ClpX)", "ATP-dep protease ClpX", infos, fixed=T)  
  infos<- gsub("Endopeptidase Clp regulatory ClpX", "ATP-dep protease ClpX", infos, fixed=T)  
  infos<- gsub(", clpA", " ClpA", infos, fixed=T)  
  infos<- gsub("Carbon starvation CstA.*", "Carbon starvation CstA", infos)  
  infos<- gsub("Carbon starvation A.*", "Carbon starvation CstA", infos) 
  infos<- gsub("Two-component system P regulon response regulator PhoB","P regulon transcript regulator PhoB",infos)
  infos<- gsub("P regulon transcript regulatory prot PhoB","P regulon transcript regulator PhoB",infos)
  infos<- gsub("P regulon sensor prot PhoR","P regulon sensor histidine kinase PhoR",infos)
  infos<- gsub("Bifunctional oligoribonuclease and PAP phosphatase NrnA", "Bifunctional oligoribonuclease/PAP phosphatase NrnA", infos)  
  infos<- gsub("ATP-dep ATP-dep protease Clp N-term", "ATP-dep protease Clp N-term", infos)  
  infos<- gsub(", MazF antagonist", "", infos)  
  infos<- gsub("GroES L", "GroES", infos)
  infos<- gsub("[C-c]old-shock", "Cold shock", infos)  
  infos<- gsub("[H-h]eat-shock", "Heat shock", infos) 
  infos<- gsub("[C-c]old shock prot, DNA bind$", "Cold shock family", infos)  
  infos<- gsub("[C-c]old-shock prot, DNA bind$", "Cold shock family", infos)  
  infos<- gsub("[C-c]old shock DNA-bind", "Cold shock", infos)  
  infos<- gsub("[C-c]old-shock DNA-bind", "Cold shock", infos)  
  infos<- gsub("Universal stress prot$", "Universal stress family", infos)  
  infos<- gsub("ATP-dep protease $", "ATP-dep protease Clp", infos)  
  infos<- gsub("ATP-dep protease$", "ATP-dep protease Clp", infos)  
  infos<- gsub("H-NS family nucleoid-associated regulatory$", "H-NS histone family", infos)  
  infos<- gsub("Nucleoid-structuring H-NS$", "H-NS histone family", infos)  
  infos<- gsub("H-NS histone$", "H-NS histone family", infos)  
  infos<- gsub("DNA-bind H-NS$", "H-NS histone family", infos)
  infos<- gsub("Bifunctional (p)ppGpp synthetase/guanosine-3,5-bis(diP) 3-pyrophosphohydrolase", "(P)ppGpp bifunctional synthase/hydrolase", infos, fixed = T)
  infos<- gsub("Bifunctional (p)ppGpp synthetase/guanosine-3',5'-bis(diP) 3'-pyrophosphohydrolase", "(P)ppGpp bifunctional synthase/hydrolase", infos, fixed = T)
  infos<- gsub("Bifunctional (P)ppGpp synthetase/guanosine-3,5-bis(DiP) 3-pyrophosphohydrolase", "(P)ppGpp bifunctional synthase/hydrolase", infos, fixed = T)
  infos<- gsub("GTP pyrophosphokinase, (P)ppGpp synthetase$", "(P)ppGpp synthetase", infos)
  infos<- gsub("GTP pyrophosphokinase rsh$", "(P)ppGpp synthetase", infos) # pyrophosphokinase = synthetase, see https://www.uniprot.org/uniprotkb/P0AG20/entry, 
  infos<- gsub("GTP pyrophosphokinase$", "(P)ppGpp synthetase", infos)  # pyrophosphokinase = synthetase, see https://www.uniprot.org/uniprotkb/P0AG20/entry
  infos<- gsub("TraR/DksA family transcript regulator", "RNA pol-bind transcription factor DksA", infos)
  infos<- gsub("RNA pol-bind DksA", "RNA pol-bind transcription factor DksA", infos)
  infos<- gsub("DNA-bind transcript regulator, Lrp family", "DNA-bind Lrp family transcript regulator", infos)
  infos<- gsub("Transcript regulator, AsnC/Lrp family", "DNA-bind Lrp family transcript regulator", infos)
  infos<- gsub(" prot$", "", infos)  
  infos<- gsub("^[P-p]robable ", "", infos)  
  infos<- gsub("^[P-p]utative ", "", infos)  
  infos<- gsub(", putative", "", infos, fixed=T)  
  infos<- gsub(" $", "", infos)  
  infos<- Hmisc::capitalize(infos)
  infos_superassign <<- infos
}
# the descriptions are repeated along the dictionary, they have to be "glommed"
starvation_genes_unique <- starvation_genes
# NB: the starvation_genes obj was already created before, see the paragraph above
modification_Stress(infos = starvation_genes_unique$Description)
starvation_genes_unique$Description <- infos_superassign
### few descriptions are still not unique for the same gene, thus I'm solving them manually according to their Description and/or gene name..
{
  # CLP
  starvation_genes_unique[ starvation_genes_unique$Entry %in% c("A0A037ZLN2","A0A1F3NAW7","UPI00048CEB86","UPI0009627D53") , "Description" ] <- "ATP-dep protease ClpP"
  starvation_genes_unique[ starvation_genes_unique$Entry %in% c("A0A7S6PV13","A0A4P5XDK1","A0A0D3MKN6","A0A077KDY8") , "Description" ] <- "ATP-dep protease ClpC"
  starvation_genes_unique[ grepl("ClpB",starvation_genes_unique$Description) , "Description" ] <- "ATP-dep protease ClpB"
  starvation_genes_unique[ grepl("^clpB",starvation_genes_unique$Gene) , "Description" ] <- "ATP-dep protease ClpB"
  starvation_genes_unique[ grepl("ClpA",starvation_genes_unique$Description) , "Description" ] <- "ATP-dep protease ClpA"
  starvation_genes_unique[ grepl("^clpA",starvation_genes_unique$Gene) , "Description" ] <- "ATP-dep protease ClpA"
  starvation_genes_unique[ grepl("ClpS",starvation_genes_unique$Description) , "Description" ] <- "ATP-dep protease adaptor ClpS"
  starvation_genes_unique[ grepl("^clpS",starvation_genes_unique$Description) , "Description" ] <- "ATP-dep protease adaptor ClpS"
  starvation_genes_unique[ grepl("ClpXP protease specificity-enhancing factor SspB",starvation_genes_unique$Description) , "Description" ] <- "ClpXP protease specificity-enhancing factor"
  starvation_genes_unique[ grepl("^clpP",starvation_genes_unique$Gene) , "Description" ] <- "ATP-dep protease ClpP"
  starvation_genes_unique[ grepl("^clpX",starvation_genes_unique$Gene) , "Description" ] <- "ATP-dep protease ClpX"
  starvation_genes_unique[ grepl("CplX",starvation_genes_unique$Description) , "Description" ] <- "ATP-dep protease ClpX"
  # HSP 
  starvation_genes_unique[ starvation_genes_unique$Entry %in% c("A0A0F2RJC5","UPI00041397C6") , "Description"] <- "Hsp60-system chaperonin GroEL"
  starvation_genes_unique[ starvation_genes_unique$Entry %in% c("A0A0B1RY98","A0A0K6HNJ5","A0A7J5V3I3","UPI00036A7D0A") , "Description"] <- "Hsp70-system chaperone DnaK"
  starvation_genes_unique[ starvation_genes_unique$Gene %in% c("dnak","dnaK") , "Description"] <- "Hsp70-system chaperone DnaK"
  starvation_genes_unique[ grepl("Hsp70",starvation_genes_unique$Description) , "Description"] <- "Hsp70-system chaperone DnaK"
  starvation_genes_unique[ grepl("^hsp70",starvation_genes_unique$Gene) , "Description"] <- "Hsp70-system chaperone DnaK"
  starvation_genes_unique[ grepl("^HSP70",starvation_genes_unique$Gene) , "Description"] <- "Hsp70-system chaperone DnaK"
  starvation_genes_unique[ grepl("Dna[J-j]",starvation_genes_unique$Description) , "Description"] <- "Hsp chaperone DnaJ"
  starvation_genes_unique[ grepl("^dnaJ",starvation_genes_unique$Description) , "Description"] <- "Hsp chaperone DnaJ"
  starvation_genes_unique[ grepl("dnaK",starvation_genes_unique$Gene) , "Description"] <- "Hsp70-system chaperone DnaK"
  starvation_genes_unique[ grepl("IbpA",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp20 family molecular chaper IbpA"
  starvation_genes_unique[ grepl("^ibpA",starvation_genes_unique$Gene) , "Description"] <- "Heat shock Hsp20 family molecular chaper IbpA"
  starvation_genes_unique[ grepl("Hsp20",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp20 family molecular chaper IbpA"
  starvation_genes_unique[ grepl(".*Hsp33.*",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp33 HslO"
  starvation_genes_unique[ grepl(".*Hsp90*",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp90"
  starvation_genes_unique[ grepl(".*hsp-90*",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp90"
  starvation_genes_unique[ grepl(".*hsp-90*",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp90"
  starvation_genes_unique[ grepl("Heat shock 90",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp90"
  starvation_genes_unique[ grepl("Heat shock 90",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp90"
  starvation_genes_unique[ grepl("Co-chaper GroES",starvation_genes_unique$Description) , "Description"] <- "Heat shock Hsp90"
  starvation_genes_unique[ grepl("Alcohol dehydrog GroES",starvation_genes_unique$Description) , "Description"] <- "Hsp90 GroES containing alcohol dehydrogenase"
  starvation_genes_unique[ grepl("Chaperonin GroEL",starvation_genes_unique$Description) , "Description"] <- "Hsp60-system chaperonin GroEL"
  starvation_genes_unique[ grepl("Molecular chaper GroEL",starvation_genes_unique$Description) , "Description"] <- "Hsp60-system chaperonin GroEL"
  starvation_genes_unique[ grepl("^hsp82",starvation_genes_unique$Gene) , "Description"] <- "Heat shock Hsp82"
  # c di GMP
  which_ones <- grepl("[D-d]iguanylate|diG-|DiG-cyclase|[D-d]iG cyclase|c-di-GMP|[C-c]yclic-di-GMP-|Cyclic di-GMP",starvation_genes_unique$Description) & grepl("cyclase",starvation_genes_unique$Description) & grepl("PDE|phosphodiester",starvation_genes_unique$Description)
  starvation_genes_unique[ which_ones , "Description"] <- "Diguanylate (c-di-GMP) cyclase/PDE"
  which_ones <- grepl("[D-d]iguanylate|diG-|DiG-cyclase|[D-d]iG cyclase|c-di-GMP|[C-c]yclic-di-GMP-|Cyclic di-GMP",starvation_genes_unique$Description) & grepl("cyclase",starvation_genes_unique$Description) & !grepl("PDE|phosphodiester",starvation_genes_unique$Description)
  starvation_genes_unique[ which_ones , "Description"] <- "Diguanylate (c-di-GMP) cyclase"
  which_ones <- grepl("[D-d]iguanylate|c-di-GMP|[C-c]yclic-di-GMP-|Cyclic di-GMP",starvation_genes_unique$Description) & grepl("PDE|phosphodiester",starvation_genes_unique$Description) & !grepl("cyclase",starvation_genes_unique$Description)
  starvation_genes_unique[ which_ones , "Description"] <- "Diguanylate (c-di-GMP) PDE"
  # Cold Shock
  starvation_genes_unique[ grepl("CspA",starvation_genes_unique$Description) , "Description"] <- "Cold shock CspA"
  starvation_genes_unique[ grepl("^cspA",starvation_genes_unique$Gene) , "Description"] <- "Cold shock CspA"
  starvation_genes_unique[ grepl("^cspC",starvation_genes_unique$Gene) , "Description"] <- "Cold shock CspC"
  # DPS
  starvation_genes_unique[ grepl("Dps|DPS",starvation_genes_unique$Description) & grepl("bind",starvation_genes_unique$Description) , "Description"] <- "DPS (DNA protection during starvation)"
  starvation_genes_unique[ grepl("^dps",starvation_genes_unique$Gene) , "Description"] <- "DPS (DNA protection during starvation)"
  starvation_genes_unique[ grepl("DNA starvation/stationary phase protection",starvation_genes_unique$Description) , "Description"] <- "DPS (DNA protection during starvation)"
  # Usp
  starvation_genes_unique[ grepl("UspA",starvation_genes_unique$Description) , "Description"] <- "Universal stress protein UspA"
  starvation_genes_unique[ grepl("^uspF",starvation_genes_unique$Gene) , "Description"] <- "Universal stress protein UspF"
  starvation_genes_unique[ grepl("UspE",starvation_genes_unique$Description) , "Description"] <- "Universal stress protein UspE"
  # P or C specific starvation
  starvation_genes_unique[ grepl("CstA",starvation_genes_unique$Description) , "Description"] <- "Carbon starvation CstA fam"
  starvation_genes_unique[ grepl("^Carbon starvation prot A$",starvation_genes_unique$Description) , "Description"] <- "Carbon starvation CstA fam"
  starvation_genes_unique[ grepl("PhoH|PsiH",starvation_genes_unique$Description) , "Description"] <- "P-starvation-inducible PsiH family" # They are the same https://pubmed.ncbi.nlm.nih.gov/8444794/
  starvation_genes_unique[ grepl("P-starvation-inducible E",starvation_genes_unique$Description) , "Description"] <- "P-starvation-inducible PsiE"
  starvation_genes_unique[ grepl("P-starvation-inducible PsiE",starvation_genes_unique$Description) , "Description"] <- "P-starvation-inducible PsiE family"
  # ppGpp
  starvation_genes_unique[ grepl("SpoT",starvation_genes_unique$Description) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase"
  starvation_genes_unique[ grepl("^spoT",starvation_genes_unique$Gene) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase"
  starvation_genes_unique[ grepl("RelA| rsh",starvation_genes_unique$Description) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase" # yes, double function also for relA
  starvation_genes_unique[ grepl("^relA",starvation_genes_unique$Gene) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase" # yes, double function also for relA
  starvation_genes_unique[ grepl("GTP pyrophosphokinase",starvation_genes_unique$Gene) , "Description"] <- "(P)ppGpp synthetase" # yes, double function also for relA
  # DksA
  starvation_genes_unique[ grepl("^dksA",starvation_genes_unique$Gene) , "Description"] <- "RNA pol-bind transcription factor DksA"
}
# starvation_genes_unique <- starvation_genes_unique[!duplicated(starvation_genes_unique$Description), ]

table_descr<-cbind.data.frame(table_abs[starvation_genes_unique$Entry, ], descr=starvation_genes_unique$Description)
table_descr <- table_descr[ ! table_descr$descr=="" , ]  # fragments and proteins without name
table_descr <- table_descr[ ! is.na(table_descr$descr), ]

table_descr<-aggregate(.~descr, table_descr, FUN=sum)
# head( table_descr[order(table_descr$descr), ] )
row.names(table_descr)<-table_descr$descr
table_descr$descr<-NULL

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
table_to_plot$Protein <- gsub("Universal stress family", "Universal stress protein family", table_to_plot$Protein)
table_to_plot$Protein <- gsub("domain$", "domain-containing", table_to_plot$Protein)
table_to_plot$Protein <- gsub("N-term domain-containing$", "N-term domain-contain", table_to_plot$Protein)
table_to_plot$Protein <- gsub("Cold shock$", "Cold shock family", table_to_plot$Protein)
table_to_plot$Protein <- gsub("Clp$", "Clp family", table_to_plot$Protein)
table_to_plot$Protein <- gsub("transp$", "transporter", table_to_plot$Protein)
# unique(table_to_plot$Protein)

# table_to_plot$Protein <- gsub("Hsp20/alpha family (uncult Planktosalinus sp.)", "Hsp20/alpha family (uncultured Planktosalinus sp.)", table_to_plot$Protein, fixed=T)
colors_type <- c("violet","pink",
                 "chartreuse","green3","chartreuse3","greenyellow","forestgreen","yellowgreen",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "cyan","dodgerblue","royalblue2","cadetblue2",
                 "yellow","gold2",
                 "chocolate","black",
                 "darkred","firebrick1","red","firebrick3","tomato3","orangered","indianred1",
                 "gray90","plum2",
                 "navyblue","gray30","maroon1",
                 "darkviolet","darkmagenta"
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
        strip.text = element_text(size=10),
        plot.title = element_text(size=6.5),
        legend.key.height = unit(0.095, "cm"),
        legend.key.width = unit(0.325, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text ( size = 9 ),
        legend.text = element_text ( size = 8.1 ),
        legend.margin = margin(-14.5,0,2,-6),
        legend.position="bottom",
        plot.margin = margin(2,1,1,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       # title = "Most abundant mRNA of starvation related genes (transcript level)"
  )
ggsave("RNA_PHA_Results/GlommingThroughDescriptions_StarvationRelated_Most_abundant_PROPORTIONS_ON_THEMSELVES.png", width= 5.45, height = 5, dpi=300)




################ STARVATION  GENES *** GLOMMING THROUGH DESCRIPTION + GENUS *** ###########################

# NB: this is the same dictionary built in the previous paragraph
table_descr<-cbind.data.frame(table_abs[starvation_genes_unique$Entry, ], descr=starvation_genes_unique$Description)
# The abs counts are used in these as the abundance is not crucial "selection" rows, then the proportion will be re-computed afterwards
table_descr <- table_descr[ ! table_descr$descr=="" , ]  # fragments and proteins without name
table_descr <- table_descr[ ! is.na(table_descr$descr), ]

# UspA counts are here lost otherwise, albeit it belongs to this family (which is conversely featured in the final plot regardless UspA)
table_descr$descr <- gsub("Universal stress protein UspA","Universal stress family", table_descr$descr)

# Genus
{ table_descr$Taxa <- starvation_genes_unique$Taxon # NB: the table is already ordered according to this dictionary
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

table_descr$descr <-paste0( table_descr$descr," (", table_descr$Taxa,")" )
# Now "descr" is Description + Genus
table_descr$Taxa <- NULL
table_descr<-aggregate(.~descr, table_descr, FUN=sum)
# head( table_descr[order(table_descr$descr), ] )
row.names(table_descr)<-table_descr$descr
table_descr$descr<-NULL

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
table_to_plot$Protein <- gsub("(P)ppGpp bifunctional ", "(P)ppGpp ", table_to_plot$Protein, fixed=T)
table_to_plot$Protein <- gsub("MFS transp ", "MFS transporters ", table_to_plot$Protein, fixed=T)
table_to_plot$Protein <- gsub("ATP-dep", "ATP-dependent", table_to_plot$Protein)
table_to_plot$Protein <- gsub("Heat shock Hsp20 family molecular chaper", "Hsp20 family chaperon", table_to_plot$Protein)
table_to_plot$Protein <- gsub("Universal stress family", "Universal stress protein family", table_to_plot$Protein)
table_to_plot$Protein <- gsub("domain$", "domain-containing", table_to_plot$Protein)
table_to_plot$Protein <- gsub("N-term domain-containing$", "N-term domain-contain", table_to_plot$Protein)
table_to_plot$Protein <- gsub("Cold shock CspA", "Cold shock protein CspA", table_to_plot$Protein)

colors_type <- c("deeppink",
                 "greenyellow","chartreuse1","darkgreen","green3",
                 # "cyan2","deepskyblue2","blue1","cadetblue2","darkblue","cyan4","lightsteelblue1","skyblue","deepskyblue","royalblue3","dodgerblue","mediumblue","lightskyblue","steelblue","cyan3","royalblue1","blue3","dodgerblue3","deepskyblue3",
                 "cyan","lightskyblue","cyan2","deepskyblue1","royalblue2",
                 "black",
                 "yellow","gold",
                 "darkred","tomato1","firebrick1","red","firebrick3","tomato3","orangered","indianred1",
                 "gray80","gray65","wheat",
                 "deeppink1",
                 "yellow4","darkblue","darkmagenta",
                 "violet","maroon1"
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
        legend.margin = margin(-16,0,2,-4),
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
ggsave("RNA_PHA_Results/GlommingThroughDescriptions_StarvationTAXON_Most_abundant.png", width= 5.45, height = 4.8, dpi=300)
