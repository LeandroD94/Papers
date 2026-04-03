# ssh -X matteo@150.217.62.209
# ssh -X matteo@150.217.62.222
# cd /home/matteo/Desktop/Works_in_progress/
# R


##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggvenn")
  library("ggh4x")
  library("ggbreak")
  # analysis packages
  library("vegan")
  # utilities
  library("xlsx")
  library("reshape")
}

{ 
  dir.create("Data_check")
  dir.create("Results_DNAseq")
  dir.create("Results_DNAseq/Abundances")
  dir.create("Results_DNAseq/Alpha_diversity")
  # dir.create("Results_DNAseq/PCoA")
}

options(scipen = 100) # disable scientific annotation


### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","orange","coral2","yellow3","darkmagenta",
                 "pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red1","chartreuse3","darkslategray3")
fill_color_30<-c("chartreuse","grey65","deeppink",
                 "orange","darkmagenta", 
                 "navyblue","black","cyan", "cadetblue",
                 "slateblue3","gold2","brown4",
                 "greenyellow", "springgreen3",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue1",
                 "yellow", "chartreuse3","yellow4", "chocolate4","coral2","darkgreen",
                 "slateblue1", "red","darkslategray3")
# NB: there is always an extra color which will be the "Others" group




####################### IMPORTING DATA #####################


load(file="Phyloseq_objs_from_Kaiju_classif.RData")

# using the IDs as names
sample<-sample_names(data)
original_names<-sample
sample

#sample<-gsub("_bracken_abund.tsv","",sample)
#sample_names(data)<-sample # update

Metadata <- as.data.frame(read.table(file = "metadata.txt", header = T, sep="\t"))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)


sample_data(data)$Sample_name <- as.factor( sample_data(data)$Sample_name ) 
# sample_data(data)$Exp_Day <- factor( sample_data(data)$Exp_Day , c(81,117) )

sample_data(data.genus)<- sample_data(data)


### proportions
data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100 )
data.genus.prop<-transform_sample_counts(data.genus, function(x) (x/sum(x))*100 )



 
############################# UNCLASSIFIED CHECK #################################
 
# data_temp<-tax_glom(data, taxrank = "Domain", NArm = F) # without glomming before, some error may be seen in ggplot2 (maybe bugs due to too many observations?)
# data_temp_prop<-transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )
# # colSums(otu_table(data_temp_prop))
# #View(cbind(otu_table(data_temp),tax_table(data_temp)))
# table<-psmelt(data_temp_prop)
# table$Domain<-factor(table$Domain, levels=c( unique(table$Domain)[unique(table$Domain)!="unclassified"] , "unclassified" ) )
# 
# ggplot(data=table, aes(x=Sample_name_16S, y=Abundance, fill=as.factor(Domain))) +
#   geom_bar( stat="identity", position="stack", na.rm = F) + 
#   facet_grid(~ Group, scales = "free_x", space = "free_x") +
#   theme_classic(base_size =8.5) +
#   scale_fill_manual(values=c("unclassified"="lightgray",
#                              "Bacteria"="deepskyblue",
#                              "Archaea"="chartreuse4",
#                              #"Eukaryota_Fungi"="gold3",
#                              #"Eukaryota_Metazoa"="gold",
#                              #"Eukaryota_Protists"="yellow",
#                              "Eukaryota"="yellow",
#                              "Viruses"="violet")
#   ) +
#   theme(axis.text.x=element_text(angle=40, vjust=1, hjust=1, size= 5.5),
#         strip.text = element_text(size=6),
#         legend.key.size = unit(0.35, "cm"),
#         legend.text = element_text ( size = 11 ),
#         legend.position="bottom",
#         legend.margin = margin(-2,30,4,0),
#         plot.margin = margin(2,2,2,10),
#   ) +
#   scale_x_discrete(expand=c(0.05,1))+
#   #scale_y_sqrt(breaks=c(0,5,10,25,50,100)) + # can't be applied due to the values below 1
#   guides(fill=guide_legend(nrow=1)) +
#   labs(x="", y="Percentual abundance of clades",
#        fill="",
#        title = "Domains and unclassified proportions")
# ggsave(file="Data_check/Domains_and_unclassified_according_to_Kaiju.png",width=7.05,height=4.8, dpi=450)
# dev.off()




###################### DOMAINS PROPORTIONS #########################

data_temp<-data
tax_temp<-as.data.frame(tax_table(data_temp))
# Despite it is a "Sar" and "Alveolata", Symbodinium is a Dinoflagellate with photosynthesis, hence an algae!
tax_temp$Euk_clade[tax_temp$Genus %in% "Symbiodinium"] <- "Protist"
# unique(paste0("Eukaryota_",tax_temp$Euk_clade[tax_temp$Domain %in% "Eukaryota"]))
tax_temp$Domain[tax_temp$Domain %in% "Eukaryota"] <- tax_temp$Euk_clade[tax_temp$Domain %in% "Eukaryota"]
tax_table(data_temp)<-as.matrix(tax_temp)
data_temp<-tax_glom(data_temp, taxrank = "Domain", NArm = F) # without glomming before, some error may be seen in ggplot2 (maybe bugs due to too many observations?)

#data_temp<-subset_taxa(data_temp, Domain!="Unclassified" )
# !!! NB: due to a strange glitch of gglopt2, the "Bacteria" color in the stacked bar plot became invisible in few sample if the Domain column is modified in a factor with the Bacteria as first level (as below)
# therefore, ggplot2 will be tricked: the "Unclassified" will still be present but their abundance will be very minimal, in order to put them as "first place" (impossible to see) and the bacteria as second level (which actually seems the first one)
otu_table(data_temp)[ tax_table(data_temp)[, "Domain"]=="unclassified" ] <- 0.000001

data_temp<-transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )
table<-psmelt(data_temp)

# unique(table$Domain)
table$Domain<-factor(table$Domain, levels =  c("unclassified","Bacteria","Archaea",
                                               "Nematoda","Rotifera","Protozoa",
                                               "Fungi","Algae","Amoeba",
                                               "Platyhelminthes","Protist","Viruses") ) # NB: unclassified are invisible (see above)
#table$Domain<-relevel(table$Domain,c("Unclassified","Bacteria") )

ggplot(data=table, aes(x=Sample_name, y=Abundance, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =10) +
  scale_fill_manual(values=c("Bacteria"="deepskyblue",
                             "Archaea"="blue3",
                             "Amoeba"="grey75",
                             "Fungi"="orange",
                             "Algae"="green3",
                             "Nematoda"="brown4",
                             "Platyhelminthes"="coral2",
                             "Rotifera"="red2",
                             "Protozoa"="yellow",
                             "Protist"="gold4",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5, size= 9),
        axis.text.y=element_text(size= 5.85),
        axis.title.y=element_text(size= 11, hjust = 0.65),
        axis.title.x = element_blank(),
        # axis.title.x = element_text(size= 10.5),
        plot.title =element_text(size= 7),
        plot.caption =element_text(size= 7),
        strip.text = element_text(size=10),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text ( size = 10 ),
        legend.position="bottom",
        # legend.key.spacing.x = unit(0.45,"cm"),
        legend.margin = margin(-3,20,-5,10),
        plot.margin = margin(0,5,-10,1),
  ) +
  scale_y_break(c(2.5, 97.5), 
                scales=0.035,
                space = 0.25) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    #breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5, 97.5, 100),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5 , 2, 2.5, 97.5, 100),
    limits = c(0,100)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Percent abundance of clades",
       fill="",
       # x ="Experiment day",
       title = "Domains proportions",
       caption= "NB: the Y axis is scaled to improve the readibility of lower values"
  )
ggsave(file="Data_check/Domains_according_to_Kaiju.png",width=4.5,height=4, dpi=300)
ggsave(file="Results_DNAseq/Abundances/Domains_according_to_Kaiju.png",width=4.5,height=4, dpi=300)
dev.off()




############################ RAREFACTION ANALYSIS ################################

png(file="Data_check/Rarefaction_curve.png",width=2500,height=2000, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=5000,label=T) # using the filtered (analysed) data
dev.off()
rm(r)




###################### ABUNDANCES BAR PLOT ##########################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.genus.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top_data <- subset_taxa(data_temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE)) [1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca.", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 6.8
                                 ),
        axis.title =element_text(size=6),
        axis.title.x = element_text(vjust=2.5),
        strip.text = element_text(size=10.5),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.45, "cm"),
        legend.text = element_text ( size = 9 , face="italic" ),
        legend.position="bottom",
        legend.margin = margin(-10,20,0,0),
        plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified microbial genera", 
       caption = " 'Others' includes every genus below rank 29")
ggsave(file="Results_DNAseq/Abundances/TOP_microbial_genera.png",width=5,height=4.5, dpi=300)
dev.off()


# means of TOP Genera
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average R"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R"]],1,mean)), 2),
                            "Average T"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="T"]],1,mean)), 2) #,
                            #"Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="Inoculum"]],1,mean)), 2
                            )
to_save<-to_save[order(to_save$`Average R`, decreasing=T), ]
write.xlsx(file = "Results_DNAseq/Abundances/TOP_Genera_Average_abundances.xlsx", row.names = F, to_save)




### TOP Species
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top_data <- subset_taxa(data_temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Species)){ # to mark the fungal taxa
  #   tabella_top$Species[i]<-paste0(tabella_top$Species[i],"_",tabella_top$Domain[i])
  # }
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Species<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the species
  tabella<-aggregate(.~Species+Group+Sample_name, tabella[ , c("Group","Abundance","Species","Sample_name")], FUN=sum)
  tabella$Species<-gsub ("Candidatus ","Ca.", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("Betaproteobacteria bacterium ","", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("agariperforans","agariperfor.", tabella$Species)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Species)) +
  facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 8.5
  ),
  plot.title =element_text(size=6),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.18, "cm"),
  legend.text = element_text ( size = 7.85),
  legend.position="bottom",
  legend.margin = margin(-10,32,0,0),
  plot.margin = margin(4,4,4,5)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x=" ", y="Percentual abundance of clades", 
       fill="",
       title = "Most abundant identified microbial species", 
       caption = " 'Others' includes every species below rank 29")
ggsave(file="Results_DNAseq/Abundances/TOP_microbial_Species.png",width=5.2,height=4.5, dpi=300)
dev.off()

# means of TOP Species
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]] ,
                            "Average R"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R"]],1,mean)), 2),
                            "Average T"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="T"]],1,mean)), 2) 
                            )
#,
                            # "Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="Inoculum"]],1,mean)), 2)

to_save<-to_save[order(to_save$`Average R`, decreasing=T), ]
write.xlsx(file = "Results_DNAseq/Abundances/TOP_Species_Average_abundances.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




###################### EXTRA CHECK: ABUNDANCES OF EUKARYOTA ##########################

### TOP eukar species
suppressWarnings(rm(top, others, tabella, unass_data))
{ 
  data.subset <- subset_taxa(data.genus, Domain=="Eukaryota")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Genus)){ # to mark the fungal taxa
  #   tabella_top$Genus[i]<-paste0(tabella_top$Genus[i],"_",tabella_top$Phylum[i])
  # }
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<- paste0( tabella$Euk_clade, ": ", tabella$Genus )
  tabella$Genus <- gsub(".*: Others","Others",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
levels(tabella$Genus) <- gsub("Protozoa: Symbiodinium","Protist: Symbiodinium", levels(tabella$Genus), fixed = T)
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=as.factor(Genus))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 9
  ),
  axis.title =element_text(size=10),
  plot.title =element_text(size=7),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.18, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 8.2 ),
  legend.position="bottom",
  legend.margin = margin(-12,30,1,0),
  plot.margin = margin(2,4,2,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Clades percentages (only eukaryotes)", 
       fill="",
       title = "Most abundant identified eukaryota genera (percentages on eukaryotes only)", 
       caption = " 'Others' includes every eukaryota below rank 29, no unclassified")
ggsave(file="Results_DNAseq/Abundances/TOP_Eukaryota_genera_ONLY_EUKAR.png",width=4.85,height=4.5, dpi=300)
dev.off()

# means of TOP Genera
to_save<- cbind.data.frame( "Euk_clade"= as.data.frame(tax_table(prune.dat_top))[["Euk_clade"]] ,
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average R"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R"]],1,mean)), 2),
                            "Average T"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="T"]],1,mean)), 2)
                            ) #"Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="Inoculum"]],1,mean)), 2)

to_save<-to_save[order(to_save$`Average R`, decreasing=T), ]
write.csv(to_save, file = "Results_DNAseq/Abundances/TOP_Eukaryota_genera_ONLY_EUKAR_average_abundances.csv", row.names = F, quote=F)


### Again, less genera this time...
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data.genus, Domain=="Eukaryota")
  data.subset <- subset_taxa(data.subset, Genus!="Symbiodinium")  # this genus makes no sense biologically!
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Genus)){ # to mark the fungal taxa
  #   tabella_top$Genus[i]<-paste0(tabella_top$Genus[i],"_",tabella_top$Phylum[i])
  # }
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<- paste0( tabella$Euk_clade, ": ", tabella$Genus )
  tabella$Genus <- gsub(".*: Others","Others",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
fill_color_customised<-c("navyblue","gray85","brown4",
                         "darkcyan", "darkgreen",
                         "gold3","orange",
                         "indianred",  "coral2",
                         "slateblue1", "grey15", "violet", 
                         "yellow","blue1", "magenta","cyan", "chocolate4",
                         "red", "green2","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_customised) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 7
  ),
  axis.title =element_text(size=9),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.12, "cm"),
  legend.key.width = unit(0.16, "cm"),
  legend.text = element_text ( size =7.1 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-12,42,0,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results_DNAseq/Abundances/TOP23_Eukaryota_genera_ONLY_EUKAR.png",width=4.8,height=4, dpi=300)
dev.off()




###################### EXTRA CHECK: ABUNDANCES OF ARCHAEA ##########################

### TOP archaea species
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data, Domain=="Archaea")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Species)){ # to mark the fungal taxa
  #   tabella_top$Species[i]<-paste0(tabella_top$Species[i],"_",tabella_top$Phylum[i])
  # }
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Species<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-gsub ("Candidatus ","Ca.", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("uncultured","uncult.", tabella$Species)
  tabella$Species<-gsub ("^archaeon$","unknown Archaea", tabella$Species)
  tabella$Species<-gsub (" archaeon "," ", tabella$Species)
  tabella<-aggregate(.~Species+Group+Sample_name, tabella[ , c("Group","Abundance","Species","Sample_name")], FUN=sum)
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=as.factor(Species))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 7.5
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.15, "cm"),
  legend.text = element_text ( size = 7.25 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-10,42,0,0),
  plot.margin = margin(2,4,2,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Clades percentages(only archaea)", 
       fill="",
       title = "Most abundant identified archaea species", 
       caption = " 'Others' includes every archaea below rank 29, no unclassified")
ggsave(file="Results_DNAseq/Abundances/TOP_Archaea_species_ONLY_ARCHAEA.png",width=5.2,height=4.5, dpi=300)
dev.off()

# means of TOP Species
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]] ,
                            "Average R"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R"]],1,mean)), 2),
                            "Average T"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="T"]],1,mean)), 2) ) # ,
                            #"Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="Inoculum"]],1,mean)), 2)

to_save<-to_save[order(to_save$`Average R`, decreasing=T), ]
write.xlsx(file = "Results_DNAseq/Abundances/TOP_Archaea_Average_abundances.xlsx", row.names = F, to_save)




###################### EXTRA CHECK: ABUNDANCES OF VIRUSES ##########################
# 
# ### TOP Viruses (species, focus on ONLY viruses)
# suppressWarnings(rm(top, others, tabella, unass_data))
# {
#   data.subset <- subset_taxa(data, Domain=="Viruses")
#   # top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
#   # prune.dat_top <- prune_taxa(top,data.subset)
#   # View(tax_table(prune.dat_top))  # check row names
#   # View(Tax_IDs) # check row name --> take the corresponding tax ID
#   # solving duplicates
#   tax_table(data.subset)["sp1570", "Domain_species"]<- "Viruses_Caudoviricetes 2788787"
#   tax_table(data.subset)["sp9099", "Domain_species"]<- "Viruses_Caudoviricetes 2731619"
#   data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
#   top <- names(sort(taxa_sums(data.subset), decreasing=TRUE))[1:29]
#   prune.dat_top <- prune_taxa(top,data.subset)
#   others<-taxa_names(data.subset)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.subset) # only viruses
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Species<-"Others"
#   tabella_others$Domain_species<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Domain_species <- gsub ("Viruses_","", tabella$Domain_species)
#   tabella$Domain_species<-factor(tabella$Domain_species, levels = c(unique(tabella$Domain_species)[! unique(tabella$Domain_species) %in% "Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample_name_16S, y=Abundance, fill=Domain_species)) +
#   facet_grid(~ Group, scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") +
#   theme_classic(base_size =8.5) +
#   scale_fill_manual(values=fill_color_30) +
#   scale_y_continuous(expand=c(0,1) ) +
#   theme(axis.text.x=element_text(angle=40,
#                                  vjust=1,
#                                  hjust=1,
#                                  size= 7
#   ),
#   axis.title =element_text(size=10),
#   axis.title.x = element_text(vjust=2.5),
#   # strip.text = element_text(size=6.7),
#   legend.key.height = unit(0.2, "cm"),
#   legend.key.width = unit(0.25, "cm"),
#   legend.text = element_text ( size = 9 ),
#   legend.position="bottom",
#   legend.margin = margin(-5,42,0,0),
#   plot.margin = margin(4,4,4,15)) +
#   guides(fill=guide_legend(nrow=10)) +
#   labs(x="", y="Percentual abundance of clades",
#        fill="",
#        title = "Most abundant identified viruses species (focus on viruses only)",
#        caption = " 'Others' includes every virus below rank 29")
# ggsave(file="Results_DNAseq/Abundances/TOP_DNAVirus_species_ONLY_VIRUS.png",
#        width=6.5,height=4.8, dpi=300)
# dev.off()
# 
# 
# # means of TOP virus
# to_save<- cbind.data.frame( "Species"= as.data.frame(tax_table(prune.dat_top))[["Domain_species"]] ,
#                             "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R1"]],1,mean)), 2),
#                             "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R2"]],1,mean)), 2),
#                             "Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="Inoculum"]],1,mean)), 2)
# )
# to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
# write.xlsx(file = "Results_DNAseq/Abundances/TOP_Viruses_Averages.xlsx", row.names = F, to_save)
# 




######################## ALPHA DIVERSITY ##############################

### PROKARYOTES (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Bacteria") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Group")
# pAlpha
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
  geom_boxplot(data=pAlpha$data, aes(x=Group, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Sample_name), color="red2" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Eubacteria genera)", 
       y="Alpha Diversity Measure",
       x="Group scale") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Genera_Eubacteria.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### PROKARYOTES (SPECIES)
data.temp<-subset_taxa(data, Domain %in% c("Bacteria") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Group")
# pAlpha
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
  geom_boxplot(data=pAlpha$data, aes(x=Group, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Sample_name), color="red2" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Eubacteria species)", 
       y="Alpha Diversity Measure",
       x="Group scale") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Species_Eubacteria.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### Archaea (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Archaea") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Group")
# pAlpha
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
  geom_boxplot(data=pAlpha$data, aes(x=Group, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Sample_name), color="red4" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Archaea genera)", 
       y="Alpha Diversity Measure",
       x="Group scale") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Genera_Archaea.png", width = 5.5,height =4.5, dpi=300)




### Eukaryota (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Eukaryota") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Group")
# pAlpha
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
  geom_boxplot(data=pAlpha$data, aes(x=Group, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Sample_name), color="chocolate3" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Eukaryota genera)", 
       y="Alpha Diversity Measure",
       x="Group scale") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Genera_Eukaryota.png", width = 5.5,height =4.5, dpi=300)




########################## FOCUS ON AOB AOA and NOB ###########################

# AOBs derive from from MIDAS, also in other articles I can't find other names...
AOB_list <- c("Nitrosomonas", "Nitrosococcus", "Nitrosomonadales",   # synonymous (see MIDAS) or higher levels
              "Nitrosospira",  "Nitrosolobus", "Nitrosovibrio",   # synonymous, see MIDAS
              "966-1",     # https://www.sciencedirect.com/science/article/pii/S0045653523014649
              "Ellin6067" )   # https://www.sciencedirect.com/science/article/pii/S0960852423008970
# adding also AOA, from PMID: 24559743
AOA_list <- c("Nitrososphaera", "Nitrososphaerota", 
              "Nitrosocaldus", 
              "Nitrosotalea",
              "Nitrocosmicus")  
AOA_list <- c(AOA_list,
              "Nitrosarchaeum", # DOI: 10.1099/ijsem.0.002926
              "Nitrosacidococcus", # https://doi.org/10.1016/j.wroa.2022.100157
              "Nitrosocosmicus",  # https://www.sciencedirect.com/science/article/pii/S0038071720300225
              "Nitrosocaldus", "Nitrosothermus" , "Nitrosocaldaceae", # https://doi.org/10.3389/fmicb.2020.608832
              "Nitrosoabyssus", "Nitrosymbiomonas", # doi: 10.1186/s40168-025-02146-2
              "Nitrosopelagicus", # doi.org/10.1002/9781118960608.gbm01969
              "Nitrosopolaris", # https://doi.org/10.1093/femsmc/xtac019
              "Nitrosopumilales", "Nitrosopumilaceae",  # https://doi.org/10.1002/9781118960608.obm00122
              "Nitrososphaeraceae",   # yes, also as "Genus", ref from https://www.sciencedirect.com/science/article/pii/S0038071720300225
              "Nitrosotenuis" # https://journals.asm.org/doi/10.1128/aem.01430-18
)

# NOBs
NOB_list <- c("Nitrospira","Nitrospirae", "Nitrospirae", "Nitrospirales","Nitrospiraceae", "Nitrospirota",  # alcuni da MIDAS, ma tutti elencati in DOI: 10.1016/S0076-6879(11)86005-2
              # few Nitrospira names are quite similar, these are synonims or, at least, they still are of Nitrospirales order which is still known for AO activity, see https://www.nature.com/articles/s41396-023-01518-6
              "Nitromaritima",
              "Nitrospina","Nitrospinae","Nitrospinaceae","Nitrospinota",
              "Nitrobacter",
              "Nitrococcus",
              "Nitronereus" , # DOI: 10.1038/s41396-023-01518-6
              "Nitrohelix", "Nitronauta", "Nitrocaldera", "Nitrotheca", # doi: 10.1021/acs.est.3c00636
              "Nitrolancea", #doi: 10.1099/ijs.0.062232-0.
              "Nitrotoga")

data.genus.temp<-subset_taxa( data.genus.prop, Genus!="unclassified")   # NB: removing unclassified!
data.genus.temp<- subset_taxa(data.genus.temp, !Domain %in% c("Viruses"))  # otherwise viruses as "Flavobacterium phage" will be included!
# searching
vector_without_Ca <- gsub("Candidatus_","", tax_table( data.genus.temp )[ ,"Genus"])  # No "Candidatus" but same order, to a easier use of the logical below
data.genus.temp<-subset_taxa( data.genus.temp, vector_without_Ca %in% c(AOB_list, AOA_list, NOB_list) )   # NB: removing unclassified!
# solving few synonimies (see notes above)
tax_table(data.genus.temp)[,"Genus"] <- gsub("Nitrospirae","Nitrospira", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"] <- gsub("Nitrospinae","Nitrospina", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"] <- gsub("Nitrosovibrio","Nitrosospira", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"] <- gsub("Nitrososphaerota","Nitrososphaera", tax_table(data.genus.temp)[,"Genus"])
data.genus.temp <- tax_glom(data.genus.temp, taxrank = "Genus", NArm = F)
# TOP
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.genus.temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:20]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Other AOOs or NOOs"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca.", tabella$Genus)
}
# adding infos on who is what ...
tabella$Genus2 <- gsub("Ca.","", tabella$Genus, fixed = T)
tabella$Genus[tabella$Genus2%in%AOB_list] <- paste( "AOB:", tabella$Genus[tabella$Genus2%in%AOB_list])
tabella$Genus[tabella$Genus2%in%AOA_list] <- paste( "AOA:", tabella$Genus[tabella$Genus2%in%AOA_list])
tabella$Genus[tabella$Genus2%in%NOB_list] <- paste( "NOB:", tabella$Genus[tabella$Genus2%in%NOB_list])
# ordering levels
tabella$Genus<-factor(tabella$Genus, levels = c( unique(tabella$Genus)[! unique(tabella$Genus) %in% "Other AOOs or NOOs"] ,
                                                 "Other AOOs or NOOs"
                                                 )
                      )

base_plot <- ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_y_continuous(expand=c(0,0.03) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 7
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=9.5),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.35, "cm"),
  legend.text = element_text ( size = 8.9 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-12,30,0,0),
  plot.margin = margin(4,4,4,10)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified AOB,AOA or NOB organisms", 
       caption = " 'Other' includes every accumulator below rank 20")

custom_colors<-c("yellow","springgreen2","wheat","royalblue1","green2","yellow3","indianred1",
                 "pink3", "gray70", "red2","red4","magenta1","darkgreen","darkmagenta",
                 "deepskyblue2","firebrick3","lightskyblue1","navyblue","red1","darkviolet",
                 "black")
base_plot +
  scale_fill_manual(values=custom_colors)
ggsave(file="Results_DNAseq/Abundances/AOB_AOA_NOB_organisms_OneColorForEach.png",width=4.85,height=4.5, dpi=300)
dev.off()

# averages
to_save<- cbind.data.frame( "Euk_clade"= as.data.frame(tax_table(prune.dat_top))[["Euk_clade"]] ,
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average R"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R"]],1,mean)), 2),
                            "Average T"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="T"]],1,mean)), 2)
)
to_save<-to_save[order(to_save$`Average R`, decreasing=T), ]
write.xlsx(file = "Results_DNAseq/Abundances/AOB_AOA_NOB_organisms_Averages.xlsx", row.names = F, to_save)




##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results_DNAseq/R_version_and_packages.txt")
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
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))


