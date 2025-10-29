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
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Alpha_diversity")
  dir.create("Results/PCoA")
  # dir.create("Results/Correlations")
}

options(scipen = 100) # disable scientific annotation


### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
fill_color_30<-c("wheat3","deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "aquamarine2",
                 "bisque2","cyan","yellow3","brown",
                 "springgreen4", "firebrick3",
                 "grey15","deepskyblue2","lightgreen","orange2",
                 "darkorchid", "lightblue1" ,"blue1", "violet", 
                 "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                 "lightgrey", "red","darkslategray3")
# NB: there is always an extra color which will be the "Others" group




####################### IMPORTING DATA #####################


load(file="Phyloseq_objs_from_Kaiju_classif.RData")

# using the IDs as names
sample<-sample_names(data)
original_names<-sample
sample

#sample<-gsub("_bracken_abund.tsv","",sample)
#sample_names(data)<-sample # update

Metadata <- as.data.frame(read.table(file = "metadata.csv", header = T, sep=";"))
Metadata$FASTQ_Code<-gsub("_S.*","",Metadata$FASTQ_Code)
row.names(Metadata)<-Metadata$FASTQ_Code # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)


sample_data(data)$Reactor <-factor( sample_data(data)$Reactor , 
                                           levels= c("Inoculum","R1","R2") 
                                           )
sample_data(data)$Time_point <- as.factor( sample_data(data)$Time_point ) 
sample_data(data)$Sample_name_16S <- factor( sample_data(data)$Sample_name_16S , sample_data(data)$Sample_name_16S  )
sample_data(data)$Exp_Day <- factor( sample_data(data)$Exp_Day , c(81,117) )

sample_data(data.genus)<- sample_data(data)


### proportions
data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100 )
data.genus.prop<-transform_sample_counts(data.genus, function(x) (x/sum(x))*100 )




############################# UNCLASSIFIED CHECK #################################

data_temp<-tax_glom(data, taxrank = "Domain", NArm = F) # without glomming before, some error may be seen in ggplot2 (maybe bugs due to too many observations?)
data_temp_prop<-transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )
# colSums(otu_table(data_temp_prop))
#View(cbind(otu_table(data_temp),tax_table(data_temp)))
table<-psmelt(data_temp_prop)
table$Domain<-factor(table$Domain, levels=c( unique(table$Domain)[unique(table$Domain)!="unclassified"] , "unclassified" ) )

ggplot(data=table, aes(x=Sample_name_16S, y=Abundance, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Reactor, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("unclassified"="lightgray",
                             "Bacteria"="deepskyblue",
                             "Archaea"="chartreuse4",
                             #"Eukaryota_Fungi"="gold3",
                             #"Eukaryota_Metazoa"="gold",
                             #"Eukaryota_Protists"="yellow",
                             "Eukaryota"="yellow",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust=1, size= 5.5),
        strip.text = element_text(size=6),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom",
        legend.margin = margin(-2,30,4,0),
        plot.margin = margin(2,2,2,10),
  ) +
  scale_x_discrete(expand=c(0.05,1))+
  #scale_y_sqrt(breaks=c(0,5,10,25,50,100)) + # can't be applied due to the values below 1
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Percentual abundance of clades",
       fill="",
       title = "Domains and unclassified proportions")
ggsave(file="Data_check/Domains_and_unclassified_according_to_Kaiju.png",width=7.05,height=4.8, dpi=450)
dev.off()




###################### DOMAINS PROPORTIONS #########################

data_temp<-data
tax_temp<-as.data.frame(tax_table(data_temp))
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

ggplot(data=table, aes(x=Exp_Day, y=Abundance, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  facet_grid(~ Reactor, scales = "free_x", space = "free_x") +
  theme_classic(base_size =10) +
  scale_fill_manual(values=c("Bacteria"="deepskyblue",
                             "Archaea"="blue3",
                             "Amoeba"="grey75",
                             "Fungi"="orange",
                             "Algae"="green3",
                             "Nematoda"="brown4",
                             "Platyhelminthes"="coral3",
                             "Rotifera"="red3",
                             "Protozoa"="yellow",
                             "Protist"="gold4",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1, size= 9),
        axis.text.y=element_text(size= 6),
        axis.title.y=element_text(size= 10.5),
        axis.title.x = element_blank(),
        # axis.title.x = element_text(size= 10.5),
        plot.title =element_text(size= 7),
        plot.caption =element_text(size= 7),
        strip.text = element_text(size=9.5),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom",
        # legend.key.spacing.x = unit(0.45,"cm"),
        legend.margin = margin(-3,10,-5,10),
        plot.margin = margin(0,5,-10,1),
  ) +
  scale_y_break(c(0.5, 97.5), 
                scales=0.035,
                space = 0.25) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    #breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5, 97.5, 100),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5, 97.5, 100),
    limits = c(0,100)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Percent abundance of clades",
       fill="",
       # x ="Experiment day",
       title = "Domains proportions",
       caption= "NB: the Y axis is scaled to improve the readibility of lower values"
  )
ggsave(file="Data_check/Domains_according_to_Kaiju.png",width=5.2,height=4, dpi=300)
dev.off()




# ############################ RAREFACTION ANALYSIS ################################
# 
# 
# png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
# r<-rarecurve(t(as(otu_table(data),"matrix")), step=300,label=F) # using the filtered (analysed) data
# dev.off()
# rm(r)
# 
# 
# 



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
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.8
                                 ),
        axis.title =element_text(size=6),
        axis.title.x = element_text(vjust=2.5),
        strip.text = element_text(size=10.5),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.45, "cm"),
        legend.text = element_text ( size = 10 ),
        legend.position="bottom",
        legend.margin = margin(-10,20,0,0),
        plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified microbial genera", 
       caption = " 'Others' includes every genus below rank 29")
ggsave(file="Results/Abundances/TOP_microbial_genera.png",width=6,height=4.8, dpi=300)
dev.off()


# means of TOP Genera
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R1"]],1,mean)), 2),
                            "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R2"]],1,mean)), 2) #,
                            #"Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="Inoculum"]],1,mean)), 2
                            )
to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_Genera_Average_abundances.xlsx", row.names = F, to_save)



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
  tabella<-aggregate(.~Species+Reactor+Exp_Day, tabella[ , c("Reactor","Abundance","Species","Exp_Day")], FUN=sum)
  tabella$Species<-gsub ("Candidatus ","Ca.", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("Betaproteobacteria bacterium ","", tabella$Species)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=Species)) +
  facet_grid(cols= vars(Reactor),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 9
  ),
  plot.title =element_text(size=6),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=10.5),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.42, "cm"),
  legend.text = element_text ( size = 8.5),
  legend.position="bottom",
  legend.margin = margin(-8,20,0,0),
  plot.margin = margin(4,4,4,5)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x=" ", y="Percentual abundance of clades", 
       fill="",
       title = "Most abundant identified microbial species", 
       caption = " 'Others' includes every species below rank 29")
ggsave(file="Results/Abundances/TOP_microbial_Species.png",width=6,height=4.8, dpi=300)
dev.off()

# means of TOP Species
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]] ,
                            "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R1"]],1,mean)), 2),
                            "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R2"]],1,mean)), 2) 
                            )
#,
                            # "Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="Inoculum"]],1,mean)), 2)

to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_Species_Average_abundances.xlsx", row.names = F, to_save)

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
  tabella$Genus<- paste0( tabella$Euk_clade, "_", tabella$Genus )
  tabella$Genus <- gsub(".*_Others","Others",tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=as.factor(Genus))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Reactor, scales = "free_x", space = "free_x") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 9
  ),
  axis.title =element_text(size=10),
  plot.title =element_text(size=7),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.20, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(-12,30,1,0),
  plot.margin = margin(2,4,2,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Clades percentages (only eukaryotes)", 
       fill="",
       title = "Most abundant identified eukaryota genera (percentages on eukaryotes only)", 
       caption = " 'Others' includes every eukaryota below rank 29, no unclassified")
ggsave(file="Results/Abundances/TOP_Eukaryota_genera_ONLY_EUKAR.png",width=6,height=4.8, dpi=300)
dev.off()

# means of TOP Genera
to_save<- cbind.data.frame( "Euk_clade"= as.data.frame(tax_table(prune.dat_top))[["Euk_clade"]] ,
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R1"]],1,mean)), 2),
                            "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R2"]],1,mean)), 2)
                            ) #"Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="Inoculum"]],1,mean)), 2)

to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
write.csv(to_save, file = "Results/Abundances/TOP_Eukaryota_genera_ONLY_EUKAR_average_abundances.csv", row.names = F, quote=F)




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
  tabella$Species<-gsub ("Candidatus ","Ca. ", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("uncultured","uncult.", tabella$Species)
  tabella$Species<-gsub ("^archaeon$","unknown Archaea", tabella$Species)
  tabella$Species<-gsub ("Methanobrevibacter sp. TLL-48-HuF1","Methanobrevibacter TLL-48-HuF1", tabella$Species)
  tabella<-aggregate(.~Species+Reactor+Sample_name_16S, tabella[ , c("Reactor","Abundance","Species","Sample_name_16S")], FUN=sum)
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name_16S, y=Abundance, fill=as.factor(Species))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Reactor, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.8
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=6.7),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.35, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(-10,42,0,0),
  plot.margin = margin(2,4,2,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Clades percentages(only archaea)", 
       fill="",
       title = "Most abundant identified archaea species", 
       caption = " 'Others' includes every archaea below rank 29, no unclassified")
ggsave(file="Results/Abundances/TOP_Archaea_species_ONLY_ARCHAEA.png",width=6.5,height=4.8, dpi=300)
dev.off()

# means of TOP Species
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]] ,
                            "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R1"]],1,mean)), 2),
                            "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R2"]],1,mean)), 2) ) # ,
                            #"Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="Inoculum"]],1,mean)), 2)

to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_Archaea_Average_abundances.xlsx", row.names = F, to_save)





# ###################### EXTRA CHECK: ABUNDANCES OF VIRUSES ##########################
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
#   facet_grid(~ Reactor, scales = "free_x", space = "free_x") +
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
# ggsave(file="Results/Abundances/TOP_DNAVirus_species_ONLY_VIRUS.png",
#        width=6.5,height=4.8, dpi=300)
# dev.off()
# 
# 
# # means of TOP virus
# to_save<- cbind.data.frame( "Species"= as.data.frame(tax_table(prune.dat_top))[["Domain_species"]] ,
#                             "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R1"]],1,mean)), 2),
#                             "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R2"]],1,mean)), 2),
#                             "Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="Inoculum"]],1,mean)), 2)
# )
# to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
# write.xlsx(file = "Results/Abundances/TOP_Viruses_Averages.xlsx", row.names = F, to_save)
# 




######################## ALPHA DIVERSITY ##############################


### PROKARYOTES (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Bacteria", "Archaea") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor")
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Sample_name_16S), color="red2" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (prokaryote genera)", 
       y="Alpha Diversity Measure",
       x="Reactor scale") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_PROKARYOTES.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### PROKARYOTES (SPECIES)
data.temp<-subset_taxa(data, Domain %in% c("Bacteria", "Archaea") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor")
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Sample_name_16S), color="red2" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (prokaryote species)", 
       y="Alpha Diversity Measure",
       x="Reactor scale") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_PROKARYOTES.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )




########################### PCoA AND BETADIVERSITY ########################

data.prop_bacteria_Genus<-subset_taxa(data.genus.prop, Domain %in% c("Bacteria","Archea"))
data.prop_bacteria_Species<-subset_taxa(data.prop, Domain %in% c("Bacteria","Archea"))
data.prop_eukar_Genus<-subset_taxa(data.genus.prop, Domain=="Eukaryota") # to focus only on fungi
data.prop_eukar_Species<-subset_taxa(data.prop, Domain=="Eukaryota")
# data.prop_virus_Genus<-subset_taxa(data.genus, Domain=="Viruses")
# data.prop_virus_Species<-subset_taxa(data.prop, Domain=="Viruses")

for( d in list(data.genus.prop, data.prop, 
               data.prop_bacteria_Genus, data.prop_bacteria_Species, #, 
               data.prop_eukar_Genus, data.prop_eukar_Species
               #data.prop_virus_Genus, data.prop_virus_Species
               ) 
     ) {
  
  # setting the automatic labels and file names
  if(identical(d,data.genus.prop) | 
     identical(d,data.prop_bacteria_Genus) | 
     identical(d,data.prop_eukar_Genus) # |
     # identical(d,data.prop_virus_Genus) 
     ){  # ){
    levels<-"genera"
  } else { levels <-"species"}
if(identical(d,data.genus.prop) | identical(d,data.prop)){
  #subset<-"PCoA_with_Every_microbes"
  subset<-"PCoA"
  file_name<-"microbial"
  plot_title<-"identified microbes"
}
if(identical(d,data.prop_bacteria_Genus) | identical(d,data.prop_bacteria_Species)){
  #subset<-"PCoA_with_Prokaryotes_only"
  subset<-"PCoA"
  file_name<-"prokaryotic"
  plot_title<-"identified prokaryotes"
}
if(identical(d,data.prop_eukar_Genus) | identical(d,data.prop_eukar_Species)){
  #subset<-"PCoA_with_Eukar_only"
  subset<-"PCoA"
  file_name<-"Eukar"
  plot_title<-"identified eukariotes"
}
# if(identical(d,data.prop_virus_Genus) | identical(d,data.prop_virus_Species)){
#   #subset<-"PCoA_with_virus_only"
#   subset<-"PCoA"
#   file_name<-"Viruses"
#   plot_title<-"identified viruses"
# }
  
# sample_names(d)<-gsub("","",sample_names(d))  # to change names in plot


# ### computing PERMANOVA analysis on groups
# sqrt_prop_table<-as.data.frame(t(sqrt(as.data.frame(otu_table(d)))))
# perm<- adonis2(sqrt_prop_table ~Project,
#                data=as(sample_data(d),"data.frame"),
#                permutations = 9999,
#                method="euclidean") # --> hellinger
# perm_result<-round( perm$`Pr(>F)`[1] , 3)

# ### computing BETADISPERSION analysis on groups
# dist<-vegan::vegdist(sqrt_prop_table, distance="euclidean")
# disper<-vegan::betadisper(dist, as(sample_data(d),"data.frame")[["Project"]] )
# disp<-vegan::permutest(disper, permutations=9999)
# disp_result<-round( disp$tab$`Pr(>F)`[1] , 3)


data.prop.labels<-d


{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor" ) +
  # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
  geom_point(size=1.5) + #  color="deepskyblue2") +
  geom_point(size=4.5, alpha= 0.6)+ # , color="deepskyblue2") +
  theme_classic(base_size = 9) +
  #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_name_16S), color="black", size=2, 
            hjust=0.7, vjust=-0.1, show.legend = FALSE) +
  theme(legend.margin = margin(10,0,10,0),
        legend.text = element_text(size=7.8)) +
  guides(shape="none", color="none") + 
  labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
       shape="",
       caption=" ",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file=paste0("Results/",subset,"/PCoA_with_ID_Hellinger_on_",file_name,"_",levels,"_Description.png"), width = 5, height = 4, dpi=300)


# # with ellipses
# plot_ordination(data.sqrt_prop, ordBC , color="Reactor2") +
#   # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
#   geom_point(size=2.5) +
#   geom_point(size=4.4, alpha= 0.6) +
#   theme_classic(base_size = 10.2) +
#   # stat_ellipse(linewidth=0.18) + 
#   labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
#        # subtitle = paste("  PERMANOVA Pr(>F):", perm_result, "\n  Betadispersion Pr(>F):", disp_result),
#        color="Project", 
#        shape="Reactor scale",
#        x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,".png"), width = 5.5, height = 4, dpi=300)
# 

}





########################## FOCUS ON PHA ACCUMULATORS ###########################

table_PHA<-read.table(file="PHA_Accumulators_maggio2025.tsv", fill = T, header = T, sep="\t")

data.genus.temp<- subset_taxa(data.genus.prop, Genus %in% table_PHA[,1] )
data.genus.temp<- subset_taxa(data.genus.temp, !Domain %in% c("Viruses"))  # otherwise viruses as "Flavobacterium phage" will be included!
# Solving synonyms (see LPSN)...
tax_table(data.genus.temp)[,"Genus"]<-gsub("Fuscovulum","Rhodobacter", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Phaeovulum","Rhodobacter", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Paucibacter","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Mitsuaria","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Kinneretia","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Roseateles","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("^Competibacter","Candidatus_Competibacter", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("^Accumulibacter","Candidatus_Accumulibacter", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Paramicrobacterium","Microbacterium", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Aureobacterium","Microbacterium", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Prescottella","Rhodococcus", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Ensifer Casida","Sinorhizobium", tax_table(data.genus.temp)[,"Genus"])
tax_table(data.genus.temp)[,"Genus"]<-gsub("Calidifontimicrobium","Azohydromonas", tax_table(data.genus.temp)[,"Genus"])
data.genus.temp<-tax_glom(data.genus.temp, taxrank = "Genus", NArm = F)  # re-glomming
data.genus.temp<- prune_taxa(taxa_sums(data.genus.temp)>0, data.genus.temp)

suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.genus.temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Other accumulators"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Other accumulators"],"Other accumulators"))
}
ggplot(data=tabella, aes(x=Exp_Day, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.8
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=6.7),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9.2 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-9,25,0,0),
  plot.margin = margin(4,4,4,10)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified accumulating PHA genera", 
       caption = " 'Other accumulators' includes every accumulator below rank 29")
ggsave(file="Results/Abundances/TOP_accum_PHA_genera.png",width=6.5,height=4.8, dpi=300)
dev.off()

# means of TOP Genera
to_save<- cbind.data.frame( "Genus"= as.data.frame(tax_table(data.genus.temp))[["Genus"]] ,
                            "Average R1"= round( as.numeric(apply(otu_table(data.genus.temp)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R1"]],1,mean)), 3),
                            "Average R2"= round( as.numeric(apply(otu_table(data.genus.temp)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="R2"]],1,mean)), 3) )# ,
                            #"Average Inoculum"= round( as.numeric(apply(otu_table(data.genus.temp)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_Code[Metadata$Reactor=="Inoculum"]],1,mean)), 3)

to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_accum_PHA_genera_Averages.xlsx", row.names = F, to_save)




##################### R AND PACKAGES VERSION #########################

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
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))


