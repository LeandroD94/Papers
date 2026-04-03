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
fill_color_19<-c("darkblue","indianred","springgreen1","wheat","black","coral1","yellow","darkmagenta","cyan","gray","firebrick3", "royalblue3","gold3","darkgreen","violet","chocolate3", "deepskyblue1","red","chartreuse3","darkslategray3")
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

# Subsetting ...
# This sample is not included in the current experiment...
Metadata <- Metadata[Metadata$Sample_name!="PN11" , ]
data <- subset_samples(data, Sample_name!="PN11")

# re-ordering the factors
sample_data(data)$Sample_name <- as.factor( sample_data(data)$Sample_name ) 
sample_data(data)$Experiment_Day <- factor( sample_data(data)$Experiment_Day , c(192,287) )

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
# ggplot(data=table, aes(x=Experiment_Day_16S, y=Abundance, fill=as.factor(Domain))) +
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
#   labs(x="Experiment Day", y="Percentual abundance of clades",
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

color_domain <- c("Bacteria"="deepskyblue",
                  "Archaea"="blue2",
                  "Amoeba"="yellow",
                  "Fungi"="darkgreen",
                  "Algae"="springgreen2",
                  "Nematoda"="gray75",
                  "Platyhelminthes"="coral1",
                  "Rotifera"="red1",
                  "Protozoa"="orange",
                  "Protist"="black",
                  "Viruses"="violet")

plot_domain <- ggplot(data=table, aes(x=Experiment_Day, y=Abundance, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  # facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =10.5) +
  scale_fill_manual(values=color_domain ) +
  theme(axis.text.x=element_text(angle=0, vjust=1, hjust=0.5, size= 9.5),
        axis.text.y=element_text(size= 7.5),
        axis.title.y=element_text(size= 13.5, hjust = 0.55, vjust=-2.2),
        # axis.title.x = element_blank(),
        axis.title.x = element_text(size= 12.5, vjust=4, hjust=0.315),
        plot.caption =element_text(size= 4),
        # strip.text = element_text(size=10),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.15, "cm"),
        legend.text = element_text ( size = 9.5, hjust=0 ),
        # legend.position="bottom",
        legend.position="right",
        legend.key.spacing.y = unit(0.25,"cm"),
        # legend.margin = margin(10,25,-15,10),
        legend.margin = margin(-14,-10,0,0),
        # plot.margin = margin(0,5,0,1),
  ) +
  scale_y_break(c(0.16, 97), 
                expand=c(0,0.001),
                scales=0.1,
                space = 0.4) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    #breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5, 97.5, 100),
    breaks = c(0, seq(0.01, 0.05, 0.01), seq(0.06, 0.18, 0.02), 97.5, 100),
    limits = c(0,100)
  ) +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Domain percent abundances",
       fill="",
       x ="Experiment day",
       # title = "Domains proportions",
       caption= "NB: the Y axis is scaled to improve the readibility of lower values"
  )
plot_domain
ggsave(file="Data_check/Domains_according_to_Kaiju.png",width=4,height=3.8, dpi=300)
ggsave(file="Results_DNAseq/Abundances/Domains_according_to_Kaiju.png",width=4.5,height=4, dpi=300)
dev.off()
  
# stesso, ma legenda contrario (per ritagli)
plot_domain +
  theme( legend.text = element_text(hjust=1))
ggsave(file="Results_DNAseq/Abundances/Domains_according_to_Kaiju2.png",width=4.5,height=4, dpi=300)
dev.off()

rm(plot_domain)


# stesso, ma inversione delle proporzioni
table$Experiment_Day2 <- factor( table$Experiment_Day, levels = c("287","192") )
ggplot(data=table, aes(x=Abundance, y=Experiment_Day2, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  # facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =10.5) +
  scale_fill_manual(values=color_domain ) +
  theme(axis.text.y=element_text(angle=0, size= 9.5),
        axis.text.x=element_text(size= 7.5, angle=25, vjust=1, hjust = 1),
        axis.title.y=element_text(size= 10.5, hjust = 0.65, vjust=-1),
        # axis.title.x = element_blank(),
        axis.title.x = element_text(size= 11.5, vjust=20.5, hjust=0.55),
        plot.caption =element_text(size= 4),
        # strip.text = element_text(size=10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.text = element_text ( size = 8 ),
        legend.position="bottom",
        legend.key.spacing.x = unit(0.2,"cm"),
        legend.margin = margin(12,32,-15,0),
        panel.grid.major.x = element_line(linewidth=0.15, color="darkgray"),
        panel.grid.minor.x = element_line(linewidth=0.1, color="gray80")
        ) +
  scale_x_break(c(0.152, 97), 
                expand=c(0,0.001),
                scales=0.1,
                space = 0.4) +
  scale_x_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    #breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5, 97.5, 100),
    # breaks = c(0, seq(0.01, 0.05, 0.01), seq(0.06, 0.18, 0.02), 97.5, 100),
    breaks = c(0, seq(0.01, 0.05, 0.01), 0.08, 0.12, 0.15, 97.5, 100),
    limits = c(0,100)
  ) +
  guides(fill=guide_legend(nrow=2, reverse=T )) +
  labs(x="Microbes proportions",
       fill="",
       y ="Experiment Day",
       # title = "Domains proportions"
  )
ggsave(file="Results_DNAseq/Abundances/Domains_according_to_Kaiju3.png",width=4.5,height=2.9, dpi=300)




############################ RAREFACTION ANALYSIS ################################

png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=5000,label=T) # using the filtered (analysed) data
dev.off()
rm(r)

png(file="Data_check/Rarefaction_curve_NoLabels.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=5000,label=F) # using the filtered (analysed) data
dev.off()
rm(r)



###################### ABUNDANCES BAR PLOT ##########################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.genus.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top_data <- subset_taxa(data_temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE)) [1:19]
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
fill_color_modified<-c("gray25","orange","aquamarine4","darkmagenta","violet","black","red3","gray70","cyan","slateblue1","yellow","gold3","springgreen1", "darkgreen","orangered1", "blue2","deeppink", "chartreuse3", "royalblue1","darkslategray3")
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Genus)) +
  # facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_modified) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 8.5
                                 ),
        axis.title.y =element_text(size=9),
        axis.title.x =element_text(size=10),
        # axis.title.x = element_text(vjust=2.5),
        # strip.text = element_text(size=10.5),
        legend.key.height = unit(0.125, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text ( size = 8.25 , face="italic" ),
        legend.position="bottom",
        legend.margin = margin(-1,34,0,0),
        plot.margin = margin(2,4,2,4)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment Day", y="Percentual abundance of clades",
       fill="",
       # title = "Most abundant identified microbial genera", 
       caption = " 'Others' includes every genus below rank 29"
       )
ggsave(file="Results_DNAseq/Abundances/TOP_microbial_genera.png",width=3.5,height=3.8, dpi=300)
dev.off()

# means of TOP Genera
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Average"= round( as.numeric(apply(otu_table(prune.dat_top) ,1,mean)), 2)
                            )
to_save<-to_save[order(to_save$Average, decreasing=T), ]
write.xlsx(file = "Results_DNAseq/Abundances/TOP_Genera_Average_abundances.xlsx", row.names = F, to_save)



### TOP Species
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top_data <- subset_taxa(data_temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:20]
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
  tabella<-aggregate(.~Species+Experiment_Day, tabella[ , c("Abundance","Species","Experiment_Day")], FUN=sum)
  tabella$Species<-gsub ("Candidatus ","Ca.", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("Betaproteobacteria bacterium ","", tabella$Species)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Species)) +
  # facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=c("gray30",fill_color_19)) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 7
  ),
  axis.title.y =element_text(size=8),
  axis.title.x =element_text(size=8),
  # strip.text = element_text(size=10),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.1, "cm"),
  legend.text = element_text ( size = 5.5 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-2,30,0,0),
  plot.margin = margin(2,4,2,4)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment Day", y="Percentual abundance of clades", 
       fill="",
       # title = "Most abundant identified microbial species", 
       caption = " 'Others' includes every species below rank 29")
ggsave(file="Results_DNAseq/Abundances/TOP_microbial_Species.png",width=3.5,height=3.8, dpi=300)
dev.off()

# means of TOP Species
to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]] ,
                            "Day192"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="192"))), 4),
                            "Day287"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="287"))), 4)
                            
)

to_save<-to_save[order(to_save$Species, decreasing=F), ]
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
levels(tabella$Genus) <- gsub("Protozoa: Symbiodinium","Symbiodinium (*)", levels(tabella$Genus), fixed = T)
fill_color_modified<-c("chocolate4","brown","chocolate","coral3","green1","red4","green3","lightgrey","forestgreen","olivedrab1","gray50","springgreen2", "olivedrab", "gray68","coral","black","gray80", "orangered", "darkgreen","darkslategray3")
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=as.factor(Genus))) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9.5) +
  scale_fill_manual(values=fill_color_modified) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 8
  ),
  axis.title.y =element_text(size=8.5),
  axis.title.x =element_text(size=9),
  # axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=10.5),
  legend.key.height = unit(0.24, "cm"),
  legend.key.width = unit(0.19, "cm"),
  legend.text = element_text ( size = 6.4 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-2,33,0,0),
  plot.margin = margin(2,4,2,4)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment Day", y="Percentual abundance of eukaryotes",
       fill="",
       # title = "Most abundant identified microbial genera", 
       caption = " 'Others' includes every eukaryote below rank 19"
  )
ggsave(file="Results_DNAseq/Abundances/TOP_Eukaryota_genera_ONLY_EUKAR.png",width=3.5,height=3.8, dpi=300)
dev.off()

# means of TOP Genera
to_save<- cbind.data.frame( "Euk_clade"= as.data.frame(tax_table(prune.dat_top))[["Euk_clade"]] ,
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Day192"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="192"))), 4),
                            "Day287"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="287"))), 4)
                            
)
to_save<-to_save[order(to_save$Euk_clade, decreasing=F), ]
write.csv(to_save, file = "Results_DNAseq/Abundances/TOP_Eukaryota_genera_ONLY_EUKAR_average_abundances.csv", row.names = F, quote=F)




###################### EXTRA CHECK: ABUNDANCES OF ARCHAEA ##########################

### TOP archaea species
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data.genus, Domain=="Archaea")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus ","Ca.", tabella$Genus)  # to save horizontal space in the plot!
  tabella$Genus<-gsub ("uncultured","uncult.", tabella$Genus)
  tabella$Genus<-gsub ("^archaeon$","unknown Archaea", tabella$Genus)
  tabella$Genus<-gsub (" archaeon "," ", tabella$Genus)
  tabella<-aggregate(.~Genus+Experiment_Day, tabella[ , c("Abundance","Genus","Experiment_Day")], FUN=sum)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=as.factor(Genus))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =10) +
  scale_fill_manual(values=fill_color_19) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 9
  ),
  axis.title.y =element_text(size=8.5),
  axis.title.x =element_text(size=9),
  # axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=10.5),
  legend.key.height = unit(0.24, "cm"),
  legend.key.width = unit(0.19, "cm"),
  legend.text = element_text ( size = 6.9 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-2,33,0,0),
  plot.margin = margin(2,4,2,4)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment Day", y="Clades percentages (only archaea)", 
       fill="",
       caption = " 'Others' includes every archaea below rank 19, no unclassified")
ggsave(file="Results_DNAseq/Abundances/TOP_Archaea_GENERA_ONLY_ARCHAEA.png",width=3.5,height=3.8, dpi=300)
dev.off()

to_save<- cbind.data.frame( "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]] , 
                            "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Day192"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="192"))), 4),
                            "Day287"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="287"))), 4)
)
to_save<-to_save[order(to_save$Day192, decreasing=T), ]
write.xlsx(file = "Results_DNAseq/Abundances/TOP_Archaea_Average_abundances.xlsx", row.names = F, to_save)




###################### EXTRA CHECK: ABUNDANCES OF VIRUSES ##########################

### TOP Viruses (species, focus on ONLY viruses)
suppressWarnings(rm(top, others, tabella, unass_data))
{
  data.subset <- subset_taxa(data, Domain=="Viruses")
  # top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  # prune.dat_top <- prune_taxa(top,data.subset)
  # View(tax_table(prune.dat_top))  # check row names
  # View(Tax_IDs) # check row name --> take the corresponding tax ID
  # solving duplicates
  # tax_table(data.subset)["sp1570", "Domain_species"]<- "Viruses_Caudoviricetes 2788787"
  # tax_table(data.subset)["sp9099", "Domain_species"]<- "Viruses_Caudoviricetes 2731619"
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  top <- names(sort(taxa_sums(data.subset), decreasing=TRUE))[1:30]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset) # only viruses
  tabella_top<-psmelt(prune.dat_top)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Species<-"Others"
  tabella_others$Domain_species<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Domain_species <- gsub ("Viruses_","", tabella$Domain_species)
  tabella$Domain_species <- gsub ("340016","Uncultured virus 340016", tabella$Domain_species)
  tabella$Domain_species<-factor(tabella$Domain_species, levels = c(unique(tabella$Domain_species)[! unique(tabella$Domain_species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Domain_species)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 9
  ),
  axis.title =element_text(size=9.5),
  plot.title =element_text(size=7),
  # axis.title.x = element_text(vjust=1),
  # strip.text = element_text(size=11),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.22, "cm"),
  legend.text = element_text ( size = 7.5 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-2,40,0,0),
  plot.margin = margin(2,4,2,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="Experiment Day", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified viruses species (focus on viruses only)",
       caption = " 'Others' includes every virus below rank 29")
ggsave(file="Results_DNAseq/Abundances/TOP_DNAVirus_species_ONLY_VIRUS.png",
       width=5.2,height=4.5, dpi=300)
dev.off()


# # means of TOP virus
# to_save<- cbind.data.frame( "Species"= as.data.frame(tax_table(prune.dat_top))[["Domain_species"]] ,
#                             "Average R1"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R1"]],1,mean)), 2),
#                             "Average R2"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="R2"]],1,mean)), 2),
#                             "Average Inoculum"= round( as.numeric(apply(otu_table(prune.dat_top)[ ,colnames(otu_table(prune.dat_top))%in%Metadata$FASTQ_ID[Metadata$Group=="Inoculum"]],1,mean)), 2)
# )
# to_save<-to_save[order(to_save$`Average R1`, decreasing=T), ]
# write.xlsx(file = "Results_DNAseq/Abundances/TOP_Viruses_Averages.xlsx", row.names = F, to_save)





######################## ALPHA DIVERSITY ##############################

### PROKARYOTES (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Bacteria") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Experiment_Day")
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
  geom_boxplot(data=pAlpha$data, aes(x=Experiment_Day, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Experiment_Day), color="red2" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Eubacteria genera)", 
       y="Alpha Diversity Measure",
       x="Experiment Day") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Genera_Eubacteria.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### PROKARYOTES (SPECIES)
data.temp<-subset_taxa(data, Domain %in% c("Bacteria") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Experiment_Day")
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
  geom_boxplot(data=pAlpha$data, aes(x=Experiment_Day, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Experiment_Day), color="red2" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Eubacteria species)", 
       y="Alpha Diversity Measure",
       x="Experiment Day") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Species_Eubacteria.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### Archaea (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Archaea") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Experiment_Day")
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
  geom_boxplot(data=pAlpha$data, aes(x=Experiment_Day, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Experiment_Day), color="red4" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Archaea genera)", 
       y="Alpha Diversity Measure",
       x="Experiment Day") +
  guides(fill="none", color="none", shape="none") + 
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results_DNAseq/Alpha_diversity/Alfa_diversity_on_Genera_Archaea.png", width = 5.5,height =4.5, dpi=300)




### Eukaryota (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Eukaryota") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Experiment_Day")
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
  geom_boxplot(data=pAlpha$data, aes(x=Experiment_Day, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  geom_text(aes(label=Experiment_Day), color="chocolate3" , size=4)+
  # scale_color_manual(values = c("BBBB"="chartreuse", "AAAA"="coral")) +
  labs(title="Alpha diversity (Eukaryota genera)", 
       y="Alpha Diversity Measure",
       x="Experiment Day") +
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
              "Nitrosopumilus", # https://doi.org/10.1002/9781118960608.gbm01290
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
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:19]
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
tabella$Genus[tabella$Genus2%in%AOB_list] <- paste(  tabella$Genus[tabella$Genus2%in%AOB_list] , "(AOB)")
tabella$Genus[tabella$Genus2%in%AOA_list] <- paste(  tabella$Genus[tabella$Genus2%in%AOA_list], "(AOA)")
tabella$Genus[tabella$Genus2%in%NOB_list] <- paste(  tabella$Genus[tabella$Genus2%in%NOB_list], "(NOB)")
# ordering levels
tabella$Genus<-factor(tabella$Genus, levels = c( unique(tabella$Genus)[! unique(tabella$Genus) %in% "Other AOOs or NOOs"] ,
                                                 "Other AOOs or NOOs"
                                                 )
                      )

base_plot <- ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Genus)) +
  # facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_y_continuous(expand=c(0,0.3) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 8.5
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2),
  strip.text = element_text(size=9.5),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.35, "cm"),
  legend.text = element_text ( size = 9 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-3,30,0,0),
  plot.margin = margin(4,4,4,10)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment Day", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified AOB,AOA or NOB organisms", 
       caption = " 'Other' includes every accumulator below rank 20")

custom_colors<-c("navyblue","gray85","green","green3","forestgreen","darkmagenta",
                 "pink3","gray35","orange","tomato1","firebrick","red","indianred","blue1","black",
                 "red3","yellow","orangered1","royalblue2",
                 "darkslategray3")
base_plot +
  scale_fill_manual(values=custom_colors)
ggsave(file="Results_DNAseq/Abundances/AOB_AOA_NOB_organisms_OneColorForEach.png",width=4.85,height=4.5, dpi=300)



### Again, but only NOBs in proportion
tabella2 <- tabella[grep("(NOB)", tabella$Genus, fixed=T)|grep("Nitrosomonas", tabella$Genus, fixed=T), ]
levels(tabella2$Genus) <- gsub(" (NOB)", "", levels(tabella2$Genus), fixed=T)
levels(tabella2$Genus) <- gsub(" (AOB)", "", levels(tabella2$Genus), fixed=T)
levels(tabella2$Genus)[ !levels(tabella2$Genus) %in% c("Nitrospira","Ca.Nitrotoga","Nitrobacter","Nitrosomonas") ] <- "Other Nitrifiers"
base_plot <- ggplot(data=tabella2, aes(x=Experiment_Day, y=Abundance, fill=Genus)) +
  # facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =10) +
  scale_y_continuous(expand=c(0,0.05) , breaks = seq(0,14,2) , limits=c(0,14) ) +
  theme(axis.text.x=element_text(angle=0,
                                 vjust=1,
                                 hjust=0.5,
                                 size= 9
  ),
  axis.title.y =element_text(size=9),
  axis.title.x = element_text(vjust=0.9, size=10),
  # strip.text = element_text(size=9.5),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.35, "cm"),
  legend.text = element_text ( size = 10 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-2,32,0,0),
  panel.grid.major.y = element_line(linewidth=0.15, color="darkgray"),
  panel.grid.minor.y = element_line(linewidth=0.1, color="gray80"),
  plot.margin = margin(4,1,4,4)) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment Day", y="Percent abundances of Nitrifiers\nin total DNA-seq",
       fill="")
custom_colors<-c("navyblue", "red3","yellow2","orangered1","royalblue")
base_plot +
  scale_fill_manual(values=custom_colors)
ggsave(file="Results_DNAseq/Abundances/AOB_NOB_Focus_OneColorForEach_V1.png",width=3.5,height=3.8, dpi=300)
custom_colors<-c("navyblue", "firebrick3","red","orangered1","chartreuse3")
base_plot +
  scale_fill_manual(values=custom_colors)
ggsave(file="Results_DNAseq/Abundances/AOB_NOB_Focus_OneColorForEach_V2.png",width=3.5,height=3.8, dpi=300)


# Flipped!
tabella2$Experiment_Day2 <- factor( tabella2$Experiment_Day, levels = c("287","192") )
ggplot(data=tabella2, aes(y=Experiment_Day2, x=Abundance, fill=Genus)) +
  # facet_grid(cols= vars(Group),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =10) +
  scale_x_continuous(expand=c(0,0.05) , breaks = seq(0,16,2) , limits=c(0,16) ) +
  theme(axis.text.y=element_text(angle=0,
                                 size= 9.2
                                 ),
        axis.text.x=element_text(size= 9, angle=30, vjust=1, hjust = 1),
  axis.title.x =element_text(size=10),
  axis.title.y = element_text(vjust=4, size=10.1),
  # strip.text = element_text(size=9.5),
  legend.key.height = unit(0.5, "cm"),
  legend.key.width = unit(0.2, "cm"),
  legend.text = element_text ( size = 8.6 ),
  legend.position="bottom",
  legend.key.spacing.x = unit(0.20,"cm"),
  legend.margin = margin(0,35,0,0),
  panel.grid.major.x = element_line(linewidth=0.15, color="darkgray"),
  panel.grid.minor.x = element_line(linewidth=0.1, color="gray80"),
  plot.margin = margin(4,5,4,10)) +
  guides(fill=guide_legend(nrow=1 , reverse=T),  ) +
  labs(y="Experiment Day", x="Percent abundances of Nitrifiers in total DNA-seq",
       fill="") +
  scale_fill_manual(values=c("navyblue", "red3","yellow2","orangered1","royalblue"))
ggsave(file="Results_DNAseq/Abundances/AOB_NOB_Focus_OneColorForEach_V3.png",width=4.5,height=2.3, dpi=300)



# averages
to_save<- cbind.data.frame( "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ,
                            "Day192"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="192"))), 4),
                            "Day287"= round( as.numeric( otu_table(subset_samples(prune.dat_top,Experiment_Day=="287"))), 4)
                            
)
to_save<-to_save[order(to_save$Day192, decreasing=T), ]
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


