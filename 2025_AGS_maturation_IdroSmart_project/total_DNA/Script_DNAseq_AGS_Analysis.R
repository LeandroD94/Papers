# ssh -X matteo@150.217.62.209
# ssh -X matteo@150.217.62.222
# cd /home/matteo/Desktop/Works_in_progress/TUDelft_DNAShotgun_check/
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
  dir.create("Results/Correlations")
}

options(scipen = 100) # disable scientific annotation


### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
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

load(file = "Phyloseq_object_for_R.RData") # from converting Kraken2 output (Species names in a sample list) to phyloseq output (taxonomy, feature table)
# using the IDs as names
sample<-sample_names(data)
original_names<-sample
sample

sample<-gsub("_bracken_abund.tsv","",sample)
sample_names(data)<-sample # update

Metadata <- as.data.frame(read.table(file = "metadata.txt", header = T, sep="\t"))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)


sample_data(data)$Experiment_day <-factor( sample_data(data)$Experiment_day , 
                                           levels= c("1","57","115","164","other_reactor") 
                                           )
sample_names(data)<-sample_data(data)$Experiment_day

if(identical(original_names,sample_names(data_NO_bracken))){
  sample_names(data_NO_bracken)<-sample_names(data)
  sample_data(data_NO_bracken)<-sample_data(data)
} else { stop("\n\nWait! The sample names of the two objects do not match\n")}




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

write.csv2(cbind(otu_table(data),tax_table(data)), file="Data_check/Raw_feature_Table_pre_filtering.csv", row.names = T)


# data_NO_bracken<- subset_taxa(data_NO_bracken, Genus!="Unclassified")  # if NA

## using the original counts from kraken2, being the TRUE sequencing counts (--> depths)
# data.temp<-transform_sample_counts(data_NO_bracken, function(x) (x/sum(x))*100)
# filtered<-taxa_names(filter_taxa(data.temp, function(x) max(x) <= 0.0001, TRUE))
## !!! NB: the percentage filter does not work here, it discard many biologically logical taxa, then using a more raw filter...

filtered<-taxa_names(filter_taxa(data_NO_bracken, function(x) max(x) < 3, TRUE))
write.csv( cbind(as.data.frame(tax_table(data_NO_bracken))[filtered, c("Domain","Phylum","Species")], as.data.frame(otu_table(data_NO_bracken))[filtered, ] ), 
           file="Data_check/Filtered_species_UNDER_max_3reads_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data_NO_bracken, function(x) max(x) < 5, TRUE)))[["Species"]]
filtered
data<-subset_taxa(data, ! Species %in% filtered) # filtering the object using in analyses
data<-prune_taxa(taxa_sums(data) > 0, data) 
suppressWarnings( rm(filtered, data.temp) )


# some Viridiplant may be present in trace due to homonimy with bacteria (see the Converter building script) then removing them from the objects...
data<-subset_taxa(data, Phylum!="Streptophyta")
data_NO_bracken<-subset_taxa(data_NO_bracken, Phylum!="Streptophyta")
# same for insects, Arthropoda and other Metazoa...
metazoa_phyla<- unique( as.data.frame(tax_table(data_NO_bracken))[as.data.frame(tax_table(data_NO_bracken))[["Kingdom"]]=="Metazoa","Phylum"] )
if(length(metazoa_phyla)>1){ # this is to avoid errors caused by re-launching portions of the script
  data<-subset_taxa(data, ! Phylum %in% metazoa_phyla )
  data_NO_bracken<-subset_taxa(data_NO_bracken,  ! Phylum %in% metazoa_phyla)
}
# if there are plasmid, they have to be removed from the classification (no taxonomic informations)
data<-subset_taxa(data, ! grepl("Plasmid", Species))
data_NO_bracken<-subset_taxa(data_NO_bracken, ! grepl("Plasmid", Species))
# few reads may be assigned to "Metagenome" labeled entries with no taxonomy at all --> removing them
data<-subset_taxa(data, ! grepl("metagenome", Species) )
# in case of data NO bracken, these metagenome assigned reads will be treated as unclassified
temp<-as.data.frame(tax_table(data_NO_bracken))
temp[grepl("metagenome",temp$Species) , "Domain"] <- "Unclassified"
tax_table(data_NO_bracken)<-as.matrix(temp)
rm(temp)


proof1<- "Marker of the filtering, it is required for the script"




############################# UNCLASSIFIED CHECK #################################

data_temp<-tax_glom(data_NO_bracken, taxrank = "Domain", NArm = F) # without glomming before, some error may be seen in ggplot2 (maybe bugs due to too many observations?)
data_NO_b_prop<-transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )
# colSums(otu_table(data_NO_b_prop))
#View(cbind(otu_table(data_temp),tax_table(data_temp)))
table<-psmelt(data_NO_b_prop)
table$Domain<-factor(table$Domain, levels=c( unique(table$Domain)[unique(table$Domain)!="Unclassified"] , "Unclassified" ) )

ggplot(data=table, aes(x=Sample, y=Abundance, fill=Domain)) +
  # facet_grid(cols= vars(Project),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =9) +
  scale_fill_manual(values=c("Unclassified"="lightgray",
                             "Bacteria"="deepskyblue",
                             "Archaea"="chartreuse4",
                             #"Eukaryota_Fungi"="gold3",
                             #"Eukaryota_Metazoa"="gold",
                             #"Eukaryota_Protists"="yellow",
                             "Eukaryota"="yellow",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust=1, size= 6),
        strip.text = element_text(size=10.2),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom",
        legend.margin = margin(-2,30,4,0),
        plot.margin = margin(4,4,4,15),
        ) +
  #scale_y_sqrt(breaks=c(0,5,10,25,50,100)) + # can't be applied due to the values below 1
  guides(fill=guide_legend(nrow=1)) +
  labs(x="Samples", y="Percentual abundance of clades",
       fill="",
       title = "Domains and unclassified proportions \n(raw Kraken2 data, after host decontamination)")
ggsave(file="Data_check/Domains_and_unclassified_according_to_rawKraken2.png",width=6,height=4.8, dpi=450)
dev.off()




###################### DOMAINS PROPORTIONS #########################

data_temp<-tax_glom(data_NO_bracken, taxrank = "Kingdom", NArm = F) # without glomming before, some error may be seen in ggplot2 (maybe bugs due to too many observations?)

#data_temp<-subset_taxa(data_temp, Domain!="Unclassified" )
# !!! NB: due to a strange glitch of gglopt2, the "Bacteria" color in the stacked bar plot became invisible in few sample if the Domain column is modified in a factor with the Bacteria as first level (as below)
# therefore, ggplot2 will be tricked: the "Unclassified" will still be present but their abundance will be very minimal, in order to put them as "first place" (impossible to see) and the bacteria as second level (which actually seems the first one)
otu_table(data_temp)[ tax_table(data_temp)[, "Domain"]=="Unclassified" ] <- 0.00001

data_temp<-transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )
table<-psmelt(data_temp)
which_eukar<- paste( table$Domain[table$Domain=="Eukaryota"], table$Kingdom[table$Domain=="Eukaryota"] , sep="_")
table$Domain[table$Domain=="Eukaryota"]<-which_eukar # this new column can't be updated during its own creation, otherwise NAs will be obtained
table$Domain[table$Domain=="Eukaryota_none"]<-"Eukaryota_Protists" # if the kingdom is "none" then they are almost always protists 

table$Domain<-factor(table$Domain, levels =  c("Unclassified","Bacteria","Archaea","Eukaryota_Metazoa","Eukaryota_Protists","Viruses") ) # NB: unclassified are invisible (see above)
#table$Domain<-relevel(table$Domain,c("Unclassified","Bacteria") )

ggplot(data=table, aes(x=Experiment_day, y=Abundance, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid(~ Project_description, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("Bacteria"="deepskyblue",
                             "Archaea"="chartreuse4",
                             "Eukaryota_Fungi"="gold3",
                             "Eukaryota_Metazoa"="gold",
                             "Eukaryota_Protists"="yellow",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=32, vjust=1, hjust=1, size= 6),
        strip.text = element_text(size=6.7),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom",
        legend.margin = margin(-2,10,0,0),
        axis.title.x = element_blank(),
        # plot.margin = margin(0,5,-10,2),
        ) +
  scale_y_break(c(1.5, 97.5), scales=0.035, space = 0.3) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5, 97.5, 100),
    limits = c(0,100)
  ) +
  guides(fill=guide_legend(nrow=1)) +
  labs(y="Percentual abundance of clades",
       fill="",
       title = "Domains proportions (raw Kraken2 data)",
       caption= "NB: the Y axis is scaled to improve the readibility of lower values"
       )
ggsave(file="Data_check/Domains_according_to_rawKraken2.png",width=7,height=5.4, dpi=450)
dev.off()

# NB: if the archaea are not displayed in the plot... just repeat this part of the script.


### AGAIN, WITHOUT BACTERIA
table2<-table[table$Domain!="Bacteria", ]
ggplot(data=table2, aes(x=Experiment_day, y=Abundance, fill=as.factor(Domain))) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid(~ Project_description, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c(
                             "Archaea"="chartreuse4",
                             "Eukaryota_Fungi"="gold3",
                             "Eukaryota_Metazoa"="gold",
                             "Eukaryota_Protists"="yellow",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust=1, size= 6),
        # strip.text = element_text(size=6.7),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom",
        legend.margin = margin(-2,10,0,0),
        axis.title.x = element_blank(),
        # plot.margin = margin(0,5,-10,2),
  ) +
  # scale_y_break(c(0.22, 0.4), scales=0.035, space = 0.3) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    # breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5,2, 97.5, 100),
  ) +
  guides(fill=guide_legend(nrow=1)) +
  labs(y="Percentual abundance of clades",
       fill="",
       title = "Domains proportions (raw Kraken2 data, excluding bacteria)",
       caption= "NB: the bacteria have not been taken into account for this plot"
  )
ggsave(file="Data_check/Domains_according_to_rawKraken2_NO_BACTERIA.png",width=7,height=5.4, dpi=450)
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

####################### PREPARATION OF THE DATA #######################

if(! "proof1" %in% ls()){
  stop("\n Wait! Did you perform the filtering steps??? \n\n")
}

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  #data.class = tax_glom(data, taxrank = "Class", NArm = F)
  #data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
  #data.class.prop <- transform_sample_counts(data.class, function(x) x/sum(x)*100)
  #data.order.prop <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(x) x/sum(x)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  #Taxa.class<-as.data.frame(tax_table(data.class))
  #Taxa.order<-as.data.frame(tax_table(data.order))
}




########################### COUNTS EXPORT ##########################################

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/counts_species.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Relative_abundances/counts_species_prop.csv",quote=F)
  #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
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
  for(i in 1:length(tabella_top$Phylum)){ # to mark the fungal taxa
    tabella_top$Phylum[i]<-paste0(tabella_top$Phylum[i],"_",tabella_top$Domain[i])
  }
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-gsub ("Candidatus","Ca.", tabella$Phylum)
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
# tabella$Project<-gsub("AAAA"," ",tabella$Project)              # if it is needed to rename
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid(~ Project_description, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
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
        legend.key.width = unit(0.42, "cm"),
        legend.text = element_text ( size = 10 ),
        legend.position="bottom",
        legend.margin = margin(-1,31,5,0),
        plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified microbial phyla",
       fill="",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Abundances/TOP_microbial_phyla.png",width=7,height=5.4, dpi=300)
dev.off()

# means of TOP phyla
write.xlsx(file = "Results/Abundances/TOP_microb_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)




### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.genus.prop)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  for(i in 1:length(tabella_top$Genus)){ # to mark the fungal taxa
    tabella_top$Genus[i]<-paste0(tabella_top$Genus[i],"_",tabella_top$Domain[i])
  }
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus","Ca.", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  # facet_grid(cols= vars(Project_description),scales = "free_x", space = "free_x") +
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
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text ( size = 9.1 ),
        legend.position="bottom",
        legend.margin = margin(-6,41,0,0),
        plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       fill="",
       title = "Most abundant identified microbial genera", 
       caption = " 'Others' includes every genus below rank 29")
ggsave(file="Results/Abundances/TOP_microbial_genera.png",width=7,height=5.4, dpi=300)
dev.off()

# means of TOP genera
write.xlsx(file = "Results/Abundances/TOP_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



### TOP Species
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.prop)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Species)){ # to mark the fungal taxa
  #   tabella_top$Species[i]<-paste0(tabella_top$Species[i],"_",tabella_top$Domain[i])
  # }
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Species<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the species
  tabella<-aggregate(.~Species+Experiment_day, tabella[ , c("Experiment_day","Abundance","Species")], FUN=sum)
  # tabella$Species<-gsub ("Candidatus","Ca.", tabella$Species)
  tabella$Species<-gsub ("Candidatus ","", tabella$Species)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
levels(tabella$Experiment_day)<-gsub("other_reactor","other_AGS",levels(tabella$Experiment_day))

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Species)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=6.7),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.text = element_text ( size = 8.8 ),
  legend.position="bottom",
  legend.margin = margin(-4,41,0,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="Experiment day", y="Percentual abundance of clades", 
       fill="",
       title = "Most abundant identified microbial species", 
       caption = " 'Others' includes every species below rank 29")
ggsave(file="Results/Abundances/TOP_microbial_Species.png",width=7,height=5.4, dpi=300)
dev.off()

# means of TOP Species
to_save<- cbind.data.frame( "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]] , 
                            "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]] ,
                            "Average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2),
                            "Average in mature pilot"= round( as.numeric(apply( otu_table(subset_samples(prune.dat_top,Experiment_day%in%c("115","164"))) ,1,mean)) ,2)
)
to_save<-to_save[order(to_save$Average, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_Species_Average_abundances.xlsx", row.names = F, to_save)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




#################### EXTRA CHECK: FOCUS ON CA. ACCUMULIBACTER #########################

colors_for_accumuli<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","darkslategray3")

suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.prop, Genus=="Candidatus Accumulibacter")
  prune.data.others<-subset_taxa(data.prop, Genus!="Candidatus Accumulibacter")
  tabella_top<-psmelt(top_data)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Species<-"Others (Not Accumulibacter)"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the species
  tabella<-aggregate(.~Species+Experiment_day, tabella[ , c("Experiment_day","Abundance","Species")], FUN=sum)
  tabella$Species<-gsub ("Candidatus","Ca.", tabella$Species)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others (Not Accumulibacter)"],"Others (Not Accumulibacter)"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Species)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_y_continuous(expand=c(0,1) ) +
  scale_fill_manual(values=colors_for_accumuli) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=3),
  # strip.text = element_text(size=6.7),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.42, "cm"),
  legend.text = element_text ( size = 9.8 ),
  legend.position="bottom",
  legend.margin = margin(-1,31,6,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="Experiment day", y="Percentual abundance of clades", 
       fill="",
       title = "Focus on Candidatus Accumulibacter species", 
       caption = " 'Others' includes every other species besides C. Accumulibacter ones")
ggsave(file="Results/Abundances/Ca_Accumulibacter_Species.png",width=7,height=5.4, dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




###################### EXTRA CHECK: ABUNDANCES OF EUKARYOTA ##########################


### TOP eukar genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ 
  data.subset <- subset_taxa(data, Domain=="Eukaryota")
  data.subset <- tax_glom(data.subset, taxrank = "Genus") # aggregation to sum to genus level (too many repetitions otherwise)
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  top_data <- subset_taxa(data.subset, Genus!="none") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset)
  tabella_top<-psmelt(prune.dat_top)
  for(i in 1:length(tabella_top$Genus)){ # to mark the fungal taxa
    tabella_top$Genus[i]<-paste0(tabella_top$Genus[i],"_",tabella_top$Phylum[i])
  }
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("_none","_eukar", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
levels(tabella$Experiment_day)<-gsub("other_reactor","other_AGS",levels(tabella$Experiment_day))
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
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
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.text = element_text ( size = 8.7 ),
  legend.position="bottom",
  legend.margin = margin(-6,42,0,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="Experiment day", y="Percentual abundance of clades (only eukaryotes)", 
       fill="",
       title = "Most abundant identified eukaryota genera", 
       caption = " 'Others' includes every eukaryota below rank 29")
ggsave(file="Results/Abundances/TOP_Eukaryota_genra_ONLY_EUKAR.png",width=7,height=5.4, dpi=300)
dev.off()




###################### EXTRA CHECK: ABUNDANCES OF ARCHAEA ##########################


### TOP eukar species
suppressWarnings(rm(top, others, tabella, unass_data))
{ 
  data.subset <- subset_taxa(data, Domain=="Archaea")
  # data.subset <- tax_glom(data.subset, taxrank = "Species") # aggregation to sum to Species level (too many repetitions otherwise)
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  top_data <- subset_taxa(data.subset, Species!="none") 
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
  tabella$Species<-gsub ("uncultured","uncult.", tabella$Species)
  tabella$Species<-gsub ("Methanobrevibacter sp. TLL-48-HuF1","Methanobrevibacter TLL-48-HuF1", tabella$Species)
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Species)) +
  geom_bar(stat="identity", position="stack") +
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
  legend.key.width = unit(0.3, "cm"),
  legend.text = element_text ( size = 9.2 ),
  legend.position="bottom",
  legend.margin = margin(-6,42,0,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="Experiment day", y="Percentual abundance of clades (only archaea)", 
       fill="",
       title = "Most abundant identified archaea species", 
       caption = " 'Others' includes every archaea below rank 29")
ggsave(file="Results/Abundances/TOP_Archaea_species_ONLY_ARCHAEA.png",width=7,height=5.4, dpi=300)
dev.off()




###################### EXTRA CHECK: ABUNDANCES OF VIRUSES ##########################

dir.create("Results/Abundances/ExtraCheck_using_only_Viruses")


### TOP Viruses (species, focus on ONLY viruses)
suppressWarnings(rm(top, others, tabella, unass_data))
{ 
  data.subset <- subset_taxa(data, Domain=="Viruses")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  top_data <- subset_taxa(data.subset, Species!="none") # those "none" are viruses without assigned Family in NCBI, then they are aggregate of different "NAs" from the species level (which contributes to the total in proportions anyways)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset) # only viruses
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Species<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species <- gsub ("phage","", tabella$Species)
  tabella$Species <- gsub ("Stenotrophomonas  vB_SM_ytsc_ply2008005","Stenotrophomonas vB_SM_ytsc", tabella$Species)
  tabella$Species <- gsub ("vB_SM","vB", tabella$Species)
  tabella$Species <- gsub ("vB_Sm","vB", tabella$Species)
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
levels(tabella$Experiment_day)<-gsub("other_reactor","other_AGS",levels(tabella$Experiment_day))
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Species)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=40,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=6.7),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.text = element_text ( size = 8.5 ),
  legend.position="bottom",
  legend.margin = margin(-5,42,0,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="Experiment day", y="Percentual abundance of clades", 
       fill="",
       title = "Most abundant identified viruses species (focus on viruses only)", 
       caption = " 'Others' includes every virus below rank 29")
ggsave(file="Results/Abundances/ExtraCheck_using_only_Viruses/TOP_DNAVirus_species_ONLY_VIRUS.png",
       width=7,height=5.4, dpi=300)
dev.off()



### TOP Viruses (Families, focus on ONLY viruses )
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data.fam, Domain=="Viruses")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  top_data <- subset_taxa(data.subset, Family!="none") # those "none" are viruses without assigned Family in NCBI, then they are aggregate of different "NAs" from the species level (which contributes to the total in proportions anyways)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.fam)
  others<-others[!(others %in% top)] 
  prune.data.others<-prune_taxa(others,data.subset)  # focus on viruses only!
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Family<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Family<-factor(tabella$Family, levels = c(unique(tabella$Family)[! unique(tabella$Family) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Family)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
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
  legend.key.width = unit(0.42, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin(-1,31,5,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="Experiment day", y="Percentual abundance of clades", 
       fill="",
       title = "Most abundant identified viruses families (focus on viruses only)", 
       caption = " 'Others' includes every virus below rank 29")
ggsave(file="Results/Abundances/ExtraCheck_using_only_Viruses/TOP_DNAVirus_Families_ONLY_viruses.png",width=8.2,height=5.4, dpi=300)
dev.off()



### TOP Viruses (Phyla, focus on ONLY viruses )
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data.phy, Domain=="Viruses")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  top_data <- subset_taxa(data.subset, Phylum!="none")
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.subset)
  if(length(taxa_names(data.subset))>=5) {
    others<-taxa_names(data.phy)
    others<-others[!(others %in% top)] 
    prune.data.others<-prune_taxa(others,data.subset)  # focus on viruses only!
    tabella_top<-psmelt(prune.dat_top)
    tabella_others<-psmelt(prune.data.others)
    tabella_others$Phylum<-"Others"
    tabella_others$Phylum<-"No assigned phylum"
    tabella<-rbind.data.frame(tabella_top,tabella_others)
    plot_title<-"Most abundant identified viruses phyla (focus on viruses only)"
    plot_caption<-" 'Others' includes every virus below rank 5 "
  } else {
    tabella_top<-psmelt(prune.dat_top)
    tabella<-tabella_top
    plot_title<-"Every identified viruses phyla (focus on viruses only)"
    plot_caption<-"There are no other viral phyla"
  }
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
  # tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "No assigned phylum"],"No assigned phylum"))
}
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  # facet_grid(cols= vars(Project),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
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
  legend.key.width = unit(0.42, "cm"),
  legend.text = element_text ( size = 10 ),
  legend.position="bottom",
  legend.margin = margin(-1,31,5,0),
  plot.margin = margin(4,4,4,15)) +
  guides(fill=guide_legend(nrow=1)) +
  labs(x="Experiment day",
       y="Percentual abundance of clades",
       fill= "",
       title = plot_title 
       #caption = plot_caption
       )
ggsave(file="Results/Abundances/ExtraCheck_using_only_Viruses/TOP_DNAVirus_Phyla_ONLY_viruses.png",width=8.2,height=5.4, dpi=300)
dev.off()




######################## ALPHA DIVERSITY ##############################


### PROKARYOTES (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Bacteria", "Archaea") )
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Experiment_day
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
mix_alpha$Sample_name<- gsub("X","",mix_alpha$Sample_name)  # added by the melt function
mix_alpha$Sample_name<-factor(mix_alpha$Sample_name, levels=levels(sample_data(data)$Experiment_day) )

mix_alpha$ID <- "this_reactor"
mix_alpha$ID[grepl("other_",mix_alpha$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 

ggplot(data=mix_alpha, aes(y=value, x=Sample_name, color=Time_color)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  geom_point(size =2, aes(shape=ID), color="deepskyblue2") +
  geom_point(size =3.8, alpha=0.4,  color="deepskyblue2", aes(shape=ID)) +
  #geom_line(aes(group=ID), size= 0.4, alpha= 0.4) + 
  geom_line(aes(group=as.character(ID)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 9) +
  labs(title="Alpha diversity (prokaryote genera)", y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  geom_text(aes(label=Sample_name), color="black", size=1.8, vjust=-0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_blank(),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_PROKARYOTES.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### PROKARYOTES (SPECIES)
data.temp<-subset_taxa(data, Domain %in% c("Bacteria", "Archaea") )
sample_names(data.temp)<-sample_data(data.temp)$Experiment_day
mix_alpha<- estimate_richness(data.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
mix_alpha$Sample_name<- gsub("X","",mix_alpha$Sample_name)  # added by the melt function
mix_alpha$Sample_name<-factor(mix_alpha$Sample_name, levels=levels(sample_data(data)$Experiment_day) )

mix_alpha$ID <- "this_reactor"
mix_alpha$ID[grepl("other_",mix_alpha$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 

ggplot(data=mix_alpha, aes(y=value, x=Sample_name, color=Time_color)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  geom_point(size =2, aes(shape=ID), color="deepskyblue2") +
  geom_point(size =3.8, alpha=0.4,  color="deepskyblue2", aes(shape=ID)) +
  #geom_line(aes(group=ID), size= 0.4, alpha= 0.4) + 
  geom_line(aes(group=as.character(ID)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 9) +
  labs(title="Alpha diversity (prokaryote species)", y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  geom_text(aes(label=Sample_name), color="black", size=1.8, vjust=-0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_blank(),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_PROKARYOTES.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### EUKARYOTES
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Eukaryota") )
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Experiment_day
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
mix_alpha$Sample_name<- gsub("X","",mix_alpha$Sample_name)  # added by the melt function
mix_alpha$Sample_name<-factor(mix_alpha$Sample_name, levels=levels(sample_data(data)$Experiment_day) )

mix_alpha$ID <- "this_reactor"
mix_alpha$ID[grepl("other_",mix_alpha$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 

ggplot(data=mix_alpha, aes(y=value, x=Sample_name, color=Time_color)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  geom_point(size =2, aes(shape=ID), color="deepskyblue2") +
  geom_point(size =3.8, alpha=0.4,  color="deepskyblue2", aes(shape=ID)) +
  #geom_line(aes(group=ID), size= 0.4, alpha= 0.4) + 
  geom_line(aes(group=as.character(ID)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 9) +
  labs(title="Alpha diversity (Eukaryote genera)", y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  geom_text(aes(label=Sample_name), color="black", size=1.8, vjust=-0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_blank(),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_EUKARYOTA.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### VIRUSES
data.temp<-subset_taxa(data, Domain %in% c("Viruses") )
sample_names(data.temp)<-sample_data(data.temp)$Experiment_day
mix_alpha<- estimate_richness(data.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<-(mix_alpha$Observed)/log((mix_alpha$Shannon))
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
mix_alpha$Sample_name<- gsub("X","",mix_alpha$Sample_name)  # added by the melt function
mix_alpha$Sample_name<-factor(mix_alpha$Sample_name, levels=levels(sample_data(data)$Experiment_day) )

mix_alpha$ID <- "this_reactor"
mix_alpha$ID[grepl("other_",mix_alpha$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 

ggplot(data=mix_alpha, aes(y=value, x=Sample_name, color=Time_color)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  geom_point(size =2, aes(shape=ID), color="deepskyblue2") +
  geom_point(size =3.8, alpha=0.4,  color="deepskyblue2", aes(shape=ID)) +
  #geom_line(aes(group=ID), size= 0.4, alpha= 0.4) + 
  geom_line(aes(group=as.character(ID)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 9) +
  labs(title="Alpha diversity (virus species)", y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  geom_text(aes(label=Sample_name), color="black", size=1.8, vjust=-0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_blank(),
        axis.text.y= element_text(angle=0, size=5.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_Viruses.png", width = 5.5,height =4.5, dpi=300)

suppressWarnings(rm(mix_alpha) )




########################### PCoA AND BETADIVERSITY ########################

data.prop_bacteria_Genus<-subset_taxa(data.genus.prop, Domain %in% c("Bacteria","Archea"))
data.prop_bacteria_Species<-subset_taxa(data.prop, Domain %in% c("Bacteria","Archea"))
data.prop_eukar_Genus<-subset_taxa(data.genus.prop, Domain=="Eukaryota") # to focus only on fungi
data.prop_eukar_Species<-subset_taxa(data.prop, Domain=="Eukaryota")
data.prop_virus_Genus<-subset_taxa(data.genus, Domain=="Viruses")
data.prop_virus_Species<-subset_taxa(data.prop, Domain=="Viruses")

for( d in list(data.genus.prop, data.prop, 
               data.prop_bacteria_Genus, data.prop_bacteria_Species, #, 
               data.prop_eukar_Genus, data.prop_eukar_Species,
               data.prop_virus_Genus, data.prop_virus_Species
               ) 
     ) {
  
  # setting the automatic labels and file names
  if(identical(d,data.genus.prop) | 
     identical(d,data.prop_bacteria_Genus) | 
     identical(d,data.prop_eukar_Genus) |
     identical(d,data.prop_virus_Genus) 
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
if(identical(d,data.prop_virus_Genus) | identical(d,data.prop_virus_Species)){
  #subset<-"PCoA_with_virus_only"
  subset<-"PCoA"
  file_name<-"Viruses"
  plot_title<-"identified viruses"
}
  
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
sample_data(data.prop.labels)$ID <- "this_reactor"
sample_data(data.prop.labels)$ID[grepl("other_",sample_data(data.prop.labels)$Sample_name)] <- "z_other" # the z is to allow ggplot to automatically treat it as the second level 


{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC , shape= "ID") +
  # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
  geom_point(size=1.5, color="deepskyblue2") +
  geom_point(size=4.5, alpha= 0.6, color="deepskyblue2") +
  theme_classic(base_size = 9) +
  #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
  geom_text(aes(label=sample_data(data.prop.labels)$Experiment_day), color="black", size=2, 
            hjust=0.7, vjust=-0.1, show.legend = FALSE) +
  theme(legend.margin = margin(10,0,10,0),
        legend.text = element_text(size=7.8)) +
  guides(shape="none") + 
  labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
       shape="",
       caption=" ",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file=paste0("Results/",subset,"/PCoA_with_ID_Hellinger_on_",file_name,"_",levels,"_Description.png"), width = 5.55, height = 4.2, dpi=300)


# # with ellipses
# plot_ordination(data.sqrt_prop, ordBC , color="Project_description2") +
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


#################### CORRELATING ABUNDANCES vs TIME (CENTERED LOG RATIO) ########################

# selecting only genera with at least 0.01% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
data.genus.temp <- subset_samples(data.genus.temp, Sample_name != "other_reactor") 
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.01, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
#data.genus.temp <- subset_samples(data.genus, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
data.genus.temp<- data.genus
data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
abundances<-abundances[ , c("1","57","115","164")] # these are the exp days, aka the ordered sample names here


Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)),
                 method = "kendall") # correlated with the flow of time (the samples have been ordered accordingly)
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
# Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")

write.csv2(Corr_results, "Results/Correlations/Kendall_correlation_CLR_bacteria_abundances_with_time.csv", row.names = T, quote = F)

# nothing here...
con <- file("Results/Correlations/Nothing_here.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("No significant correlations regarding microbes with at least 0.01% abund in 2 samples", fill=TRUE)
sink() # restore STR OUTPUT to R console
close(con)
suppressWarnings(rm(con))


######################### FOCUS ON SESSILE CILIATES ########################

data.genus.temp<-data.genus.prop
data.genus.temp <- subset_samples(data.genus.temp, Sample_name != "other_reactor") 
data.genus.temp<-subset_taxa(data.genus.temp, Order %in% "Sessilida")
# extracting the abundances
{ abundances<-otu_table(data.genus.temp)
  row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])
  # ordering the sample names according to the time flow 
  abundances<-abundances[ , c("1","57","115","164")]
  abundances<-rbind.data.frame(abundances, Total=apply(abundances,MARGIN=2,FUN=sum))  # NB: with margin 1 it would have summed the colums!
  abundances$Genus<-row.names(abundances)
}

table_to_plot<-melt(abundances, id.vars = "Genus")
table_to_plot$Genus<-factor( table_to_plot$Genus, levels = c(abundances$Genus[abundances$Genus!="Total"] , "Total") )

ggplot(data=table_to_plot, aes(y=value, x=variable, color=Genus)) +
  theme_classic(base_size =12) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), linewidth= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_5 ) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=25, vjust=1, hjust=1, size=8),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom",
        legend.margin = margin(1,42,0,0)
  ) +
  guides(color=guide_legend(nrow=1)) + 
  #scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  # scale_y_continuous(breaks = seq( 0, max(table_to_plot$value)+0.01, 0.01) ) +
  labs(x="Experiment day", y="ABS counts %",
       color="",
       title = "Sessile ciliates (Sessilida order) in pilot scale reactor")
ggsave(file="Results/Abundances/Sessile_ciliates_over_time.png",width=6.2,height=5, dpi=300) 
dev.off()




########################## FOCUS ON PAO AND GAO ###########################

### common PAO GAO list (PMID:38524765 ,  PMID: 22827168,  PMID: 15774634,  PMID: 33404187, MIDAS website)
PAO_list<-c("Candidatus Phosphoribacter","Ca Phosphoribacter","Phosphoribacter",
            "Tetrasphaera",
            "Candidatus Accumulimonas","Accumulimonas","Ca Accumulimonas",
            "Candidatus Accumulibacter","Accumulibacter","Ca Accumulibacter",
            "Candidatus Microthrix","Microthrix","Ca Microthrix",
            "Candidatus Dechloromonas","Dechloromonas",
            "Malikia",
            "Quatrionicoccus",
            "Beggiatoa",
            "Gemmatimonas",
            "Friedmaniella",
            "Tessaracoccus",
            "Azonexus",
            # "Pseudomonas",   # too aspecific!
            "Candidatus Accumulimonas","Accumulimonas","Ca Accumulimonas",
            "Microlunatis","Microlunatus",
            "Candidatus Lutibacillus","Lutibacillus","Ca Lutibacillus")
GAO_list<-c("Candidatus Competibacter","Candidatus Contendobacter","Candidatus Proximibacter",
            "Ca Competibacter","Ca Contendobacter","Ca Proximibacter",
            "Competibacter","Contendobacter",
            # "Proximibacter", # may also be a PAO!
            "Defluviicoccus","Micropruina","Propionivibrio")

data.genus.temp<-data.genus.prop
PAO_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% PAO_list )
GAO_identified_phylo <- subset_taxa(data.genus.temp, Genus %in% GAO_list)
tabella_PAO<-psmelt(PAO_identified_phylo)
tabella_GAO<-psmelt(GAO_identified_phylo)
tabella_PAO$GroupPG <- "PAO"
tabella_GAO$GroupPG <- "GAO"
tabella<- rbind.data.frame(tabella_PAO, tabella_GAO)
tabella$GroupPG<-factor(tabella$GroupPG, levels=c("PAO","GAO"))

tabella$GroupPG2<- "temp" # extra line to allow a quick reset of this part (the factors hate the g_subs )
tabella$GroupPG2[ tabella$Genus %in% c("Candidatus Accumulibacter","Candidatus Competibacter","Accumulibacter","Competibacter") ] <- tabella$Genus [tabella$Genus %in% c("Candidatus Accumulibacter","Candidatus Competibacter","Accumulibacter","Competibacter") ]
tabella$GroupPG2[ ! tabella$Genus %in% c("Candidatus Accumulibacter") &  tabella$Genus %in% PAO_list ] <- "Other PAOs"
tabella$GroupPG2[ ! tabella$Genus %in% c("Candidatus Competibacter") &  tabella$Genus %in% GAO_list ] <- "Other GAOs"
tabella$GroupPG2<- gsub("Candidatus ","Ca. ", tabella$GroupPG2)
tabella$GroupPG2<- factor( tabella$GroupPG2 , levels= c("Ca. Accumulibacter", "Other PAOs", "Ca. Competibacter", "Other GAOs") )
ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=GroupPG2)) +
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
        legend.margin =  margin(-10,2,0,0)
  ) +
  scale_y_continuous(lim=c(0,100), 
                     # breaks=seq(0,50,5)
  ) +
  labs(x="", y="Percent abundance in the sample", 
       fill="",
       title = paste("PAO / GAO abundances in the samples"))
ggsave(file="Results/Abundances/PAO_GAO_comparison.png",width=5.6,height=4.2,dpi=300)
dev.off()  


# exporting the abundances
PAO_identified_phylo<- phyloseq(otu_table(PAO_identified_phylo),tax_table(PAO_identified_phylo)) # removing the two trees
GAO_identified_phylo<- phyloseq(otu_table(GAO_identified_phylo),tax_table(GAO_identified_phylo))
tot <- merge_phyloseq(PAO_identified_phylo, GAO_identified_phylo)
# means of every PAO and GAO
tot2<- cbind.data.frame("Averages"=as.numeric(apply(otu_table(tot),1,mean)),
                        "Inoculum"=as.numeric(otu_table( prune_samples(sample_names(tot)%in%c("1"),tot) )),
                       "Averages_mature_pilot"=as.numeric(apply(otu_table( prune_samples(sample_names(tot)%in%c("57","115","164"),tot) ),1,mean)),
                       "Other_AGS"=as.numeric(otu_table( prune_samples(sample_names(tot)%in%c("other_reactor"),tot) )),
                       "genus"= as.data.frame(tax_table(tot))[["Genus"]])
tot2$Type <- ifelse( tot2$genus %in% PAO_list, "PAO", "GAO")
write.xlsx(file = "Results/Abundances/PAO_GAO_names_and_averages.xlsx", row.names = F,
           tot2
)


suppressWarnings(rm(tot, tabella, tabella_GAO, tabella_PAO) )




####################### FOCUS ON AOB, NOB and N reducer ########################


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
              "Ca Competibacter","Candidatus Competibacter", "Competibacter",
              # "Nitrosomonas",   # PMID: 11948173 and PMID: 17298375 ... it looks like it can also reduce NO2- ... but including it in this list would cause problems in the plot with AOB
              # "Nitrosospira",   # PMID: 17298375 , same as above
              # "Nitrospira", https://doi.org/10.1016/j.tim.2016.05.004 ... this is a NOB, yet same as above
              "Ca Accumulibacter","Candidatus Accumulibacter","Accumulibacter",
              "Haliangium",
              "Ca Promineofilum","Candidatus Promineofilum","Promineofilum",
              "Rhodoplanes", "Bradyrhizobium", "Iamia", "Thiothrix", "Sulfuritalea",
              "Zoogloea", "Thermogutta", "Diaphorobacter" , "Simplicispira",
              "Ca Phosphitivorax", "Candidatus Phosphitivorax", "Phosphitivorax",
              "Steroidobacter" ,"Sterolibacterium", "Azospira", "Acidovorax",
              "Pseudomonas",
              "Uruburuella", "Anaerolinea", "Corynebacterium",
              "Luteimonas", "Hyphomonas", "Halioglobus", "Azoarcus",
              "Ca Solibacter", "Candidatus Solibacter", "Solibacter",
              "Ca Brocadia", "Candidatus Brocadia", "Brocadia",
              "Skermania" , "Pannonibacter", "Sulfurimonas" , "Desulfitibacter" ,
              "Thiopseudomonas" , "Azonexus" ,
              "Ca Proximibacter", "Candidatus Proximibacter", "Proximibacter",
              #the following denitrifiers are from  PMID: 24559743
              "Enterobacter", "Micrococcus", "Spirillus","Proteus", "Aerobacter","Flavobacterium",
              # the following are from the FISH hand book  (Vedi anche tabella 3.1 del relativo capitolo)
              "Curvibacter", "Paracoccus", "Rhizobium", "Delftia", "Simplicispira","Dechloromonas", "Halomonas", "Thermomonas", "Aminomonas", "Methylophilaceae",
              "Azovibrio","Ralstonia", "Shinella","Ochrobactrum","Hyphomicrobium","Blastobacter", "Pseudoxanthomonas", "Stenotrophomonas", "Castellaniella",
              "Herbaspirillum","Citrobacter","Methylophilaceae","Alicycliphilus","Ottowia","Diaphorobacter", "Rhodobacter", "Microvirgula","Aquaspirillum", "Vogesella"
)

data.genus.temp<-data.genus.prop
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

ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=GroupPG)) + 
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
  # scale_y_continuous(lim=c(0,40), breaks=seq(0,50,5)) +
  labs(x="", y="Percent abundance in the sample", 
       fill="",
       title = paste("AOB / NOB / denitrifiers abundances in the samples"))
ggsave(file="Results/Abundances/AOB_NOB_denitrif_comparison.png",width=5.6,height=4.2,dpi=300)
dev.off()


# exporting the abundances
AOB_identified_phylo<- phyloseq(otu_table(AOB_identified_phylo),tax_table(AOB_identified_phylo)) # removing the two trees
NOB_identified_phylo<- phyloseq(otu_table(NOB_identified_phylo),tax_table(NOB_identified_phylo))
Nred_identified_phylo<- phyloseq(otu_table(Nred_identified_phylo),tax_table(Nred_identified_phylo))
tot <- merge_phyloseq(AOB_identified_phylo, NOB_identified_phylo, Nred_identified_phylo)
# means of every AOB, NOB and denitrif
tot<- cbind.data.frame("Averages"=as.numeric(apply(otu_table(tot),1,mean)),
                       # otu_table(tot), 
                       "genus"= as.data.frame(tax_table(tot))[["Genus"]]
)
tot$Type <- tot$genus
tot$Type[tot$Type %in% AOB_list] <-"AOB"
tot$Type[tot$Type %in% NOB_list] <-"NOB"
tot$Type[tot$Type %in% NO2_r_list] <-"N_reducer"
write.xlsx(file = "Results/Abundances/AOB_NOB_names_and_averages.xlsx", row.names = F,
           tot 
)


suppressWarnings(rm(tot, tabella, tabella_NOB, tabella_AOB) )




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

