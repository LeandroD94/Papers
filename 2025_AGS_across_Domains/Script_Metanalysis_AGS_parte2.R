# 26/07/2025

### Script part 2: main analysis



############################# IMPORTING THE OBJ ##################################

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
  library("reshape")
}
options(scipen = 100)


load("Data_prepared_for_analysis.RData")  # see script part 1


### these reactors are not SBR AGS ...
if ( "PRJNA1082061" %in% sample_data(data)$Project ){
  data_complete<-data
  data.genus_complete<- data.genus
}
data<-subset_samples(data, ! Project %in% c("PRJNA1082061","PRJNA606844","PRJNA1015110"))
data.genus<-subset_samples(data.genus, ! Project %in% c("PRJNA1082061","PRJNA606844","PRJNA1015110"))


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
levels(table$Project_grid)<-gsub("_","\n",levels(table$Project_grid))
table$Domain <- factor(table$Domain, levels = c("Bacteria","Archaea","Eukaryota","Viruses","unclassified"))
ggplot(data=table, aes(x=Sample, y=Abundance, fill=Domain)) +
  facet_grid(cols= vars(Project_grid),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("unclassified"="lightgray",
                             "Bacteria"="deepskyblue",
                             "Archaea"="blue4",
                             #"Eukaryota_Fungi"="gold3",
                             #"Eukaryota_Metazoa"="gold",
                             #"Eukaryota_Protists"="yellow",
                             "Eukaryota"="yellow",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust=1, size= 6),
        strip.text = element_text(size=6),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom",
        legend.margin = margin(-5,20,4,0),
        plot.margin = margin(2,1,2,1),
        panel.spacing.x = unit(1.5,"pt")
  ) +
  scale_x_discrete(expand=c(-0.1,1))+
  scale_y_continuous(expand=c(0,1))+
  #scale_y_sqrt(breaks=c(0,5,10,25,50,100)) + # can't be applied due to the values below 1
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Percentual abundance of clades",
       fill="",
       title = "Domains and unclassified proportions")
ggsave(file="Data_check/Domains_and_unclassified_according_to_Kaiju.png",width=6.8,height=4.2, dpi=450)
dev.off()




###################### DOMAINS PROPORTIONS #########################

data_temp<-data
tax_temp<-as.data.frame(tax_table(data_temp))
# Despite it is a "Sar" and "Alveolata", Symbodinium is a Dinoflagellate with photosynthesis, hence an algae!
tax_temp$Euk_clade[tax_temp$Genus %in% "Symbiodinium"] <- "Algae"
tax_temp$Domain[tax_temp$Domain %in% "Eukaryota"] <- tax_temp$Euk_clade[tax_temp$Domain %in% "Eukaryota"]
tax_table(data_temp)<-as.matrix(tax_temp)
data_temp<-tax_glom(data_temp, taxrank = "Domain", NArm = F) # without glomming before, some error may be seen in ggplot2 (maybe bugs due to too many observations?)

#data_temp<-subset_taxa(data_temp, Domain!="Unclassified" )
# !!! NB: due to a strange glitch of gglopt2, the "Bacteria" color in the stacked bar plot became invisible in few sample if the Domain column is modified in a factor with the Bacteria as first level (as below)
# therefore, ggplot2 will be tricked: the "Unclassified" will still be present but their abundance will be very minimal, in order to put them as "first place" (impossible to see) and the bacteria as second level (which actually seems the first one)
otu_table(data_temp)[ tax_table(data_temp)[, "Domain"]=="unclassified" ] <- 0.0000000001

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
  facet_grid(~ Reactor_scale, scales = "free", space = "free") +
  # facet_nested(~ Reactor_scale + Project_number, scales = "free", space = "free") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=c("Bacteria"="deepskyblue",
                             "Archaea"="blue3",
                             "Amoeba"="grey75",
                             "Fungi"="gold3",
                             "Algae"="green3",
                             "Nematoda"="brown4",
                             "Platyhelminthes"="coral3",
                             "Rotifera"="red3",
                             "Protozoa"="yellow",
                             "Protist"="grey30",
                             "Viruses"="violet")
  ) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size= 5.9),
        axis.title.y = element_text(vjust=0.85, size=10.5, hjust=0.65),
        strip.text = element_text(size=11.5),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.45, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom",
        legend.key.spacing.x = unit(0.3,"cm"),
        legend.margin = margin(-1,15,0,0),
        axis.title.x = element_blank() ,
        plot.margin = margin(-5,5,-10,1)
        # plot.title= element_text( size= 8), 
  ) +
  scale_x_discrete(c(-0.05,1)) +
  scale_y_break(c(3.5, 97.5), 
                scales=0.035,
                space = 0.3) +
  scale_y_continuous(
    #breaks = c(0,10,25,50,seq(50,1200,50)) ,
    #breaks = c(0,0.1,0.2,0.3,0.4,0.5, 1, 1.5, 97.5, 100),
    breaks = c(0,0.2,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 97.5, 100),
    limits = c(0,100.1)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(y="Percentual abundance of clades",
       fill="",
       caption= "NB: the Y axis is scaled to improve the readibility of lower values"
  )
ggsave(file="Data_check/Domains_according_to_Kaiju.png",width=6.6,height=4.5, dpi=300)
dev.off()

# means
data_temp<-subset_taxa(data_temp, ! Domain %in% "unclassified" )
to_save<- cbind.data.frame("Overall_average"=as.numeric(apply(otu_table(data_temp),1,mean)), 
                           "Averages_Pilot_scales"=as.numeric(apply (otu_table(subset_samples(data_temp, Reactor_scale=="Pilot")),1,mean)),
                           "Averages_Full_scales"=as.numeric(apply (otu_table(subset_samples(data_temp, Reactor_scale=="Full")),1,mean)),
                           "Averages_Lab_scales"=as.numeric(apply (otu_table(subset_samples(data_temp, Reactor_scale=="Lab")),1,mean)),
                           "Domain_Kingdom"= as.data.frame(tax_table(data_temp))[["Domain"]]
)
write.csv(to_save, file = "Data_check/Domains_averages_according_to_Kaiju.csv", row.names = F, quote=F)





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
{ #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  # write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Relative_abundances/counts_species_prop.csv",quote=F)
  #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
}




###################### ABUNDANCES BAR PLOT ##########################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.genus.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top <- names(sort(taxa_sums(data_temp), decreasing=TRUE))[1:24]
  prune.dat_top <- prune_taxa(top,data_temp)
  others<-taxa_names(data_temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the species
  tabella<-aggregate(.~Genus+Euk_clade+Reactor_scale+Influent+Project_number+Sample_name, tabella[ , c("Reactor_scale","Project_number","Euk_clade","Influent","Sample_name","Abundance","Genus")], FUN=sum)
  tabella$Genus<-gsub ("Candidatus_","Ca.", tabella$Genus)
  tabella$Genus<-gsub ("bacterium_","", tabella$Genus)
  tabella$Genus<-paste0(tabella$Euk_clade,": ",tabella$Genus)
  tabella$Genus<-gsub("No_euk","Prokaryota", tabella$Genus)
  tabella$Genus<-gsub(".*Others","Others", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
fill_color_25_BIS <- fill_color_25
fill_color_25_BIS[fill_color_25=="yellow3"]<-"darkcyan"
to_plot <- ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  # scale_fill_manual(values=fill_color_30_BIS) +
  scale_fill_manual(values=fill_color_25_BIS) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=6.7),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.55, "cm"),
  legend.text = element_text ( size = 8.5 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-15,15,0,0),
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  panel.spacing.x = unit(1.45,"pt"),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(2,2,2,2)
  ) +
  guides(fill=guide_legend(nrow=9)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
# ID
to_plot +
  facet_nested(~ Reactor_scale + Project_number, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(size=8.5))
ggsave(file="Results/Abundances/TOP_microbial_genera.png",width=6.58,height=4.5, dpi=300)
dev.off()
# # Influent
# to_plot + facet_nested(~ Influent + Reactor_scale, scales = "free_x", space = "free_x") +
#   theme(strip.text = element_text(size=9.2))
# ggsave(file="Results/Abundances/TOP_microbial_Genus_Influent.png",width=6.6,height=4.5, dpi=300)
# dev.off()

# means of TOP genera
to_save<- cbind.data.frame("Overall_average"=round( as.numeric(apply(otu_table(prune.dat_top),1,mean)) ,3 ),
                           # "Overall_median"=round( as.numeric(apply(otu_table(prune.dat_top),1,median)) ,3 ), 
                           "Averages_Full scales"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)), 3),
                           "Averages Pilot scales"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)), 3),
                           "Averages_Lab scales"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)), 3),
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]],
                           "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]]
)
to_save$Genus<- gsub("_"," ", to_save$Genus)
to_save$Genus<- gsub("bacterium ","", to_save$Genus)
to_save$Main_role<-""
to_save$Main_role[to_save$Genus%in% PAO_list] <- "PAO"
to_save$Main_role[to_save$Genus%in% GAO_list] <- "GAO"
to_save$Main_role[to_save$Genus%in% AOB_list] <- "AOO"
to_save$Main_role[to_save$Genus%in% NOB_list] <- "NOO"
# to_save$Main_role[to_save$Genus%in% table_PHA[,1]] <- "PHA producer"
to_save <- to_save[order(to_save$Overall_average, decreasing = T), ]
                    
write.csv(to_save, file = "Results/Abundances/TOP_genera_Average_abundances_no_unclassified.csv", row.names = F, quote=F)

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



### TOP Species
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top_data <- subset_taxa(data_temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:24]
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
  tabella<-aggregate(.~Species+Reactor_scale+Influent+Project_number+Sample_name, tabella[ , c("Reactor_scale","Project_number","Influent","Sample_name","Abundance","Species")], FUN=sum)
  tabella$Species<-gsub ("Candidatus ","Ca.", tabella$Species)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
to_plot <- ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Species)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6,
  ),
  axis.title =element_text(size=10),
  # axis.title.x = element_text(vjust=2.5),
  # strip.text = element_text(size=6.7),
  # legend.spacing.x  = unit(2, "cm"),
  legend.key.height = unit(0.15, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 8.6 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-13,25,0,0),
  panel.spacing.x = unit(1.45,"pt"),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(1,2,2,2)
  ) +
  guides(fill=guide_legend(nrow=9)) +
  labs(x="", y="Percentual abundance of clades", 
       fill="")
# ID
to_plot +
  facet_nested(~ Reactor_scale + Project_number, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(size=8.5))
ggsave(file="Results/Abundances/TOP_microbial_Species.png",width=6.6,height=4.5, dpi=300)
dev.off()
# # Influent
# to_plot + facet_nested(~ Influent + Reactor_scale, scales = "free_x", space = "free_x") +
#   theme(strip.text = element_text(size=9.25))
# ggsave(file="Results/Abundances/TOP_microbial_Species_Influent.png",width=6.6,height=4.5, dpi=300)
# dev.off()

# means of TOP Species
to_save<- cbind.data.frame("Overall_average"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), 
                           "Averages_Pilot_scales"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)),
                           "Averages_Full_scales"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)),
                           "Averages_Lab_scales"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)),
                           "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]],
                           "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]]
)
write.csv(to_save, file = "Results/Abundances/TOP_Species_Average_abundances_no_unclassified.csv", row.names = F, quote=F)


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




##### AGAIN, BUT FOCUSING ON FULL SCALES ONLY (SPECIES) ...

suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.prop, Genus!="unclassified")
data_temp<- subset_samples(data_temp, Reactor_scale=="Full")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top_data <- subset_taxa(data_temp)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:24]
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
  tabella<-aggregate(.~Species+Reactor_scale+Influent+Project_number+Sample_name, tabella[ , c("Reactor_scale","Project_number","Influent","Sample_name","Abundance","Species")], FUN=sum)
  tabella$Species<-gsub ("Candidatus ","Ca. ", tabella$Species)  # to save horizontal space in the plot!
  tabella$Species<-gsub ("bacterium ","", tabella$Species)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Species<-factor(tabella$Species, levels = c(unique(tabella$Species)[! unique(tabella$Species) %in% "Others"],"Others"))
}
to_plot <- ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Species)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  #scale_fill_manual(values=fill_color_30_BIS) +
  scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 6.5
  ),
  axis.title =element_text(size=10),
  # axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=11),
  # legend.spacing.x  = unit(3, "cm"),
  legend.key.height = unit(0.15, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 8.6 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-12,25,0,0),
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  panel.spacing.x = unit(1.45,"pt"),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(2,2,2,2)) +
  guides(fill=guide_legend(nrow=9)) +
  labs(x="", y="Percentual abundance of clades", 
       fill="")
# ID
to_plot + 
  facet_nested(~ Reactor_scale + Project_number, scales = "free_x", space = "free_x")
ggsave(file="Results/Abundances/TOP_microbial_Species_ONLY_FULL_SCALES_ID.png",width=6,height=4.5, dpi=300)
dev.off()

# means
to_save<- cbind.data.frame("Overall_average"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), 
                           #"Averages_Pilot_scales"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)),
                           "Averages_Full_scales"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)),
                           #"Averages_Lab_scales"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)),
                           "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]],
                           "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]]
)
write.csv(to_save, file = "Results/Abundances/TOP_Species_Average_abundances_no_unclassified_FULLscales.csv", row.names = F, quote=F)




########################## FOCUS ON AGS KEY MEMBERS ###########################

### PREPARING THE MAIN OBJ ...
data.genus.temp<-subset_taxa( data.genus.prop, Genus!="unclassified")   # NB: removing unclassified!
data.genus.temp<- subset_taxa(data.genus.temp, !Domain %in% c("Viruses"))  # otherwise viruses as "Flavobacterium phage" will be included!
tax_table(data.genus.temp)[ ,"Genus"] <- gsub("Candidatus_","Candidatus ", tax_table( data.genus.temp )[ ,"Genus"])
data.genus.temp<-subset_taxa( data.genus.temp, Genus %in% c(NOB_list,AOB_list,PAO_list,GAO_list, table_PHA[,1]) )   # NB: removing unclassified!

# Solving eventual synonyms (see LPSN) regarding PHA accumulators...
{
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Fuscovulum","Rhodobacter", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Phaeovulum","Rhodobacter", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Paucibacter","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Mitsuaria","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Kinneretia","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Roseateles","Pelomonas", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("^Competibacter","Candidatus Competibacter", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("^Accumulibacter","Candidatus Accumulibacter", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Paramicrobacterium","Microbacterium", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Aureobacterium","Microbacterium", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Prescottella","Rhodococcus", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Ensifer Casida","Sinorhizobium", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("Calidifontimicrobium","Azohydromonas", tax_table(data.genus.temp)[,"Genus"])
}

data.genus.temp<-tax_glom(data.genus.temp, taxrank = "Genus", NArm = F)  # re-glomming

# Labelling the functional groups
colnames(tax_table(data.genus.temp))[colnames(tax_table(data.genus.temp))=="Euk_clade"] <- "Group_Function"   # overwriting this column (empty at genus level)
easy_subset <- as.character(tax_table(data.genus.temp)[ ,"Genus"])
tax_table(data.genus.temp)[  easy_subset %in% AOB_list ,"Group_Function"] <- "AOO"
tax_table(data.genus.temp)[  easy_subset %in% NOB_list ,"Group_Function"] <- "NOO"
tax_table(data.genus.temp)[  easy_subset %in% PAO_list ,"Group_Function"] <- "PAO"
tax_table(data.genus.temp)[  easy_subset %in% GAO_list ,"Group_Function"] <- "GAO"
tax_table(data.genus.temp)[  easy_subset %in% table_PHA[,1] ,"Group_Function"] <- "Others PHA"
data.genus.temp<- prune_taxa(taxa_sums(data.genus.temp)>0, data.genus.temp)

# Obj to plot
tabella <- psmelt(data.genus.temp)
tabella$GroupFunction<-factor(tabella$Group_Function, levels=c("PAO","GAO","AOO","NOO","Others PHA"))
# extra focus on Accumuli and Competi for the plots ...
tabella$GroupFunction2<- as.character(tabella$GroupFunction) # extra line to allow a quick repeating of this section if needed...
tabella$GroupFunction2[ tabella$Genus %in% c("Candidatus Accumulibacter","Candidatus Competibacter","Accumulibacter","Competibacter") ] <- tabella$Genus [tabella$Genus %in% c("Candidatus Accumulibacter","Candidatus Competibacter","Accumulibacter","Competibacter") ]
tabella$GroupFunction2[ ! tabella$Genus %in% c("Candidatus Accumulibacter") &  tabella$Genus %in% PAO_list ] <- "Other PAOs"
tabella$GroupFunction2[ ! tabella$Genus %in% c("Candidatus Competibacter") &  tabella$Genus %in% GAO_list ] <- "Other GAOs"
tabella$GroupFunction2<- gsub("Candidatus ","Ca. ", tabella$GroupFunction2)
tabella$GroupFunction2<- factor( tabella$GroupFunction2 , levels= c("Ca. Accumulibacter", "Other PAOs", "Ca. Competibacter", "Other GAOs","AOO","NOO","Others PHA") )

# PLOTTING ...
ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=GroupFunction2)) +
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  facet_nested(~ Reactor_scale + Project_number , scales = "free_x", space = "free_x") +  
  scale_fill_manual(values=c("Ca. Accumulibacter"="coral", "Ca. Competibacter"="chartreuse2",
                             "Other PAOs"="red2","Other GAOs"="chartreuse4",
                             "AOO"="deepskyblue","NOO"="blue3",
                             "Others PHA"="darkred"
  )
  ) +
  theme(axis.text.x=element_text(angle=43, hjust=1,vjust=1, size=6),
        axis.title.y = element_text(size=10.5), 
        axis.text.y = element_text(size=8), 
        strip.text.x = element_text(size=11.2), 
        legend.key.height = unit(0.17, "cm"),
        legend.key.width = unit(0.55, "cm"),
        legend.text = element_text ( size = 12.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,2,2,2),
        panel.spacing.x = unit(1.45,"pt"),
        legend.margin =  margin(-10,14,0,0)
  ) +
  scale_y_continuous(lim=c(0,35),
                     breaks=seq(0,100,5),
                     expand = c(0,0)
  ) +
  scale_x_discrete(expand=c(-0.1,1)) +
  labs(x="", y="Percentual abundance of groups", 
       fill="")
ggsave(file="Results/Abundances/Functional_groups_averages.png",width=6.5,height=4.5,dpi=300)
dev.off()  


# means of every PAO and GAO
tot2<- cbind.data.frame("Overall_Averages"=as.numeric(apply(otu_table(data.genus.temp),1,mean)),
                        "Averages_FullScale"=as.numeric(apply (otu_table(subset_samples(data.genus.temp, Reactor_scale=="Full")),1,mean)) ,
                        "Averages_PilotScale"=as.numeric(apply (otu_table(subset_samples(data.genus.temp, Reactor_scale=="Pilot")),1,mean)),
                        "Averages_LabScale"=as.numeric(apply (otu_table(subset_samples(data.genus.temp, Reactor_scale=="Lab")),1,mean)),
                        "Genus"= as.data.frame(tax_table(data.genus.temp))[["Genus"]],
                        "Functional_Group"= as.data.frame(tax_table(data.genus.temp))[["Group_Function"]]
)
tot2 <- tot2[order(tot2$Overall_Averages, decreasing = T), ]

write.csv(file = "Results/Abundances/Functional_groups_averages_WITHOUT_UNCLASSIFIED.csv", row.names = F, quote = F, 
          tot2
)


suppressWarnings(rm(tot, tot2, tabella, data.genus.temp) )




################# FOCUS ON  ON SESSILE CILIATES ###########################

data.genus.temp<-subset_taxa( data.genus.prop, Domain=="Eukaryota")   # NB: also removing unclassified!
data.genus.temp <- transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100 )
data.genus.temp<- subset_taxa(data.genus.temp, Genus %in%  sessilida_list )
tabella<-psmelt(data.genus.temp)
ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  #facet_nested(~ Influent + Reactor_scale, scales = "free_x", space = "free_x") +  
  facet_nested(~  Reactor_scale + Project_number, scales = "free_x", space = "free_x") +  
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1, size=6),
        axis.title.y = element_text(size=10.5), 
        axis.text.y = element_text(size=8), 
        strip.text.x = element_text(size=11.5), 
        legend.key.height = unit(0.18, "cm"),
        legend.key.width = unit(0.58, "cm"),
        legend.text = element_text ( size = 12.5 , face="italic"),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,2,2,2),
        panel.spacing.x = unit(1.5,"pt"),
        legend.margin =  margin(-10,8,0,0)
  ) +
  
  scale_y_continuous(lim=c(0,14.3),
                     breaks=seq(0,14,2),
                     expand = c(0,0)
  ) +
  scale_x_discrete(expand=c(-0.025,1)) +
  labs(x="", y="Percentual abundance of clades", 
       fill="")
ggsave(file="Results/Abundances/FOCUS_ON_SESSILIDA.png",width=6.5,height=4.5,dpi=300)
dev.off()  




############ EXTRA CHECK: ABUNDANCES OF EUKARYOTA (TOTAL ON EUKARYOTA) ####################

# TOP eukar genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ 
  data.subset <- subset_taxa(data.genus, Domain=="Eukaryota")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:24]
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
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_nested(~ Reactor_scale + Project_number, scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  #scale_fill_manual(values=fill_color_30_BIS) +
  scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=9),
  legend.key.height = unit(0.18, "cm"),
  legend.key.width = unit(0.42, "cm"),
  legend.text = element_text ( size = 8.7 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-14,30,0,0),
  panel.spacing.x = unit(1.45,"pt"),
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(1,1,2,1)) +
  guides(fill=guide_legend(nrow=9)
  ) +
  labs(x="", y="Percentual abundance of clades", fill="")
ggsave(file="Results/Abundances/TOP_Eukaryota_ONLY_EUKAR.png",width=6.6,height=4.5, dpi=300)
dev.off()

to_save<- cbind.data.frame("Overall_average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2 ),
                           "Overall_average_on_full_dataset" =  round( as.numeric(apply ( otu_table( prune_taxa(top,data.genus.prop) ) , 1,mean) ), 4) ,
                           "Averages_Pilot_scales"= round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)), 2),
                           "Averages_Full_scales"= round(as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)), 2),
                           "Averages_Lab_scales"= round(as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)), 2),
                           "Genera"= as.data.frame(tax_table(prune.dat_top))[["Genus"]],
                           "Euk_clade"= as.data.frame(tax_table(prune.dat_top))[["Euk_clade"]]
)
to_save <- to_save[order(to_save$Overall_average , decreasing = T) , ]

round( as.numeric(apply ( otu_table( prune_taxa(top,data.genus.prop) ) , 1,mean) ), 4)

write.csv(to_save, file = "Results/Abundances/TOP_Eukaryota_ONLY_EUKAR_averages.csv", row.names = F, quote=F)




###################### EXTRA CHECK: ABUNDANCES OF ARCHAEA ##########################

### TOP archaea
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data.genus, Domain=="Archaea")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified") 
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:24]
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
  tabella$Genus<-gsub ("Candidatus ","Ca.", tabella$Genus)  # to save horizontal space in the plot!
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)  # to save horizontal space in the plot!
  tabella$Genus<-gsub ("uncultured","uncult.", tabella$Genus)
  tabella$Genus<-gsub ("^archaeon$","unknown Archaea", tabella$Genus)
  tabella$Genus<-gsub ("Methanobrevibacter sp. TLL-48-HuF1","Methanobrevibacter TLL-48-HuF1", tabella$Genus)
  tabella<-aggregate(.~Genus+Reactor_scale+Influent+Project_number+Sample_name, tabella[ , c("Reactor_scale","Project_number","Influent","Sample_name","Abundance","Genus")], FUN=sum)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_nested(~ Reactor_scale + Project_number, scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  #scale_fill_manual(values=fill_color_30_BIS) +
  scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=9.2),
  legend.key.height = unit(0.18, "cm"),
  legend.key.width = unit(0.34, "cm"),
  legend.text = element_text ( size = 8.7, face="italic" ),
  legend.position="bottom",
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  panel.spacing.x = unit(1.45,"pt"),
  legend.margin = margin(-14,30,0,0),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(2,2,2,2)) +
  guides(fill=guide_legend(nrow=9)) +
  labs(x="", y="Percentual abundance of clades", 
       fill="")
ggsave(file="Results/Abundances/TOP_Archaea_species_ONLY_ARCHAEA.png",width=6.6,height=4.5, dpi=300)
dev.off()

to_save<- cbind.data.frame("Overall_average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2),
                           "Overall_average_on_full_dataset" =  round( as.numeric(apply ( otu_table( prune_taxa(top,data.genus.prop) ) , 1,mean) ), 4) ,
                           "Averages_Pilot_scales"= round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)), 2) ,
                           "Averages_Full_scales"= round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)), 2),
                           "Averages_Lab_scales"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)), 2),
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]],
                           "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]]
)
to_save<-to_save[order(to_save$Overall_average, decreasing = T), ]
to_save$Genus<-gsub ("^archaeon$","unknown Archaea", to_save$Genus)


write.csv(to_save, file = "Results/Abundances/TOP_Archaea_species_ONLY_ARCHAEA_averages.csv", row.names = F, quote=F)




###################### EXTRA CHECK: ABUNDANCES OF VIRUSES ##########################

### TOP viruses species
suppressWarnings(rm(top, others, tabella, unass_data))
{ data.subset <- subset_taxa(data, Domain=="Viruses")
  data.subset <- transform_sample_counts(data.subset, function(x) (x/sum(x))*100 )
  # two obs are omonymous! The tax ID will be add (seen in the raw Kaiju output)
  # View(tax_table(prune.dat_top))
  tax_table(data.subset)["sp3177","Domain_species"] <- "Viruses_Caudoviricetes 2788787"
  tax_table(data.subset)["sp9747","Domain_species"] <- "Viruses_Caudoviricetes 2731619"
  otu_table(data.subset)[is.nan(otu_table(data.subset))] <- 0
  top_data <- subset_taxa(data.subset, Genus!="unclassified")
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:24]
  prune.dat_top <- prune_taxa(top,data.subset)
  others<-taxa_names(data.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.subset)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Domain_species)){ # to mark the fungal taxa
  #   tabella_top$Domain_species[i]<-paste0(tabella_top$Domain_species[i],"_",tabella_top$Phylum[i])
  # }
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Domain_species<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Domain_species<-gsub ("Viruses_","", tabella$Domain_species)  # to save horizontal space in the plot!
  tabella$Domain_species<-gsub ("_ply2008005c","", tabella$Domain_species)  # to save horizontal space in the plot!
  tabella<-aggregate(.~Domain_species+Reactor_scale+Influent+Project_number+Sample_name, tabella[ , c("Reactor_scale","Project_number","Influent","Sample_name","Abundance","Domain_species")], FUN=sum)
  tabella$Domain_species<-factor(tabella$Domain_species, levels = c(unique(tabella$Domain_species)[! unique(tabella$Domain_species) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Domain_species)) +
  facet_nested(~ Reactor_scale + Project_number, scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  #scale_fill_manual(values=fill_color_30_BIS) +
  scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  strip.text = element_text(size=9.2),
  legend.key.height = unit(0.18, "cm"),
  legend.key.width = unit(0.34, "cm"),
  legend.text = element_text ( size = 8.5, face="italic" ),
  legend.position="bottom",
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  legend.margin = margin(-10,30,0,0),
  scale_x_discrete(expand=c(-0.1,1)),
  panel.spacing.x = unit(1.45,"pt"),
  plot.margin = margin(2,2,2,2)) +
  guides(fill=guide_legend(nrow=9)) +
  labs(x="", y="Percentual abundance of clades", 
       fill="")
ggsave(file="Results/Abundances/TOP_Viruses_species_ONLY_VIRUSES.png",width=6.6,height=4.5, dpi=300)
dev.off()

to_save<- cbind.data.frame("Overall_average"= round( as.numeric(apply(otu_table(prune.dat_top),1,mean)), 2),
                           # "Overall_average_on_full_dataset" =  round( as.numeric(apply ( otu_table( prune_taxa(top,data.prop) ) , 1,mean) ), 4) ,
                           "Averages_Pilot_scales"= round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)), 2),
                           "Averages_Full_scales"= round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)), 2) ,
                           "Averages_Lab_scales"= round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)), 2) ,
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]],
                           "Species"= as.data.frame(tax_table(prune.dat_top))[["Species"]],
                           "Domain_species"= as.data.frame(tax_table(prune.dat_top))[["Domain"]]
)
to_save <- to_save[ order(to_save$Overall_average, decreasing = T) , ]
write.csv(to_save, file = "Results/Abundances/TOP_Viruses_species_ONLY_VIRUSES_averages.csv", row.names = F, quote=F)





######################## ALPHA DIVERSITY ##############################


### PROKARYOTES (GENERA)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Bacteria", "Archaea") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_scale", 
                      color="Project_number",
                      #color="Reactor_scale",
                      #shape="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_PROKARYOTES.png", width = 5.2,height =4, dpi=300)

suppressWarnings(rm(mix_alpha) )



### PROKARYOTES (GENERA) WITH ABUND FILTER
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Bacteria", "Archaea") )
who_over <- tax_table(data.genus)[ rowMeans( otu_table(data.genus.prop) )>0.0001 , "Genus"]
data.genus.temp<-subset_taxa(data.genus.temp, Genus %in% who_over )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_scale", 
                      color="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
    # caption="filtered"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth =0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_PROKARYOTES_IF_FILTER_OVER_00001MEAN.png", width = 5.2,height =4, dpi=300)

suppressWarnings(rm(mix_alpha) )




### PROKARYOTES (SPECIES)
data.temp<-subset_taxa(data, Domain %in% c("Bacteria", "Archaea") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor_scale", 
                      color="Project_number",
                      #color="Reactor_scale",
                      #shape="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1),
        legend.margin = margin(-7,0,2,0)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_PROKARYOTES.png", width = 5,height =3.5, dpi=300)

suppressWarnings(rm(mix_alpha) )




### PROKARYOTES (Species) WITH ABUND FILTER
data.temp<-subset_taxa(data, Domain %in% c("Bacteria", "Archaea") )
who_over <- tax_table(data.prop)[ rowMeans( otu_table(data.prop) )>0.0001 , "Species"]
data.temp<-subset_taxa(data.temp, Species %in% who_over )

pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor_scale", 
                      color="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9),
        axis.text.y= element_text(size=7.8 )
        ) +
  theme(panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1),
        legend.margin = margin(-7,0,2,0)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_PROKARYOTES_IF_FILTER_OVER_00001MEAN.png", width = 5,height =3.5, dpi=300)




### EUKARYOTES
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Eukaryota") )
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_scale",
                      color="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) +
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1),
        legend.margin = margin(-7,0,2,0)
        )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_EUKARYOTA.png", width = 5,height =3.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### EUKARIOTA (WITH ABUND FILTER)
data.genus.temp<-subset_taxa(data.genus, Domain %in% c("Eukaryota") )
data.genus.temp.prop<- transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100 )
who_over <- tax_table(data.genus.temp.prop)[ rowMeans( otu_table(data.genus.temp.prop) )>0.0001 , "Genus"]
data.genus.temp<-subset_taxa(data.genus.temp, Genus %in% who_over )

pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_scale",
                      color="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) +
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1),
        legend.margin = margin(-7,0,2,0)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_GENUS_EUKARYOTA_IF_FILTER_OVER_00001MEAN.png", width = 5,height =3.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### VIRUSES
data.temp<-subset_taxa(data, Domain %in% c("Viruses") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor_scale",
                      color="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) +
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1),
        legend.margin = margin(-7,0,2,0)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_Viruses.png", width = 5,height =3.5, dpi=300)

suppressWarnings(rm(mix_alpha) )



### VIRUSES (WITH ABUND FILTER)
data.temp<-subset_taxa(data, Domain %in% c("Viruses") )
data.temp.prop <- transform_sample_counts(data.temp, function(x) (x/sum(x))*100 )
who_over <- tax_table(data.temp.prop)[ rowMeans( otu_table(data.temp.prop) )>0.0001 , "Species"]
data.temp<-subset_taxa(data.temp, Species %in% who_over )

pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor_scale",
                      color="Project_number"
)
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
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) +
  theme_classic() +
  scale_color_manual(values = fill_color_10 ) +
  labs( 
    y="Alpha Diversity Measure",
    x="",
    color="Dataset"
  ) +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9)) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1),
        legend.margin = margin(-7,0,2,0)
  )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_SPECIES_Viruses_IF_FILTER_OVER_00001MEAN.png", width = 5,height =3.5, dpi=300)

suppressWarnings(rm(mix_alpha) )




### Prokaryotes species (DEPTH)
data.temp<-subset_taxa(data, Domain %in% c("Bacteria", "Archaea") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor_scale", 
                      color="Seq_depth",
                      #color="Reactor_scale",
                      #shape="Project_number"
)
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
pAlpha$data$Seq_depth<-factor(pAlpha$data$Seq_depth, levels=c("<5","<10","10","15","20"))
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) + 
  theme_classic() +
  scale_color_manual(values = c("<5"="red", "<10"="coral",
                                "10"="gold2", "15"="darkgreen",
                                "20"="chartreuse")) +
  labs(title="Alpha diversity (prokaryotes species)", 
       y="Alpha Diversity Measure",
       x="Reactor scale") +
  guides(fill="none", shape="none") + 
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(linewidth=0.4), 
        panel.grid.minor.y = element_line(linewidth=0.25) , 
        axis.title.x = element_text(vjust = -1),
        legend.position = "bottom",
  )
ggsave(file="Data_check/Alpha_div_Prokaryotes_DepthSeq.png", width = 5.2, height =4, dpi=300)



### VIRUSES (DEPTH)
data.temp<-subset_taxa(data, Domain %in% c("Viruses") )
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Reactor_scale",
                      color="Seq_depth",
                      #color="Reactor_scale",
                      #shape="Project_number"
)
#pAlpha
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
pAlpha$data$Seq_depth<-factor(pAlpha$data$Seq_depth, levels=c("<5","<10","10","15","20"))
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_scale, y=value, color=NULL, shape=NULL), alpha=0) +
  theme_classic() +
  # geom_text(aes(label=Seq_depth), color="black" , size=2)+
  scale_color_manual(values = c("<5"="red", "<10"="coral",
                                "10"="gold2", "15"="darkgreen",
                                "20"="chartreuse")) +
  labs(title="Alpha diversity (virus species)",
       y="Alpha Diversity Measure",
       x="Reactor scale") +
  guides(fill="none", shape="none") +
  theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.25) ,
        axis.title.x = element_text(vjust = -1),
        legend.position = "bottom",
  )
ggsave(file="Data_check/Alpha_div_Viruses_DepthSeq.png", width = 5.2, height =4, dpi=300)




########################### PCoA AND BETA DIVERSITY #############################

data.prop_bacteria_Genus<-subset_taxa(data.genus.prop, Domain %in% c("Bacteria","Archea"))
data.prop_bacteria_Species<-subset_taxa(data.prop, Domain %in% c("Bacteria","Archea"))
data.prop_eukar_Genus<-subset_taxa(data.genus.prop, Domain=="Eukaryota")
data.prop_eukar_Species<-subset_taxa(data.prop, Domain=="Eukaryota")
data.prop_virus_Genus<-subset_taxa(data.genus, Domain=="Viruses")
data.prop_virus_Species<-subset_taxa(data.prop, Domain=="Viruses")


for( d in list(data.genus.prop, data.prop, 
               data.prop_bacteria_Genus, data.prop_bacteria_Species,
               data.prop_eukar_Genus, data.prop_eukar_Species,
               data.prop_virus_Genus, data.prop_virus_Species
) 
) {
  
  # setting the automatic labels and file names
  if(identical(d,data.genus.prop) | identical(d,data.prop_bacteria_Genus) ){  # identical(d,data.prop_eukar_Genus)){
    levels<-"genera"
  } else { levels <-"species"}
  if(identical(d,data.genus.prop) | identical(d,data.prop)){
    subset<-"PCoA_with_Every_microbes"
    file_name<-"microbial"
    plot_title<-"identified microbes"
  }
  if(identical(d,data.prop_bacteria_Genus) | identical(d,data.prop_bacteria_Species)){
    subset<-"PCoA_with_Prokaryotes_only"
    file_name<-"prokaryotic"
    plot_title<-"identified prokaryotes"
  }
  if(identical(d,data.prop_eukar_Genus) | identical(d,data.prop_eukar_Species)){
    subset<-"PCoA_with_Eukar_only"
    file_name<-"Eukar"
    plot_title<-"identified eukariotes"
  }
  if(identical(d,data.prop_virus_Genus) | identical(d,data.prop_virus_Species)){
    #subset<-"PCoA_with_virus_only"
    subset<-"PCoA_with_Viruses_only"
    file_name<-"Viruses"
    plot_title<-"identified viruses"
  }
  
  
  data.prop.labels<-d
  sample_data(data.prop.labels)$Influent_extra<- gsub("hypersaline,pharmaceutical","hypersaline\npharmaceutical",sample_data(data.prop.labels)$Influent_extra)
  {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
    DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
    ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
    eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
    eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
  }
  
  
  
  # Reactor number 
  plot_ordination(data.sqrt_prop, ordBC , color="Reactor_number", shape= "Reactor_scale") +
    # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
    geom_text(aes(label=sample_data(data.prop.labels)$Reactor_number), color="grey20", size=2.8, show.legend = FALSE) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    guides(color="none")+
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Reactor number",
         shape="Reactor scale",
         caption="The labels on the points are the reactor IDs",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_with_ID_Hellinger_on_",file_name,"_",levels,"_ReactorIDs.png"), width = 5.55, height = 4.2, dpi=300)
  
  
  # Project number 
  plot_ordination(data.sqrt_prop, ordBC , color="Project_number", shape="Reactor_scale") +
    scale_color_manual(values=fill_color_15) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
    geom_text(aes(label=sample_data(data.prop.labels)$Reactor_number), color="grey20", size=2.8, show.legend = FALSE) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    guides(color="none")+
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Reactor number",
         shape="Reactor scale",
         caption="The labels on the points are the reactor IDs",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_Hellinger_on_",file_name,"_",levels,"_projectIDs.png"), width = 5.55, height = 4.2, dpi=300)
  
  
  # with ellipses (Reactor scale)
  plot_ordination(data.sqrt_prop, ordBC , color="Reactor_scale",
                  shape= "Influent_type") +
    scale_color_manual(values=c("chartreuse","coral","deepskyblue")) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    stat_ellipse(linewidth=0.18, aes(color=Reactor_scale, shape=NULL)) + 
    #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
    #geom_text(aes(label=sample_data(data.prop.labels)$Reactor_number), color="grey20", size=2.8, show.legend = FALSE) +
    geom_text(aes(label=sample_data(data.prop.labels)$Project_number), color="grey20", size=2.8, show.legend = FALSE) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    # guides(color="none")+
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Reactor scale",
         shape="Influent type",
         caption="The labels on the points are the project numbers",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_ScaleEllipse.png"), width = 5.1, height = 3.8, dpi=300)
  
  
  # with ellipses (Influent)
  plot_ordination(data.sqrt_prop, ordBC , color="Influent", 
                  shape= "Reactor_scale") +
    # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    stat_ellipse(linewidth=0.18, aes(color=Influent)) + 
    #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
    geom_text(aes(label=sample_data(data.prop.labels)$Reactor_number), color="grey20", size=2.8, show.legend = FALSE) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    # guides(color="none")+
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Influent type",
         shape="Reactor scale",
         caption="The labels on the points are the reactor IDs",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_InfluentEllipse.png"), width = 5.5, height = 4, dpi=300)
  
  
  # with ellipses (Influent + details)
  plot_ordination(data.sqrt_prop, ordBC , color="Influent", 
                  shape= "Reactor_scale") +
    # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    stat_ellipse(linewidth=0.18, aes(color=Influent)) + 
    geom_text(aes(label=sample_data(data.prop.labels)$Influent_extra), 
              color="grey20", size=2, show.legend = FALSE,
              lineheight = 0.8) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    # guides(color="none")+
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Influent type",
         shape="Reactor scale",
         caption="The labels on the points are the reactor details",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_InfluentEllipse_v2.png"), width = 5.5, height = 4, dpi=300)
  
  
  # EXTRA CHECK: DEPTH
  plot_ordination(data.sqrt_prop, ordBC , color="Seq_depth", shape= "Reactor_scale") +
    scale_color_manual(values = c("<5"="red", "<10"="coral",
                                  "10"="gold2", "15"="darkgreen",
                                  "20"="chartreuse")) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
    geom_text(aes(label=sample_data(data.prop.labels)$Seq_depth), color="grey20", size=2.8, show.legend = FALSE) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Reactor number",
         shape="Reactor scale",
         caption="The labels on the points are sequencing depths (GB)",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_SeqDEPTH.png"), width = 5.5, height = 4, dpi=300)
  
  
  # EXTRA: Removing the reactor 15 in case of the virus profiling, it is too anomalous
  if( identical(d,data.prop_virus_Species) | identical(d, data.prop_virus_Genus)){
    file_name2<-"Viruses_NO_REACT15"
    d2<-subset_samples(data.prop.labels, Reactor_number!="15")
    data.sqrt_prop2<-transform_sample_counts(d2, sqrt) # square root of proportion
    DistBC2 = phyloseq::distance(data.sqrt_prop2, method = "euclidean")
    ordBC2 = ordinate(data.sqrt_prop2, method = "PCoA", distance = DistBC2)
    eigval2<-ordBC2$values$Eigenvalues # to get eigen values of every PC
    eigval2<- round((eigval2/sum(eigval2))*100, 1) # to get variation explained by every PC
    
    # plot with reactor IDs
    plot_ordination(data.sqrt_prop2, ordBC2 , color="Influent",
                    shape= "Reactor_scale") +
      geom_point(size=1.5) +
      geom_point(size=4.5, alpha= 0.6) +
      theme_classic(base_size = 9) +
      stat_ellipse(linewidth=0.18, aes(color=Influent)) +
      geom_text(aes(label=sample_data(data.sqrt_prop2)$Influent_extra),
                color="grey20", size=2, show.legend = FALSE,
                lineheight = 0.8) +
      theme(legend.margin = margin(10,0,10,0),
            legend.text = element_text(size=7.8)) +
      # guides(color="none")+
      labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,"), without reactor 15"),
           color="Influent type",
           shape="Reactor scale",
           caption="The labels on the points are the reactor IDs",
           x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
    ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_InfluentEllipse_NO15_IDs.png"), width = 5.5, height = 4, dpi=300)
    
    # plot with reactors details
    plot_ordination(data.sqrt_prop2, ordBC2 , color="Influent",
                    shape= "Reactor_scale") +
      geom_point(size=1.5) +
      geom_point(size=4.5, alpha= 0.6) +
      theme_classic(base_size = 9) +
      stat_ellipse(linewidth=0.18, aes(color=Influent)) +
      #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE)
      geom_text(aes(label=sample_data(d2)$Reactor_number), color="grey20", size=2.8, show.legend = FALSE) +
      theme(legend.margin = margin(10,0,10,0),
            legend.text = element_text(size=7.8)) +
      # guides(color="none")+
      labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,"), without reactor 15"),
           color="Influent type",
           shape="Reactor scale",
           caption="The labels on the points are the reactor details",
           x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
    ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_InfluentEllipse_NO15_details.png"), width = 5.5, height = 4, dpi=300)
    
  } # end of the "if" related to reactor 15
  
  
  
  # EXTRA CHECK: ONLY FULL SCALES AND PILOT OF PROJ 1
  d2 <- subset_samples(d, Reactor_scale=="Full" | Reactor_number=="1")
  data.prop.labels<-d2
  sample_data(data.prop.labels)$Influent_extra<- gsub("small granules","Denmark\n(small)",sample_data(data.prop.labels)$Influent_extra)
  sample_data(data.prop.labels)$Influent_extra<- gsub("large granules","Denmark\n(large)",sample_data(data.prop.labels)$Influent_extra)
  sample_data(data.prop.labels)$Influent_extra[sample_data(data.prop.labels)$FASTQ_ID=="1623908"] <- "Tuscany\n(mature)"
  {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
    DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
    ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
    eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
    eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
  }
  plot_ordination(data.sqrt_prop, ordBC , color="Project_number", 
                  shape= "Reactor_scale") +
    # scale_color_manual(values=c("AAAA"="coral","BBBB"="deepskyblue")) +
    geom_point(size=1.5) +
    geom_point(size=4.5, alpha= 0.6) +
    theme_classic(base_size = 9) +
    geom_text(aes(label=sample_data(data.prop.labels)$Influent_extra), 
              color="grey20", size=2, show.legend = FALSE,
              lineheight = 0.8) +
    theme(legend.margin = margin(10,0,10,0),
          legend.text = element_text(size=7.8)) +
    # guides(color="none")+
    labs(title=paste0("PCoA on ",plot_title," using Hellinger distance\n (euclidean on Hellinger transformed ",levels,")"), 
         color="Project",
         shape="Scale",
         x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
  ggsave(file=paste0("Results/",subset,"/PCoA_PERMANOVA_Hellinger_on_",file_name,"_",levels,"_FULLSCALES.png"), width = 5.5, height = 4, dpi=300)
  
  
  
}  # end of the loop

rm(d, d2, data.sqrt_prop, data.sqrt_prop2, data.prop.labels, DistBC,ordBC, eigval)





################# VEEN DIAGRAMS FULL vs PILOT vs SCALE (EVERY MICROBIAL GENUS, no virus) ##########################

data.venn<-subset_taxa(data.genus.prop, !Genus %in% c("unclassified") )
data.venn<-subset_taxa(data.venn, !Domain %in% c("Viruses") )
threshold <- 0.0025
always<-taxa_names(filter_taxa(data.venn, function(x) max(x) >= threshold, TRUE))
data.venn<-prune_taxa(always, data.venn)
sample_percent <- 75  # used below

Full<-subset_samples(data.venn, Reactor_scale=="Full")
Full<- prune_taxa(taxa_sums(Full) > 0, Full) 
who<-apply(otu_table(Full), MARGIN=2, function(x) ifelse(x>0, 1, 0)) # if more than 0 --> "a point"
how_many_samples <- round( length(sample_names(Full))* sample_percent / 100 , 0 )
who<- who[rowSums(who)>=how_many_samples, ] # more than X "points" --> at least in the X% of samples
Full<-as.vector(tax_table(Full)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

Lab<-subset_samples(data.venn, Reactor_scale=="Lab")
Lab<- prune_taxa(taxa_sums(Lab) > 0, Lab)
#always<-taxa_names(filter_taxa(Lab, function(x) min(x) > threshold, TRUE))
#Lab<-prune_taxa(always, Lab)
who<-apply(otu_table(Lab), MARGIN=2, function(x) ifelse(x>0, 1, 0))
how_many_samples <- round( length(sample_names(Lab))* sample_percent / 100 , 0 )
who<- who[rowSums(who)>=how_many_samples, ]
Lab<-as.vector(tax_table(Lab)[row.names(who),"Genus"])

Pilot<-subset_samples(data.venn, Reactor_scale=="Pilot")
Pilot<- prune_taxa(taxa_sums(Pilot) > 0, Pilot)
#always<-taxa_names(filter_taxa(Pilot, function(x) min(x) > threshold, TRUE))
#Pilot<-prune_taxa(always, Pilot)
#Pilot<-as.character(tax_table(prune_taxa(taxa_sums(Pilot)>threshold, Pilot))[,"Genus"])
who<-apply(otu_table(Pilot), MARGIN=2, function(x) ifelse(x>0, 1, 0))
how_many_samples <- round( length(sample_names(Pilot))* sample_percent / 100 , 0 )
who<- who[rowSums(who)>=how_many_samples, ] 
Pilot<-as.vector(tax_table(Pilot)[row.names(who),"Genus"])


ONLY_IN_Lab<- Lab[! Lab %in% Pilot & ! Lab %in% Full]
ONLY_IN_Lab<- paste(ONLY_IN_Lab, collapse = ", ")
head(ONLY_IN_Lab)

ONLY_IN_Full<- Full[! Full %in% Pilot & ! Full %in% Lab]
ONLY_IN_Full<- paste(ONLY_IN_Full, collapse = ", ")
head(ONLY_IN_Full)

ONLY_IN_Pilot<- Pilot[! Pilot %in% Full & ! Pilot %in% Lab]
ONLY_IN_Pilot<- paste(ONLY_IN_Pilot, collapse = ", ")
head(ONLY_IN_Pilot)


x<-list(Full=Full,Lab=Lab,Pilot=Pilot)
ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse2","coral","deepskyblue")) +
  labs(
    #title="Microbial genera featured in every sample\n"
    ) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) )
ggsave(filename = paste0("Results/Venn_Diagramm_MICROBIAL_GENERA_Reactor_scales_over",threshold,"min.png"), width = 4, height = 4, dpi=300, bg = "white")
dev.off()


common <- Full[ Full %in% Pilot & Full %in% Lab ]
# length(common)
data_Genus_common<- subset_taxa(data.venn, Genus %in% common)
write.csv2(tax_table(data_Genus_common), file = "Results/EVERY_COMMON_MICROBE_TAXONOMY.csv" , row.names = F, quote = F)



# Exporting common obs beyond bacteria
data.temp <- subset_taxa(data_Genus_common, !Domain %in% "Bacteria" )
data.temp <- subset_taxa(data.temp, !Genus %in% "archaeon" )
to_save<- cbind.data.frame("Overall_averages"= round( as.numeric(apply(otu_table(data.temp),1,mean)), 5),
                           "Averages_Full_scales"= round( as.numeric(apply (otu_table(subset_samples(data.temp, Reactor_scale=="Full")),1,mean)), 5),
                           "Averages_Pilot_scales"= round( as.numeric(apply (otu_table(subset_samples(data.temp, Reactor_scale=="Pilot")),1,mean)), 5),
                           "Averages_Lab_scales"= round( as.numeric(apply (otu_table(subset_samples(data.temp, Reactor_scale=="Lab")),1,mean)), 5),
                           "Genus"= as.data.frame(tax_table(data.temp))[["Genus"]],
                           "Euk_clade"= as.data.frame(tax_table(data.temp))[["Euk_clade"]],
                           "Domain"= as.data.frame(tax_table(data.temp))[["Domain"]]
)
to_save <- to_save[order(to_save$Overall_average, decreasing = T), ]
write.csv2( to_save ,
            file= "Results/Abundances/COMMON_EUKA_OR_ARCHAEA.csv",
            quote=F, row.names = F
            )




# ANNOTATING THE BACTERIA FEATURED ONLY IN FULL SCALES AND PILOT SCALES
FULL_and_PILOT_only <- Full[ !Full%in%Lab ]
FULL_and_PILOT_only <- FULL_and_PILOT_only[ !FULL_and_PILOT_only %in% Pilot[Pilot%in%Lab] ]
data.temp <- subset_taxa(data.venn, Genus %in% FULL_and_PILOT_only )
# View(otu_table(data.temp))
temp <- cbind.data.frame( "Overall_mean"=round(apply(otu_table(data.temp),1,mean),4) , tax_table(data.temp)[,c("Euk_clade","Genus")] )
temp <- temp[order(temp$Overall_mean, decreasing = T) , ]
write.csv2( temp ,
            file= "Results/Abundances/ONLY_IN_FULL_OR_PILOT_SBR_TAXA.csv",
            quote=F, row.names = F)



suppressWarnings(rm(ONLY_IN_Lab, ONLY_IN_Full, ONLY_IN_Pilot,
                    x, con, Lab, Full, Pilot, data.venn, who))




############### ABUNDANCES OF THE MOST ABUNDANT *COMMON* GENUS  (FULL vs PILOT vs SCALE ) ###################

### TOP COMMON Genus
suppressWarnings(rm(top, others, tabella, unass_data))
# unclassified already removed during the venn paragraph...

data_temp<-tax_glom(data_Genus_common, taxrank = "Genus", NArm = F) # to avoid a graphical glitch of ggplot2
data_temp_prop<- transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )
# MAX_THRESH<- 0.1
# data_temp_prop <- filter_taxa(data_temp_prop, function(x) max(x) >= MAX_THRESH, TRUE)  # ONLY WITH ABUND HIGHER THAN X
# data_temp_prop<- transform_sample_counts(data_temp, function(x) (x/sum(x))*100 )   # re-computing the proportions

{ top_data <- data_temp_prop  # using the data with the always featured genus to search for the most abundant ones
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:30]   # NB: common...
  prune.dat_top <- prune_taxa(top,data.genus.prop) # ... but using the COMPLETE object to show the (true) percentages of the observations in the whole samples
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  # for(i in 1:length(tabella_top$Genus)){ # to mark the fungal taxa
  #   tabella_top$Genus[i]<-paste0(tabella_top$Genus[i],"_",tabella_top$Domain[i])
  # }
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the genus
  tabella<-aggregate(.~Genus+Reactor_scale+Reactor_number+Influent+Sample_name, tabella[ , c("Reactor_scale","Reactor_number","Influent", "Sample_name","Abundance","Genus")], FUN=sum)
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)  # to save horizontal space in the plot!
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_grid(~ Reactor_scale, scales = "free", space="free") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =9.5) +
  scale_fill_manual(values= c("coral4", fill_color_30_BIS)) +
  scale_y_break(c(10, 65),
                scales=7,
                space = 0.35,
                expand = c(0,-1)
  ) +
  scale_y_continuous(
    expand = c(0,-1),
    breaks = c(0,5,10,seq(65,100,5)) ,
    limits = c(0,101)
  ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.title.y = element_text(vjust=0.85, size=10.5, hjust=0.65),
  strip.text = element_text(size=11),
  # legend.key.spacing.y = unit(0.1, "cm"),
  # legend.key.spacing.x = unit(1.5, "cm"),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9 ),
  legend.position="bottom",
  legend.margin = margin(-5,25,-20,0),
  plot.margin = margin(1,2,1,2),
  panel.spacing.x = unit(1.45,"pt")
  ) +
  scale_x_discrete(expand=c(-0.02,1)) +
  # scale_y_continuous(expand=c(0,1) ) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance of clades",
       fill="" #,
       #title = "Most abundant COMMON microbial genus",
       #caption = " 'Others' includes every other common genus across samples"
  )
ggsave(file="Results/Abundances/TOP_COMMON_microbial_genera_WITH_MAX_OVER_THRESHOLD.png",width=6.6,height=4.5, dpi=300)
dev.off()


to_save<- cbind.data.frame("Overall_average_across_common"=round( as.numeric(apply(otu_table(prune.dat_top),1,mean)) ,3 ),
                           # "Overall_median"=round( as.numeric(apply(otu_table(prune.dat_top),1,median)) ,3 ), 
                           "Averages_Full_scales"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Full")),1,mean)), 3),
                           "Averages Pilot scale"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Pilot")),1,mean)), 3),
                           "Averages_Lab_scales"=round( as.numeric(apply (otu_table(subset_samples(prune.dat_top, Reactor_scale=="Lab")),1,mean)), 3),
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]],
                           "Domain"= as.data.frame(tax_table(prune.dat_top))[["Domain"]]
)
to_save$Genus<- gsub("_"," ", to_save$Genus)
to_save$Genus<- gsub("bacterium ","", to_save$Genus)
to_save$Main_role<-""
to_save$Main_role[to_save$Genus%in% PAO_list] <- "PAO"
to_save$Main_role[to_save$Genus%in% GAO_list] <- "GAO"
to_save$Main_role[to_save$Genus%in% AOB_list] <- "AOO"
to_save$Main_role[to_save$Genus%in% NOB_list] <- "NOO"
# to_save$Main_role[to_save$Genus%in% table_PHA[,1]] <- "PHA producer"
to_save <- to_save[order(to_save$Overall_average_across_common, decreasing = T), ]

write.csv(to_save, 
          file = "Results/Abundances/TOP_COMMON_microbial_genera_Averages_no_unclassified.csv", 
          row.names = F, quote=F
          )



# ### Common Genera INCLUDING also every other genus (also not common)
# suppressWarnings(rm(top, others, tabella, unass_data))
# data_temp<- subset_taxa(data.genus.prop, Genus!="unclassified")
# data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
# tabella<- psmelt(data_temp)
# original_name_sequence<- tabella$Genus
# tabella$Genus[original_name_sequence %in% common]  <- "Common genera"
# tabella$Genus[! original_name_sequence %in% common]  <- "Other genera"
# tabella <- aggregate( Abundance ~ Genus+ Sample_name + Genus + Reactor_scale, tabella, FUN= sum)
# ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
#   # facet_nested(~ Reactor_scale + Reactor_number, scales = "free_x", space = "free_x") +
#   facet_nested(~ Reactor_scale , scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") +
#   theme_classic(base_size =8.5) +
#   scale_fill_manual(values=c("deepskyblue2","gray45")) +
#   scale_y_continuous(expand=c(0,1) ) +
#   theme(axis.text.x=element_text(angle=45,
#                                  vjust=1,
#                                  hjust=1,
#                                  size= 6
#   ),
#   axis.title =element_text(size=10),
#   axis.title.x = element_text(vjust=2.5),
#   # strip.text = element_text(size=6.7),
#   legend.spacing.y  = unit(0.5, "cm"),
#   legend.spacing.x  = unit(1.5, "cm"),
#   legend.key.height = unit(0.2, "cm"),
#   legend.key.width = unit(0.4, "cm"),
#   legend.text = element_text ( size = 10.8 ),
#   legend.position="bottom",
#   legend.margin = margin(-12,41,0,0),
#   plot.margin = margin(2,4,2,15)) +
#   guides(fill=guide_legend(nrow=1)) +
#   labs(x="", y="Percentual abundance",
#        fill="",
#        title = "Common vs uncommon genera")
# ggsave(file="Results/Abundances/COMMON_microbial_genera_vs_uncommon_BUT_OVER_THRESHOLD.png",width=6.6,height=5, dpi=300)
# dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))

