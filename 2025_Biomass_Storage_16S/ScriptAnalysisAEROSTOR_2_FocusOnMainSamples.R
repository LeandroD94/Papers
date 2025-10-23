##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("DESeq2")
  library("ALDEx2")
  library("vegan")
  library("ecodist")
  # utilities
  library("reshape")
  library("xlsx")  
  library("qiime2R")
}

options(scipen = 100) # disable scientific annotation

load("data_prepared_for_analysis.RData")   # see script part 1




########################### FOCUS ON THE *AS* SAMPLE WITH FRIDGE-FROZEN ##########################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_samples(data.genus.prop, Reactor_ID=="A 11/23" )
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
    tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  } else {
    stop("Something went wrong when updating the taxonomy")
  }
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("aceae$","ac.", tabella$Genus)
  tabella$Genus<-gsub ("Chitinophagales","Chitinophag.", tabella$Genus)
  tabella$Genus<-gsub ("f Xanthobacteraceae","Xanthobacter.", tabella$Genus)
  tabella$Genus<-gsub ("marinales","marin.", tabella$Genus)
  tabella$Genus<-gsub ("Microtrichales","Microtrich.", tabella$Genus)
  tabella$Genus<-gsub ("uncultured_ ","uncult. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
  levels(tabella$Storage)<-gsub ("_"," ", levels(tabella$Storage))
}
levels(tabella$Storage) <- gsub("C"," °C +Gly", levels(tabella$Storage) )
levels(tabella$Storage)  <- gsub("+Gly noGly","no Gly", levels(tabella$Storage) , fixed=T)
levels(tabella$Storage)  <- gsub("4 °C +Gly","4 °C", levels(tabella$Storage) , fixed=T)

fill_color_30_customised<-c("darkblue","deeppink", "firebrick3",
                     "darkcyan","darkmagenta", 
                     "brown4","grey55", "aquamarine3",
                     "bisque","cyan","darkgreen","yellow3",
                     "springgreen4", "grey20", 
                     "deepskyblue4","lightgreen","orange2",
                     "darkorchid", "lightblue1" , "violet", "blue",
                     "yellow", "pink3","yellow4", "chocolate4","coral2","red",
                     "lightgrey", "green2","darkslategray3")
ggplot(data=tabella, aes(x=Time_before_extr2, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_30_customised) +
  facet_grid( . ~ Storage ,
                scales = "free_x", space="free"
                #strip = strip_nested(size = "variable")
            ) +
  theme_classic(base_size =11) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1.25,
                                 hjust=1,
                                 size= 7.5
  ),
  axis.text.y=element_text(size=7),
  axis.ticks.x = element_blank(),
  strip.text.x = element_text(size=8.75, angle=0),
  legend.key.height = unit(0.14, "cm"),
  legend.key.width = unit(0.28, "cm"),
  legend.text = element_text ( size = 8 , face="italic"  ),
  legend.position="bottom",
  legend.margin = margin(-14 ,34, 1 ,0),
  plot.margin = margin(2,1,2,1),
  panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results/Abundances/AS_FRIDGEvsFROZEN_TOP_genera.png",width=5.2,height=4.25, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



#### PCoA on genera
data.prop.labels<-subset_samples(data.genus.prop, Reactor_ID=="A 11/23" )
sample_data(data.prop.labels)$Glycerol <- gsub(" ","no Gly", sample_data(data.prop.labels)$Glycerol)
sample_data(data.prop.labels)$Glycerol <- gsub("noGly","no Gly", sample_data(data.prop.labels)$Glycerol)
levels(sample_data(data.prop.labels)$Storage) <- gsub("_noGly","", levels(sample_data(data.prop.labels)$Storage) )
Time<- gsub("m","", as.character(sample_data(data.prop.labels)$Time_before_extr2))
Time_black <- Time
Time_black [ as.character(sample_data(data.prop.labels)$Storage)=="-80C" ] <- ""
Time_white <- Time
Time_white [ as.character(sample_data(data.prop.labels)$Storage)!="-80C" ] <- ""
# levels(sample_data(data.prop.labels)$Time_before_extr2) <- gsub("m","", levels(sample_data(data.prop.labels)$Time_before_extr2))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral2",
                "4C" = "yellow2",
                "-20C" = "deepskyblue",
                "-80C" = "darkblue",
                "No_storage"="red3")
plot_ordination(data.sqrt_prop, ordBC, color="Storage" , shape = "Glycerol" )+
  scale_shape_manual( values = c(15,16)) +
  geom_point(size=4.52, alpha=0.4,  show.legend = F) +
  theme_classic(base_size = 11) +
  theme(legend.text =element_text(size=10.8)) +
  geom_text(aes(label=Time_black),
            color="black",
            size=4.25, hjust=0.45,show.legend = FALSE) +
  geom_text(aes(label=Time_white),
            color="gray95",
            size=4.25, hjust=0.55,show.legend = FALSE) +
  labs(
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="", shape="") +
  scale_color_manual(values= values_color)
ggsave(file="Results/Beta_diversity_PCoA/AS_FridgeFrozen_PCoA_Hellinger_Genus.png", width = 4.2, height = 3.85, dpi=300)




###################### FOCUS ON THE *PN-AGS2* SAMPLE WITH FRIDGE-FROZEN #######################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_samples(data.genus.prop, Reactor=="PN_AGS2" )
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
    tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  } else {
    stop("Something went wrong when updating the taxonomy")
  }
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("Chitinophagales","Chitinophag.", tabella$Genus)
  tabella$Genus<-gsub ("aceae$",".", tabella$Genus)
  tabella$Genus<-gsub ("f Xanthobacteraceae","Xanthobacter.", tabella$Genus)
  tabella$Genus<-gsub ("marinales","marin.", tabella$Genus)
  tabella$Genus<-gsub ("Microtrichales","Microtrich.", tabella$Genus)
  tabella$Genus<-gsub ("uncultured_ ","uncult. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
  levels(tabella$Storage)<-gsub ("_"," ", levels(tabella$Storage))
}
levels(tabella$Storage) <- gsub("C"," °C +Gly", levels(tabella$Storage) )
levels(tabella$Storage)  <- gsub("+Gly noGly","no Gly", levels(tabella$Storage) , fixed=T)
levels(tabella$Storage)  <- gsub("4 °C +Gly","4 °C", levels(tabella$Storage) , fixed=T)

fill_color_30_customised<-c("chartreuse","deeppink", "firebrick3",
                     "darkcyan","darkmagenta", 
                     "grey55", "aquamarine3",
                     "bisque","chocolate4","cyan","yellow3","brown4",
                     "springgreen4", "grey20", 
                     "deepskyblue4","lightgreen","orange2",
                     "darkorchid", "lightblue1" , "violet", "deepskyblue",
                     "yellow", "pink3","yellow4", "red","darkgreen","coral2",
                     "lightgrey", "darkblue","darkslategray3")
ggplot(data=tabella, aes(x=Time_before_extr2, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_30_customised) +
  facet_grid( . ~ Storage ,
              scales = "free_x", space="free"
              #strip = strip_nested(size = "variable")
  ) +
  theme_classic(base_size =11) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1.25,
                                 hjust=1,
                                 size= 7.5
  ),
  axis.text.y=element_text(size=7),
  axis.ticks.x = element_blank(),
  strip.text.x = element_text(size=8.75, angle=0),
  legend.key.height = unit(0.14, "cm"),
  legend.key.width = unit(0.2, "cm"),
  legend.text = element_text ( size = 8.05 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-14.5 ,33, 1 ,0),
  plot.margin = margin(2,1,2,1),
  panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results/Abundances/PNAGS2_TOP_genera_FRIDGEvsFROZEN.png",width=5.2,height=4.25, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



#### PCoA on genera
data.prop.labels<-subset_samples(data.genus.prop, Reactor=="PN_AGS2" )
sample_data(data.prop.labels)$Glycerol <- gsub(" ","no Gly", sample_data(data.prop.labels)$Glycerol)
sample_data(data.prop.labels)$Glycerol <- gsub("noGly","no Gly", sample_data(data.prop.labels)$Glycerol)
levels(sample_data(data.prop.labels)$Storage) <- gsub("_noGly","", levels(sample_data(data.prop.labels)$Storage) )
levels(sample_data(data.prop.labels)$Storage) <- gsub("C"," °C", levels(sample_data(data.prop.labels)$Storage) )
Time<- gsub("m","", as.character(sample_data(data.prop.labels)$Time_before_extr2))
Time_black <- Time
Time_black [ as.character(sample_data(data.prop.labels)$Storage)=="-80 °C" ] <- ""
Time_white <- Time
Time_white [ as.character(sample_data(data.prop.labels)$Storage)!="-80 °C" ] <- ""
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral2",
                "4 °C" = "yellow2",
                "-20 °C" = "deepskyblue",
                #"-20C_noGly" = "deepskyblue",
                "-80 °C" = "darkblue",
                #"-80C_noGly" = "royalblue3",
                "No_storage"="red3")

plot_ordination(data.sqrt_prop, ordBC, color="Storage" , shape = "Glycerol" )+
  scale_shape_manual( values = c(15,16)) +
  geom_point(size=4.6, alpha=0.4,  show.legend = F) +
  theme_classic(base_size = 11) +
  theme(legend.text =element_text(size=10.8)) +
  geom_text(aes(label=Time_black),
            color="black",
            size=4.2, hjust=0.5,show.legend = FALSE) +
  geom_text(aes(label=Time_white),
            color="gray92",
            size=4.25, hjust=0.51,show.legend = FALSE) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="", shape="") +
  scale_color_manual(values= values_color)
ggsave(file="Results/Beta_diversity_PCoA/PNAGS2_FridgeFrozen_PCoA_Hellinger_Genus.png", width = 4.2, height = 3.85, dpi=300)




############# FOCUS ON THE CUOIODEPUR SAMPLE WITH "no STORAGE" and 3 Y #############

##### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_samples(data.genus.prop, Reactor_ID=="A 07/22" )
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,top_data)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
    tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  } else {
    stop("Something went wrong when updating the taxonomy")
  }
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("Flavobacteriaceae","Flavobacter.", tabella$Genus)
  tabella$Genus<-gsub ("Chitinophagales","Chitinophag.", tabella$Genus)
  tabella$Genus<-gsub ("Steroidobacteraceae","Steroidobact.", tabella$Genus)
  tabella$Genus<-gsub ("Peptostreptococcales-Tissierellales","Peptostr-Tiss.", tabella$Genus)
  tabella$Genus<-gsub ("Microtrichales","Microtrich.", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
levels(tabella$Storage)<-gsub ("No_storage","Room T", levels(tabella$Storage))
levels(tabella$Storage) <- gsub("C"," °C +Gly", levels(tabella$Storage) )
levels(tabella$Storage)  <- gsub("+Gly noGly","no Gly", levels(tabella$Storage) , fixed=T)
levels(tabella$Storage)  <- gsub("4 °C +Gly","4 °C", levels(tabella$Storage) , fixed=T)

fill_color_30_customised<-c("gray15","deeppink", "firebrick3",
                     "darkcyan","darkmagenta", 
                     "darkblue","grey85", "aquamarine3",
                     "bisque","springgreen4","yellow3","brown4",
                     "cyan", "grey50", 
                     "deepskyblue4","lightgreen","orange2",
                     "darkorchid", "lightblue1" , "violet", "blue",
                     "yellow", "pink3","yellow4", "chocolate4","coral2",
                     "lightgrey", "red","green2","darkslategray3")
ggplot(data=tabella, aes(x=Time_before_extr2, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_30_customised) +
  facet_nested( . ~ Storage ,
                scales = "free_x", space="free",
                strip = strip_nested(size = "variable")
  ) +
  theme_classic(base_size =11) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1.35,
                                 hjust=1,
                                 size= 7.65),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=7),
        strip.text.x = element_text(size=9, angle=0),
        legend.key.height = unit(0.14, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.text = element_text ( size = 8.05 , face="italic" ),
        legend.position="bottom",
        legend.margin = margin(-15 ,33, 1 ,0),
        plot.margin = margin(2,1,2,1),
        panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="")
ggsave(file="Results/Abundances/ASmain_TOP_genera_SAMPLE_WITH_ROOMT.png",width=5.2,height=4.25, dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



##### PCoA

data.prop.labels<-subset_samples(data.genus.prop, Reactor_ID=="A 07/22" )

# updating labels
sample_data(data.prop.labels)$Storage <- gsub("No_storage","Room T", sample_data(data.prop.labels)$Storage)
sample_data(data.prop.labels)$Storage <- gsub("C"," °C", sample_data(data.prop.labels)$Storage)
levels(sample_data(data.prop.labels)$Time_before_extr2) <- gsub("m","",levels(sample_data(data.prop.labels)$Time_before_extr2))
levels(sample_data(data.prop.labels)$Time_before_extr2)[ levels(sample_data(data.prop.labels)$Time_before_extr2)=="6" ] <- "\n  6"  # to avoid overlapping
levels(sample_data(data.prop.labels)$Time_before_extr2)[ levels(sample_data(data.prop.labels)$Time_before_extr2)=="1y" ] <- "1y    "  # to avoid overlapping
levels(sample_data(data.prop.labels)$Time_before_extr2)[ levels(sample_data(data.prop.labels)$Time_before_extr2)=="3y" ] <- "\n3y"  # to avoid overlapping
sample_data(data.prop.labels)$Collection_data2 <- as.character(sample_data(data.prop.labels)$Collection_data)
sample_data(data.prop.labels)$Collection_data2[ sample_data(data.prop.labels)$Storage!="Direct" ] <- ""
# new groups to separate the samples of the same reactor in geom path...
new_line1 <- 1:length(sample_data(data.prop.labels)$Reactor_ID)  # numeri causali, tutti diversi, non saranno collegati
new_line1[sample_data(data.prop.labels)$Storage %in% "Room T" | sample_data(data.prop.labels)$Sample_name=="A1" ] <- "Stesso_gruppo_stessa_linea"
new_line2 <- as.character(sample_data(data.prop.labels)$Reactor_ID)
new_line2[sample_data(data.prop.labels)$Storage %in% "Room T"] <- c("B","A")

#computation
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "4 °C" = "lightblue2",
                "-20 °C" = "deepskyblue3",
                "Room T"="red2")
plot_AS <-plot_ordination(data.sqrt_prop, ordBC, color="Storage" )+
  guides(shape="none") +
  geom_point(size=4.5, alpha=0.5,  show.legend = F) +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9.65),
        legend.margin = margin(1,1,1,-10)
  ) +
  geom_text(aes(label=Time_before_extr2),
            color="black", size=4.45, show.legend = FALSE,
            #hjust=0.1,
            lineheight=0.45
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Collection_data2), color="gray25", 
            size=3, vjust=-0.25, hjust=0.52,
            show.legend = FALSE, lineheight=0.87
  ) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="") +
  scale_color_manual(values= values_color) +
  geom_path(aes(group=new_line1), color="darkgray" , linewidth = 0.25, alpha=0.8,
            arrow=arrow(length =unit(0.285,"cm"), type = "closed")
  ) +
  geom_path(aes(group=new_line2), color="darkgray" , linewidth = 0.25, alpha=0.8,
            arrow=arrow(length =unit(0.285,"cm"), type = "closed"),
  )

new_order <- c("1408525","1408527","1408531","1408529","1744973","1744980","1778921")  # to ensure the proper sequence of the arrows through the rownames (FASTQ ID)
plot_AS$data <- plot_AS$data[ new_order ,]

plot_AS
ggsave(file="Results/Beta_diversity_PCoA/ASmain_PCoA_Hellinger_Genus_WithNoStorage.png", width = 4.2, height = 3.65, dpi=300)


