# Part 3: Analysing the differences between directly extracted samples and stored samples


##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggvenn")
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




####################### EVERY SAMPLE PCoA ##########################

### PCOA on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Reactor <- gsub("CuoioDepur_Aero","AS tank", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Reactor <- gsub("Pilot_AGS","Pilot AGS", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Reactor <- gsub("PN_AGS","Lab PN AGS", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Store_or_not <- sample_data(data.prop.labels)$Storage
levels(sample_data(data.prop.labels)$Store_or_not) <- gsub("No_storage","Room T", levels(sample_data(data.prop.labels)$Store_or_not) )
levels(sample_data(data.prop.labels)$Store_or_not) <- gsub("C"," °C", levels(sample_data(data.prop.labels)$Store_or_not) )
levels(sample_data(data.prop.labels)$Store_or_not) <- gsub("_noGly","", levels(sample_data(data.prop.labels)$Store_or_not) )
# levels(sample_data(data.prop.labels)$Store_or_not) <- gsub("-80C_noGly","Frozen", levels(sample_data(data.prop.labels)$Store_or_not) )
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "4 °C" = "yellow2",
                "-20 °C" = "deepskyblue3",
                #"-20 °C (no Gly)" = "cyan",
                "-80 °C" = "darkblue",
                #"-80 °C (no Gly)" = "darkblue",
                "Room T"="red4")
# plot_ordination(data.sqrt_prop, ordBC, color="Store_or_not" , shape="Reactor")+
plot_ordination(data.sqrt_prop, ordBC, color="Store_or_not" , shape="Glycerol")+
  geom_point(size=3, alpha=0.45,  show.legend = F) +
  theme_classic(base_size = 11) +
  theme(legend.text =element_text(size=10),
        legend.key.size = unit(0.5,"cm"),
        legend.position = "bottom"
        ) +
  scale_shape_manual(values=c(1,16,17)) +
  # geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="black", size=3.85, show.legend = FALSE) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="", shape="") +
  scale_color_manual(values= values_color)
# ggsave(file="Results/Beta_diversity_PCoA/EVERY_SAMPLE_PCoA_Hellinger_Genus.png", width = 4.25, height = 4, dpi=300)
ggsave(file="Results/Beta_diversity_PCoA/EVERY_SAMPLE_PCoA_Hellinger_Genus_LAAARGE.png", width = 6, height = 4, dpi=300)




### PCOA ONLY AGS
data.prop.labels<- subset_samples(data.genus.prop, Reactor!="CuoioDepur_Aero")
sample_data(data.prop.labels)$Reactor <- gsub("Pilot_AGS","Pilot AGS", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Reactor <- gsub("PN_AGS","Lab PN AGS", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Biomass[!sample_data(data.prop.labels)$Sample_name %in% c("PG6","L3","L6","L10")] <- "" # to avoid text overlapping in the plot
sample_data(data.prop.labels)$Store_or_not <- sample_data(data.prop.labels)$Storage
levels(sample_data(data.prop.labels)$Store_or_not) <- gsub("-20C","Frozen", levels(sample_data(data.prop.labels)$Store_or_not) )
levels(sample_data(data.prop.labels)$Store_or_not) <- gsub("-80C","Frozen", levels(sample_data(data.prop.labels)$Store_or_not) )
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "Frozen" = "deepskyblue2",
                "Frozen_noGly" = "cyan"
)
plot_ordination(data.sqrt_prop, ordBC, color="Store_or_not" , shape="Reactor")+
  geom_point(size=3.5, alpha=0.4,  show.legend = F) +
  theme_classic(base_size = 11) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Biomass), color="black", size=3.65, 
            hjust=0.5, vjust= 0,
            show.legend = FALSE) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="", shape="") +
  scale_color_manual(values= values_color)
ggsave(file="Results/Beta_diversity_PCoA/PCoA_Hellinger_Genus_AGS_only.png", width = 4.2, height = 4, dpi=300)




######################## ALPHA DIVERSITY (-20 gly) ##############################

data.genus.temp<-subset_samples(data.genus, Reactor%in%c("CuoioDepur_Aero","PN_AGS2") )
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Sample_name
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<- mix_alpha$Shannon /log(mix_alpha$Observed)
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
{
  mix_alpha$Reactor_ID <- metadata[mix_alpha$Sample_name, "Reactor_ID" ]  # NB: the metadata has to be ordered according to time!
  mix_alpha$Storage <- metadata[mix_alpha$Sample_name, "Storage" ]
  mix_alpha$Reactor <- as.character(metadata[mix_alpha$Sample_name, "Reactor" ])
  mix_alpha$Sample_type <- metadata[mix_alpha$Sample_name, "Sample_type" ]
  mix_alpha$Collection_data <- metadata[mix_alpha$Sample_name, "Collection_data" ]
  mix_alpha$Collection_data <- factor(mix_alpha$Collection_data, levels = unique(metadata$Collection_data))
  mix_alpha$Glycerol <- metadata[mix_alpha$Sample_name, "Glycerol" ]
  mix_alpha$Time_before_extr <- metadata[mix_alpha$Sample_name, "Time_before_extr" ]
  mix_alpha$Time_before_extr <- gsub("None","",mix_alpha$Time_before_extr)
  mix_alpha$Time_before_extr <- gsub("month","",mix_alpha$Time_before_extr)
  mix_alpha$Time_before_extr <- gsub("year","y",mix_alpha$Time_before_extr)
  mix_alpha$Time_before_extr <- gsub("day","d",mix_alpha$Time_before_extr)
  mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
  mix_alpha$ID_storage <- paste0(mix_alpha$Reactor_ID, mix_alpha$Storage)
  levels(mix_alpha$Storage) <- gsub("-20C","-20 °C\n(Gly)", levels(mix_alpha$Storage) )
}

mix_alpha2 <- mix_alpha[ mix_alpha$Storage %in% c("Direct","-20 °C\n(Gly)") , ]  
mix_alpha2$Time_before_extr[ ! mix_alpha2$Reactor_ID %in% c("A 11/23","PN_AGS_4")] <- ""

mix_alpha2$Collection_data <- as.character((mix_alpha2$Collection_data))
mix_alpha2$Collection_data[ as.character(mix_alpha2$Reactor) %in% c("PN_AGS2") ] <- ""
mix_alpha2$Collection_data[ as.character(mix_alpha2$Storage)!= "Direct" ] <- ""
mix_alpha2$Collection_data <- gsub("/20","/", mix_alpha2$Collection_data , fixed = T)

ggplot(data=mix_alpha2, aes(y=value, x=Storage, color=Storage, fill=Storage)) +
  geom_boxplot( width=0.6,  alpha=0.15 , color="black", linewidth=0.15, outliers =F) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  scale_color_manual(values= c("Direct"="coral",
                               "-20C" = "deepskyblue3",
                               "-20 °C\n(Gly)" = "deepskyblue3" )
                     )+
  geom_path(aes(group=Reactor_ID), color="darkgrey", linewidth = 0.15, alpha=0.75,
  ) +
  geom_point(size =0.35 , aes(shape=Reactor, color=Storage, fill=Storage) ) +
  geom_point(size =0.95, alpha=0.5, aes(shape=Reactor, color=Storage, fill=Storage)) +
  theme_classic2(base_size = 8) +
  scale_x_discrete(expand = c(0.5,0)) + # to have more space on the borders
  geom_text(aes(label= Time_before_extr), color="black", size=1.65,
            hjust=-0.85, 
            # position=position_jitter(width=0.005,seed=1994),
            show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  # geom_text(aes(label= Collection_data), color="gray20", size=0.95,
  #           hjust=1.25,  show.legend = FALSE, lineheight =0.7) +
  theme(axis.text.x = element_text(size=6, angle=20, hjust=1, vjust=1),
        axis.text.y= element_text(angle=0, size=4),
        strip.text = element_text(size=5),
        #axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.3),
        panel.grid.minor.y = element_line(linewidth=0.15),
        plot.margin = margin(5,1,1,1)
  ) +
  labs(x="", y="Alpha Diversity Measure") +  
  guides(fill="none", color="none", shape="none")
ggsave(file="Results/Alfa_diversity_GENUS_level_BoxPlot_20Cgly.png", width = 2,height =2.4, dpi=300)

suppressWarnings(rm(mix_alpha) )




######################## ALPHA DIVERSITY (-80 C) ##############################

data.genus.temp<-subset_samples(data.genus, Reactor%in%c("Pilot_AGS","PN_AGS2") | Reactor_ID=="A 11/23" )
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$Sample_name
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<- mix_alpha$Shannon /log(mix_alpha$Observed)
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
{
  mix_alpha$Reactor_ID <- metadata[mix_alpha$Sample_name, "Reactor_ID" ]  # NB: the metadata has to be ordered according to time!
  mix_alpha$Storage <- metadata[mix_alpha$Sample_name, "Storage" ]
  mix_alpha$Reactor <- as.character(metadata[mix_alpha$Sample_name, "Reactor" ])
  mix_alpha$Sample_type <- metadata[mix_alpha$Sample_name, "Sample_type" ]
  mix_alpha$Collection_data <- metadata[mix_alpha$Sample_name, "Collection_data" ]
  mix_alpha$Collection_data <- factor(mix_alpha$Collection_data, levels = unique(metadata$Collection_data))
  mix_alpha$Glycerol <- metadata[mix_alpha$Sample_name, "Glycerol" ]
  mix_alpha$Time_before_extr <- metadata[mix_alpha$Sample_name, "Time_before_extr" ]
  mix_alpha$Time_before_extr <- gsub("None","",mix_alpha$Time_before_extr)
  mix_alpha$Time_before_extr <- gsub("month","",mix_alpha$Time_before_extr)
  mix_alpha$Time_before_extr <- gsub("year","y",mix_alpha$Time_before_extr)
  mix_alpha$Time_before_extr <- gsub("day","d",mix_alpha$Time_before_extr)
  mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
  mix_alpha$ID_storage <- paste0(mix_alpha$Reactor_ID, mix_alpha$Storage)
  levels(mix_alpha$Storage) <- gsub("C_noGly"," °C\n(no Gly)", levels(mix_alpha$Storage) )
}

mix_alpha2 <- mix_alpha[ mix_alpha$Storage %in% c("Direct","-80 °C\n(no Gly)") , ]
mix_alpha2$Time_before_extr[ mix_alpha2$Reactor %in% c("Pilot_AGS")] <- ""
mix_alpha2$LINE2 <- "a line missing otherwise"
mix_alpha2$LINE2[ ! mix_alpha2$Sample_name %in% c("A25","A40") ] <- ""

mix_alpha2$Collection_data <- as.character((mix_alpha2$Collection_data))
mix_alpha2$Collection_data[ as.character(mix_alpha2$Reactor) %in% c("PN_AGS2","CuoioDepur_Aero") ] <- ""
mix_alpha2$Collection_data[ as.character(mix_alpha2$Storage)!= "Direct" ] <- ""
mix_alpha2$Collection_data <- gsub("/20","/", mix_alpha2$Collection_data , fixed = T)

ggplot(data=mix_alpha2, aes(y=value, x=Storage, color=Storage, fill=Storage)) +
  geom_boxplot( width=0.6,  alpha=0.15 , color="black", linewidth=0.15, outliers =F) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  scale_shape_manual(values= c("CuoioDepur_Aero"=16,
                               "PN_AGS2" = 17,
                               "Pilot_AGS" = 15
                               )
  )+
  scale_color_manual(values= c("Direct"="coral",
                               "-80 °C\n(no Gly)" = "royalblue4"
                               )
  )+
  geom_path(aes(group=Reactor_ID), color="darkgrey", linewidth = 0.20, alpha=0.75 ) +
  geom_path(aes(group=LINE2), color="darkgrey", linewidth = 0.20, alpha=0.75 ) +
  geom_point(size =0.20 , aes(shape=Reactor, color=Storage, fill=Storage) ) +
  geom_point(size =0.95, alpha=0.5, aes(shape=Reactor, color=Storage, fill=Storage)) +
  theme_classic2(base_size = 8) +
  scale_x_discrete(expand = c(0.5,0)) + # to have more space on the borders
  geom_text(aes(label= Time_before_extr), color="black", size=1.65,
            hjust=-0.85, 
            # position=position_jitter(width=0.005,seed=1994),
            show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  # geom_text(aes(label= Collection_data), color="gray25", size=1.15,
  #              hjust=1.25,  show.legend = FALSE, lineheight =0.7) +
  theme(axis.text.x = element_text(size=6, angle=20, hjust=1, vjust=1),
        axis.text.y= element_text(angle=0, size=4),
        strip.text = element_text(size=5),
        panel.grid.major.y = element_line(linewidth=0.3),
        panel.grid.minor.y = element_line(linewidth=0.15),
        plot.margin = margin(5,1,1,1)
  ) +
  labs(x="", y="Alpha Diversity Measure") +  
  guides(fill="none", color="none", shape="none")
ggsave(file="Results/Alfa_diversity_GENUS_level_BoxPlot_80C.png", width = 2,height =2.4, dpi=300)

suppressWarnings(rm(mix_alpha) )




########################### PCoA BETA DIV (CUOIODEPUR STORAGE ONLY) ##########################

# on genera
data.prop.labels<-subset_samples(data.genus.prop, Reactor=="CuoioDepur_Aero" )

# new group to separate the samples of the same reactor in geom path...
new_line1 <- 1:length(sample_data(data.prop.labels)$Reactor_ID)  # numeri causali, tutti diversi, non saranno collegati
new_line1[sample_data(data.prop.labels)$Storage %in% "No_storage" | sample_data(data.prop.labels)$Sample_name=="A1" ] <- "Stesso_gruppo_stessa_linea"
new_line2 <- as.character(sample_data(data.prop.labels)$Reactor_ID)
new_line2[sample_data(data.prop.labels)$Storage %in% "No_storage"] <- c("A","B")
direct_extr_line <- 1:length(sample_data(data.prop.labels)$Reactor_ID)
direct_extr_line[sample_data(data.prop.labels)$Storage %in% "Direct" ] <- "Diretto_stessa_linea"


# updating labels
sample_data(data.prop.labels)$Storage<-gsub("No_storage","Room T", sample_data(data.prop.labels)$Storage)
sample_data(data.prop.labels)$Seq_batch<-as.character(sample_data(data.prop.labels)$Seq_batch)

# computation
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "4C" = "lightblue2",
                "-20C" = "deepskyblue3",
                "-20C_noGly" = "deepskyblue3",
                "-80C" = "darkblue",
                "-80C_noGly" = "darkblue",
                "Room T"="red4")
# plotting...
base_plot <- plot_ordination(data.sqrt_prop, ordBC, color="Storage" , shape = "Glycerol" )+
  geom_point(size=3.7, alpha=0.5,  show.legend = F) +
  guides(shape="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="black", size=2.5, show.legend = FALSE) +
  labs(
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="")
#base_plot$data <- base_plot$data[order(base_plot$data$Sample_order), ]   # aggiungere exp day per ordinare la direzione freccia
base_plot + scale_color_manual(values= values_color)+
  geom_path(aes(group=new_line1), color="darkgray" , linewidth = 0.18, alpha=0.9,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) + 
  geom_path(aes(group=new_line2), color="darkgray" , linewidth = 0.18, alpha=0.9,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  )
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_focus_on_TimeStorageS.png", width = 5, height = 4.5, dpi=300)
# # again, but with ellipses
# base_plot + stat_ellipse(linewidth=0.15) + scale_color_manual(values= values_color)
# ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_WITH_ELLIPSES.png", width = 5, height = 4.5, dpi=300)
# again, but each sample is a color ...
plot_ordination(data.sqrt_prop, ordBC, color="Reactor_ID") + # shape = "Gly or criostorage" +
  scale_color_manual(values= fill_color_10)+
  geom_path(aes(group=new_line1), color="darkgray" , linewidth = 0.18, alpha=0.9,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) + 
  geom_path(aes(group=new_line2), color="darkgray" , linewidth = 0.18, alpha=0.9,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) + 
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(shape="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="black", size=2.5, show.legend = FALSE) +
  labs(
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="")
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_EACH_react_a_color.png", width = 5, height = 4.5, dpi=300)
# with sample names
plot_ordination(data.sqrt_prop, ordBC, color="Reactor_ID") + # shape = "Gly" +
  scale_color_manual(values= fill_color_10)+
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(shape="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), color="black", size=3.2, show.legend = FALSE) +
  labs(
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="")
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_SAMPLE_NAME.png", width = 5, height = 4.5, dpi=300)
# again, checking the batches
plot_ordination(data.sqrt_prop, ordBC, color="Seq_batch", shape = "Seq_batch") +
  scale_color_manual(values= fill_color_15)+
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(shape="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="black", size=2.8, show.legend = FALSE) +
  labs(
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       color="Batches")
ggsave(file="Data_check/Beta_div_Hellinger_according_to_batch.png", width = 5, height = 4.5, dpi=300)



####### MAIN PCOA: CUOIODEPUR FROZEN ONLY, AVOIDING NO GLY

data.prop.labels<-subset_samples(data.genus.prop, Reactor=="CuoioDepur_Aero" )
data.prop.labels<-subset_samples(data.prop.labels, ! Storage %in% c("4C","-20C_noGly","-80C_noGly","-80C"))
# new groups to separate the samples of the same reactor in geom path...
new_line1 <- 1:length(sample_data(data.prop.labels)$Reactor_ID)  # numeri causali, tutti diversi, non saranno collegati
new_line1[sample_data(data.prop.labels)$Storage %in% "No_storage" | sample_data(data.prop.labels)$Sample_name=="A1" ] <- "Stesso_gruppo_stessa_linea"
new_line2 <- as.character(sample_data(data.prop.labels)$Reactor_ID)
new_line2[sample_data(data.prop.labels)$Storage %in% "No_storage"] <- c("B","A")
# reformatting the collection data
sample_data(data.prop.labels)$Collection_data<-gsub("^[0-9][0-9]/","", sample_data(data.prop.labels)$Collection_data)
sample_data(data.prop.labels)$Collection_data<-gsub("/202","/2", sample_data(data.prop.labels)$Collection_data, fixed=T)
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Storage!="Direct"] <-""
# updating labels
sample_data(data.prop.labels)$Storage <- gsub("No_storage","Room T", sample_data(data.prop.labels)$Storage)
#sample_data(data.prop.labels)$Storage <- gsub("-20C","Froozen", sample_data(data.prop.labels)$Storage)
#computation
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "4C" = "lightblue2",
                "-20C" = "deepskyblue3",
                "Room T"="red4")
plot_this<-plot_ordination(data.sqrt_prop, ordBC, color="Storage" )+
  guides(shape="none") +
  geom_point(size=3.1, alpha=0.5,  show.legend = F) +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=9)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="black", size=2.5, show.legend = FALSE) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="") +
  scale_color_manual(values= values_color) +
  geom_path(aes(group=new_line2), color="darkgray" , linewidth = 0.19, alpha=0.6,
            # arrow=arrow(length =unit(0.15,"cm"), type = "closed"),
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Collection_data), color="gray25", 
            size=2.1, vjust=0.6, hjust=0.55,
            show.legend = FALSE)
plot_this
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_OnlyGLY_20C.png", width = 4, height = 3.5, dpi=300)
# plot_this + theme(legend.position = "none")
# ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_OnlyGLY_20C_LAAARGE.png", width = 6.5, height = 3.5, dpi=300)



####### MAIN PCOA: CUOIODEPUR FROZEN ONLY, AVOIDING NO GLY, NO "no STORAGE"

data.prop.labels<-subset_samples(data.genus.prop, Reactor=="CuoioDepur_Aero" )
data.prop.labels<-subset_samples(data.prop.labels, ! Storage %in% c("4C","-20C_noGly","-80C_noGly","-80C","No_storage"))
data.prop.labels<-subset_samples(data.prop.labels, ! Sample_name == "A30" )  # to avoid overlap (the other time points are more than enough)
# new groups to separate the samples of the same reactor in geom path...
new_line2 <- as.character(sample_data(data.prop.labels)$Reactor_ID)
new_line2[sample_data(data.prop.labels)$Storage %in% "No_storage"] <- c("B","A")
# reformatting the collection data
sample_data(data.prop.labels)$Collection_data<-gsub("^[0-9][0-9]/","", sample_data(data.prop.labels)$Collection_data)
sample_data(data.prop.labels)$Collection_data<-gsub("/202","/2", sample_data(data.prop.labels)$Collection_data, fixed=T)
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Storage!="Direct"] <-""
sample_data(data.prop.labels)$Collection_data [ sample_data(data.prop.labels)$Collection_data =="02/24"] <- "\n02/24"  # to avoid overlapping
sample_data(data.prop.labels)$Collection_data [ sample_data(data.prop.labels)$Collection_data =="01/24"] <- "01/24  "

levels(sample_data(data.prop.labels)$Time_before_extr2) <- gsub("m","",levels(sample_data(data.prop.labels)$Time_before_extr2))
# updating labels
sample_data(data.prop.labels)$Storage <- gsub("No_storage","Room T", sample_data(data.prop.labels)$Storage)
#computation
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "-20C" = "deepskyblue3")
plot_this <- plot_ordination(data.sqrt_prop, ordBC, color="Storage" )+
  guides(shape="none") +
  geom_point(size=4, alpha=0.5,  show.legend = F) +
  theme_classic(base_size = 10.25) +
  theme(legend.text =element_text(size=8.5),
        plot.margin = margin(1,1,1,1),
        legend.margin = margin(1,1,1,-15)
        ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="black", 
            size=3.8, show.legend = FALSE) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="") +
  scale_color_manual(values= values_color) +
  geom_path(aes(group=new_line2), color="darkgray" , linewidth = 0.22, alpha=0.9,
            # arrow=arrow(length =unit(0.15,"cm"), type = "closed"),
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Collection_data), color="gray25", 
            size=2.25, vjust=0.65, hjust=0.52,
            show.legend = FALSE, lineheight=0.85
            )
plot_this
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_OnlyGLY_20C_withouthNOstorage.png", width = 3.65, height = 3.5, dpi=300)
plot_this + theme(legend.position = "none")
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_OnlyGLY_20C_withouthNOstorage_LAAARGE.png", width = 5.5, height = 3.6, dpi=300)


### Again, but focus on the samples which almost clustered
data.prop.labels2 <- subset_samples( data.prop.labels , Reactor_ID %in% c("A 04/10/22", "A 08/23", "A 02/05/23", "A 31/05/23", "A 01/24", "A 02/24") )
sample_data(data.prop.labels2)$Time_before_extr2 <- gsub("6","6\n", sample_data(data.prop.labels2)$Time_before_extr2)
sample_data(data.prop.labels2)$Collection_data <- gsub("10/22","  10/22", sample_data(data.prop.labels2)$Collection_data)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels2, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
base_plot <- plot_ordination(data.sqrt_prop, ordBC, color="Storage" )+
  guides(shape="none") +
  theme_classic(base_size = 10) +
  theme(legend.text =element_text(size=8.5),
        plot.margin = margin(1,1,1,1),
        legend.margin = margin(1,1,1,-15)
  ) +
  geom_point(size=3.52, alpha=0.5,  show.legend = F) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Collection_data), color="gray25", 
            size=3, vjust=-0.25, hjust=0.52,
            show.legend = FALSE, lineheight=0.85
  ) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="") +
  scale_color_manual(values= values_color) +
  geom_path(aes(group=Reactor_ID), color="darkgray" , linewidth = 0.25, alpha=0.85,
            arrow=arrow(length =unit(0.28,"cm"), type = "closed")
  )
# to match it with the PCoA with every sample...
base_plot + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), 
            color="black", size=3.2, show.legend = FALSE)
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_OnlyGLY_20C_ZOOM_on_CLUSTERED.png", width = 3.65, height = 3.5, dpi=300)
# to match it with the PCoA with 07/22...
base_plot +
  theme_classic(base_size = 10) +
  geom_text(aes(label=Time_before_extr2),
            color="black", size=4.35,
            show.legend = FALSE, lineheight=0.75 ) +
theme(legend.text =element_text(size=9),
      legend.margin = margin(1,1,1,-10)
      )
ggsave(file="Results/Beta_diversity_PCoA/CUOID_PCoA_Hellinger_Genus_OnlyGLY_20C_ZOOM_on_CLUSTERED2.png", width = 4.2, height = 3.65, dpi=300)




###################### PCoA BETA DIV (AGS Pilot Scale, vs -80C no Gly) #######################

# on genera
data.prop.labels<-subset_samples(data.genus.prop, Reactor=="Pilot_AGS" )
# updating labels
levels(sample_data(data.prop.labels)$Time_before_extr2) <- gsub("m","",levels(sample_data(data.prop.labels)$Time_before_extr2))
levels(sample_data(data.prop.labels)$Storage) <- gsub("-80C_noGly","-80 °C (no Gly)",levels(sample_data(data.prop.labels)$Storage))
sample_data(data.prop.labels)$Collection_data<-gsub("/202","/2", sample_data(data.prop.labels)$Collection_data, fixed=T)
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Storage !="Direct" ] <-""
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Collection_data %in% c("02/11/23") ] <-"    02/11/23"
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Collection_data %in% c("07/11/23") ] <-"    07/11/23"
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Collection_data %in% c("22/03/24") ] <-"22/03/24   "
sample_data(data.prop.labels)$Collection_data[ sample_data(data.prop.labels)$Collection_data %in% c("27/03/24") ] <-"27/03/24  "

# computation
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "-80 °C (no Gly)" = "darkblue")
plot_this <- plot_ordination(data.sqrt_prop, ordBC, color="Storage" )+
  guides(shape="none") +
  geom_point(size=3.35, alpha=0.5,  show.legend = F) +
  theme_classic(base_size = 10.25) +
  theme(legend.text =element_text(size=8.5),
        plot.margin = margin(1,1,1,1),
        legend.margin = margin(1,1,1,-15)
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Time_before_extr2), color="gray85", size=3.2, show.legend = FALSE) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="") +
  scale_color_manual(values= values_color) +
  geom_path(aes(group=Reactor_ID), color="darkgray" , linewidth = 0.19, alpha=0.86,
            # arrow=arrow(length =unit(0.15,"cm"), type = "closed"),
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Collection_data), color="gray25", 
            size=2.25, vjust=0.45, hjust=0.45,
            show.legend = FALSE, lineheight=0.85
  )
plot_this
ggsave(file="Results/Beta_diversity_PCoA/AGS_PCoA_Hellinger_Genus_80C_withouthNOstorage.png", width = 3.65, height = 3.5, dpi=300)
plot_this + theme(legend.position = "none")
ggsave(file="Results/Beta_diversity_PCoA/AGS_PCoA_Hellinger_Genus_80C_withouthNOstorage_LAAARGE.png", width = 5, height = 3.5, dpi=300)




############### TOP ABUNDANCES OF AGS PILOT SCALE (DIRECT vs -80) ####################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_samples(data.genus.prop, Reactor=="Pilot_AGS" )
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
  tabella$Genus<-gsub ("vibrionales","vibrion.", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
  levels(tabella$Storage)<-gsub ("_"," ", levels(tabella$Storage))
}
levels(tabella$Storage) <- gsub("C noGly"," °C (no Gly)", levels(tabella$Storage) )
levels(tabella$Collection_data) <- gsub("/202","/2", levels(tabella$Collection_data) )

fill_color_30_customised<-c("firebrick3","darkblue","lightblue1","black",
                            "darkcyan","darkmagenta", 
                            "brown4","grey50", "aquamarine3",
                            "bisque","cyan","darkgreen","yellow3",
                            "springgreen4", "grey80", 
                            "deepskyblue4","lightgreen","orange2",
                            "green", "violet", "blue",
                            "yellow", "pink3","yellow4", "chocolate4","darkorchid","red",
                            "lightgrey", "orange2","darkslategray3")
ggplot(data=tabella, aes(x= Storage, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_30_customised) +
  facet_grid( . ~ Collection_data ,
              scales = "free_x", space="free"
              #strip = strip_nested(size = "variable")
  ) +
  theme_classic(base_size =11) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1.1,
                                 hjust=1,
                                 size= 7.2
  ),
  axis.text.y=element_text(size=7),
  axis.ticks.x = element_blank(),
  strip.text.x = element_text(size=8.75, angle=0),
  legend.key.height = unit(0.14, "cm"),
  legend.key.width = unit(0.2, "cm"),
  legend.text = element_text ( size = 8 , face="italic"  ),
  legend.position="bottom",
  legend.margin = margin(-14 ,34, 1 ,0),
  plot.margin = margin(2,1,2,1),
  panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results/Abundances/AGSpilot_TOP_genera.png",width=5.2,height=4.25, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




################# VEEN DIAGRAMS Direct vs -20no GLY ( CUOIODEPUR ) ##########################

dir.create("Results/Searching_unique/")

data.venn<-subset_taxa(data.genus.prop, !is.na(Genus))
data.venn<-subset_taxa(data.venn, Genus != "uncultured" )
data.venn<-subset_samples(data.venn,  Reactor%in%c("CuoioDepur_Aero","PN_AGS2") )
data.venn<-subset_samples(data.venn, Storage %in% c( "Direct", "-20C"))
# data.venn<-subset_samples(data.venn, ! Reactor_ID %in% c( "A 11/23 ")) # it has only -20Gly and -20
data.venn<- filter_taxa(data.venn, function(x) min(x) < 0.01, TRUE) 
# otu_table(data.venn)

Direct<-subset_samples(data.venn, Storage=="Direct")
Direct<- prune_taxa(taxa_sums(Direct) > 0, Direct) 
Direct<-as.vector(tax_table(Direct)[ ,"Genus"]) # to translate in ASV code (of genus collapsed object)

data20C<-subset_samples(data.venn, Storage=="-20C")
data20C<- prune_taxa(taxa_sums(data20C) > 0, data20C)
data20C<-as.vector(tax_table(data20C)[ ,"Genus"])

ONLY_IN_Direct<- Direct[! Direct %in% data20C ]
vectorDirect <- ONLY_IN_Direct
ONLY_IN_Direct<- paste(ONLY_IN_Direct, collapse = ", ")
head(ONLY_IN_Direct)

ONLY_IN_data20C<- data20C[! data20C %in% Direct ]
vector20 <- ONLY_IN_data20C
ONLY_IN_data20C<- paste(ONLY_IN_data20C, collapse = ", ")
head(ONLY_IN_data20C)  # Methanosarcina


x<-list(Direct=Direct, data20C=data20C)
ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse2","coral","deepskyblue") ) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) )
ggsave(filename = paste0("Results/Searching_unique/Unique_bacteria_Storages_20CwithGly_OVER001.png"), width = 4, height = 4, dpi=300, bg = "white")
dev.off()


### saving the abundances ...
selected_Direct_20 <- c(vectorDirect , vector20 )
data_Genus_venn <- subset_taxa(data.venn, Genus %in% selected_Direct_20 )
sample_names(data_Genus_venn) <- sample_data(data_Genus_venn)$Sample_code
table_venn <- cbind.data.frame( round(otu_table(data_Genus_venn),4), tax_table(data_Genus_venn)[ ,c("Phylum","Genus")] )
table_venn$Which_storage_unique <- ifelse( table_venn$Genus %in% vectorDirect, "Direct" , "-20C")
table_venn <- table_venn[ order(rowMeans(otu_table(data_Genus_venn))) ,]
write.csv2( table_venn , 
            file = "Results/Searching_unique/Unique_bacteria_Storages_20CwithGly_OVER001_feature_table.csv" , row.names = F, quote = F)


### plotting the abundances
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.venn, Genus %in% selected_Direct_20 )
  tabella<-psmelt(top_data)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
  levels(tabella$Storage)<-gsub ("_"," ", levels(tabella$Storage))
}
levels(tabella$Storage) <- gsub("C"," °C (Gly)", levels(tabella$Storage) )
levels(tabella$Collection_data) <- gsub("/202","/2", levels(tabella$Collection_data) )

fill_color_custom<-c("brown4","orange","violet","gray15",
                     "darkmagenta", "blue","gray75","gold","darkgreen",
                     "deepskyblue2","red","chartreuse2")

ggplot(data=tabella, aes(x= Storage, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_custom) +
  facet_grid2( . ~ Reactor + Collection_data ,
               scales = "free_x", space="free",
               strip = strip_nested(size = "variable")
  ) +
  theme_classic(base_size =11) +
  # scale_y_continuous(expand=c(0,0.5) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1.1,
                                 hjust=1,
                                 size= 7.2
  ),
  axis.text.y=element_text(size=7),
  axis.ticks.x = element_blank(),
  strip.text.x = element_text(size=7, angle=0),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.2, "cm"),
  legend.text = element_text ( size = 9 , face="italic"  ),
  legend.position="bottom",
  legend.margin = margin(-14 ,25, 1 ,0),
  plot.margin = margin(2,1,2,1),
  panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results/Searching_unique/Unique_bacteria_Storages_AGSs_20C_withGly_OVER001_BARPLOT.png",width=6,height=4.25, dpi=300)




################# VEEN DIAGRAMS Direct vs -80 with no GLY (AGSs) ##########################

suppressWarnings(dir.create("Results/Searching_unique/"))

data.venn<-subset_taxa(data.genus.prop, !is.na(Genus))
data.venn<-subset_taxa(data.venn, Genus != "uncultured" )
data.venn<-subset_samples(data.venn,  Reactor%in%c("Pilot_AGS","PN_AGS2") )
data.venn<-subset_samples(data.venn, Storage %in% c( "Direct", "-80C_noGly"))
data.venn<- filter_taxa(data.venn, function(x) min(x) < 0.01, TRUE) 
# otu_table(data.venn)

Direct<-subset_samples(data.venn, Storage=="Direct")
Direct<- prune_taxa(taxa_sums(Direct) > 0, Direct) 
Direct<-as.vector(tax_table(Direct)[ ,"Genus"]) # to translate in ASV code (of genus collapsed object)

data80C<-subset_samples(data.venn, Storage=="-80C_noGly")
data80C<- prune_taxa(taxa_sums(data80C) > 0, data80C)
data80C<-as.vector(tax_table(data80C)[ ,"Genus"])

# data20C<-subset_samples(data.venn, Storage=="-20C")
# data20C<- prune_taxa(taxa_sums(data20C) > 0, data20C)
# data20C<-as.vector(tax_table(data20C)[ ,"Genus"])

ONLY_IN_data80C<- data80C[! data80C %in% Direct]
vector80C <- ONLY_IN_data80C
ONLY_IN_data80C<- paste(ONLY_IN_data80C, collapse = ", ")
head(ONLY_IN_data80C)

ONLY_IN_Direct<- Direct[! Direct %in% data80C]
vectorDirect <- ONLY_IN_Direct
ONLY_IN_Direct<- paste(ONLY_IN_Direct, collapse = ", ")
head(ONLY_IN_Direct)

# ONLY_IN_data20C<- data20C[! data20C %in% Direct & ! data20C %in% data80C]
# vector20 <- ONLY_IN_data20C
# ONLY_IN_data20C<- paste(ONLY_IN_data20C, collapse = ", ")
# head(ONLY_IN_data20C)  # Methanosarcina


x<-list(Direct=Direct, data80C=data80C)
ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse2","coral","deepskyblue") ) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) )
ggsave(filename = paste0("Results/Searching_unique/Unique_bacteria_Storages_AGSs_80C_noGly.png"), width = 4, height = 4, dpi=300, bg = "white")
dev.off()


### saving the abundances ...
selected_Direct_80 <- c(vectorDirect , vector80C )
data_Genus_venn <- subset_taxa(data.venn, Genus %in% selected_Direct_80 )
sample_names(data_Genus_venn) <- sample_data(data_Genus_venn)$Sample_code
table_venn <- cbind.data.frame( round(otu_table(data_Genus_venn),4), tax_table(data_Genus_venn)[ ,c("Phylum","Genus")] )
table_venn$Which_storage_unique <- ifelse( table_venn$Genus %in% vectorDirect, "Direct" , "-80C")
table_venn <- table_venn[ order(rowMeans(otu_table(data_Genus_venn))) ,]
write.csv2( table_venn , 
            file = "Results/Searching_unique/Unique_bacteria_AGSs_Direct_vs_80C_NOgly_feature_table.csv" , row.names = F, quote = F)


### plotting the abundances
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.venn, Genus %in% selected_Direct_80 )
  tabella<-psmelt(top_data)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
  levels(tabella$Storage)<-gsub ("_"," ", levels(tabella$Storage))
}
levels(tabella$Storage) <- gsub("C noGly"," °C (no Gly)", levels(tabella$Storage) )
levels(tabella$Collection_data) <- gsub("/202","/2", levels(tabella$Collection_data) )

fill_color_custom<-c("brown4","springgreen2","wheat","lightcoral","orange",
                     "yellow3","darkmagenta", "blue","gray15","gold","darkgreen",
                     "violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")

ggplot(data=tabella, aes(x= Storage, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_custom) +
  facet_grid( . ~ Reactor + Collection_data ,
              scales = "free_x", space="free"
              #strip = strip_nested(size = "variable")
  ) +
  theme_classic(base_size =11) +
  # scale_y_continuous(expand=c(0,0.5) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1.1,
                                 hjust=1,
                                 size= 7.2
  ),
  axis.text.y=element_text(size=7),
  axis.ticks.x = element_blank(),
  strip.text.x = element_text(size=8.75, angle=0),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.2, "cm"),
  legend.text = element_text ( size = 9 , face="italic"  ),
  legend.position="bottom",
  legend.margin = margin(-14 ,25, 1 ,0),
  plot.margin = margin(2,1,2,1),
  panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results/Searching_unique/Unique_bacteria_Storages_AGSs_80C_noGly_OVER001_BARPLOT.png",width=4.8,height=4.25, dpi=300)




############# DA WITH DESEQ2 (Direct vs OLDEST STORAGE, only CUOIODEPUR at -20gly) ##################

### building the list of samples to test
direct_samples <- as.character(metadata$Sample_name[metadata$Reactor=="CuoioDepur_Aero" & metadata$Storage=="Direct"])
direct_ID <- as.character(metadata$Collection_data[metadata$Sample_name%in%direct_samples])
# picking the oldest sample (stored for a longer period) for each collection data (ID)
already_picked_ID <- NULL
selected_samples <- direct_samples
for( time in c("3y","1y","6m","3m")){
  selected_table <- metadata[metadata$Collection_data%in%direct_ID & metadata$Time_before_extr2==time & metadata$Storage=="-20C", ]
  selected_samples <- c(selected_samples, as.character(selected_table$Sample_name[!selected_table$Sample_name%in%selected_samples & ! selected_table$Collection_data%in%already_picked_ID]) )
  already_picked_ID <- c(already_picked_ID, as.character(selected_table$Collection_data))
}
suppressWarnings(rm(data_pruned, data.genus_pruned, subset_target))
subset_target <- subset_samples(data.genus, Sample_name %in% selected_samples ) # Sample_name%in%c("PGB10","PGB1"))
# View(sample_data(subset_target))   # to check --> OK

# the lines below format the according to the other analysis and thus and thus using the same code (the time will be displayed correctly in the plot due to the other column)
levels(sample_data(subset_target)$Time_before_extr) <- gsub("[1-9]month","6month", levels(sample_data(subset_target)$Time_before_extr))
levels(sample_data(subset_target)$Time_before_extr) <- gsub("1year","6month", levels(sample_data(subset_target)$Time_before_extr))
levels(sample_data(subset_target)$Time_before_extr) <- gsub("3year","6month", levels(sample_data(subset_target)$Time_before_extr))


##### STARTING THE DIFFERENTIAL ANALYSIS
data_pruned<- prune_taxa(taxa_sums(subset_target) > 10, subset_target) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

sample_data(data_pruned)$Reactor_ID  <- gsub("/", "_", sample_data(data_pruned)$Reactor_ID , fixed=T)
sample_data(data_pruned)$Reactor_ID  <- gsub(" ", "_", sample_data(data_pruned)$Reactor_ID )
sample_data(data_pruned)$Reactor_ID  <- gsub("-", "_", sample_data(data_pruned)$Reactor_ID , fixed=T)

Table_tot<-NULL
Res_tot<-NULL
for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  if(t=="Genus"){ # updating missing names (NA and uncultured) but only for genus level
    taxa_temp<-as.data.frame(tax_table(d))
    for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
      taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
    for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
    for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
      taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
    for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
    for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
    tax_table(d)<-as.matrix(taxa_temp)
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Reactor_ID + Time_before_extr)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time_before_extr", "None", "6month"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>0 ,]
  #res<-res[(res$padj < 0.05) ,]
  res<-res[res$baseMean > 50, ] # arbitrary threshold to avoid the most noisy result
  res
  if(length(res$log2FoldChange)>0){ # if there are results...
    cat(paste(length(res$log2FoldChange),"results for the",t,"level\n"))
    Sys.sleep(1)
    r<-as.data.frame(res)
    r$ASV<-row.names(r)
    Taxa.d<-as.data.frame(tax_table(d))
    Taxa.d$ASV<-row.names(Taxa.d)
    r<-dplyr::left_join(r, Taxa.d, by="ASV")
    r$Kingdom<-NULL
    r$Species<-NULL
    assign(paste(t,"results",sep="_"), r)
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Direct_extract_vs_Mixed_liquor.csv"), row.names = F, quote=F, na = "")
    r_level<-r
    r_level[, "Taxon"]<- rep(t)
    Res_tot<-rbind.data.frame(Res_tot,r_level)
    ### single box plots
    target<-r[[t]]
    colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
    target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
    Table_DE<-psmelt(target)
    colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
    Table_DE$ASV<-NULL
    # Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
    assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
    ### appending to unique box plot
    index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t)
    index<- index[-length(index)] # removing the last index, regarding the taxa of interest
    Table_DE[,index]<-NULL
    Table_DE$Taxa<-t
    colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
    Table_tot<-rbind.data.frame(Table_tot, Table_DE)
  } else {
    cat("Any results for the",t,"level\n")
    Sys.sleep(1)
  }
}

# View(Res_tot)
write.csv2(Res_tot, file="Results/Differently_abundant_bact/Storage_20gly_Differently_abundant.csv", row.names = F)

Res_tot <- Res_tot[!is.na(Res_tot$Order), ]
Table_tot <- Table_tot[!grepl("NA_ o NA", Table_tot$Bacteria), ]
Table_tot <- Table_tot[!grepl("uncultured_ o uncultured", Table_tot$Bacteria), ]



############ PLOTTING THESE RESULTS ...

Table_tot$Time_before_extr<-factor(Table_tot$Time_before_extr, levels = unique(Table_tot$Time_before_extr))
Table_tot<-Table_tot[order(match(Table_tot$Time_before_extr, levels(Table_tot$Time_before_extr))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Time_before_extr)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Time_before_extr)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

## to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Time_before_extr)


#unique(tabella_g$Bacteria)
tabella_g$Time_before_extr <- gsub("_"," ",tabella_g$Time_before_extr)
tabella_g$Bacteria<-gsub("_marine","\nmarine",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("Candidatus_","Ca. ",tabella_g$Bacteria)


#### removing redundant results
{ 
  # 1) if redundant then same abundances
  auto_Max<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, max )
  auto_Min<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, min )
  auto_Median<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, median )
  # reordering the names (lower ranks first --> the redundants are those at higher levels)
  ordered_names <- c( unique(Table_tot[Table_tot$Taxa=="Genus", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Family", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Order", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Class", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Phylum", "Bacteria" ])
  )
  ordered_names[!is.na(ordered_names)]
  auto_Max<-auto_Max[ ordered_names  ]
  auto_Min<-auto_Min[ ordered_names  ]
  auto_Median<-auto_Median[ ordered_names ] # the table it self follows the correct taxonomic order because it is built in this way
  # same abundances
  auto_redund<-names(auto_Max)[ duplicated(auto_Min) & duplicated(auto_Max) & duplicated(auto_Median)  ]
  
  # 2) if redundant then same ASV
  auto_ASV<-Res_tot$ASV[duplicated(Res_tot$ASV)]
  auto_ASV_Names<-unique(Table_tot$Bacteria[Table_tot$OTU %in% auto_ASV])
  auto_redund<-auto_redund[auto_redund %in% auto_ASV_Names]
}

Redund<-auto_redund

tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)

# tabella3 <-rbind.data.frame(tabella_g, tabella_2)   # few results, it's better to re-join them ...
tabella3 <- tabella_g  

tabella3$Bacteria<-gsub("_ ", "\n", tabella3$Bacteria)
tabella3 <- tabella3[ ! is.na(tabella3$Bacteria) , ]
tabella3$Time_before_extr <- gsub("None", "Direct" , tabella3$Time_before_extr )
tabella3$Time_before_extr <- gsub("6month", "-20C Gly" , tabella3$Time_before_extr )
tabella3$Time_before_extr <- factor( tabella3$Time_before_extr , levels = c("Direct","-20C Gly") )
levels(tabella3$Time_before_extr2) <- gsub("m","", levels(tabella3$Time_before_extr2) )
levels(tabella3$Time_before_extr2) <- gsub("1y"," 1y", levels(tabella3$Time_before_extr2) )
levels(tabella3$Time_before_extr2) <- gsub("3y"," 3y", levels(tabella3$Time_before_extr2) )
tabella3$Seq_batch <- paste("Seq batch", tabella3$Seq_batch)
these_colors <- c( "Direct"="coral3", "-20C Gly"="deepskyblue" )
plot_DESEQ <- ggplot(tabella3, aes(x= Xaxis, y=Abundance, fill=Time_before_extr)) +
  theme_classic2(base_size = 12) +
  scale_fill_manual(values = these_colors ) +
  scale_color_manual(values = these_colors ) +
  #facet_wrap2(nrow=2,factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus"))~Bacteria, 
  facet_wrap2(nrow=2, . ~Bacteria, 
              #labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", 
              strip=strip_nested(size = "variable", bleed = T),
              drop = TRUE) +
  geom_line(aes(group=Reactor_ID), linewidth=0.15, color="gray30", alpha= 0.5) +
  scale_x_discrete(labels=unique(levels(tabella_g$Time_before_extr)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=6.3,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=10.5),
        axis.text.y = element_text(size=7),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        panel.grid.minor.y= element_blank(),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt")
  ) +
  # guides(fill="none") +
  labs(y="Percentual abundance", fill="", color="", x="", shape="") +
  theme(legend.margin=margin(-15, 0, 0, 0),
        legend.position = "bottom")

# plot_DESEQ + geom_point(aes(color=Time_before_extr, shape=Seq_batch), size=0.85, alpha=1) +
#   geom_point(aes(color=Time_before_extr, shape=Seq_batch), size=2.5, alpha=0.5) +
#   geom_text(aes(label= Sample_name  ), size=2.1, color="gray25", show.legend = F)
# ggsave(filename = "Results/Differently_abundant_bact/StorageCUOIOD_Differently_abundant_sampleName_BATCH.png", width = 6, height = 3.65, dpi=300)

plot_DESEQ + geom_point(aes(color=Time_before_extr, shape=Seq_batch), size=0.85, alpha=1) +
  geom_point(aes(color=Time_before_extr, shape=Seq_batch), size=2.5, alpha=0.5)
ggsave(filename = "Data_check/DA_StorageCUOIOD_Differently_abundant_BATCH.png", width = 6, height = 3.65, dpi=300)

# plot_DESEQ + geom_point(aes(color=Time_before_extr), size=0.85, alpha=1) +
#   geom_point(aes(color=Time_before_extr), size=2.5, alpha=0.5)
# ggsave(filename = "Results/Differently_abundant_bact/StorageCUOIOD_Differently_abundant.png", width = 6, height = 3.65, dpi=300)

plot_DESEQ + geom_point(aes(color=Time_before_extr), size=0.85, alpha=1) +
  geom_point(aes(color=Time_before_extr), size=2.5, alpha=0.3) +
  geom_text(aes(label= Time_before_extr2  ), size=1.8, hjust=-0.4 , color="gray25", show.legend = F)
ggsave(filename = "Results/Differently_abundant_bact/StorageCUOIOD_Differently_abundant_Time.png", width = 5.6, height = 3.65, dpi=300)




############# DA WITH DESEQ2 (Direct vs stored, ONLY AGS and PN-AGS2 at -80 no gly) ##################

suppressWarnings(rm(data_pruned, data.genus_pruned, subset_target))

subset_target <- subset_samples(data, Reactor=="Pilot_AGS" | Sample_name%in%c("PGB11","PGB1") )
# the line below formats the table according to the other analysis and thus using the same code (the time will be displayed correctly in the plot due to the other column)
levels(sample_data(subset_target)$ Time_before_extr) <- gsub("[1-9]month","6month", levels(sample_data(subset_target)$Time_before_extr))
data_pruned<- prune_taxa(taxa_sums(subset_target) > 10, subset_target)
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

sample_data(data_pruned)$Reactor_ID  <- gsub("/", "_", sample_data(data_pruned)$Reactor_ID )
sample_data(data_pruned)$Reactor_ID  <- gsub(" ", "_", sample_data(data_pruned)$Reactor_ID )

Table_tot<-NULL
Res_tot<-NULL
for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  if(t=="Genus"){ # updating missing names (NA and uncultured) but only for genus level
    taxa_temp<-as.data.frame(tax_table(d))
    for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
      taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
    for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
    for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
      taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
    for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
    for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
      taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
    tax_table(d)<-as.matrix(taxa_temp)
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Reactor_ID + Time_before_extr)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time_before_extr", "None", "6month"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.5) & abs(res$log2FoldChange)>0,]
  #res<-res[(res$padj < 0.05), ]
  res<-res[res$baseMean > 50, ]
  res
  if(length(res$log2FoldChange)>1){ # if there are results...
    cat(paste(length(res$log2FoldChange),"results for the",t,"level\n"))
    Sys.sleep(1)
    r<-as.data.frame(res)
    r$ASV<-row.names(r)
    Taxa.d<-as.data.frame(tax_table(d))
    Taxa.d$ASV<-row.names(Taxa.d)
    r<-dplyr::left_join(r, Taxa.d, by="ASV")
    r$Kingdom<-NULL
    r$Species<-NULL
    assign(paste(t,"results",sep="_"), r)
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Direct_extract_vs_Mixed_liquor.csv"), row.names = F, quote=F, na = "")
    r_level<-r
    r_level[, "Taxon"]<- rep(t)
    Res_tot<-rbind.data.frame(Res_tot,r_level)
    ### single box plots
    target<-r[[t]]
    colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
    target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
    Table_DE<-psmelt(target)
    colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
    Table_DE$ASV<-NULL
    # Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
    assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
    ### appending to unique box plot
    index<- which(colnames(Table_DE)=="Kingdom") : which(colnames(Table_DE)==t)
    index<- index[-length(index)] # removing the last index, regarding the taxa of interest
    Table_DE[,index]<-NULL
    Table_DE$Taxa<-t
    colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
    Table_tot<-rbind.data.frame(Table_tot, Table_DE)
  } else {
    cat("Any results for the",t,"level\n")
    Sys.sleep(1)
  }
}

# View(Res_tot)

write.csv2(Res_tot, file="Results/Differently_abundant_bact/StorageAGS_Differently_abundant.csv", row.names = F)



############ PLOTTING THESE RESULTS ...

Table_tot$Time_before_extr<-factor(Table_tot$Time_before_extr, levels = c("None","6month"))
Table_tot<-Table_tot[order(match(Table_tot$Time_before_extr, levels(Table_tot$Time_before_extr))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Time_before_extr)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Time_before_extr)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

## to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Time_before_extr)


#unique(tabella_g$Bacteria)
tabella_g$Time_before_extr <- gsub("_"," ",tabella_g$Time_before_extr)
tabella_g$Bacteria<-gsub("aceae$","ac.",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_marine","\nmarine",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("Candidatus_","Candidatus\n",tabella_g$Bacteria)

##### removing redundant results
{ 
  # 1) if redundant then same abundances
  auto_Max<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, max )
  auto_Min<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, min )
  auto_Median<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, median )
  # reordering the names (lower ranks first --> the redundants are those at higher levels)
  ordered_names <- c( unique(Table_tot[Table_tot$Taxa=="Genus", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Family", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Order", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Class", "Bacteria" ]),
                      unique(Table_tot[Table_tot$Taxa=="Phylum", "Bacteria" ])
  )
  ordered_names[!is.na(ordered_names)]
  auto_Max<-auto_Max[ ordered_names  ]
  auto_Min<-auto_Min[ ordered_names  ]
  auto_Median<-auto_Median[ ordered_names ] # the table it self follows the correct taxonomic order because it is built in this way
  # same abundances
  auto_redund<-names(auto_Max)[ duplicated(auto_Min) & duplicated(auto_Max) & duplicated(auto_Median)  ]
  
  # 2) if redundant then same ASV
  auto_ASV<-Res_tot$ASV[duplicated(Res_tot$ASV)]
  auto_ASV_Names<-unique(Table_tot$Bacteria[Table_tot$OTU %in% auto_ASV])
  auto_redund<-auto_redund[auto_redund %in% auto_ASV_Names]
}

Redund<-auto_redund
Redund<-c(Redund,"Bacteroidota","Chloroflexi") # to select redundants (same ASV, same results at different levels)
tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)

# tabella3 <-rbind.data.frame(tabella_g, tabella_2)   # few results, it's better to re-join them ...
tabella3 <- tabella_g

tabella3$Bacteria<-gsub("_ ", "\n", tabella3$Bacteria)
tabella3 <- tabella3[ ! is.na(tabella3$Bacteria) , ]
tabella3$Time_before_extr <- gsub("None", "Direct" , tabella3$Time_before_extr )
tabella3$Time_before_extr <- gsub("6month", "-80C" , tabella3$Time_before_extr )
tabella3$Time_before_extr <- factor( tabella3$Time_before_extr , levels = c("Direct","-80C") )
levels(tabella3$Time_before_extr2) <- gsub("m","", levels(tabella3$Time_before_extr2))
these_colors <- c( "Direct"="coral3", "-80C"="blue" )
plot_DESEQ <- ggplot(tabella3, aes(x= Xaxis, y=Abundance, fill=Time_before_extr, shape=Reactor)) +
  theme_classic2(base_size = 12) +
  scale_fill_manual(values = these_colors ) +
  scale_color_manual(values = these_colors ) +
  #facet_wrap2(nrow=2,factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus"))~Bacteria, 
  facet_wrap2(nrow=2, . ~Bacteria, 
              #labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", 
              # strip=strip_nested(size = "variable", bleed = T),
              drop = TRUE) +
  geom_line(aes(group=Reactor_ID), linewidth=0.15, color="gray30", alpha= 0.5) +
  scale_x_discrete(labels=unique(levels(tabella_g$Time_before_extr)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=5.6,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=11.8),
        axis.text.y = element_text(size=7),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        panel.grid.minor.y= element_blank(),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt")
  ) +
  guides(shape="none") +
  # guides(fill="none") +
  labs(y="Percentual abundance", fill="", color="", x="", shape="") +
  theme(legend.margin=margin(-15, 0, 0, 0),
        legend.position = "bottom")

plot_DESEQ + geom_point(aes(color=Time_before_extr), size=0.75, alpha=1) +
  geom_point(aes(color=Time_before_extr), size=2.5, alpha=0.3) +
  geom_text(aes(label= Time_before_extr2  ), size=1.8, hjust=-0.9,
            color="gray25", show.legend = F)
ggsave(filename = "Results/Differently_abundant_bact/StorageAGS_Differently_abundant_TRIANGLE_IS_PN.png", width = 5.75, height = 3.65, dpi=300)


### plotting again but without the PN-AGS2
tabella4 <- tabella3[ tabella3$Reactor!="PN_AGS2" , ]
plot_DESEQ <- ggplot(tabella4, aes(x= Xaxis, y=Abundance, fill=Time_before_extr)) +
  theme_classic2(base_size = 9.65) +
  scale_fill_manual(values = these_colors ) +
  scale_color_manual(values = these_colors ) +
  #facet_wrap2(nrow=2,factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus"))~Bacteria, 
  facet_wrap2(nrow=2, . ~Bacteria, 
              #labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", 
              # strip=strip_nested(size = "variable", bleed = T),
              drop = TRUE) +
  geom_line(aes(group=Reactor_ID), linewidth=0.15, color="gray30", alpha= 0.5) +
  scale_x_discrete(labels=unique(levels(tabella_g$Time_before_extr)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=5.8,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=11.8),
        axis.text.y = element_text(size=7),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        panel.grid.minor.y= element_blank(),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt")
  ) +
  guides(shape="none") +
  # guides(fill="none") +
  labs(y="Percentual abundance", fill="", color="", x="", shape="") +
  theme(legend.margin=margin(-12.5, 0, 0, 0),
        legend.position = "bottom")

plot_DESEQ + geom_point(aes(color=Time_before_extr), size=0.75, alpha=1) +
  geom_point(aes(color=Time_before_extr), size=2.5, alpha=0.3) +
  geom_text(aes(label= Time_before_extr2  ), size=1.8, hjust=-0.9,
            color="gray25", show.legend = F)
ggsave(filename = "Results/Differently_abundant_bact/StorageAGS_Differently_abundant_REMOVING_PNAGS2.png", width = 5.75, height = 3.65, dpi=300)




####################### OBSERVING MAIN THRENDS THROUGH TIME ####################### 

data_selected_all <- subset_taxa(data.genus.prop, Genus %in% c("Nitrosomonas","Candidatus_Kaiserbacteria","Leptospira",
                                                               "Arcobacter","Sphaerotilus","JG30-KF-CM45",
                                                               "NS11-12_marine_group","NS9_marine_group",
                                                               "env.OPS_17",  # Sphingobacteriales
                                                               "Candidatus_Competibacter", "Flavobacterium",  # add because they are important!
                                                               "Lysobacter","Prosthecobacter", # gly results
                                                               "Defluviimonas", "Rhodocyclaceae","Thauera"))
colors_taxa <- c("Nitrosomonas"="blue3",
                 "Ca. Kaiserbacteria"="darkgreen",   # Candidatus will be modified in Ca. below
                 "Ca. Competibacter"="pink3",
                 "Flavobacterium"="lightblue",
                 "Leptospira"="red2",
                 "Arcobacter"="brown4",
                 "Sphaerotilus"="yellow2",
                 "JG30-KF-CM45"="gold4",
                 "NS11-12_marine"="coral2",
                 "NS9_marine"="orange",
                 "env.OPS_17"="green",
                 "Lysobacter"="deepskyblue4",
                 "Prosthecobacter"="darkmagenta",
                 "Defluviimonas"="cyan2",
                 "Rhodocyclaceae"="gray25",
                 "Thauera"="violet"
)


#### -20 with GLY
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-20C"))

for(x in c("A 07/22","PN_AGS_4","A 11/23", "A 31/05/23", "A 04/10/22") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  table_to_plot$Genus <- gsub("marine_group","marine", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=8),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-20C gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/20C_gly_",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}



#### -20 no GLY
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-20C_noGly"))

for(x in c("PN_AGS_4","A 11/23") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  table_to_plot$Genus <- gsub("marine_group","marine", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=8),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-20C NO gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/20C__NOgly__",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}



#### -80 with GLY
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-80C"))

for(x in c("PN_AGS_4","A 11/23") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  table_to_plot$Genus <- gsub("marine_group","marine", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=8),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-80C gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/80C_gly_",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}



#### -80 no GLY
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-80C_noGly"))

for(x in c("PN_AGS_4","A 11/23") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  table_to_plot$Genus <- gsub("marine_group","marine", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=8),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-80C NO gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/80C__NOgly__",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}




################# EXTRA: OBSERVING KEY TAXA THRENDS ##############

dir.create("Results/Abundances/Abund_over_time/Extra_check_on_certain_taxa/")

data_selected_all <- subset_taxa(data.genus.prop, Genus %in% c("Candidatus_Accumulibacter","Propionivibrio",
                                                               "Tetrasphaera", "Candidatus_Brocadia",
                                                               #"Nitrospira","Nitrobacter",
                                                               "Candidatus_Competibacter",  # add because they are important!
                                                               "Thauera","Neomegalonema")
)
colors_taxa <- c(
  "Ca. Accumulibacter"="red2",   # Candidatus will be modified in Ca. below
  "Ca. Competibacter"="green2",
  "Propionivibrio"="darkgreen",
  "Tetrasphaera"="brown4",
  "Ca. Brocadia"="gray20",
  # "Nitrosomonas"="darkblue",
  "Neomegalonema"="orange",
  "Thauera"="violet"
)



#### -20 no GLY
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-20C_noGly"))

for(x in c("PN_AGS_4","A 11/23") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=9.5),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-20C NO gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/Extra_check_on_certain_taxa/EXTRA__20C__NOgly__focusOnKeyTaxa",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}




#### -80 no GLY
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-80C_noGly"))

for(x in c("PN_AGS_4","A 11/23") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=9.5),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-80C NO gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/Extra_check_on_certain_taxa/EXTRA__80C__NOgly__focusOnKeyTaxa",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}



#### -20 C Gly
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-20C"))

for(x in c("PN_AGS_4","A 11/23", "A 07/22") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=9.5),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-20C gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/Extra_check_on_certain_taxa/EXTRA_20C_gly__focusOnKeyTaxa",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}



#### -80 C Gly
data_selected_storage <- subset_samples(data_selected_all, Storage %in% c("Direct","-80C"))

for(x in c("PN_AGS_4","A 11/23") ){
  
  data_selected <- subset_samples(data_selected_storage, Reactor_ID == x)
  table_to_plot <- psmelt(data_selected)
  table_to_plot$Genus <- gsub("Candidatus_","Ca. ", table_to_plot$Genus)
  file_name<- gsub("/","", x)
  file_name<- gsub(" ","", file_name)
  file_name<- gsub("PN_AGS_4","PN_AGS_2", file_name)
  title_here<- gsub("A ","AS ", x)
  title_here<- gsub("PN_AGS_4","PN-AGS 2", title_here)
  
  ggplot( data=table_to_plot, aes(x=Time_before_extr2, y=Abundance, color=Genus)) +
    theme_classic(base_size = 9) +
    geom_point() +
    scale_color_manual(values=colors_taxa) +
    geom_line(aes(group=Genus) , show.legend = F, linewidth=0.15) +
    scale_x_discrete(expand = c(0,0.1)) +
    theme( axis.title.x = element_blank(),
           axis.text.x = element_text(size=6.85, angle=25),
           axis.ticks.x = element_blank(),
           panel.grid.major.y = element_line(linewidth = 0.25, color="gray"),
           panel.grid.minor.y = element_line(linewidth = 0.1, color="darkgray"),
           legend.text = element_text(size=9.5),
           legend.position = "bottom",
           legend.key.spacing.y =  unit(0.2,"cm"),
           legend.key.size = unit(0.25, "cm"),
           legend.margin = margin(-1,24,0,1),
           title= element_text(size=9)
    ) +
    guides(color=guide_legend(nrow=4,byrow=T)) +
    labs(color="",group="", y= "Percent abundance in the sample",
         title= paste("Abundances in",title_here,"over time (-80C gly)")
    )
  ggsave(paste0("Results/Abundances/Abund_over_time/Extra_check_on_certain_taxa/EXTRA_80C_gly__focusOnKeyTaxa",file_name,".png"),
         dpi=300, width=4.25, height = 3.5)
}

