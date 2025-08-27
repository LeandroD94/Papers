# 26/07/2025

### Script part 3: comparison with other AGSs (not SBR)



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


### Definying the main sample groups of this part
temp <- rep("SBR AGS", length(sample_names(data) ) )
temp[ sample_data(data)[,"Project"][[1]] %in% c("PRJNA606844","PRJNA1015110","PRJNA1082061")  ] <- rep("CFR")
# temp[ sample_data(data)[,"Project"]=="PRJNA1082061" ] <- rep("CFR-PN")
temp <- factor( temp , levels=c("SBR AGS","CFR") )
sample_data(data)[,"AGS"] <- temp
sample_data(data.genus)<- sample_data(data)


### proportions
data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100 )
data.genus.prop<-transform_sample_counts(data.genus, function(x) (x/sum(x))*100 )




###################### ABUNDANCES BAR PLOT (SBR vs CFR) ##########################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.genus.prop, Genus!="unclassified")
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top <- names(sort(taxa_sums(data_temp), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_temp)
  others<-taxa_names(data_temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the species
  tabella<-aggregate(.~Genus+Euk_clade+Reactor_scale+Influent+Project_number+Sample_name+AGS, tabella[ , c("Reactor_scale","Project_number","Euk_clade","Influent","Sample_name","Abundance","Genus","AGS")], FUN=sum)
  tabella$Genus<-gsub ("Candidatus_","Ca.", tabella$Genus)
  tabella$Genus<-gsub ("bacterium_","", tabella$Genus)
  tabella$Genus<-paste0(tabella$Euk_clade,": ",tabella$Genus)
  tabella$Genus<-gsub("No_euk","Prokaryota", tabella$Genus)
  tabella$Genus<-gsub(".*Others","Others", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
to_plot <- ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30_BIS) +
  #scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 8.5 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-14,30,0,0),
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  panel.spacing.x = unit(1.41,"pt"),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(2,1,2,1)
  ) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
to_plot + 
  facet_nested(~ AGS + Project_number, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(size=8.4))
ggsave(file="Results/Comparison_with_particular_AGS/TOP_microbial_genera_comparison.png",width=6.8,height=4.5, dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




###################### ABUNDANCES BAR PLOT (Only CFR) ##########################

### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
data_temp<- subset_taxa(data.genus.prop, Genus!="unclassified")
data_temp<- subset_samples(data_temp,AGS %in% "CFR" )
data_temp<- transform_sample_counts(data_temp, fun= function(x) x/sum(x)*100 )
{ top <- names(sort(taxa_sums(data_temp), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_temp)
  others<-taxa_names(data_temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  # the aggregation below solves a graphical glitch with ggplot2 ... and also re-orders in alphabetic order the species
  tabella<-aggregate(.~Genus+Euk_clade+Reactor_scale+Influent+Project_number+Sample_name+AGS, tabella[ , c("Reactor_scale","Project_number","Euk_clade","Influent","Sample_name","Abundance","Genus","AGS")], FUN=sum)
  tabella$Genus<-gsub ("Candidatus_","Ca.", tabella$Genus)
  tabella$Genus<-gsub ("bacterium_","", tabella$Genus)
  tabella$Genus<-paste0(tabella$Euk_clade,": ",tabella$Genus)
  tabella$Genus<-gsub("No_euk","Prokaryota", tabella$Genus)
  tabella$Genus<-gsub(".*Others","Others", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[! unique(tabella$Genus) %in% "Others"],"Others"))
}
to_plot <- ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30_BIS) +
  #scale_fill_manual(values=fill_color_25) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.title =element_text(size=10),
  axis.title.x = element_text(vjust=2.5),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 8.5 , face="italic"),
  legend.position="bottom",
  legend.margin = margin(-12.5,30,0,0),
  legend.spacing.y  = unit(0.1, "cm"),
  legend.spacing.x  = unit(3, "cm"),
  panel.spacing.x = unit(1.45,"pt"),
  scale_x_discrete(expand=c(-0.1,1)),
  plot.margin = margin(2,2,2,2)
  ) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
to_plot + 
  facet_nested(~ AGS + Reactor_scale, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(size=11.5))
ggsave(file="Results/Comparison_with_particular_AGS/TOP_microbial_genera_ONLY_CFR.png",width=5.95,height=4.5, dpi=300)
dev.off()




########################## PCoA (WITH CFR) ############################################

data.prop.labels<-subset_taxa(data.prop, Domain %in% c("Bacteria","Archea"))
levels(sample_data(data.prop.labels)$AGS)<-gsub("SBR AGS", "SBR", levels(sample_data(data.prop.labels)$AGS) )
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# with ellipses (Reactor scale)
plot_ordination(data.sqrt_prop, ordBC , color="AGS", shape="Reactor_scale") +
  scale_color_manual(values=c("deepskyblue3","orange")) +
  scale_shape_manual(values=c(16,18,15)) +
  geom_point(size=1.5) +
  geom_point(size=4.65, alpha= 0.4) +
  theme_classic(base_size = 9) +
  #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2, show.legend = FALSE) 
  #geom_text(aes(label=sample_data(data.prop.labels)$Reactor_number), color="grey20", size=2.8, show.legend = FALSE) +
  geom_text(aes(label=sample_data(data.prop.labels)$Project_number), color="grey25", size=2.9, show.legend = FALSE) +
  theme(legend.margin = margin(10,0,10,0),
        legend.text = element_text(size=7.8)) +
  # guides(color="none")+
  labs(title=paste0("PCoA on prokaryota using Hellinger distance\n (euclidean on Hellinger transformed genera"), 
       color="Type",
       shape="Scale",
       #caption="The labels on the points are the project numbers",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")
       )
aggsave(file=paste0("Results/Comparison_with_particular_AGS/PCoA_Hellinger_prokaryota_genera.png"), width = 5, height = 3.8, dpi=300)




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
  tax_table(data.genus.temp)[,"Genus"]<-gsub("^Competibacter","Candidatus_Competibacter", tax_table(data.genus.temp)[,"Genus"])
  tax_table(data.genus.temp)[,"Genus"]<-gsub("^Accumulibacter","Candidatus_Accumulibacter", tax_table(data.genus.temp)[,"Genus"])
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
tax_table(data.genus.temp)[  easy_subset %in% AOB_list ,"Group_Function"] <- "AOB"
tax_table(data.genus.temp)[  easy_subset %in% NOB_list ,"Group_Function"] <- "NOB"
tax_table(data.genus.temp)[  easy_subset %in% PAO_list ,"Group_Function"] <- "PAO"
tax_table(data.genus.temp)[  easy_subset %in% GAO_list ,"Group_Function"] <- "GAO"
tax_table(data.genus.temp)[  easy_subset %in% table_PHA[,1] ,"Group_Function"] <- "Others PHA"
data.genus.temp<- prune_taxa(taxa_sums(data.genus.temp)>0, data.genus.temp)


# Obj to plot
tabella <- psmelt(data.genus.temp)
tabella$GroupFunction<-factor(tabella$Group_Function, levels=c("PAO","GAO","AOB","NOB","Others PHA"))
# extra focus on Accumuli and Competi for the plots ...
tabella$GroupFunction2<- as.character(tabella$GroupFunction) # extra line to allow a quick repeating of this section if needed...
tabella$GroupFunction2[ tabella$Genus %in% c("Candidatus Accumulibacter","Candidatus Competibacter","Accumulibacter","Competibacter") ] <- tabella$Genus [tabella$Genus %in% c("Candidatus Accumulibacter","Candidatus Competibacter","Accumulibacter","Competibacter") ]
tabella$GroupFunction2[ ! tabella$Genus %in% c("Candidatus Accumulibacter") &  tabella$Genus %in% PAO_list ] <- "Other PAOs"
tabella$GroupFunction2[ ! tabella$Genus %in% c("Candidatus Competibacter") &  tabella$Genus %in% GAO_list ] <- "Other GAOs"
tabella$GroupFunction2<- gsub("Candidatus ","Ca. ", tabella$GroupFunction2)
tabella$GroupFunction2<- factor( tabella$GroupFunction2 , levels= c("Ca. Accumulibacter", "Other PAOs", "Ca. Competibacter", "Other GAOs","AOB","NOB","Others PHA") )

# PLOTTING ...
ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=GroupFunction2)) +
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  facet_nested(~ AGS + Project_number , scales = "free_x", space = "free_x") +  
  scale_fill_manual(values=c("Ca. Accumulibacter"="coral", "Ca. Competibacter"="chartreuse2",
                             "Other PAOs"="red2","Other GAOs"="chartreuse4",
                             "AOB"="deepskyblue","NOB"="blue3",
                             "Others PHA"="darkred"
  )
  ) +
  theme(axis.text.x=element_text(angle=43, hjust=1,vjust=1, size=6),
        axis.title.y = element_text(size=11), 
        axis.text.y = element_text(size=8), 
        strip.text.x = element_text(size=11), 
        legend.key.height = unit(0.18, "cm"),
        legend.key.width = unit(0.52, "cm"),
        legend.text = element_text ( size = 12.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,2,2,2),
        panel.spacing.x = unit(1.45,"pt"),
        legend.margin =  margin(-10,1,0,0)
  ) +
  scale_y_continuous(lim=c(0,35), 
                     breaks=seq(0,100,5),
                     expand = c(0,1)
  ) +
  scale_x_discrete(expand=c(-0.1,1)) +
  labs(x="", y="Percentual abundance of groups", 
       fill="")
ggsave(file="Results/Comparison_with_particular_AGS/Functional_groups_averages_comparisons.png",width=6.75,height=4.5,dpi=300)
dev.off()  


# means of every PAO and GAO
tot2<- cbind.data.frame("Averages_without_unclassif"=as.numeric(apply(otu_table(data.genus.temp),1,mean)),
                        "Averages_in_FullScale"=as.numeric(apply (otu_table(subset_samples(data.genus.temp, Reactor_scale=="Full")),1,mean)) ,
                        "Averages_in_PilotScale"=as.numeric(apply (otu_table(subset_samples(data.genus.temp, Reactor_scale=="Pilot")),1,mean)),
                        "Averages_in_LabScale"=as.numeric(apply (otu_table(subset_samples(data.genus.temp, Reactor_scale=="Lab")),1,mean)),
                        "Genus"= as.data.frame(tax_table(data.genus.temp))[["Genus"]],
                        "Functional_Group"= as.data.frame(tax_table(data.genus.temp))[["Group_Function"]]
)
tot2 <- tot2[order(tot2$Averages_without_unclassif, decreasing = T), ]

write.csv(file = "Results/Comparison_with_particular_AGS/Functional_groups_averages_comparisons.csv", row.names = F, quote = F, 
          tot2
)


suppressWarnings(rm(tot, tot2, tabella, data.genus.temp) )




################# FOCUS ON  ON SESSILE CILIATES ###########################

data.genus.temp<-subset_taxa( data.genus.prop, Domain=="Eukaryota")   # NB: removing unclassified!
data.genus.temp<- subset_taxa(data.genus.temp, Genus %in%  sessilida_list )
tabella<-psmelt(data.genus.temp)
ggplot(data=tabella, aes( x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.7) +
  geom_bar(stat="identity", position="stack", width = 0.7, alpha= 1) +
  theme_classic(base_size =8.5) +
  #facet_nested(~ Influent + Reactor_scale, scales = "free_x", space = "free_x") +  
  facet_nested(~  AGS + Project_number, scales = "free_x", space = "free_x") +  
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1, size=6),
        axis.title.y = element_text(size=11), 
        axis.text.y = element_text(size=8), 
        strip.text.x = element_text(size=11), 
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.45, "cm"),
        legend.text = element_text ( size = 12.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,2,2,2),
        panel.spacing.x = unit(1.5,"pt"),
        legend.margin =  margin(-10,2,0,0)
  ) +
  # scale_y_continuous(lim=c(0,0.0001),
  #                    breaks=seq(0,0.0001,0.000025),
  #                    expand = c(0,1)
  # ) +
  scale_x_discrete(expand=c(-0.025,1)) +
  labs(x="", y="Percentual abundance (%) of clades", 
       fill="")
ggsave(file="Results/Comparison_with_particular_AGS/FOCUS_ON_SESSILIDA_comparison.png",width=6.75,height=4.5,dpi=300)
dev.off()  




################# VEEN DIAGRAMS AGSs (EVERY MICROBIAL GENUS, no virus) ##########################

data.venn<-subset_taxa(data.genus.prop, !Genus %in% c("unclassified") )
data.venn<-subset_taxa(data.venn, !Domain %in% c("Viruses") )
threshold <- 0.0001   # already filtered in SBR, let's bee more lenient here
always<-taxa_names(filter_taxa(data.venn, function(x) max(x) >= threshold, TRUE)) 
data.venn<-prune_taxa(always, data.venn)
sample_percent <- 50   # there are much less samples and with different settings for each project, thus let's keep this more lenient than in SBR ... 

SBR<-subset_samples(data.venn, AGS=="SBR AGS" )
SBR<- prune_taxa(taxa_sums(SBR) > 0, SBR)
#always<-taxa_names(filter_taxa(Lab, function(x) min(x) > threshold, TRUE))
#Lab<-prune_taxa(always, Lab)
who<-apply(otu_table(SBR), MARGIN=2, function(x) ifelse(x>0, 1, 0))
how_many_samples <- round( length(sample_names(SBR))* sample_percent / 100 , 0 )
who<- who[rowSums(who)>= how_many_samples, ]
SBR<-as.vector(tax_table(SBR)[row.names(who),"Genus"])

PN<-subset_samples(data.venn, AGS=="CFR" & Project=="PRJNA1082061" )
PN<- prune_taxa(taxa_sums(PN) > 0, PN)
who<-apply(otu_table(PN), MARGIN=2, function(x) ifelse(x>0, 1, 0))
how_many_samples <- round( length(sample_names(PN))* sample_percent / 100 , 0 )
who<- who[rowSums(who)>= how_many_samples, ]
PN<-as.vector(tax_table(PN)[row.names(who),"Genus"])

CFR<-subset_samples(data.venn, AGS=="CFR" & Project!="PRJNA1082061")
CFR<- prune_taxa(taxa_sums(CFR) > 0, CFR)
who<-apply(otu_table(CFR), MARGIN=2, function(x) ifelse(x>0, 1, 0))
how_many_samples <- round( length(sample_names(CFR))* sample_percent / 100 , 0 )
who<- who[rowSums(who)>= how_many_samples, ]
CFR<-as.vector(tax_table(CFR)[row.names(who),"Genus"])


ONLY_IN_PN<- PN[! PN %in% CFR & ! PN %in% SBR]
ONLY_IN_PN<- paste(ONLY_IN_PN, collapse = ", ")
head(ONLY_IN_PN)

ONLY_IN_SBR<- SBR[! SBR %in% CFR & ! SBR %in% PN]
ONLY_IN_SBR<- paste(ONLY_IN_SBR, collapse = ", ")
head(ONLY_IN_SBR)

ONLY_IN_CFR<- CFR[! CFR %in% SBR & ! CFR %in% PN]
ONLY_IN_CFR<- paste(ONLY_IN_CFR, collapse = ", ")
head(ONLY_IN_CFR)


x<-list(SBR=SBR,"CFR PN"=PN,CFR=CFR)
ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse2","deepskyblue","coral")) +
  labs(title="Microbial genera featured in every sample\n") +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) )
# ggsave(filename = paste0("Results/Comparison_with_particular_AGS/Venn_Diagramm_MICROBIAL_GENERA_AGSs_over",threshold,"min.png"), width = 4, height = 4, dpi=300, bg = "white")
# dev.off()

### NB: HIGHER "COMMON" GENERA NUMBER, AS THIS FILTER IS LESS STRINGENT WHEN APPLIED ON *EVERY* SBR AT THE SAME TIME (SEE PART 2)


common <- SBR[ SBR %in% CFR & SBR %in% PN ]
common_CFR <- CFR [ CFR %in% PN ]
# length(common)
data_Genus_common<- subset_taxa(data.venn, Genus %in% common)

tot2<- cbind.data.frame("Overall_averages"=as.numeric(apply(otu_table(data_Genus_common),1,mean)),
                        "Averages_SBRAGS"=as.numeric(apply (otu_table(subset_samples(data_Genus_common, AGS=="SBR AGS")),1,mean)) ,
                        "Averages_CFR"=as.numeric(apply (otu_table(subset_samples(data_Genus_common, AGS=="CFR" & Project!="PRJNA1082061")),1,mean)),
                        "Averages_PNCFR"=as.numeric(apply (otu_table(subset_samples(data_Genus_common, AGS=="CFR" & Project=="PRJNA1082061")),1,mean)),
                        "Genus"= as.data.frame(tax_table(data_Genus_common))[["Genus"]]
)
tot2 <- tot2[order(tot2$Overall_averages, decreasing = T), ]
write.csv(file = paste0("Results/Comparison_with_particular_AGS/Venn_Diagr_MICROBIAL_GENERA_AGSs_averages___over",threshold,"min.csv"), row.names = F, quote = F, 
          tot2
)



### Adding a symbol in the SBR-only tables

temp <- read.csv("Results/Abundances/TOP_COMMON_microbial_genera_Averages_no_unclassified.csv")
temp$Also_in_every_CFR <- rep("")
temp$Also_in_CFR <- rep("")
temp$Also_in_PN_CFR <- rep("")
temp$Also_in_every_CFR [ temp$Genus %in% common_CFR ] <- "*"
temp$Also_in_CFR [ temp$Genus %in% CFR ] <- "*"
temp$Also_in_PN_CFR [ temp$Genus %in% PN ] <- "*"
write.csv(temp, 
          file = "Results/Abundances/TOP_COMMON_microbial_genera_Averages_no_unclassified.csv", 
          row.names = F, quote=F
          )

temp <- read.csv("Results/Abundances/COMMON_EUKA_OR_ARCHAEA.csv")
temp$Also_in_every_CFR <- rep("")
temp$Also_in_CFR <- rep("")
temp$Also_in_PN_CFR <- rep("")
temp$Also_in_every_CFR [ temp$Genus %in% common_CFR ] <- "*"
temp$Also_in_CFR [ temp$Genus %in% CFR ] <- "*"
temp$Also_in_PN_CFR [ temp$Genus %in% PN ] <- "*"
write.csv(temp, 
          file= "Results/Abundances/COMMON_EUKA_OR_ARCHAEA.csv",
          quote=F, row.names = F
)



suppressWarnings(rm(ONLY_IN_PN, ONLY_IN_SBR, ONLY_IN_CFR,
                    x, con, PN, SBR, CFR, data.venn, who))



