##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggvenn")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("DESeq2")
  library("vegan")
  # utilities
  library("xlsx")  
  library("Hmisc")
}

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one



####################### IMPORTING DATA #####################

load(file = "Data_prepared_by_1st_part_of_the_script.RData")

table_PHA<-read.table(file="Lista_accumulatori_PHA_da_articoli_5aprile2023.tsv", fill = T, header = F, sep="\t")



##################### PCoA OF R1 and R2 (including the period of equal feast/famine) ##################

if("data_complete" %in% ls() & "proof1" %in% ls() & "proof2" %in% ls() & "proof3" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete",
                                          "colors","fill_color_5","fill_color_10",
                                          "table_PHA",
                                          "unfiltered_data")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Sample_name !="F13" ) # the related S sample was unsatured --> discarded during the paired analysis
data<-subset_samples(data, ! Sample_Type %in% c("CuoioDepur_Aero","Inoculum")  )
data<-subset_samples(data, Experiment_state !="Adherent"  ) 
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day ))

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
### with the color shift
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2", shape="Sample_Type") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(315,409), midpoint = 350,
                        breaks=c(315, 335, 365, 390, 409)
  ) +
  #scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.28,
            arrow=arrow(length =unit(0.23,"cm"), type = "closed") 
  ) +
  #geom_text(aes(label=sample_names(data.prop.labels)), color="black")+
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.title = element_text(size=9),
        legend.margin = margin(0,0,0,-10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none",
  ) +
  labs(title="",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_div_Hellinger_on_Genera_shades_R1R2_INCLUDING_DAYS_SAME_RATIO_FF.png", width = 4.2, height = 3.8, dpi=300)
# with the color shift and labels
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2", shape="Sample_Type") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(315,409), midpoint = 350,
                        breaks=c(315, 335, 365, 390, 409)
  ) +
  #scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.35) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.28,
            arrow=arrow(length =unit(0.23,"cm"), type = "closed") 
  ) +
  #geom_text(aes(label=sample_names(data.prop.labels)), color="black")+
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.title = element_text(size=9),
        legend.margin = margin(0,0,0,-10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none",
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day2), color="black", vjust= 0.08, hjust=0.8, size=2.05, show.legend = FALSE) +
  labs(title="",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_div_Hellinger_on_Genera_shadesANDtxt_R1R2_INCLUDING_DAYS_SAME_RATIO_FF.png", width = 4.2, height = 3.8, dpi=300)




#######################  !!!!!  STARTING THE ANALYSIS OF R1 vs R2 !!!!!   #############################

if("data_complete" %in% ls() & "proof1" %in% ls() & "proof2" %in% ls() & "proof3" %in% ls() ){
to_remove<-ls()
to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete",
                                        "colors","fill_color_5","fill_color_10",
                                        "table_PHA",
                                        "unfiltered_data")]
rm(list=to_remove)
}

data<-subset_samples(data_complete, Experiment_state %in% c("Started","Almost_stabilized") )
data<-subset_samples(data, Sample_name !="F13" ) # the related S sample was unsatured --> discarded during the paired analysis
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
  data.class.prop <- transform_sample_counts(data.class, function(x) x/sum(x)*100)
  data.order.prop <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(x) x/sum(x)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

# adding informations to missing names
taxa_temp<-Taxa.genus
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
  Taxa.genus.update<-taxa_temp
}

rm(taxa_temp)




###################### ABUNDANCES BAR PLOT _ R1 versus R2 ##########################


# TOP 5 Phylum con stacked bar plot
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.phy.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
# tabella$Sampling_date<-gsub("_22","",tabella$Sampling_date)
unique(tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
tabella$Sampling_date<-gsub("_","/",tabella$Sampling_date, fixed=T)
# order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
order_time<-tabella$Sampling_date[c(grep("/04",tabella$Sampling_date),grep("/05",tabella$Sampling_date),grep("/06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Phylum)) +
  theme_classic(base_size =14) + 
  #facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  facet_grid2( . ~ Experiment_day+Sample_Type, scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
        strip.text = element_text(size=10.5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 14.1 ),
        legend.position="bottom",legend.margin = margin(2,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Experiment day", y="Percent abundance",
       fill="",
       # title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
#ggsave(file="Results/P_lim_reactor/R1_versus_R2/Abundances/TOP5_phyla_abundances.png",width=7,height=10.3, dpi=300) 
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Abundances/TOP5_phyla_abundances.png",height=6.5,width=8.2, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/P_lim_reactor/R1_versus_R2/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))



# TOP 10 Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:10]
prune.dat_top <- prune_taxa(top,data.genus.prop)
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_selected<-Taxa.genus.update[row.names(tax_selected),]
tax_table(prune.dat_top)<-as.matrix(tax_selected)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
tabella$Sampling_date<-gsub("_","/",tabella$Sampling_date, fixed=T)
# order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
order_time<-tabella$Sampling_date[c(grep("/04",tabella$Sampling_date),grep("/05",tabella$Sampling_date),grep("/06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Genus)) + 
  theme_classic(base_size =14) + 
  #facet_grid2(Sampling_date+Sample_Type~.,
   #           scales = "free", space="free",
   #           strip = strip_nested(size="constant"))+
  facet_grid2( . ~ Experiment_day+Sample_Type, scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
        strip.text = element_text(size=10.5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12.5 ),
        legend.position="bottom",
        legend.margin = margin(0,0,0,0)) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Experiment day", y="Percent abundance",
       fill="",
       #title = "Ten most abundant genera",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Abundances/TOP10_genera_abundances.png",height =6.5,width = 8.2, dpi=300) 
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/P_lim_reactor/R1_versus_R2/Abundances/TOP_10_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



####################### ABUNDANCES OF PHA BACTERIA _ R1 versus R2 ###################

lista_PHA<-table_PHA[,1]     # they derive from papers

prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))

tabella$Experiment_state<-gsub("_"," ",tabella$Experiment_state)
# factoring also the sample names to ensure a certain order in the plot (=time) ; in the sample_names they are already ordered
tabella$Sample<-factor(tabella$Sample, levels=sample_names(data))
tabella$Genus<-gsub("Candidatus_","",tabella$Genus)

fill_color_14<-c("bisque2","deeppink","darkgreen","darkmagenta","cyan","yellow2","brown","firebrick3","wheat3","violet","darkslategray3","springgreen3","blue","deepskyblue2") # , "gray", "pink3","yellow4","red","springgreen3","coral") 

ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=Genus)) + 
  facet_grid2( Sample_Type~ . , scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", alpha= 0.8) +
  theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_14)) +
  theme(axis.text.x=element_text(angle=35, hjust=1,vjust=1, size=9.5),
        axis.text.y = element_text(size=9), 
        axis.title.y = element_text(size=10.8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.7 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        legend.margin =  margin(3,22,0,0)) +
  guides(fill=guide_legend(nrow=3)) + 
  scale_y_continuous(lim=c(0,100)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = paste("PHA accumulating genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"))
ggsave(file="Results/P_lim_reactor/R1_versus_R2/PHA_accumulators_along_the_time.png",width=6.8,height=5.8,dpi=300)
dev.off()



########################## ALFA DIVERSITY _ R1 versus R2 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Sample_Type", color="Sample_Type")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H <- pAlpha$data[pAlpha$data$variable=="Shannon", ]
obs <- pAlpha$data[pAlpha$data$variable=="Observed", ]
# adding evenness
{identical(H$Sampling_date, obs$Sampling_date) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Sample_Type, New_data$Sampling_date),]
pAlpha$data<-New_data
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
# with points only
pAlpha +
  theme_classic2() +
  scale_color_manual(values = c("R1"="coral", "R2"="red3")) +
  geom_point(size=3.1, alpha= 0.4) +
  geom_line(aes(group = pAlpha$data$Sampling_date),col="grey",size=0.15) +
  labs(x="", title="Alpha diversity between R1 and R2 reactors") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Sample_Type), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.47) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Alfa_diversity_GENUS_paired_wilcoxon.png", width = 6,height =5.2, dpi=300)
# with ID
pAlpha + theme_classic2() +
  scale_color_manual(values = c("R1"="coral", "R2"="red3")) +
  geom_point(size=3.1, alpha= 0.4) +
  geom_line(aes(group = pAlpha$data$Sampling_date),col="grey",size=0.15) +
  geom_text(aes(label=Sample_name), size= 2.5, color= "black") +
  labs(x="", title="Alpha diversity between R1 and R2 reactors") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Sample_Type), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.47) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Alfa_diversity_GENUS_paired_wilcoxon_with_ID.png", width = 6,height =5.2, dpi=300)
# with experiment day label
pAlpha + theme_classic2() +
  scale_color_manual(values = c("R1"="coral", "R2"="red3")) +
  geom_point(size=3.1, alpha= 0.4) +
  geom_line(aes(group = pAlpha$data$Sampling_date),col="grey",size=0.15) +
  geom_text(aes(label=Experiment_day), size= 2.65,
            vjust= 0.75,
            hjust= ifelse(pAlpha$data$Sample_Type=="R1",1,0),
            color= "black", alpha = 0.9) +
  labs(x="") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Sample_Type), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.47) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Alfa_diversity_GENUS_paired_wilcoxon_with_day.png", width = 5.2,height =4.8, dpi=300)



##################### BETA DIVERSITY _ R1 versus R2 #######################

suppressWarnings(rm(ASV.prop))

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

###### PERMANOVA

metadata<-as(sample_data(data.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Sampling_date + Sample_Type, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Sampling_date + Sample_Type, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Sampling_date + Sample_Type, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Sampling_date + Sample_Type, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Sampling_date + Sample_Type, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Sampling_date + Sample_Type, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/P_lim_reactor/R1_versus_R2/Beta_div/Beta_divers_permanova_Helling.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Sample_Type)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/P_lim_reactor/R1_versus_R2/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)


######################## PCoA BETA DIV _ R1 versus R2 #########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day))

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  geom_line(aes(group=Sampling_date),col="grey", size=0.15)+
  scale_color_manual(values=c("R2"="red3","R1"="coral")) +
  geom_point(size=3, alpha=0.3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sampling_date), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), 
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H),
       caption="lines connect paired samples")
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 7, height = 5, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  geom_line(aes(group=Sampling_date),col="grey", size=0.15)+
  scale_color_manual(values=c("R2"="red3","R1"="coral")) +
  geom_point(size=3, alpha=0.3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples")
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_div_Hellinger_genera_no_ellipse.png", width = 7, height = 5, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  geom_line(aes(group=Sampling_date),col="grey", size=0.15)+
  scale_color_manual(values=c("R2"="red3","R1"="coral")) +
  geom_point(size=3, alpha=0.3) + theme_classic(base_size = 14) +
  labs(title="", color="", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples")
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_points.png", width = 7, height = 5, dpi=300)


### again but with the color shift
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2", shape="Sample_Type") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(339,409), midpoint = 370,
                        breaks=c(339, 365, 390, 409)
  ) +
  #scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.23,"cm"), type = "closed") 
  ) +
  #geom_text(aes(label=sample_names(data.prop.labels)), color="black")+
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.title = element_text(size=9),
        legend.margin = margin(0,0,0,-10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none",
  ) +
  labs(title="",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_dive_Hellinger_on_Genera_shades_R1R2.png", width = 4.2, height = 3.8, dpi=300)
# with shades and also day text
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2", shape="Sample_Type") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(339,409), midpoint = 370,
                        breaks=c(339, 365, 390, 409)
  ) +
  #scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.23,"cm"), type = "closed") 
  ) +
  #geom_text(aes(label=sample_names(data.prop.labels)), color="black")+
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.title = element_text(size=9),
        legend.margin = margin(0,0,0,-10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none",
  ) +
  labs(title="",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/R1_versus_R2/Beta_div/PCoA_Beta_dive_Hellinger_on_Genera_shades_and_text_R1R2.png", width = 4.2, height = 3.8, dpi=300)



############## DA WITH DESEQ2 _ R1 versus R2 #############

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering steps yet?\n", fill = T)
}


suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DESeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
# removing the unpaired sample
system(" echo 'In case of pair analysis, the log2foldchange displayes the change from R1 to R2, see https://support.bioconductor.org/p/105981/ \n' > Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/NB.txt ")  
system(" echo 'Every result with a baseMean value lower than 50 has been arbitrarly filtered.' >> Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/NB.txt ")  

Table_tot<-NULL
Res_tot<-NULL

# correcting also the time (with more time, the difference increases)
sample_data(data_pruned)$time<-as.factor(gsub("^[0-9][0-9]_","",sample_data(data_pruned)$Sampling_date))
### --> error if included in the model: too much auto-correlations among the variables (the time itself is the paired factor already!)
# then, the april samples (low difference) will be discarded to avoid this bias
data_pruned2<- subset_samples(data_pruned, ! time %in% c("03_23","04_23") )


for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned2, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # if the paired samples are seen as numbers (it depends from on levels name)
  sample_data(d)[["Sampling_date"]]<-as.character(sample_data(d)[["Sampling_date"]])
  
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sampling_date +Sample_Type) # interaction: more time passed = even more differences!
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sample_Type", "R1", "R2"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/DA_",t,"_ratio_R1_vs_R2.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Sample_Type<-factor(Table_tot$Sample_Type, levels = unique(Table_tot$Sample_Type))
Table_tot<-Table_tot[order(match(Table_tot$Sample_Type, levels(Table_tot$Sample_Type))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Sample_Type)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Sample_Type)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

### to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Sample_Type)


# marking the PHA accumulators
lista_PHA<-table_PHA[,1]     # they derive from papers


# plotting
unique(tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_ f ","\nf_",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_ o ","\no_",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_marine_group","\nmarine_group",tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Sample_Type)) +
  theme_classic(base_size = 9.1) +
  scale_color_manual(values = c("R2"="red3","R1"="coral")) +
  facet_wrap2(nrow=4,~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Sample_Type), size=0.5, alpha=1) +
  geom_point(aes(color=Sample_Type), size=2.8, alpha=0.3) +
  theme(strip.text.x=element_text(size=8.2,colour="black"),
        strip.switch.pad.wrap = unit(3,"line"),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=9.5)
        )+
  scale_x_discrete(labels=unique(levels(tabella_g$Sample_Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
plot_1 +
  geom_line(aes(group=Sampling_date), size=0.2) +
  labs(title= "Differently abundant genera",
       y="Percent abundance", fill=" ", x="")
ggsave(filename = "Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/Plot_genera_DeSeq2_Type.png", width = 10, height = 7, dpi=300)
# again but with name
plot_1 +
  geom_text(aes(label=Experiment_day), size=2.3, color="black", alpha= 0.8,
            hjust=ifelse(tabella_g$Sample_Type=="R1",1,0)
            ) +
  geom_line(aes(group=Sampling_date), size=0.2, col="darkgray") +
  labs(title= "Differently abundant genera",
    y="Percent abundance", fill="Sample_Type", x="")
ggsave(filename = "Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/Plot_genera_DeSeq2_Type_with_name.png", width = 10, height = 7, dpi=300)
dev.off()


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
#Redund<-c("") # to select redundants (same ASV, same results at different levels)
tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)

tabella_3<-subset(tabella_2, Taxa %in% c("Phylum","Class"))
plot_3<-ggplot(tabella_3, aes(x= Xaxis, y=Abundance, fill=Sample_Type)) +
  theme_classic(base_size = 9) +
  scale_color_manual(values = c("R2"="red3","R1"="coral")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Sample_Type), size=0.7, alpha=1) +
  geom_point(aes(color=Sample_Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=8.4,colour="black"), 
        strip.switch.pad.wrap = unit(2,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_3$Sample_Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  labs(title= "Differently abundant families and orders",
    y="Percent Abundance", fill=" ", x="") +
  geom_line(aes(group=Sampling_date), size=0.2) 

tabella_4<-subset(tabella_2, ! Taxa %in% c("Phylum","Class"))
tabella_4$Bacteria<-gsub("-","\n",tabella_4$Bacteria, fixed = F)
# tabella_4$Bacteria<-gsub("Peptostreptococcales","Peptostreptococc.",tabella_4$Bacteria)
plot_4<-ggplot(tabella_4, aes(x= Xaxis, y=Abundance, fill=Sample_Type)) +
  theme_classic(base_size = 9) +
  scale_color_manual(values = c("R2"="red3","R1"="coral")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Sample_Type), size=0.7, alpha=1) +
  geom_point(aes(color=Sample_Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=7.8,colour="black"), 
        strip.switch.pad.wrap = unit(2,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_4$Sample_Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  labs(title= "Differently abundant families and orders",
       y="Percent Abundance", fill=" ", x="") +
  geom_line(aes(group=Sampling_date), size=0.2) 

# now an unique plot (no dates)
head_plot<-plot_4 + # refining first plot with the title
  labs(title= "Differently abundant families, orders, class and phyla", y="Percent abundance", x="") + 
  #theme(plot.title = element_text(size=12.5))
png(filename = "Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants.png",
    width = 2750, height = 1500, res=300)
grid.arrange(head_plot,
             plot_3 + labs(title= "", y="Percent abundance", x=""))
dev.off()

# again, an unique plot (with dates)
head_plot<-plot_4 + # refining first plot with the title
  labs(title= "Differently abundant families, orders, class and phyla", y="Percent abundance", x="") + 
  theme(plot.title = element_text(size=12.5)) + 
  geom_text(aes(label=Experiment_day), size=2.3, color="black", alpha= 0.8,
            hjust=ifelse(plot_4$data$Sample_Type=="R1",1,0)
  )
plot3_bis<- plot_3 +
  geom_text(aes(label=Experiment_day), size=2.3, color="black", alpha= 0.8,
                                 hjust=ifelse(plot_3$data$Sample_Type=="R1",1,0)
            )
png(filename = "Results/P_lim_reactor/R1_versus_R2/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants_WITH_DATE.png", width = 2750, height = 1500, res=300)
grid.arrange(head_plot,
             plot3_bis + labs(title= "", y="Percent abundance", x=""))
dev.off()


# for comparisons...
Genera.DESEQ2_R1_vs_R2<-unique(tabella_g[tabella_g$Taxa=="Genus","Bacteria"])
rm(tabella_2, plot_2, plot_1, head_plot)



#######################   !!!!!  STARTING THE ANALYSIS OF TIME SERIES in R1 !!!!!   #############################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete",
                                          "colors","fill_color_5","fill_color_10",
                                          "table_PHA",
                                          "unfiltered_data")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Sample_Type == "R1" )
data<-subset_samples(data, Experiment_state != "Adherent" )
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
  data.class.prop <- transform_sample_counts(data.class, function(x) x/sum(x)*100)
  data.order.prop <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(x) x/sum(x)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

# adding informations to missing names
taxa_temp<-Taxa.genus
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
  Taxa.genus.update<-taxa_temp
}

rm(taxa_temp)

sample_data(data)$Experiment_state<-factor(sample_data(data)$Experiment_state, levels=c("Still_preparing","Started","Almost_stabilized")) # decides the order in plots

days_ordered<-sample_data(data)$Sampling_date[order(sample_data(data)$Sampling_date)]
order_time<-days_ordered[c( grep("_03",days_ordered), grep("_04",days_ordered),grep("_05",days_ordered),grep("_06",days_ordered))] # months
sample_data(data)$Sampling_date<-factor(sample_data(data)$Sampling_date, levels = unique(order_time))

metadata<-metadata_complete[metadata_complete$Sample_Type=="R1" & metadata_complete$Experiment_state!="Adherent", ]


######################## ABUNDANCES BAR PLOT _ TIME SERIE R1 ##########################


# TOP 5 Phyla
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
tabella$Experiment_state<-factor(tabella$Experiment_state, levels = unique(tabella$Experiment_state))
ggplot(data=tabella, aes(x=Abundance, y=Experiment_day, fill=Phylum)) +
  theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 14 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Experiment day", x="Percent abundance",
       title = "Five most abundant phyla (R1)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Abundances/TOP_5_phyla.png",width=7,height=7.2, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/P_lim_reactor/Time_series_in_R1/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, prune.dat_top,tabella, tabella_top)


# TOP 10 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
tabella$Experiment_state<-factor(tabella$Experiment_state, levels = unique(tabella$Experiment_state))
ggplot(data=tabella, aes(y=Experiment_day, x=Abundance, fill=Genus)) +
  # facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 14 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Experiment day", x="Percent abundance",
       title = "Ten most abundant genera (R1)", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Abundances/TOP_10_Generi.png",width=7,height=7.2,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/P_lim_reactor/Time_series_in_R1/Abundances/TOP_10_genera_of_all_dataset_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



########### LINE PLOTS OF MOST ABUNDANT GENERA ##############

#fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# TOP 5 Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =14) + 
  geom_point(aes(color=Genus), size =1.8) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_continuous(breaks = seq(0, max(tabella$Abundance), 10)) +
  labs(x="Experiment day", y="Percent abundance", title = "5 most abundant genera through time (P limiting R1)")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Abundances/TOP5_Genera_along_time.png",width=7,height=5, dpi=300) 
dev.off()


########################## ALFA DIVERSITY _ TIME SERIE R1 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

dir.create("Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/")
pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Experiment_state", color="Experiment_state")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
{ H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
  obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
  # adding evenness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Experiment_state<-factor(pAlpha$data$Experiment_state, levels = levels(sample_data(data)$Experiment_state))
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Experiment_state, y=value, color=NULL), alpha=0) + theme_classic2() + 
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.2, alpha= 0.4) + 
  labs(x="Experiment state", title="Alpha diversity between time groups") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=28, vjust=1, hjust=1, size=9.8)) +
  stat_compare_means(aes(group = Experiment_state), label="p.format", method = "kruskal.test",
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.4, hjust=0.18)
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/Alfa diversity_on genera_Kruskal test.png", width = 6,height =5.2, dpi=300)
# again, but smaller (poster version)
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Experiment_state, y=value, color=NULL),
               alpha=0, size = 0.25) +
  theme_classic2(base_size = 8.6) + 
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=2.9, alpha= 0.3) + 
  labs(x="", title="") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5),
        axis.title.y = element_text(size=7)
        ) +
  stat_compare_means(aes(group = Experiment_state), label="p.format", method = "kruskal.test",
                     label.x= 1.5, size=2.3, label.y.npc = "top", vjust=-0.1, hjust=0.22)
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/Alfa diversity_on genera_Kruskal_poster.png", width = 3.7,height =2.4, dpi=300)


# Statistical analyses
alphadt<- as.data.frame(pAlpha$data)
Obser_value<-alphadt[alphadt$variable=="Observed richness", ]
factor<-Obser_value$Experiment_state
kruskal.test(Obser_value$value~factor) # just to test the plotted p value

# pairwise comparisons
con <- file("Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_of_experiment_state.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Experiment_state sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Shannon", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Evenness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)


# correlation through time (Observed r)
Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=Sample_name), size= 3.7) +
  geom_line(aes(group=Sample_Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Observed richness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 2))
       )
ggsave("Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/Dot_plot_Observed_richness_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Shannon)
Filter_value<-alphadt[alphadt$variable=="Shannon", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=Sample_name), size= 3.7) +
  geom_line(aes(group=Sample_Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Shannon",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/Dot_plot_Shannon_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Evenness)
Filter_value<-alphadt[alphadt$variable=="Evenness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=Sample_name), size= 3.7) +
  geom_line(aes(group=Sample_Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Evenness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/P_lim_reactor/Time_series_in_R1/Alpha_diversity/Dot_plot_Evenness_through_time.png", width = 5, height = 3, dpi=300)

rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor, Resulting_corr)



######################## BETA DIVERSITY _ TIME SERIE R1 #######################

suppressWarnings(rm(ASV.prop))

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_H<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_genus<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
pair_genus<- pairwise.adonis(sample_genus, factors=metadata$Experiment_state, p.adjust.m = "BH", sim.method="euclidean", perm = 9999)
pair_genus

# exporting beta diversity
suppressWarnings(rm(con))
con<-file("Results/P_lim_reactor/Time_series_in_R1/Beta_div/Beta_divers_general_and_pairwise_between_Experiment_state_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity on Hellinger \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity on Hellinger (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_genus
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
{BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,metadata$Experiment_state)
  disp_genus<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_genus$tab
  a<-as.data.frame(disp_genus$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/P_lim_reactor/Time_series_in_R1/Beta_div/Beta_dispersion_General_and_Pairwise_between_Experiment_state_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger \n")
disp_genus$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_genus,a, con)



####################### PCoA BETA DIV _ TIME SERIE R1 #########################

### on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day))
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
#   scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
#   geom_point(size=3.4, alpha=0.4) +
#   geom_point(size=1, alpha=1) +
#   theme_classic(base_size = 11) + stat_ellipse(size=0.2)  +
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment_state",
#        x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 6.5, height = 4.5, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
#   scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
#   geom_point(size=3.4, alpha=0.4) +
#   geom_point(size=1, alpha=1) +
#   theme_classic(base_size = 11) +
#   geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day), color="black", size=2.8, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="", 
#        x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_no_ellipses.png", width = 6.5, height = 4.5, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(315,409), midpoint = 380) +
  #scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_line(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.2,"cm"), type = "closed") 
            ) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  labs(title="PCoA with Hellinger distance\n on reactor R1 in P limitation",
       color="Experiment day",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
       )
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 4.2, height = 3.8, dpi=300)
# with shades and day text
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(315,409), midpoint = 380) +
  #scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_line(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.2,"cm"), type = "closed") 
  ) +
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n on reactor R1 in P limitation",
       color="Experiment day",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 4.2, height = 3.8, dpi=300)



################ CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R1 (PROPORTION) #################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "Sample_name"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "spearman") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_BH<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_BH<0.05,"*","")

write.csv2(Corr_results, "Results/P_lim_reactor/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1<-row.names(temp_for_save)


con <- file("Results/P_lim_reactor/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R1 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R1[1:6]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Experiment day", 
       y="Percent abundance",
       color="",
       title = "the 6 genera most positely correlated with time (R1) \naccording to spearman correlation on % counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R1[(length(significative_abundances_R1)-5):length(significative_abundances_R1)]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Experiment day", y="Percent abundance",
       color="",
       title = "the 6 genera most negatively correlated with time (R1)\n according to spearman correlation on % counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_5, prune.dat_top, tax_selected, tabella, corr_6_top))



########## CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R1 (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "Sample_name"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "spearman") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_BH<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_BH<0.05,"*","")

write.csv2(Corr_results, "Results/P_lim_reactor/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_R1_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )



con <- file("Results/P_lim_reactor/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R1 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))



########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R1_pos
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 9 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=6)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Experiment day", y="Percent abundance",
       color="",
       title = "genera positely correlated with time (R1)\n according to spearman correlation on CLR counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/Genera_positively_correlated_with_time_R1.png",width=7,height=5, dpi=300) 
dev.off()



# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
#corr_6_top<-significative_abundances_R1[(length(significative_abundances_R1)-5):length(significative_abundances_R1)]
corr_6_top<-significative_abundances_R1_neg
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=6)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Experiment day",
       y="Percent abundance",
       color="",
       title = "genera negatively correlated with time (R1)\n according to spearman correlation on CLR counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/Genera_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_6, prune.dat_top, tax_selected, tabella, corr_6_top))



#######################   !!!!!  STARTING THE ANALYSIS OF TIME SERIES in R2 !!!!!   #############################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete",
                                          "colors","fill_color_5","fill_color_10",
                                          "table_PHA",
                                          "unfiltered_data")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Sample_Type == "R2" )
data<-subset_samples(data, Experiment_state != "Adherent" )
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
  data.class.prop <- transform_sample_counts(data.class, function(x) x/sum(x)*100)
  data.order.prop <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(x) x/sum(x)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

# adding informations to missing names
taxa_temp<-Taxa.genus
{for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
  for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
    taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
  for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
  for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
    taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
  Taxa.genus.update<-taxa_temp
}

rm(taxa_temp)

sample_data(data)$Experiment_state<-factor(sample_data(data)$Experiment_state, levels=c("Still_preparing","Started","Almost_stabilized")) # decides the order in plots

days_ordered<-sample_data(data)$Sampling_date[order(sample_data(data)$Sampling_date)]
order_time<-days_ordered[c( grep("_03",days_ordered), grep("_04",days_ordered),grep("_05",days_ordered),grep("_06",days_ordered))] # months
sample_data(data)$Sampling_date<-factor(sample_data(data)$Sampling_date, levels = unique(order_time))

metadata<-metadata_complete[metadata_complete$Sample_Type=="R2" & metadata_complete$Experiment_state!="Adherent", ]



######################## ABUNDANCES BAR PLOT _ TIME SERIE R2 ##########################


# TOP 5 Phyla
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
tabella$Experiment_state<-factor(tabella$Experiment_state, levels = unique(tabella$Experiment_state))
ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Experiment day", y="Percent abundance", title = "Five most abundant phyla (R2)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Abundances/TOP_5_phyla.png",width=7,height=4.5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/P_lim_reactor/Time_series_in_R2/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, prune.dat_top,tabella, tabella_top)


# TOP 10 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella_top$Genus<-gsub("uncultured_ f Paludibacteraceae","uncultured_ f Paludibact.",tabella_top$Genus)
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
tabella$Experiment_state<-factor(tabella$Experiment_state, levels = unique(tabella$Experiment_state))

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 10.3 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Experiment day", y="Percent abundance", title = "Ten most abundant genera (R2)", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Abundances/TOP_10_Generi.png",width=7,height=4.5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/P_lim_reactor/Time_series_in_R2/Abundances/TOP_10_genera_of_all_dataset_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))


# TOP 10 Genera of only "ALMOST STABILIZED" group
data.genus.subset<-subset_samples(data.genus.prop, Experiment_state=="Almost_stabilized")
{top <- names(sort(taxa_sums(data.genus.subset), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.genus.subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella_top$Genus<-gsub("uncultured_ f Paludibacteraceae","uncultured_ f Paludibact.",tabella_top$Genus)
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 10.5 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percent abundance", title = "Ten most abundant genera (R2)\n considering only the 'almost stabilized' samples", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Abundances/TOP_10_Generi_ONLY_STABILIZED_SAMPLES.png",width=7,height=4.8,dpi=300)
dev.off()



########### LINE PLOTS OF MOST ABUNDANT GENERA ##############

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# TOP 5 Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =14) + 
  geom_point(aes(color=Genus), size =1.8) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_continuous(breaks = seq(0, max(tabella$Abundance), 10)) +
  labs(x="Experiment day", y="Percent abundance", title = "5 most abundant genera through time (R2)")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Abundances/TOP5_Genera_along_time.png",width=7,height=5, dpi=300) 
dev.off()



########################## ALFA DIVERSITY _ TIME SERIE R2 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

dir.create("Results/P_lim_reactor/Time_series_in_R2/Alpha_diversity/")
pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Experiment_state", color="Experiment_state")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
{ H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
  obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
  # adding evenness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Experiment_state<-factor(pAlpha$data$Experiment_state, levels = levels(sample_data(data)$Experiment_state))
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Experiment_state, y=value, color=NULL), alpha=0) + theme_classic2() + 
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.2, alpha= 0.4) + 
  labs(x="Experiment state", title="Alpha diversity between time groups") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=28, vjust=1, hjust=1, size=9.8)) +
  stat_compare_means(aes(group = Experiment_state), label="p.format", method = "kruskal.test",
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.4, hjust=0.18)
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Alpha_diversity/Alfa diversity_on genera_Kruskal test.png", width = 6,height =5.2, dpi=300)


# Statistical analyses
alphadt<- as.data.frame(pAlpha$data)
Obser_value<-alphadt[alphadt$variable=="Observed richness", ]
factor<-Obser_value$Experiment_state
kruskal.test(Obser_value$value~factor) # just to test the plotted p value

# pairwise comparisons
con <- file("Results/P_lim_reactor/Time_series_in_R2/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_of_experiment_state.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Experiment_state sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Shannon", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Evenness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)


# correlation through time (Observed r)
Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=Sample_name), size= 3.7) +
  geom_line(aes(group=Sample_Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Observed richness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 2))
  )
ggsave("Results/P_lim_reactor/Time_series_in_R2/Alpha_diversity/Dot_plot_Observed_richness_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Shannon)
Filter_value<-alphadt[alphadt$variable=="Shannon", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=Sample_name), size= 3.7) +
  geom_line(aes(group=Sample_Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Shannon",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/P_lim_reactor/Time_series_in_R2/Alpha_diversity/Dot_plot_Shannon_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Evenness)
Filter_value<-alphadt[alphadt$variable=="Evenness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=Sample_name), size= 3.7) +
  geom_line(aes(group=Sample_Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Evenness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/P_lim_reactor/Time_series_in_R2/Alpha_diversity/Dot_plot_Evenness_through_time.png", width = 5, height = 3, dpi=300)

rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor, Resulting_corr)



######################## BETA DIVERSITY _ TIME SERIE R2 #######################

suppressWarnings(rm(ASV.prop))

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_H<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_genus<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
pair_genus<- pairwise.adonis(sample_genus, factors=metadata$Experiment_state, p.adjust.m = "BH", sim.method="euclidean", perm = 9999)
pair_genus

# exporting beta diversity
suppressWarnings(rm(con))
con<-file("Results/P_lim_reactor/Time_series_in_R2/Beta_div/Beta_divers_general_and_pairwise_between_Experiment_state_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity on Hellinger \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity on Hellinger (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_genus
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
{BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,metadata$Experiment_state)
  disp_genus<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_genus$tab
  a<-as.data.frame(disp_genus$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/P_lim_reactor/Time_series_in_R2/Beta_div/Beta_dispersion_General_and_Pairwise_between_Experiment_state_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger \n")
disp_genus$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_genus,a, con)


####################### PCoA BETA DIV _ TIME SERIE R2 #########################

### on Genera
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) + stat_ellipse(size=0.2)  +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment_state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 6.5, height = 4.5, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment state", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_no_ellipses.png", width = 6.5, height = 4.5, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_line(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.2,"cm"), type = "closed") 
  ) +
  theme_classic(base_size = 8.5) +
  theme(title=element_text(size=8),
        axis.title.x = element_text(size=7.8),
        legend.margin = margin(2,-2,2,-7)
  ) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera) on reactor R2", color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation","\n            (the X axis corresponds with the passing of time, but it is NOT computed to be so)"),
       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 5, height = 3.8, dpi=300)



################ CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R2 (PROPORTION) #################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "Sample_name"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "spearman") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_BH<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_BH<0.05,"*","")

write.csv2(Corr_results, "Results/P_lim_reactor/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R2<-row.names(temp_for_save)


con <- file("Results/P_lim_reactor/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R2 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[1:6]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percent abundance",
       color="",
       title = "the 6 genera most positely correlated with time (R2) \naccording to spearman correlation on % counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[(length(significative_abundances_R2)-5):length(significative_abundances_R2)]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling",
       color="",
       y="Percent abundance", title = "the 6 genera most negatively correlated with time (R2)\n according to spearman correlation on % counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(prune.dat_top, tax_selected, tabella, corr_6_top))



########## CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R2 (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "Sample_name"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "spearman") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_BH<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_BH<0.05,"*","")

write.csv2(Corr_results, "Results/P_lim_reactor/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R2<-row.names(temp_for_save)


con <- file("Results/P_lim_reactor/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R2 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3")

# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[1:6]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella$Genus<-gsub("uncultured_","uncult. ", tabella$Genus)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12.6 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percent abundance",
       color="",
       title = "the 6 genera most positely correlated with time (R2)\n according to spearman correlation on CLR counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[(length(significative_abundances_R2)-5):length(significative_abundances_R2)]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella$Genus<-gsub("uncultured_","uncult. ", tabella$Genus)
  tabella$Genus<-gsub("Rhodobacteraceae","Rhodobacterac. ", tabella$Genus)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.6 ),
        legend.position="bottom",
        legend.margin = margin(1,25,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percent abundance",
       color="",
       title = "the 6 genera most negatively correlated with time (R2)\n according to spearman correlation on CLR counts")
ggsave(file="Results/P_lim_reactor/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_6, prune.dat_top, tax_selected, tabella, corr_6_top))




################### DIFFERENTIAL ABUNDANCES WITH DESEQ2 _ Adherent vs Suspension (R1 and R2 together) #####################


if(! "unfiltered_data" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete",
                                          "colors","fill_color_5","fill_color_10",
                                          "table_PHA",
                                          "unfiltered_data")]
  rm(list=to_remove)
}

# data<-subset_samples(data_complete, Sample_Type == "R1" )  # using   ~ Sample_Type + Sampling area   below
data<-subset_samples(data_complete, Sampling_date %in% c("06_06_23","08_06_23","09_06_23",
                                                "22_05_23","22_05_23",
                                                "15_06_23","13_06_23","19_06_23")
) # only the data around the sampling of adherent biomass (only one in pair) have been selected

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
sample_data(data_pruned)$Sampling_area<-ifelse(sample_data(data_pruned)$Experiment_state=="Adherent","Adherent","Suspension")

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sample_Type + Sampling_area)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sampling_area", "Suspension", "Adherent"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    # write.csv2(r, file=paste0("Results/P_lim_reactor/DA_DESeq2/DA_",t,"_ratio_Suspension_vs_Adherent.csv"), row.names = F, quote=F, na = "")
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

#View(Res_tot)
write.csv(Res_tot, file="Results/P_lim_reactor/Adherent_vs_Suspended_R1_AND_R2_differential_abundances/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/P_lim_reactor/Adherent_vs_Suspended_R1_AND_R2_differential_abundances/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Sampling_area)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Suspension"="chartreuse","Adherent"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", 
       fill="Sampling area", x="")
#ggsave(filename = "Results/P_lim_reactor/DA_DESeq2/DA_Sampling_area_every_result.png", width = 13, height = 7, dpi=300)
#dev.off()

Redund<- c("Rhodocyclaceae" ,"Desulfovibrionales","Burkholderiales",
           "Verrucomicrobiaceae","Acholeplasmataceae","Desulfomicrobiaceae",
           "Desulfovibronia","Gammaproteobacteria","Verrucomicrobiae")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
Table_tot2<-subset(Table_tot2, ! Taxa %in% c("Phylum")) # to remove redundant results
#Table_tot2$Sampling_date<-gsub("_23","",Table_tot2$Sampling_date)
Table_tot2$Same_day<-Table_tot2$Sample_name
Table_tot2$Same_day[Table_tot2$Sampling_date!="22_05_23" ]<-"*" # if it is not the same day then it is labeled as *
plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Sampling_area) ) + 
  theme_classic() +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.7, alpha=0.15, size=0.22, outlier.alpha = 0) +
  theme_bw(base_size = 8.5) +
  scale_fill_manual(values=c("Suspension"="chartreuse","Adherent"="coral")) +
  theme(strip.text.x=element_text(size=8,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-15, 20, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=8),
        axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=7.5), 
        axis.text.y = element_text(size=5.8), 
        plot.title = element_text(size=10),
        legend.title = element_text(size=9.5),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,2, 4,seq(8,16,4),20, seq(25,max(Table_tot$Abundance+4),5))) +
  guides(color="none")
# with points
plot +
  theme_classic() +
  geom_point(aes(color=Sampling_area, shape=Sample_Type), position = position_dodge(width=0.7, preserve = "total"), size=0.45, alpha=0.7) +
  geom_point(aes(color=Sampling_area, shape=Sample_Type), position = position_dodge(width=0.7, preserve = "total"), size=1.5, alpha=0.3) +
  scale_color_manual(values=c("Suspension"="chartreuse4","Adherent"="coral")) +
  labs(title= "Differently abundant Taxa (using both R1 and R2 observations)", y="Proportional Abundance %", 
       fill="     Sampling area", shape="Reactor", x="")
ggsave(filename = "Results/P_lim_reactor/Adherent_vs_Suspended_R1_AND_R2_differential_abundances/DA_Sampling_area_no_redundants_POINTS.png", width = 5.5, height = 4.5, dpi=300)
dev.off()
# with names
plot + geom_text(aes(label=Same_day), 
            position = position_dodge(width = 0.7, preserve = "total"),
            size=2, alpha=0.9) +
  #scale_color_manual(values=c("same_day"="red","no"="black")) +
  labs(title= "Differently abundant Taxa (using both R1 and R2 observations)", y="Percent Abundance", 
       fill="     Sampling area",  x="",
       caption = "Only the samples collected the same day (paired) are evidenced through a label ( AF3---F14 , AS3---S14 )")
ggsave(filename = "Results/P_lim_reactor/Adherent_vs_Suspended_R1_AND_R2_differential_abundances/DA_Sampling_area_no_redundants_SAME_DAY.png", width = 5.5, height = 4.5, dpi=300)
dev.off()


system(" echo 'Only the sample with the sampling data around the one of the adherent biomass have been selected (only one adherent sample is in pair!). The model design then is unpaired BUT it considers the reactor Sample_Type as cofactor (~Reactor + Sampling area). Moreover, every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/P_lim_reactor/'Adherent_vs_Suspended_R1_AND_R2_differential_abundances'/NB.txt ")

