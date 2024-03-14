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

colors<-c("CuoioDepur Aero"="deepskyblue",  # groups color in PCoA plots
          "N limiting"="lightblue3" )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one



####################### IMPORTING DATA #####################

load(file = "Data_prepared_by_1st_part_of_the_script.RData")

table_PHA<-read.table(file="Lista_accumulatori_PHA_da_articoli_5aprile2023.tsv", fill = T, header = F, sep="\t")



############# PCoA WITH CUIODEPUR BEFORE SUBSETTING (CuoioD and Nlim Reactor) ##########################

suppressWarnings(rm(data.prop.labels, data.sqrt_prop))

data.prop.labels<-subset_samples(data_complete, Reactor_Type != "P_limitation")
data.prop.labels<-transform_sample_counts(data.prop.labels, function(x) (x/sum(x))*100)
sample_data(data.prop.labels)$Reactor_Type<-gsub("CuoioD_Aero","CuoioDepur Aero",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Reactor_Type<-gsub("N_limitation","N limiting",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Reactor_Type2<- sample_data(data.prop.labels)$Reactor_Type # to get unique lines for N limiting only in the plot
sample_data(data.prop.labels)$Reactor_Type2[sample_data(data.prop.labels)$Reactor_Type2=="CuoioDepur Aero"][1]<-"A"
sample_data(data.prop.labels)$Reactor_Type2[sample_data(data.prop.labels)$Reactor_Type2=="CuoioDepur Aero"][1]<-"B"
sample_data(data.prop.labels)$Reactor_Type2[sample_data(data.prop.labels)$Reactor_Type2=="CuoioDepur Aero"][1]<-"C"
sample_data(data.prop.labels)$Reactor_Type2[sample_data(data.prop.labels)$Reactor_Type2=="CuoioDepur Aero"][1]<-"D"

{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# Type 1 (names)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.2) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_text(aes(label= Sample_name), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/N_lim_reactor/Inoculum_comparison_with_CuoioDep/PCoA_with_CuioDep_names.png", width = 6, height = 4.5, dpi=300)
# Type 1 (no names)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.2) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Reactor_Type2), color="gray", size= 0.12) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/N_lim_reactor/Inoculum_comparison_with_CuoioDep/PCoA_with_CuioDep_lines.png", width = 6, height = 4.5, dpi=300)



################ VENN DIAGRAMM BETWEEN CUOIOD AND N limiting #####################

data.temp<-prune_taxa(taxa_sums(data_complete)>2, data_complete) # avoiding singletons and doubletons (maybe sequencing noises!)
data.genus.temp<-tax_glom(data.temp, taxrank = "Genus", NArm = F)
taxa_temp<-as.data.frame(tax_table(data.genus.temp))
for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])
}
tax_table(data.genus.temp)<-as.matrix(taxa_temp)
data.venn<-data.genus.temp

Inoculum<-subset_samples(data.venn, Sample_name %in% c("P1","P2") )
Inoculum<-as.character(tax_table(prune_taxa(taxa_sums(Inoculum)>0, Inoculum))[,"Genus"])

CuoioD<-subset_samples(data.venn, Reactor_Type=="CuoioD_Aero")
CuoioD<-as.character(tax_table(prune_taxa(taxa_sums(CuoioD)>0, CuoioD))[,"Genus"])

ONLY_IN_CuoioD<- CuoioD[ ! CuoioD %in% Inoculum]
ONLY_IN_CuoioD<- paste(ONLY_IN_CuoioD, collapse = ", ")
head(ONLY_IN_CuoioD)

ONLY_IN_inoculum<- Inoculum[ ! Inoculum %in% CuoioD]
ONLY_IN_inoculum<- paste(ONLY_IN_inoculum, collapse = ", ")
head(ONLY_IN_inoculum)


con<-file("Results/N_lim_reactor/Inoculum_comparison_with_CuoioDep/Venn_Diagr_exclusive_bacteria_of_inoculum.txt")
sink(con, append=TRUE)
cat("\n\nONLY IN Inoculum", fill=TRUE)
cat(ONLY_IN_inoculum)
cat("\n\nONLY IN CuoioD", fill=TRUE)
cat(ONLY_IN_CuoioD)
sink()
close(con)

x<-list(Inoculum_N_lim_Reactor=Inoculum,CuoioDepur=CuoioD)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("lightblue","deepskyblue")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) )
ggsave(filename = "Results/N_lim_reactor/Inoculum_comparison_with_CuoioDep/Venn_Diagr_exclusive_bacteria_of_inoculum.png", width = 4, height = 4.2, dpi=300, bg = "white")
dev.off()

suppressWarnings(rm(ONLY_IN_CuoioD, ONLY_IN_inoculum, ONLY_IN_Malignant_tumor, 
                    x, con, CuoioD, Inoculum, Malignant_tumor, data.venn, who))




####################### GENERAL PCoA BEFORE SUBSETTING (Only N and P Reactors) ##########################

suppressWarnings(rm(data.prop_temp, data.sqrt_prop))

if( ! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}
data.temp<-subset_samples(data_complete, Experiment_state!="Adherent")
data.temp<-subset_samples(data.temp, Reactor_Type!="CuoioD_Aero")
data.prop_temp<-transform_sample_counts(data.temp, function(x) (x/sum(x))*100)
data.prop.labels<-data.prop_temp

sample_data(data.prop.labels)$Reactor_Type<-gsub("N_limitation","N limiting",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Reactor_Type<-gsub("P_limitation","P limiting",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Sample_Type<-gsub("Inoculum","R1",sample_data(data.prop.labels)$Sample_Type) # with "inoculum" here is meant the "original" reactor of the first part of the experiment (it's the same R1)
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day))


{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
colors<-c(
  "N limiting"="lightblue3",
  "P limiting"="coral"
)
sample_names(data.sqrt_prop)<-as.factor(sample_names(data.sqrt_prop))

# Type 1 (everything)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type", shape = "Sample_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.2) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.12) +
  stat_ellipse(size=0.15) + 
  geom_text(aes(label= Experiment_day), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(
       color="",
       shape="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_Ellypses_Names_Lines.png",width = 5.5, height = 4.2, dpi=300)
# Type 2 (no ellypses)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type", shape = "Sample_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.6) +
  geom_point(size=5, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.14) +
  geom_text(aes(label= Experiment_day), 
            color="black", size=1.6,
            alpha=0.9,
            show.legend = FALSE) +
  labs(
       color="",
       shape="",
       x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation")
       )
ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_Names_Lines.png", width = 5.5, height = 4.2, dpi=300)
# # As Type2 but with shades of colors
# plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2", shape = "Sample_Type") +
#   #scale_color_gradient2(low="lightblue", mid = "coral", high="red3", lim=c(0,409), midpoint = 340) +
#   scale_colour_multi(colours = )
#   # scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(0,270), midpoint = 150) +
#   geom_point(size=2.6) +
#   geom_point(size=5, alpha= 0.5) +
#   theme_classic(base_size = 9) +
#   geom_path(aes(group=Sample_Type), color="gray", size= 0.14) +
#   geom_text(aes(label= Experiment_day), 
#             color="black", size=1.6,
#             alpha=0.9,
#             show.legend = FALSE) +
#   labs(
#     color="",
#     shape="",
#     x=paste("PC1: ",eigval[1],"% variation"),
#     y=paste("PC2: ",eigval[2],"% variation")
#   )
# ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_Shades_of_colors.png", width = 5.5, height = 4.2, dpi=300)
# Type 3 (no ellypses, no names)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=1.8) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.12) +
  labs(title="",
       color="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_NO_Names.png", width = 5.5, height = 4.2, dpi=300)
# Type 4 (no ellypses, no lines)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type", shape = "Sample_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=3) +
  geom_point(size=4.8, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_text(aes(label= Experiment_day), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(
    color="",
    shape="",
    x=paste("PC1: ",eigval[1],"% variation"),
    y=paste("PC2: ",eigval[2],"% variation")
  )
ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_NO_Lines_NO_Ellypses.png", width = 5.5, height = 4.2, dpi=300)
# Type 5 (no lines)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type", shape = "Sample_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.2) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  stat_ellipse(size=0.15) + 
  geom_text(aes(label= Experiment_day), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(
    color="",
    shape="",
    x=paste("PC1: ",eigval[1],"% variation"),
    y=paste("PC2: ",eigval[2],"% variation")
  )
ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_NO_Lines.png", width = 5.5, height = 4.2, dpi=300)
# Type 6 (ONLY Ellypses)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type", shape = "Sample_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.5) +
  geom_point(size=4.5, alpha= 0.5) +
  theme_classic(base_size = 9) +
  stat_ellipse(size=0.15) + 
  labs(title="",
       color="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_with_Nlim_and_Plim/GENERAL_PCoA_Hellinger_on_genera_ONLY_Ellypses.png", width = 5.5, height = 4.2, dpi=300)



#################### ABUNDANCES BAR PLOTS (Only N and P Reactors _ including R2) ####################

### TOP 5 Phyla
data.target<-subset_samples(data_complete, Reactor_Type!="CuoioD_Aero")
data.target<-subset_samples(data.target, Experiment_state!="Adherent")
data.target<-tax_glom(data.target, taxrank = "Phylum")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.target)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("N_limitation","N limiting R1",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("P_limitation","P limiting",tabella$Reactor_Type)
tabella$Reactor_Type[tabella$Sample_Type=="R2"]<-"P limiting R2"
tabella$Reactor_Type[tabella$Sample_Type=="R1"]<-"P limiting R1"


ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.6), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12.5 ),
        legend.position="bottom",
        legend.margin = margin(-2,0,5,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment day", y="Percent abundance", 
       title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/General_data_analysis/TOP_5_phyla_Nlim_Plim_with_R2.png",width=8.2,height=5, dpi=300)
dev.off()

write.xlsx(file = "Results/General_data_analysis/TOP_5_phyla_ONLY_Nlim_and_Plim_Averages_with_R2.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)



### TOP 10 Genera
data.target<-subset_samples(data_complete, Reactor_Type!="CuoioD_Aero")
data.target<-subset_samples(data.target, Experiment_state!="Adherent")
data.target<-tax_glom(data.target, taxrank = "Genus")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

# adding informations to missing names
taxa_temp<-as.data.frame(tax_table(data.target))
for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])
}
tax_table(data.target)<-as.matrix(taxa_temp)

{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("N_limitation","N limiting R1",tabella$Reactor_Type)
tabella$Reactor_Type[tabella$Sample_Type=="R2"]<-"P limiting R2"
tabella$Reactor_Type[tabella$Sample_Type=="R1"]<-"P limiting R1"

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.6), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.5 ),
        legend.position="bottom", 
        legend.margin = margin(-5,22,0,0)
  ) + 
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="Experiment day", y="Percent abundance",
       title = "Ten most abundant genera",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/General_data_analysis/TOP_10_Genera_Nlim_Plim_with_R2.png",width=8.2,height=5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/General_data_analysis/TOP_10_genera_ONLY_Nlim_and_Plim_Averages_with_R2.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target))



#################### ABUNDANCES BAR PLOTS (Only N and P Reactors _ without R2) ####################

### TOP 5 Phyla
data.target<-subset_samples(data_complete, Reactor_Type!="CuoioD_Aero")
data.target<-subset_samples(data.target, Sample_Type!="R2")
data.target<-subset_samples(data.target, Experiment_state!="Adherent")
data.target<-tax_glom(data.target, taxrank = "Phylum")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.target)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("P_limitation","P limiting",tabella$Reactor_Type)

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + 
  theme_classic(base_size =12.8) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust = 1, size= 8.6), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12 ),
        plot.title = element_text ( size = 11.5 ),
        legend.position="bottom",
        #legend.margin = margin(0,0,0,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment day", y="Percent abundance", 
       title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/General_data_analysis/TOP_5_phyla_Nlim_Plim_no_R2.png",width=8,height=5, dpi=300)
dev.off()

write.xlsx(file = "Results/General_data_analysis/TOP_10_phyla_ONLY_Nlim_and_Plim_Averages_no_R2.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)



### TOP 10 Genera
data.target<-subset_samples(data_complete, Reactor_Type!="CuoioD_Aero")
data.target<-subset_samples(data.target, Experiment_state!="Adherent")
data.target<-subset_samples(data.target, Sample_Type!="R2")
data.target<-tax_glom(data.target, taxrank = "Genus")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

# adding informations to missing names
taxa_temp<-as.data.frame(tax_table(data.target))
for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
  taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])
}
tax_table(data.target)<-as.matrix(taxa_temp)

{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)
tabella$Reactor_Type[tabella$Sample_Type=="R1"]<-"P limiting"

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12.8) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=40, vjust=1, hjust = 1, size= 8.6), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.8 ),
        plot.title = element_text ( size = 11.5 ),
        legend.position="bottom", 
        legend.margin = margin(0,22,0,0)
  ) + 
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="Experiment day", y="Percent abundance",
       title = "Ten most abundant genera",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/General_data_analysis/TOP_10_Genera_Nlim_Plim_no_R2.png",width=8,height=5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/General_data_analysis/TOP_10_genera_ONLY_Nlim_and_Plim_Averages_no_R2.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target))



################# STUDING THE ABUNDANCES OF THE PHA ACCUMULATORS (N lim AND P lim WITHOUT R2) ######################

lista_PHA<-table_PHA[,1]     # they derive from papers

data.genus.temp<-subset_samples(data_complete, Reactor_Type!="CuoioD_Aero")
data.genus.temp<-subset_samples(data.genus.temp, Sample_Type!="R2")
data.genus.temp<-subset_samples(data.genus.temp, Experiment_state!="Adherent")
data.genus.temp<-tax_glom(data.genus.temp, taxrank ="Genus", NArm=F)
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-gsub("Candidatus_","",tabella$Genus)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))

# fill_color_20<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","firebrick3","springgreen4","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","darkgreen","lightgrey","coral","deepskyblue2") 
fill_color_19<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","darkgreen","firebrick3","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","springgreen3","coral","deepskyblue2") 

tabella$Sample_Type<-gsub("Inoculum","N limiting", tabella$Sample_Type)
tabella$Sample_Type<-gsub("R1","P limiting", tabella$Sample_Type)
# tabella$Sample_Type<-gsub("R2","P limiting (R2)", tabella$Sample_Type)

ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=Genus)) + 
  facet_grid2( ~Sample_Type, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.85) +
  theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_19)) +
  theme(axis.text.x=element_text(angle=40, hjust=1,vjust=1, size=9),
        axis.title.y = element_text(size=11), 
        axis.text.y = element_text(size=9.5), 
        strip.text.x = element_text(size=12.5), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.8 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,2,2,2),
        legend.margin =  margin(-2,40,0,0)
  ) +
  guides(fill=guide_legend(nrow=4)) + 
  scale_y_continuous(lim=c(0,100)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = paste("PHA accumulating genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"))
ggsave(file="Results/General_data_analysis/PHA_Accumulators_in_both_N_and_Plim_without_R2.png",width=7,height=5.8,dpi=300)
dev.off()



################### ALPHA DIVERSITIES VALUES (N lim vs both P lim reactors) #####################

# selecting only the last 4 sample for each reactor (the one that SHOULD be at least near the steady state)
data.genus.temp<-subset_samples(data_complete, Sample_name %in% c("P1","P2","P9","P10","P11","P12",
                                                                  "S19","S20","S21","S22",
                                                                  "F19","F20","F21","F22"))
data.genus.temp<-tax_glom(data.genus.temp, taxrank ="Genus", NArm=F)
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_Type", color="Reactor_Type")
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
pAlpha$data$Reactor_Type<-gsub("N_limitation","N limiting R1",pAlpha$data$Reactor_Type)
pAlpha$data$Reactor_Type[pAlpha$data$Sample_Type=="R2"]<-"P limiting R2"
pAlpha$data$Reactor_Type[pAlpha$data$Sample_Type=="R1"]<-"P limiting R1"

pAlpha$data$Reactor_Type2<-pAlpha$data$Reactor_Type
pAlpha$data$Reactor_Type2[pAlpha$data$Sample_name=="P1"]<-"Inoculum"
pAlpha$data$Reactor_Type2[pAlpha$data$Sample_name=="P2"]<-"Inoculum"

colors4<-c("Inoculum"="deepskyblue",
           "N limiting R1"="lightblue3",
           "P limiting R1"="coral",
           "P limiting R2"="red"
)

# with ID
pAlpha + theme_classic2() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL, group= Reactor_Type2),
               position = position_dodge(width = 0),
               size= 0.3, alpha= 0) +
  geom_point(aes(color=Reactor_Type2), size=3.1, alpha= 0.4) +
  geom_text(aes(label=Experiment_day), size= 2.5, color= "black") +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "kruskal.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top",
                     vjust=-0.4, hjust=0) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_diversity_Nlim_vs_BOTH_Plim_reactors_with_Kruskal_and_sampleID.png", width = 6.8,height =5.2, dpi=300)
#### again, without the inocula
pAlpha$data<-pAlpha$data[pAlpha$data$Reactor_Type2!="Inoculum", ]
pAlpha +
  theme_classic() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL, group= Reactor_Type),
               position = position_dodge(width = 0),
               size= 0.3, alpha= 0) +
  geom_point(aes(color=Reactor_Type), size=3.1, alpha= 0.4) +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "kruskal.test", paired = T,
                     label.x= 1.5, size=3.4, label.y.npc = "top",
                     vjust=-0.4, hjust=0) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_diversity_Nlim_vs_BOTH_Plim_reactors_with_Kruskal_WITHOUT_INOCULA.png", width = 6.8,height =5.2, dpi=300)


######## pairwise comparisons
alphadt<-pAlpha$data
Filter_value_obs<-alphadt[alphadt$variable=="Observed richness", ]
Filter_value_sha<-alphadt[alphadt$variable=="Shannon", ]
Filter_value_eve<-alphadt [alphadt$variable=="Evenness", ]

# with Dunnett (one factor only)
con <- file("Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_diversity_Nlim_and_BOTH_Plim_reactors_pairwise_comparisons.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on sub groups Alpha diversity  \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_obs, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_sha, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_eve, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
sink()
close(con)


################### ALPHA DIVERSITIES VALUES (N lim vs P lim , only R1) #####################

# selecting only the last 4 sample for each reactor (the one that SHOULD be at least near the steady state)
data.genus.temp<-subset_samples(data_complete, Sample_name %in% c("P1","P2","P9","P10","P11","P12",
                                                                  "S19","S20","S21","S22"
))
data.genus.temp<-tax_glom(data.genus.temp, taxrank ="Genus", NArm=F)
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_Type", color="Reactor_Type")
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
pAlpha$data$Reactor_Type<-gsub("N_limitation","N limiting",pAlpha$data$Reactor_Type)
pAlpha$data$Reactor_Type[pAlpha$data$Sample_Type=="R1"]<-"P limiting"

pAlpha$data$Reactor_Type[pAlpha$data$Sample_name=="P1"]<-"Inoculum"
pAlpha$data$Reactor_Type[pAlpha$data$Sample_name=="P2"]<-"Inoculum"

colors4<-c("Inoculum"="deepskyblue",
           "N limiting"="lightblue3",
           "P limiting"="coral"
)

# with points only
pAlpha +
  theme_classic() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL, group= Reactor_Type),
               position = position_dodge(width = 0),
               size= 0.3, alpha= 0) +
  geom_point(size=3.1, alpha= 0.4) +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "kruskal.test", paired = F,
                     label.x= 1.5, size=3.5, label.y.npc = "top",
                     vjust=-0.4, hjust=0) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_div_Nlim_vs_Plim_onlyR1_with_Kruskal.png", width = 6.8,height =5.2, dpi=300)
# with ID
pAlpha + theme_classic2() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL, group= Reactor_Type),
               position = position_dodge(width = 0),
               size= 0.3, alpha= 0) +
  geom_point(size=3.1, alpha= 0.4) +
  geom_text(aes(label=gsub(" ","",Experiment_day)),
            size= 2.5, color= "black") +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "kruskal.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top",
                     vjust=-0.4, hjust=0) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_div_Nlim_vs_Plim_onlyR1_with_Kruskal_and_EXPERIM_DAY.png", width = 6.8,height =5.2, dpi=300)


######## pairwise comparisons
alphadt<-pAlpha$data
Filter_value_obs<-alphadt[alphadt$variable=="Observed richness", ]
Filter_value_sha<-alphadt[alphadt$variable=="Shannon", ]
Filter_value_eve<-alphadt [alphadt$variable=="Evenness", ]

# with Dunnett (one factor only)
con <- file("Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_div_Nlim_and_Plim_ONLYR1_pairwise_comparisons.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on sub groups Alpha diversity  \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_obs, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_sha, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_eve, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
sink()
close(con)


#### again, but without the inocula
pAlpha$data<-pAlpha$data[pAlpha$data$Reactor_Type!="Inoculum", ]
pAlpha +
  theme_classic() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL, group= Reactor_Type),
               position = position_dodge(width = 0),
               size= 0.3, alpha= 0) +
  geom_point(size=3.1, alpha= 0.4) +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "wilcox.test", paired = F,
                     label.x= 1.2, size=3.5, label.y.npc = "top",
                     vjust=-0.4, hjust=0) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim/Alpha_div_Nlim_vs_Plim_onlyR1_with_WILCOX_test_WITHOUT_INOCULA.png", width = 6.8,height =5.2, dpi=300)



#################### \\\\\\ STARTING THE ANALYSIS OF N lim's SUBSET \\\\\ #######################

data<-subset_samples(data_complete, Reactor_Type=="N_limitation")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa

### everything ready, let's start with the analysis
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
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

metadata<-metadata_complete[metadata_complete$Reactor_Type=="N_limitation",  ]



###################### ABUNDANCES BAR PLOT _ N_lim Only ##########################

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
  labs(y="Experiment day", x="Percent abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/N_lim_reactor/Abundances/TOP5_phyla_abundances.png",width=7,height=6.5, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/N_lim_reactor/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
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
ggplot(data=tabella, aes(x=Abundance, y=Experiment_day, fill=Genus)) + theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_10) +
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom",
        legend.margin = margin(1,-0.5,0,0)) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(y="Experiment day", x="Percent abundance", title = "Ten most abundant genera", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/N_lim_reactor/Abundances/TOP10_genera_abundances.png",width=7,height=6.5, dpi=300) 
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/N_lim_reactor/Abundances/TOP_10_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



######################## PCoA BETA DIV _ N limiting #########################

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
plot_ordination(data.sqrt_prop, ordBC , color = "Sample_Type") +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  scale_color_manual(values = c("lightblue3")) +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none") +
  theme_classic(base_size = 12.5) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day), color="black", size=2.8, show.legend = FALSE) +
  labs(title="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/N_lim_reactor/Beta_divers_Hellinger_on_genera_WITH_NAMES.png", width = 6.8, height = 5, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Sample_Type") +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  scale_color_manual(values = c("lightblue3")) +
  geom_point(size=3.5, alpha=0.5, show.legend = F) +
  guides(color="none") +
  theme_classic(base_size = 12.5) + 
  labs(title="", color="", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")
  )
# with shade of color along the time
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,30,90,150,210,262),
                          midpoint = 90) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  labs(title=" ",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/N_lim_reactor/Beta_divers_Hellinger_on_genera_shade.png", width = 4.2, height = 3.8, dpi=300)
# with shade of color AND also text
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,30,90,150,210,262),
                        midpoint = 90) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(title=" ",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/N_lim_reactor/Beta_divers_Hellinger_on_genera_shade_and_ExpDayText.png", width = 4.2, height = 3.8, dpi=300)



################# STUDING THE ABUNDANCES OF THE PHA ACCUMULATORS (only N limiting) ######################

fill_color_19<-c("wheat3","deeppink","darkmagenta","coral","cyan","yellow2","brown","firebrick3","springgreen2","violet","bisque","deepskyblue2","darkslategray3","blue", "gray", "yellow4","red","darkgreen","lightgrey") # "others" will be setted as the last one

lista_PHA<-table_PHA[,1]     # they derive from papers

prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-gsub("Candidatus_","",tabella$Genus) # synonyms
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))

# tabella2<-tabella[tabella$Experiment_state!="Adherent", ]
fill_color_19<-c("wheat3","deeppink","darkmagenta","coral","cyan","yellow2","brown","firebrick3","deepskyblue2","springgreen2","violet","bisque","darkslategray3","blue", "gray", "yellow4","red","darkgreen","lightgrey") # "others" will be setted as the last one

ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.85) +
  theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_19)) +
  theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1, size=8.4),
        axis.title.y = element_text(size=10), 
        axis.text.y = element_text(size=8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom",
        plot.margin = margin(2,2,2,2),
        legend.margin =  margin(0,14,0,0)
  ) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Experiment_day", y="Percent abundance", 
       fill="",
       title = paste("PHA related genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"))
ggsave(file="Results/N_lim_reactor/PHA_Accumulators/PHA_accumulators_along_the_time_N_lim.png",width=8,height=5.8,dpi=300)
dev.off()



########## CORRELATIONS ABUNDANCES vs TIME _ N limiting (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("P2","P7") ) # removing the technical replicates
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp <- subset_samples(data.genus, ! Sample_name %in% c("P2","P7") ) # removing the technical replicate
data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sample_name
order_of_samples<-metadata[levels(sample_data(data)$Sample_name), "Sample_name"]   # the sample data has been ordered already during the preparation
order_of_samples <- order_of_samples[! order_of_samples %in% c("P2","P7") ]
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
Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")

write.csv2(Corr_results, "Results/N_lim_reactor/Time_series_in_N_lim/Correlations_of_CLR_bacteria_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1<-row.names(temp_for_save)

con <- file("Results/N_lim_reactor/Time_series_in_N_lim/LOG_RATIO_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the N_lim samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

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
tabella<-tabella[order(tabella$Sample_name) , ] # days
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 0.5)) +
  labs(x="Experiment day", y="Percent abundance", title = "the 6 genera most positely correlated with time (N limiting)\n according to spearman correlation on CLR counts")
ggsave(file="Results/N_lim_reactor/Time_series_in_N_lim/6_Genera_CLR_most_positively_correlated_with_time.png",width=7,height=5, dpi=300) 
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
tabella<-tabella[order(tabella$Sample_name) , ] # days
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 0.5)) +
  labs(x="Experiment day", y="Percent abundance", title = "the 6 genera most negatively correlated with time (N limiting)\n according to spearman correlation on CLR counts")
ggsave(file="Results/N_lim_reactor/Time_series_in_N_lim/6_Genera_CLR_most_positively_correlated_with_time.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_6, prune.dat_top, tax_selected, tabella, corr_6_top))


