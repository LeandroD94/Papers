# Comparison gly vs no gly between frozen samples


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




########################### FOCUS ON THE *PN-AGS2* SAMPLE WITH FRIDGE-FROZEN Gly vs noGly ##########################
# 
# ### TOP Genera
# suppressWarnings(rm(top, others, tabella, unass_data))
# { top_data <- subset_samples(data.genus.prop, Reactor=="PN_AGS2" )
#   top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:29]
#   prune.dat_top <- prune_taxa(top,top_data)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
#     tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   } else {
#     stop("Something went wrong when updating the taxonomy")
#   }
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(top_data)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,top_data)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
#   tabella$Genus<-gsub ("Flavobacteriaceae","Flavobacter.", tabella$Genus)
#   tabella$Genus<-gsub ("Chitinophagales","Chitinophag.", tabella$Genus)
#   tabella$Genus<-gsub ("Gemmatimonadaceae","Gemmatimonad.", tabella$Genus)
#   tabella$Genus<-gsub ("Sphingomonadaceae","Sphingomonad.", tabella$Genus)
#   tabella$Genus<-gsub ("Microbacteriaceae","Microbacteriac.", tabella$Genus)
#   tabella$Genus<-gsub ("f Xanthobacteraceae","Xanthobacter.", tabella$Genus)
#   tabella$Genus<-gsub ("marinales","marin.", tabella$Genus)
#   tabella$Genus<-gsub ("Microtrichales","Microtrich.", tabella$Genus)
#   tabella$Genus<-gsub ("uncultured_ ","uncult. ", tabella$Genus)
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
#   levels(tabella$Storage) <- gsub("C","°C (Gly)", levels(tabella$Storage) )
#   levels(tabella$Storage) <- gsub("4°C (Gly)","4°C", levels(tabella$Storage) )
#   levels(tabella$Storage) <- gsub("(Gly)_noGly","(no Gly)", levels(tabella$Storage) )
#   # levels(tabella$Storage)<-gsub ("_"," ", levels(tabella$Storage))
# }
# fill_color_30_BIS<-c("chartreuse","deeppink", "firebrick3",
#                      "darkcyan","darkmagenta", 
#                      "darkblue","grey55", "aquamarine3",
#                      "bisque","cyan","yellow3","brown4",
#                      "springgreen4", "grey20", 
#                      "deepskyblue4","lightgreen","orange2",
#                      "darkorchid", "lightblue1" , "violet", "blue",
#                      "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
#                      "lightgrey", "red","darkslategray3")
# ggplot(data=tabella, aes(x=Time_before_extr2, y=Abundance, fill=Genus)) +
#   geom_bar( stat="identity", position="stack", na.rm = F) + 
#   scale_fill_manual(values=fill_color_30_BIS) +
#   facet_grid( . ~ Storage ,
#               scales = "free_x", space="free"
#               #strip = strip_nested(size = "variable")
#   ) +
#   theme_classic(base_size =11) +
#   scale_y_continuous(expand=c(0,1) ) +
#   theme(axis.text.x=element_text(angle=30,
#                                  vjust=1.25,
#                                  hjust=1,
#                                  size= 7.5
#   ),
#   axis.text.y=element_text(size=7),
#   axis.ticks.x = element_blank(),
#   strip.text.x = element_text(size=8.75, angle=0),
#   legend.key.height = unit(0.14, "cm"),
#   legend.key.width = unit(0.2, "cm"),
#   legend.text = element_text ( size = 7.8 ),
#   legend.position="bottom",
#   legend.margin = margin(-14 ,33, 1 ,0),
#   plot.margin = margin(2,1,2,1),
#   panel.spacing.x = unit(0.05,"cm"),
#   ) +
#   guides(fill=guide_legend(nrow=8)) +
#   labs(x="", y="Percentual abundance of clades",
#        fill="")
# ggsave(file="Results/Abundances/PNAGS2_TOP_genera_FRIDGEvsFROZEN.png",width=5.2,height=4.25, dpi=300)
# 
# dev.off()
# 
# suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))
# 
# 
# 
# #### PCoA on genera
# data.prop.labels<-subset_samples(data.genus.prop, Reactor=="PN_AGS2" )
# sample_data(data.prop.labels)$Glycerol <- gsub(" ","Direct", sample_data(data.prop.labels)$Glycerol)
# levels(sample_data(data.prop.labels)$Storage) <- gsub("_noGly","", levels(sample_data(data.prop.labels)$Storage) )
# Time<- gsub("m","", as.character(sample_data(data.prop.labels)$Time_before_extr2))
# Time_black <- Time
# Time_black [ as.character(sample_data(data.prop.labels)$Storage)=="-80C" ] <- ""
# Time_white <- Time
# Time_white [ as.character(sample_data(data.prop.labels)$Storage)!="-80C" ] <- ""
# {data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# values_color<-c("Direct"="coral",
#                 "4C" = "yellow2",
#                 "-20C" = "deepskyblue3",
#                 #"-20C_noGly" = "deepskyblue",
#                 "-80C" = "darkblue",
#                 #"-80C_noGly" = "royalblue3",
#                 "No_storage"="red3")
# plot_ordination(data.sqrt_prop, ordBC, color="Storage" , shape = "Glycerol" )+
#   geom_point(size=4.5, alpha=0.4,  show.legend = F) +
#   theme_classic(base_size = 11) +
#   theme(legend.text =element_text(size=10.8)) +
#   geom_text(aes(label=Time_black),
#             color="black",
#             size=4.4, hjust=0.42,show.legend = FALSE) +
#   geom_text(aes(label=Time_white),
#             color="gray95",
#             size=4.4, hjust=0.52,show.legend = FALSE) +
#   labs(
#     x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#     color="", shape="") +
#   scale_color_manual(values= values_color)
# ggsave(file="Results/Beta_diversity_PCoA/PNAGS2_FridgeFrozen_PCoA_Hellinger_Genus.png", width = 4.2, height = 3.85, dpi=300)




###################### ABUNDANCES GLY vs no GLY ######################

### TOP GENERA
subset_target <- subset_samples(data.genus.prop, Storage %in% c("Direct","-20C_noGly","-20C") )
test_AS_samples <-c("A25","A38","A39" )  # selected the older -20 samples
selected_from_PN2 <- c("PGB1","PGB10","PGB9")  # as above
subset_target <- subset_samples(subset_target , Reactor=="PN_AGS" | Sample_name%in% c(test_AS_samples , selected_from_PN2) ) # only paired samples!

suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_target
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
  tabella$Genus<-gsub ("Comamonadaceae","Comamonad.", tabella$Genus)
  tabella$Genus<-gsub ("Xanthobacteraceae","Xanthobacter.", tabella$Genus)
  tabella$Genus<-gsub ("Sphingomonadaceae","Sphingomon.", tabella$Genus)
  tabella$Genus<-gsub ("Microbacteriaceae","Microbacter.", tabella$Genus)
  tabella$Genus<-gsub ("Microtrichales","Microtrich.", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
levels(tabella$Glycerol)<-gsub (" ","Direct", levels(tabella$Glycerol))
levels(tabella$Glycerol) <- gsub("noGly"," no Gly", levels(tabella$Glycerol) )

tabella$Reactor_ID <- gsub("A .*","AS tank\n(Main A)", tabella$Reactor_ID )
tabella$Reactor_ID <- gsub("PN_AGS_4","PN-AGS 2\n(Main B)", tabella$Reactor_ID )
tabella$Reactor_ID <- gsub("PN_AGS_3","PN-AGS 1\ntime point 3", tabella$Reactor_ID )
tabella$Reactor_ID <- gsub("PN_AGS_2","PN-AGS 1\ntime point 2", tabella$Reactor_ID )
tabella$Reactor_ID <- gsub("PN_AGS_1","PN-AGS 1\ntime point 1", tabella$Reactor_ID )

fill_color_30_BIS<-c("springgreen3","deeppink", "firebrick3",
                     "cyan", "darkmagenta", 
                     "grey55", "aquamarine3", "lightgreen",
                     "bisque","yellow3","brown4",
                     "green2", "grey20", 
                     "deepskyblue4","orange2", "lightblue4",
                     "darkorchid", "lightblue1" ,"violet", "royalblue2",
                     "yellow", "pink3","yellow4", "chocolate4","coral2","darkblue","darkgreen",
                     "lightgrey", "red","darkslategray3")
ggplot(data=tabella, aes(x=Glycerol, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=fill_color_30_BIS) +
  facet_nested( . ~ Reactor_ID ,
                scales = "free_x", space="free",
                strip = strip_nested(size = "variable")
  ) +
  theme_classic(base_size =10) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=20,
                                 vjust=1,
                                 hjust=1,
                                 size= 7.25
                                 ),
  axis.text.y=element_text(size=7),
  # axis.title.y = element_text(size=10),
  # axis.title.x = element_text(vjust=3),
  strip.text.x = element_text(size=7.8, angle=0),
  legend.key.height = unit(0.14, "cm"),
  legend.key.width = unit(0.245, "cm"),
  legend.text = element_text ( size = 8.65 , face="italic" ),
  legend.position="bottom",
  legend.margin = margin(-15 ,34, 1 ,0),
  plot.margin = margin(2,1,2,1),
  panel.spacing.x = unit(0.05,"cm"),
  ) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance of clades",
       fill="")
ggsave(file="Results/Abundances/TOP_genera_SAMPLE_GlyVSnoGly.png",width=4.9,height=4.25, dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




#################### PCoA only GLY vs no GLY ################

# subset_target <- subset_samples(data.genus.prop, Storage %in% c("Direct","-20C_noGly","-20C") )
subset_target <- subset_samples(data.genus.prop, Storage %in% c("-20C_noGly","-20C") )
test_AS_samples <-c("A25","A38","A39" )  # selected the older -20 samples
selected_from_PN2 <- c("PGB1","PGB10","PGB9")  # as above
subset_target <- subset_samples(subset_target , Reactor=="PN_AGS" | Sample_name%in% c(test_AS_samples , selected_from_PN2) ) # only paired samples!

data.prop.labels<-subset_target
sample_data(data.prop.labels)$Reactor <- gsub("CuoioDepur_Aero","AS tank", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Reactor <- gsub("Pilot_AGS","Pilot AGS", sample_data(data.prop.labels)$Reactor )
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "-20C" = "deepskyblue2",
                "-20C_noGly"="cyan")
plot_ordination(data.sqrt_prop, ordBC, color="Storage" , shape="Reactor")+
  geom_point(size=4, alpha=0.5,  show.legend = F) +
  theme_classic(base_size = 10.5) +
  geom_path(aes(group=Reactor_ID), color="darkgrey", linewidth = 0.25, alpha=0.65,
            # arrow=arrow(length =unit(0.22,"cm"), 
            #             type = "closed")
  ) +
  theme(legend.text =element_text(size=10.8)) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="", shape="") +
  scale_color_manual(values= values_color)
ggsave(file="Results/Beta_diversity_PCoA/PCoA_Hellinger_Genus_GlyVSnoGly.png", width = 4.5, height = 4.2, dpi=300)



##### ONLY SIMILAR GLY (AGSs)

# subset_target <- subset_samples(data.genus.prop, Storage %in% c("Direct","-20C_noGly","-20C") )
subset_target <- subset_samples(data.genus.prop, Storage %in% c("-20C_noGly","-20C") )
selected_from_PN2 <- c("PGB1","PGB10","PGB9")  # as above
subset_target <- subset_samples(subset_target , Reactor=="PN_AGS" | Sample_name%in% c(selected_from_PN2) ) # only paired samples!

data.prop.labels<-subset_target
levels(sample_data(data.prop.labels)$Sample_code) <- gsub("PG4.*","", levels(sample_data(data.prop.labels)$Sample_code) )
levels(sample_data(data.prop.labels)$Sample_code) <- gsub("_gly","", levels(sample_data(data.prop.labels)$Sample_code) )
levels(sample_data(data.prop.labels)$Sample_code) <- gsub("PG","", levels(sample_data(data.prop.labels)$Sample_code) )
levels(sample_data(data.prop.labels)$Storage) <- gsub("C"," °C", levels(sample_data(data.prop.labels)$Storage) )
levels(sample_data(data.prop.labels)$Storage) <- gsub("_noGly"," (no Gly)", levels(sample_data(data.prop.labels)$Storage) )
sample_data(data.prop.labels)$Reactor <- gsub("PN_AGS2","PN-AGS 2\n(Main B)", sample_data(data.prop.labels)$Reactor )
sample_data(data.prop.labels)$Reactor <- gsub("PN_AGS","PN-AGS 1", sample_data(data.prop.labels)$Reactor )
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
values_color<-c("Direct"="coral",
                "-20 °C" = "deepskyblue2",
                "-20 °C (no Gly)"="cyan")
plot_ordination(data.sqrt_prop, ordBC, color="Storage" , shape="Reactor")+
  geom_point(size=4, alpha=0.5,  show.legend = F) +
  theme_classic(base_size = 10.5) +
  geom_path(aes(group=Reactor_ID), color="darkgrey", linewidth = 0.25, alpha=0.65,
            # arrow=arrow(length =unit(0.22,"cm"), 
            #             type = "closed")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_code), color="black", size=2.5, show.legend = FALSE) +
  theme(legend.text =element_text(size=10.8)) +
  labs(
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
    color="", shape="") +
  scale_color_manual(values= values_color)
ggsave(file="Results/Beta_diversity_PCoA/PCoA_Hellinger_Genus_GlyVSnoGly_onlyAGSs.png", width = 4.5, height = 3.85, dpi=300)




############# DA WITH DESEQ2 (Gly vs noGly) ##################

suppressWarnings(rm(data_pruned, data.genus_pruned, subset_target))

subset_target <- subset_samples(data.genus.prop, Storage %in% c("-20C_noGly","-20C") )
test_AS_samples <-c("A38","A39" )  # selected the older -20 samples
selected_from_PN2 <- c("PGB10","PGB9")  # as above
subset_target <- subset_samples(subset_target , Reactor=="PN_AGS" | Sample_name%in% c(test_AS_samples , selected_from_PN2) ) # only paired samples!
data_pruned<- prune_taxa(taxa_sums(subset_target) > 10, subset_target) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

sample_data(data_pruned)$Reactor_ID  <- gsub("/", "_", sample_data(data_pruned)$Reactor_ID )
sample_data(data_pruned)$Reactor_ID  <- gsub(" ", "_", sample_data(data_pruned)$Reactor_ID )

Table_tot<-NULL
Res_tot<-NULL
for( t in c("Genus","Family") ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Reactor_ID + Glycerol)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Glycerol", "Gly", "noGly"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  #res<-res[res$baseMean > 0, ] # the environments are too much diverse, better to not use this filter here...
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
# write.csv2("No differences", file="Results/Differently_abundant_bact/Differently_abundant_GlyVSnoGLY.csv", row.names = F)



############ PLOTTING THESE RESULTS ...

Table_tot$Glycerol<-factor(Table_tot$Glycerol, levels = unique(Table_tot$Glycerol))
Table_tot<-Table_tot[order(match(Table_tot$Glycerol, levels(Table_tot$Glycerol))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Glycerol)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Glycerol)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

## to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Glycerol)

tabella_tot <- rbind.data.frame(tabella_g,tabella_2)
tabella_tot$Reactor_ID2<-gsub("PN_AGS_","PN-AGS", tabella_tot$Reactor_ID)
tabella_tot$Reactor_ID2<-gsub("AGS[1-3]","AGS", tabella_tot$Reactor_ID2)
tabella_tot$Reactor_ID2<-gsub("AGS4","AGS2", tabella_tot$Reactor_ID2)   # the other (2nd) reactor from the other exp
tabella_tot$Reactor_ID2<-gsub("A_11_23","AS11/23", tabella_tot$Reactor_ID2)

plot_1<-ggplot(tabella_tot, aes(x= Xaxis, y=Abundance, fill=Glycerol, shape=Glycerol)) +
  theme_classic2(base_size = 12) +
  scale_shape_manual(values= c(15,17))+
  scale_fill_manual(values = c("Gly"="deepskyblue","noGly"="deepskyblue3")) +
  scale_color_manual(values = c("Gly"="deepskyblue","noGly"="deepskyblue3")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Genus","Family"))~Bacteria,
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Glycerol), size=1.8, alpha=1) +
  geom_point(aes(color=Glycerol), size=3, alpha=0.5) +
  geom_text(aes(label=Reactor_ID2), size=2.6, vjust= -0.3,
            color="gray40", show.legend = F) +
  geom_line(aes(group=Reactor_ID), linewidth=0.2) +
  scale_x_discrete(labels=unique(levels(tabella_g$Glycerol)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=8.4,colour="black"),
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=13),
        axis.text.y = element_text(size=8),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt"),
        panel.grid.minor.y= element_blank()
  ) +
  guides(fill="none",color="none") +
  labs(y="Percentual abundance", fill="", color="", x="", shape="") +
  theme(legend.margin=margin(-15, 0, 0, 0), # overwriting
        legend.position = "bottom")
plot_1
ggsave(filename = "Results/Differently_abundant_bact/Gly_noGly_20C_DeSeq2.png", width = 6.8, height = 4.65, dpi=300)
dev.off()
