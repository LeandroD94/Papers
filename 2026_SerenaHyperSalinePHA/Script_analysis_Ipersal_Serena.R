############## PREPARING THE ENVIRONMENT ##############

{ 
  library("phyloseq")
  library("qiime2R")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  # analysis packages
  library("vegan")
  library("DESeq2")
  # utilities
  library("reshape")
  library("xlsx")  
}


rm(list=ls())

{dir.create("Data_check_Ipersal")
  dir.create("Results_Ipersal")
  dir.create("Results_Ipersal/Abundances")
}
options(scipen = 100) # disable scientific annotation


### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral",
                 "chartreuse","yellow3","darkmagenta","pink3", "gray25",
                 "firebrick2","gray","gold","darkgreen","violet",
                 "royalblue","wheat3","coral","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group


table_PHA<-read.table(file="PHA_Accumulators_agosto2025.tsv", fill = T, header = T, sep="\t")


# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME_IperSalino/table.qza", taxonomy="QIIME_IperSalino/taxonomy_SILVA.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

# Microaerophilia
metadata <- as.data.frame(read.table(file="Metadata_Ipers_Serena.txt", sep="\t", header = T))
row.names(metadata)<-metadata$FASTQ_ID
sample_data(data)<-metadata

sample_data(data)$Group <- factor(  sample_data(data)$Group , levels = c("Mixed liquor","Adherent"))
sample_data(data)$Experiment_Day <- as.factor(sample_data(data)$Experiment_Day )




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check_Ipersal/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check_Ipersal/Filtered_genus_under_0005_max_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured" & !is.na(filtered)] # to avoid the removal of other uncultured genera
data<-subset_taxa(data, ! Genus %in% filtered)
data<-prune_taxa(taxa_sums(data) > 0, data) 
rm(filtered, data.genus.temp)


######## Checking unassigned in prokaryote kingdom
{Unass<-tax_glom(data, taxrank = "Kingdom", NArm = F) # or domain
  Unass.prop<-transform_sample_counts(Unass, function (x) (x/sum(x)*100) )
  a<-cbind( apply(otu_table(Unass), sum, MARGIN = 1) )
  total<-sum(a)
  b<-cbind( (a/total)*100, apply(otu_table(Unass), min, MARGIN = 1) , apply(otu_table(Unass), max, MARGIN = 1 ) )
  c<-otu_table(Unass)
  c_a<-NULL
  for(x in 1:length(row.names(c))) {
    c_a<-c(c_a,colnames(c)[which.min(c[x,])]) }
  c_b<-NULL
  for(x in 1:length(row.names(c))) {
    c_b<-c(c_b,colnames(c)[which.max(c[x,])]) }
  d<-as(tax_table(Unass),"matrix")
  e<-cbind.data.frame(a,b, c_a, c_b, d[,"Kingdom"])
  colnames(e)<-c("total_ASV_count","Percent_proportion","Minor_count_among_samples","Major_count_among_samples", "FASTQ_with_less_count","FASTQ_with_more_count","Kingdom")
  e$Kingdom<-gsub("d__","",e$Kingdom)
  row.names(e)<-e$Kingdom
}
e
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check_Ipersal/Bacteria_Archaea_proportion_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)



# removing every Unassigned Phylum (even the phylum was not in SILVA? Probable contaminant!)
data<- subset_taxa(data, ! is.na(Phylum) )



# removing chloroplast and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria\n\n")
  data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mitochondria"))
}


proof1<-"This is the proof of having performed the filtering steps"





########################### % ASSIGNED IN SILVA ##################################

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}
{ Taxa.genus<-as.data.frame(tax_table(data.genus))
  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
  Taxa.class<-as.data.frame(tax_table(data.class))
  Taxa.order<-as.data.frame(tax_table(data.order))
}

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
  b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
  c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
  d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
  e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
  assigned<-rbind.data.frame(a,b,c,d,e)
  colnames(assigned)<-c("Total","Assigned","%","Taxa")
  assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
}
assigned
write.csv2(assigned,file="Data_check_Ipersal/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)




############################ RAREFACTION ANALYSIS ################################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as generated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of the original samples (and in the same order!!!)
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    check<-"red"
    if(length(x) < 10) {
      points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col="cyan",pch=16,cex=1)
    } else {
      #the slope is estimated (average) at the last b rarefaction points
      differ<- length(v)-t
      if(differ<0){differ<-0} # to avoid the error with negative values
      slope<-mean(diff(v[ (differ):length(v) ])/dx)
      if(slope<lim) {
        check<-"blue"
        sat = sat+1
      }
      cat(i,slope,check,"\n")
      #			text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],labels=slope,col=check,cex=0.5)
      #			points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col=check,pch=16,cex=1)
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  #legend("bottomright",paste(sat,"saturated samples"),bty="n")
}


png(file="Data_check_Ipersal/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data.genus),"matrix")), step=100,label=T, ylab = "Genera", xlab= "Reads amount")
# evalslopes(r,sample_names(data.genus),lim=0.001,cex=1)
dev.off()
rm(r)


proof2<-"No unsaturated sample to remove"




#################### \\\\\\ STARTING THE ANALYSIS \\\\\ #######################


if(! "proof1" %in% ls() | ! "proof2" %in% ls() ){
  stop("\nDid you perform the filtering steps yet?\n")
}


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



#save.image("Data_prepared_and_Filtered_for_Analysis.RData")



########################### COUNTS EXPORT ##########################################

dir.create("Results_Ipersal/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results_Ipersal/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_Ipersal/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results_Ipersal/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results_Ipersal/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results_Ipersal/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_Ipersal/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results_Ipersal/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results_Ipersal/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
}




###################### ABUNDANCES BAR PLOT (EVERY SAMPLE) ##########################

### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_taxa(data.phy.prop)
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
tabella$Experiment_Day<-factor(tabella$Experiment_Day, levels=unique(metadata$Experiment_Day))
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9.8 ),
  legend.position="bottom",
  legend.margin = margin(3 ,15, 8 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  scale_y_continuous(expand = c(0,1)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified phyla",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results_Ipersal/Abundances/TOP_phyla_EVERY_SAMPLE.png",width=4.2,height=3.2, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results_Ipersal/Abundances/TOP_phyla_averages_EVERY_SAMPLE.xlsx", row.names = F,to_save)

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)




### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:19]
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
  tabella$Genus<-gsub ("uncultured","uncult.", tabella$Genus)
  tabella$Genus<-gsub ("Marine_Benthic_Group_D_and_DHVEG-1","Marine_Benthic_D_DHVEG-1", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Experiment_Day<-factor(tabella$Experiment_Day, levels=unique(metadata$Experiment_Day))
fill_color_modified<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick4","gray","chartreuse","darkgreen","violet", "deepskyblue2","wheat3","chartreuse3","red2","darkslategray3")
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =7.5) +
  scale_fill_manual(values=fill_color_modified) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.15, "cm"),
  legend.key.width = unit(0.21, "cm"),
  legend.text = element_text ( size = 7 ),
  legend.position="bottom",
  legend.margin = margin(0 ,32, 4 ,5),
  plot.margin = margin(4,1,4,1)) +
  scale_y_continuous(expand = c(0,1)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified genera",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results_Ipersal/Abundances/TOP_genera_EVERY_SAMPLE.png",width=4.2,height=3.3, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "Averages_in_Adherent"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="Adherent")),1,mean)),
                           "Averages_in_MixedLiquor_without_Inoculum"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="Mixed liquor" & Sample_name!="Y1")),1,mean)),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results_Ipersal/Abundances/TOP_genera_Averages_EVERY_SAMPLE.xlsx", row.names = F, to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



###################### ABUNDANCES BAR PLOT (MIXED LIQUOR) ##########################

### TOP Phyla
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_samples( data.phy.prop , Group=="Mixed liquor")
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,top_data)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
}
{
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[! unique(tabella$Phylum) %in% "Others"],"Others"))
}
tabella$Experiment_Day<-factor(tabella$Experiment_Day, levels=unique(metadata$Experiment_Day))
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  plot.title = element_text( size= 5) , 
  # strip.text = element_text(size=8),
  legend.key.height = unit(0.35, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9.8 ),
  legend.position="bottom",
  legend.margin = margin(6 ,18, 12 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  scale_y_continuous(expand = c(0,1)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified phyla",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results_Ipersal/Abundances/TOP_phyla_MixedLiquor.png",width=4.2,height=3.6, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results_Ipersal/Abundances/TOP_phyla_averages_MixedLiquor.xlsx", row.names = F,to_save)

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)



### TOP Genera
suppressWarnings(rm(top, others, tabella, unass_data))
{ top_data <- subset_samples( data.genus.prop , Group=="Mixed liquor")
  top <- names(sort(taxa_sums(top_data), decreasing=TRUE))[1:19]
  prune.dat_top <- prune_taxa(top,top_data)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(top_data)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,top_data)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("uncultured","uncult.", tabella$Genus)
  tabella$Genus<-gsub ("Marine_Benthic_Group_D_and_DHVEG-1","Marine_Benthic_D_DHVEG-1", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Experiment_Day<-factor(tabella$Experiment_Day, levels=unique(metadata$Experiment_Day))
fill_color_modified<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick4","gray","chartreuse","darkgreen","violet", "deepskyblue2","wheat3","chartreuse3","red2","darkslategray3")

ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  # facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8) +
  scale_fill_manual(values=fill_color_modified) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  plot.title = element_text(size= 5),
  # strip.text = element_text(size=8),
  legend.key.height = unit(0.1, "cm"),
  legend.key.width = unit(0.21, "cm"),
  legend.text = element_text ( size = 7 ),
  legend.position="bottom",
  legend.margin = margin(-1 ,32, 3 ,5),
  plot.margin = margin(4,1,4,1)) +
  scale_y_continuous(expand = c(0,1)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified genera",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results_Ipersal/Abundances/TOP_genera_MixedLiquor.png",width=4.2,height=3.6, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "Averages_in_MixedLiquor_without_Inoculum"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="Mixed liquor" & Sample_name!="Y1")),1,mean)),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results_Ipersal/Abundances/TOP_genera_Averages_MixedLiquor.xlsx", row.names = F, to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




######################## ALPHA DIVERSITY ##############################

#data.genus.temp<-subset_samples(data.genus, Type=="Mixed liquor")
data.genus.temp<-data.genus
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$FASTQ_ID
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<- mix_alpha$Shannon /log(mix_alpha$Observed)
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha$Sample_name <- gsub("^X","",mix_alpha$Sample_name)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$Exp_day <- metadata[mix_alpha$Sample_name, "Experiment_Day" ]  # NB: the metadata has to be ordered according to time!
mix_alpha$Type <- metadata[mix_alpha$Sample_name, "Group" ]
mix_alpha$Type <- gsub("Suspended","Mixed liquor", mix_alpha$Type)
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
days <- unique(metadata$Experiment_Day)
mix_alpha$Time_axisX <- factor(as.character(mix_alpha$Exp_day), levels =days )

mix_alpha <- mix_alpha[ order(mix_alpha$Exp_day) ,]

ggplot(data=mix_alpha, aes(y=value, x=Time_axisX, color=Type)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  scale_color_manual(values= c("Mixed liquor"="deepskyblue2", "Adherent"="coral"  )
  )+
  geom_point(size =2, aes(shape=Type)) +
  geom_point(size =3.8, alpha=0.4, aes(shape=Type)) +
  geom_path(aes(group=Type), linewidth = 0.18, alpha=0.9,
            arrow=arrow(length =unit(0.22,"cm"), type = "closed")
  ) +
  theme_classic2(base_size = 10) +
  labs(y="Alpha Diversity Measure") +
  guides(fill="none", color="none", shape="none") +
  scale_x_discrete(expand = c(0.02,0)) + # to have more space on the borders
  # geom_text(aes(label=Exp_day), color="black", size=2, vjust=0.2, show.legend = FALSE, lineheight =0.7) +  # lineheight set the return ("\n") spacing
  theme(axis.text.x = element_text(size=8.5, angle=35, hjust=1, vjust=1),
        axis.text.y= element_text(angle=0, size=5.5),
        #axis.title.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.4),
        panel.grid.minor.y = element_line(linewidth=0.25),
        plot.margin = margin(5,1,5,1)
  ) +
  labs(x="Experiment day")
ggsave(file="Results_Ipersal/Alfa_diversity_at_GENUS_level.png", width = 4.65,height =3.6, dpi=300)




########################### PCoA BETA DIV ##########################

# on genera
data.prop.labels<-subset_samples( data.genus.prop, Group =="Mixed liquor" )

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color="Sample_name") + # , color = "Sample_Group") +
  scale_color_manual(values = c("Y1"="cyan2",
                                "YS2"="deepskyblue2",
                                "Y3"="royalblue",
                                "Y4"="blue"
                                )
  ) +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none", shape="none") +
  theme_classic(base_size = 12.5) +
  theme(title=element_text(size=9)) +
  geom_path(aes(group=as.character(Group)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  )+
  # geom_text(aes(label=Experiment_Day), color="grey60", size=3.5, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results_Ipersal/Beta_diversity_on_genera_MIXEDLIQUOR.png", width = 3.8, height = 3, dpi=300)




################# STUDING THE ABUNDANCES OF THE PHA ACCUMULATORS ######################

lista_PHA<-table_PHA[,1]     # they derive from papers

data.genus.temp<-data.genus.prop
prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-gsub("Candidatus_","",tabella$Genus)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))

tabella$Type <- gsub("Inoculum","", tabella$Type)

#fill_color_modified<-c("grey92","lightblue4","wheat","darkgreen","lightgreen","grey20","deepskyblue2","deeppink","cyan","coral","wheat4","darkmagenta","brown4","springgreen","yellow2","bisque3","pink3","darkslategray3","red2","chartreuse2", "green3","black","gold","gray","pink3","violet") 
fill_color_modified<-c("darkgreen","red4","deepskyblue4","deeppink","black","coral","gold3","wheat4","darkmagenta","springgreen","yellow2","darkslategray3","red2") 

ggplot(data=tabella, aes( x=Experiment_Day, y=Abundance, fill=Genus)) + 
  facet_grid2( ~Group, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.85) +
  theme_classic(base_size =11) +
  scale_fill_manual(values=c(fill_color_modified)) +
  theme(axis.text.x=element_text(angle=40, hjust=1,vjust=1, size=8.5),
        axis.title.y = element_text(size=10), 
        axis.text.y = element_text(size=8.5), 
        strip.text.x = element_text(size=10.5), 
        legend.key.height = unit(0.52, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text ( size = 11 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="right",
        plot.margin = margin(2,2,2,2),
        legend.margin =  margin(125,1,0,-8)
  ) +
  guides(fill=guide_legend(nrow=22)) + 
  scale_y_continuous(lim=c(0,100)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = "PHA accumulating genera found"
  )
ggsave(file="Results_Ipersal/Abundances/PHA_Accumulators.png",width=5,height=3.8, dpi=300)
dev.off( )


# WITHOUT HALOMONAS
tabella2 <- tabella[tabella$Genus!="Halomonas", ]
ggplot(data=tabella2, aes( x=Experiment_Day, y=Abundance, fill=Genus)) + 
  facet_grid2( ~Group, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.85) +
  theme_classic(base_size =11) +
  scale_fill_manual(values=c(fill_color_modified)) +
  theme(axis.text.x=element_text(angle=40, hjust=1,vjust=1, size=8.5),
        axis.title.y = element_text(size=10), 
        axis.text.y = element_text(size=8.5), 
        strip.text.x = element_text(size=10.5), 
        legend.key.height = unit(0.52, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text ( size = 11 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="right",
        plot.margin = margin(2,2,2,2),
        legend.margin =  margin(125,1,0,-8)
  ) +
  guides(fill=guide_legend(nrow=22)) + 
  scale_y_continuous(lim=c(0,1.5)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = "PHA accumulating genera found (no Halomonas)"
  )
ggsave(file="Results_Ipersal/Abundances/PHA_Accumulators_WITHOUT_HALOMONAS.png",width=5,height=3.8, dpi=300)
dev.off( )



# means of PHA acc
to_save<- cbind.data.frame(
                           "Averages_in_Inoculum"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="Mixed liquor" & Sample_name =="Y1" )),1,mean)),
                           "Averages_in_MixedLiq_without_Inoculum"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="Mixed liquor" & Sample_name !="Y1")),1,mean)),
                           "Averages_in_Adherent_Biomass"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="Adherent")),1,mean)),
                           "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]],
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
)
to_save<-to_save[order(to_save$Averages_in_MixedLiq_without_Inoculum, decreasing=T), ]
write.csv(to_save, file = "Results_Ipersal/Abundances/PHA_Accumulators_averages.csv", row.names = F, quote=F)






##################### \\\\\\ R AND PACKAGES VERSION \\\\\\ #########################


### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results_Ipersal/R_version_and_packages.txt")
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
print("ggh4x")
packageVersion("ggh4x")
cat("\n", "\n", fill=TRUE)
print("vegan")
packageVersion("vegan")
cat("\n", "\n", fill=TRUE)
print("microbiome")
packageVersion("microbiome")
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))

