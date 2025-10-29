##################### PREPARING THE R PACKAGES #################

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




############## PREPARING THE ENVIRONMENT FOR THE ANALYSIS #############

rm(list=ls())

{dir.create("Data_check")
dir.create("Results_Microaero")
dir.create("Results_Microaero/Abundances")
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


# # importing the sample 1 phy obj ...
# load("S1_sample_phyloseq.RData")
# sample_names(S1)<-sample_data(S1)$FASTQ_ID
# S1 <- phyloseq(otu_table(S1),tax_table(S1)) # removing the sample data to merge the objs afterwards

# importing the other samples ...
# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME_Microaeroph/table.qza", taxonomy="QIIME_Microaeroph/taxonomy_SILVA.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

#data <- merge_phyloseq(S1,data)

metadata <- as.data.frame(read.table(file="Metadata_Microaero.tsv", sep="\t", header = T))
row.names(metadata)<-metadata$FASTQ_ID
sample_data(data)<-metadata

# setting the order in the plots
sample_data(data)$Type<-factor(sample_data(data)$Type, levels = c("Inoculum","Control","Microaerophilic") )
sample_data(data)$Sample_name<-factor(sample_data(data)$Sample_name, levels = metadata$Sample_name )  # NB: the order of the metadata decides the order of the levels
sample_data(data)$Experiment_Day<-factor(sample_data(data)$Experiment_Day, levels = unique(metadata$Experiment_Day[order(metadata$Experiment_Day)]) )



#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_max_cutoff.csv")

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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Bacteria_Archaea_proportion_checking.csv", row.names = T, quote = F)

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
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
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


png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
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




########################### COUNTS EXPORT ##########################################

dir.create("Results_Microaero/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results_Microaero/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_Microaero/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results_Microaero/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results_Microaero/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results_Microaero/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_Microaero/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results_Microaero/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_genus.csv",quote=F)
    #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results_Microaero/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
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
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 8
                                ),
  axis.title.x = element_text(vjust=0.1),
  axis.title.y = element_text(size=8.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=9),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 11.5 ),
  legend.position="bottom",
  legend.margin = margin(5 ,5, 8 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified phyla",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results_Microaero/Abundances/TOP_phyla.png",width=5,height=4, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results_Microaero/Abundances/TOP_phyla.xlsx", row.names = F,to_save)

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
  tabella$Genus<-gsub ("f Paludibacteraceae","f Paludibacter.", tabella$Genus)
  tabella$Genus<-gsub ("o Rhodospirillales","o Rhodospiril.", tabella$Genus)
  tabella$Genus<-gsub ("f Rhodobacteraceae","f Rhodobacter.", tabella$Genus)
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Experiment_Day, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_19) +
  theme(axis.text.x=element_text(angle=42,
                                 vjust=1,
                                 hjust=1,
                                 size= 8
  ),
  axis.title.x = element_text(vjust=0.1),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=9),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.text = element_text ( size = 8.4 ),
  legend.position="bottom",
  legend.margin = margin(0,26,0,1),
  plot.margin = margin(4,2,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="Experiment day", y="Percentual abundance of clades",
       title = "Most abundant identified genera",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results_Microaero/Abundances/TOP_genera.png",width=5,height=4, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame("Overall_averages"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "Averages_in_R1"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="R1")),1,mean)),
                           "Averages_in_R2"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="R2")),1,mean)),
                           "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]],
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
                           )
to_save<-to_save[order(to_save$Overall_averages, decreasing=T), ]
write.xlsx(file = "Results_Microaero/Abundances/TOP_genera_Averages.xlsx", row.names = F, to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



######################## ALPHA DIVERSITY ##############################

data.genus.temp<-data.genus
sample_names(data.genus.temp)<-sample_data(data.genus.temp)$FASTQ_ID
mix_alpha<- estimate_richness(data.genus.temp , measures = c("Observed","Shannon"))
mix_alpha$Evenness<- mix_alpha$Shannon /log(mix_alpha$Observed)
mix_alpha$Sample_name<-row.names(mix_alpha)
mix_alpha$Sample_name <- gsub("^X","",mix_alpha$Sample_name)
mix_alpha <- melt(mix_alpha,id.vars = "Sample_name")
mix_alpha$variable<-gsub("Observed", "Observed richness", mix_alpha$variable)
mix_alpha$Exp_day <- metadata[mix_alpha$Sample_name, "Experiment_Day" ]  # NB: the metadata has to be ordered according to time!
mix_alpha$Type <- metadata[mix_alpha$Sample_name, "Type" ]
mix_alpha$variable<-factor(mix_alpha$variable, levels= c("Observed richness","Shannon","Evenness"))
days <- unique(metadata$Experiment_Day)
mix_alpha$Time_axisX <- factor(as.character(mix_alpha$Exp_day), levels =days )

mix_alpha <- mix_alpha[ order(mix_alpha$Exp_day) ,]

ggplot(data=mix_alpha, aes(y=value, x=Time_axisX, color=Type)) +
  facet_grid( variable ~ . , scales = "free", space = "free_x") +
  scale_color_manual(values= c("Control"="deepskyblue2", "Microaerophilic"="coral"  )
  )+
  geom_point(size =1.8, aes(shape=Type)) +
  geom_point(size =3.2, alpha=0.4, aes(shape=Type)) +
  geom_path(aes(group=Type, color="gray99"), linewidth = 0.15, alpha=0.8,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
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
        plot.margin = margin(4,1,4,1)
  ) +
  labs(x="Experiment day")
ggsave(file="Results_Microaero/Alfa_diversity_at_GENUS_level.png", width = 4.6,height =3.8, dpi=300)




########################### PCoA BETA DIV ##########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Group<-factor(sample_data(data.prop.labels)$Group, levels=c("R1","R2"))
# sample_data(data.prop.labels)$line1 <- as.character(sample_data(data.prop.labels)$Type) # to connect the points of the same reactor
# sample_data(data.prop.labels)$line2<- seq(1, length(sample_data(data.prop.labels)$Sample_name), 1)
# sample_data(data.prop.labels)$line2[sample_data(data.prop.labels)$Sample_name %in% c("S1","S2") ] <-"pezzo A"
# sample_data(data.prop.labels)$line3<- seq(1, length(sample_data(data.prop.labels)$Sample_name), 1)
# sample_data(data.prop.labels)$line3[sample_data(data.prop.labels)$Sample_name %in% c("S1","F2") ] <-"pezzo B"


# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color="Type") + # , color = "Sample_Type") +
  scale_color_manual(values = c("Inoculum"="deepskyblue",
                                "Control"="chartreuse3",
                                "Microaerophilic"="coral")
  ) +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none", shape="none") +
  theme_classic(base_size = 12.5) +
  theme(title=element_text(size=10)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_name), color="black", size=3, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results_Microaero/Beta_divers_Hellinger_on_genera_WITH_NAMES.png", width = 4, height = 3.2, dpi=300)


# with shade of color AND also exp day
to_plot <- plot_ordination(data.sqrt_prop, ordBC, color = "Group") +
  scale_color_manual(values = c("R1"="deepskyblue2",
                                "R2"="coral2")
  ) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       color="",
       x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
to_plot$data <- to_plot$data[order(as.numeric(to_plot$data$Experiment_Day)) , ]
to_plot$data$Group <- as.factor(to_plot$data$Group)
to_plot$data$Experiment_Day <- as.factor(to_plot$data$Experiment_Day)
to_plot + 
  geom_path(aes(group=as.character(Group)), col= "darkgray", linewidth = 0.3,
          arrow=arrow(length =unit(0.25,"cm"), type = "closed")
          )+
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  # geom_line(aes(group=as.character(line2)), col= "darkgray", linewidth = 0.3) +
  # geom_line(aes(group=as.character(line3)), col= "darkgray", linewidth = 0.3) +
  theme_classic(base_size = 10) +
  theme(title=element_text(size=8.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(-5,0,0,-10),
        legend.text = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none") +
  geom_text(aes(label=Experiment_Day), color="grey25",
            vjust= 0, hjust= 0.3, size=3.2,
            show.legend = FALSE)
ggsave(file="Results_Microaero/Beta_divers_Hellinger_on_genera_ShadeAlongTheTime.png", width = 3.8, height = 3.2, dpi=300)





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

# fill_color_20<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","firebrick3","springgreen4","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","darkgreen","lightgrey","coral","deepskyblue2") 
# fill_color_19<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","darkgreen","firebrick3","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","springgreen3","coral","deepskyblue2") 
fill_color_19_modified<-c("grey92","lightblue4","darkviolet","wheat","lightgreen","grey20","deepskyblue2", "darkblue", "deeppink","cyan","coral","wheat4","darkmagenta","brown4","lightpink","red","darkgreen", "darkslategray3", "chartreuse2","black","gold","gray","pink3","violet")  

ggplot(data=tabella, aes( x=Experiment_Day, y=Abundance, fill=Genus)) + 
  facet_grid2( ~Group, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.85) +
  theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_19_modified)) +
  theme(axis.text.x=element_text(angle=40, hjust=1,vjust=1, size=8.8),
        axis.title.y = element_text(size=11), 
        axis.text.y = element_text(size=9.5), 
        strip.text.x = element_text(size=10), 
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text ( size = 9.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="right",
        plot.margin = margin(2,2,2,2),
        legend.margin =  margin(25,1,0,-8)
  ) +
  guides(fill=guide_legend(nrow=22)) + 
  scale_y_continuous(lim=c(0,100)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = "PHA accumulating genera found"
       )
ggsave(file="Results_Microaero/Abundances/PHA_Accumulators.png",width=5.25,height=3.8, dpi=300)
dev.off( )

# means of PHA acc
to_save<- cbind.data.frame("Overall_averages"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "Averages_in_R1"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="R1")),1,mean)),
                           "Averages_in_R2"=as.numeric(apply (otu_table(subset_samples(prune.dat_top, Group=="R2")),1,mean)),
                           "Family"= as.data.frame(tax_table(prune.dat_top))[["Family"]],
                           "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
)
to_save<-to_save[order(to_save$Overall_averages, decreasing=T), ]
write.csv(to_save, file = "Results_Microaero/Abundances/PHA_Accumulators_averages.csv", row.names = F, quote=F)



##################### DA WITH DESEQ2 (R1 vs R2) #######################

suppressWarnings(rm(data_pruned, data.genus_pruned, subset_target))
subset_target <- subset_samples(data, Sample_name %in% c("S9","F9","S10","F10","S11","F11") )
data_pruned<- prune_taxa(taxa_sums(subset_target) > 10, subset_target)
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

sample_data(data_pruned)$Exp_day_character <- paste0("Day",sample_data(data_pruned)$Experiment_Day)

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
  DEseq_data<-phyloseq_to_deseq2(d, ~ Exp_day_character + Group)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Group", "R1", "R2"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1 ,]
  #res<-res[(res$padj < 0.05) ,]
  res<-res[res$baseMean > 100, ] # arbitrary threshold to avoid the most noisy result
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
write.csv2(Res_tot, file="Results_Microaero/Differently_abundant_bacteria_between_mature_R1_and_R2.csv", row.names = F)



############ PLOTTING THESE RESULTS ...

Table_tot$Group<-factor(Table_tot$Group, levels = unique(Table_tot$Group))
Table_tot<-Table_tot[order(match(Table_tot$Group, levels(Table_tot$Group))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Group)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Group)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

## to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
#levels(tabella_g$Xaxis)
#levels(tabella_g$Group)


#unique(tabella_g$Bacteria)
tabella_g$Group <- gsub("_"," ",tabella_g$Group)
tabella_g$Bacteria<-gsub("_marine","\nmarine",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("Candidatus_","Candidatus\n",tabella_g$Bacteria)


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

Redund<-c(auto_redund , "Sphingobacteriales")
tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)

tabella3 <-rbind.data.frame(tabella_g, tabella_2)   # few results, it's better to re-join them ...
tabella3$Bacteria<-gsub("_ ", "\n", tabella3$Bacteria)
tabella3$Bacteria<-gsub("_", "\n", tabella3$Bacteria)
tabella3$Bacteria<-gsub("Peptostreptococcales-", "Peptostreptococc.\n", tabella3$Bacteria, fixed=T)
tabella3 <- tabella3[ ! is.na(tabella3$Bacteria) , ]
these_colors <- c( "R2"="coral3", "R1"="deepskyblue" )
plot_DESEQ <- ggplot(tabella3, aes(x= Xaxis, y=Abundance, fill=Group)) +
  theme_classic2(base_size = 12) +
  scale_fill_manual(values = these_colors ) +
  scale_color_manual(values = these_colors ) +
  facet_wrap2(nrow=3,factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus"))~Bacteria,
              #labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free",
              strip=strip_nested(size = "variable", bleed = T),
              drop = TRUE) +
  geom_line(aes(group=Time_point), linewidth=0.15, color="gray30", alpha= 0.5) +
  scale_x_discrete(labels=unique(levels(tabella_g$Group)), expand=c(0,0.5)) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(strip.text.x=element_text(size=6,colour="black"),
        strip.switch.pad.wrap = unit(10,"line"),
        # axis.text.x = element_text(size=10, angle=40, hjust=1, vjust=1),
        axis.title.y = element_text(size=9.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=10.5),
        axis.text.y = element_text(size=6.5),
        panel.grid.major.y = element_line(linewidth=0.2, color="gray"),
        panel.grid.minor.y= element_blank(),
        plot.title= element_text(size=12),
        panel.spacing.x = unit(1, "pt"),
        plot.caption = element_text(size=6)
  ) +
  # guides(fill="R1") +1
  labs(y="Percentual abundance", fill="", color="", x="", shape="",
       caption="The numbers of each point is the experiment day"
       ) +
  theme(legend.margin=margin(-20, 0, 0, 0),
        legend.position = "bottom")

plot_DESEQ + geom_point(aes(color=Group), size=1.8, alpha=1) +
  geom_point(aes(color=Group), size=3.8, alpha=0.5) +
  geom_text(aes(label= Experiment_Day  ), size=2.4, color="gray25", show.legend = F)
ggsave(filename = "Results_Microaero/Differently_abundant_bacteria_between_mature_R1_and_R2.png", width = 5.15, height = 5, dpi=300)



#### CIRCOTAX
library(circotax)
Res_tot_CIRCO <- Res_tot[Res_tot$Taxon=="Genus" , ]
Res_tot2 <- Res_tot[ Res_tot$Taxon %in% c("Family") , ]
Res_tot2 <- Res_tot2[ !Res_tot2$Family %in% Redund , ]
Res_tot2 <- Res_tot2[ !is.na(Res_tot2$Family), ]
Res_tot_CIRCO <- rbind(Res_tot_CIRCO, Res_tot2)
Res_tot_CIRCO$Genus <- gsub("uncultured_","uncultured\n", Res_tot_CIRCO$Genus)

CircoTax(Res_tot_CIRCO,  fc_col = 2 , tax_col = 8:12,title = "")
ggsave(file="Results_Microaero/Differently_abundant_bacteria_FoldChange_representation.png" , width = 6.85, height = 6.8)




##################### \\\\\\ R AND PACKAGES VERSION \\\\\\ #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results_Microaero/R_version_and_packages.txt")
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



