##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("vegan")
  # utilities
  library("xlsx")  
  library("qiime2R")
}

dir.create("Data_check")
dir.create("Results")
dir.create("Results/Abundances")

options(scipen = 100) # disable scientific annotation



### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
# NB: there is always an extra color which will be the "Others" group



####################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample

Metadata <- as.data.frame(read.table("Metadata_Alex_Saline.tsv",sep="\t",header = T))
row.names(Metadata)<-Metadata$FASTQ_Code
head(Metadata)
original_length<-length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_Code[!is.na(Metadata$FASTQ_Code)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,sample)

# head(sample_data(data))

sample_names(data)<-as.factor(sample_data(data)$Sample_name) # NB: they are already have been labeled according to the time! Moreover, also in the Metadata their order is correct (this is important)
sample_data(data)$time_points<-as.numeric(gsub("AX","",sample_data(data)$Sample_name))
sample_data(data)$Experiment_day_char<-factor( as.character(sample_data(data)$Experiment_day), levels= sample_data(data)$Experiment_day )




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
      slope<-mean(diff(v[(length(v)-t):length(v)])/dx)
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
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

png(file="Data_check/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
evalslopes(r,sample_data(data)$Sample_name,lim=0.001,cex=1)
dev.off()
rm(r)




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




#################### \\\\\\ STARTING THE ANALYSIS \\\\\ #######################


if(! "proof1" %in% ls()  ){
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

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  #write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
}




###################### ABUNDANCES BAR PLOTS ##########################


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

ggplot(data=tabella, aes(x=Experiment_day_char, y=Abundance, fill=Phylum)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=6.5),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 9.8 ),
  legend.position="bottom",
  legend.margin = margin(-8 ,5, 5 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant identified phyla",
       fill="",
       caption = " 'Others' includes every phylum below rank 10 ")
ggsave(file="Results/Abundances/TOP_phyla.png",width=6,height=4.5, dpi=300)
dev.off()

# means of TOP phyla
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), 
                           "standard_dev"=as.numeric(apply(otu_table(prune.dat_top),1,sd)),
                           "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]])
to_save<-to_save[order(to_save$mean, decreasing=T), ]
to_save
write.xlsx(file = "Results/Abundances/TOP_phyla_averages_EverySample.xlsx", row.names = F, to_save)


rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, prune.dat_mature)



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
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

ggplot(data=tabella, aes(x=Experiment_day_char, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_19) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=6.5),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 8.1 ),
  legend.position="bottom",
  legend.margin = margin(-10,10,5,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="", y="Percentual abundance of clades",
       title = "Most abundant identified genera",
       fill="",
       caption = " 'Others' includes every genus below rank 19 ")
ggsave(file="Results/Abundances/TOP_genera.png",width=6,height=4.5, dpi=300)

dev.off()


# means of TOP genera
to_save<- cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)),
                           "standard_dev"=as.numeric(apply(otu_table(prune.dat_top),1,sd)),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
                           )
to_save<-to_save[order(to_save$mean, decreasing=T), ]
write.xlsx(file = "Results/Abundances/TOP_genera_Averages_EverySample.xlsx", row.names = F, to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, to_save))





########################### PCoA BETA DIV ##########################

# on genera
data.prop.labels<-data.genus.prop


# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC , color = "Sample_Type") +
  scale_color_manual(values="chartreuse2") + 
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none", shape="none") +
  theme_classic(base_size = 12.5) +
  theme(title=element_text(size=10)) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day), color="black", size=3, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_divers_Hellinger_on_genera_WITH_NAMES.png", width = 6.5, height = 5, dpi=300)


# with shade of color AND also exp day
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(86,274,410,555),
                        midpoint = 270) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=as.character(Sample_Type)), col= "darkgray", linewidth = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed")
  ) +
  theme_classic(base_size = 10) +
  theme(title=element_text(size=8.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9.2),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  guides(shape="none") +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day_char), color="grey35", vjust= -0.5, size=2.2, show.legend = FALSE) +
  labs(title="Beta diversity computed on Hellinger distance \nusing genera proportional abundances",
       color="Experiment day",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/Beta_divers_Hellinger_on_genera_ShadeAlongTheTime.png", width = 5, height = 4.2, dpi=300)




################### CORRELATIONS ABUNDANCES vs TIME _  (CENTERED LOG RATIO) ####################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
#data.genus.temp <- subset_samples(data.genus, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
data.genus.temp<- data.genus
data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow
abundances<-abundances[ , sample_names(data)]  # NB: the sample names are already ordered, as set in the loaded metadata tsv


Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "kendall") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA", "uncultured_ f Unknown_Family") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")

write.csv2(Corr_results, "Results/Kendall_correlation_CLR_bacteria_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )

con <- file("Results/README_Only_certain_genera_have_been_tested_in_correlation.txt")
sink(con, append = T)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))




########### LINE PLOTS OF STRONGEST CORRELATIONS


### THE GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
# corr_top<-significative_abundances_pos[1:12] # NB: the table is already ordered
corr_top<-significative_abundances_pos # each of them
suppressWarnings(rm(top, others, tabella))


data.genus.temp <- data.genus.prop
{ tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% corr_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella<-tabella[order(tabella$Sample_name) , ] # days
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day_char, color=Genus)) +
  theme_classic(base_size =12) +
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), linewidth= 0.4, alpha= 0.4) +
  # scale_color_manual(values= fill_color_19) +
  scale_color_manual(values= c("brown4","darkblue")) +   # matching the colors in the top abundances bar plot
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 13 ),
        legend.position="bottom", 
        legend.margin = margin(1,20,10,0)) +
  guides(color=guide_legend(nrow=1)) +
  # scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  scale_y_continuous(breaks = seq(0, 2.5, 0.5) , lim=c(0,2.5) ) +
  labs(x="Experiment day", y="Percent abundance",
       color="",
       title = "genera positely correlated with time\n according to Kendall correlation on CLR counts")
ggsave(file="Results/Genera_CLR_positively_correlated_with_time.png",width=7,height=5, dpi=300)
dev.off()



### THE GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
# corr_top<-significative_abundances_neg[(length(significative_abundances_neg)-11):length(significative_abundances_neg)]
corr_top<-significative_abundances_neg # all of them
suppressWarnings(rm(top, others, tabella))
data.genus.temp <- data.genus.prop
{
  tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% corr_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella<-tabella[order(tabella$Sample_name) , ] # days
tabella$Sample_name<-factor(tabella$Sample_name, levels= sample_names(data))
ggplot(data=tabella, aes(y=Abundance, x=Sample_name, color=Genus)) +
  theme_classic(base_size =12) +
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  scale_color_manual(values= fill_color_19) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=9),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12 ),
        legend.position="bottom", legend.margin = margin(1,36,5,0)) +
  guides(color=guide_legend(nrow=2)) +
  # scale_y_sqrt(breaks = c(seq(0, 8, 1), seq( 10, max(tabella$Abundance +2 ), 2) ) ) +
  labs(x="Time point", y="Percent abundance",
       color="",
       title = "genera negatively correlated with time \n according to Kendall correlation on CLR counts")
ggsave(file="Results/Genera_CLR_negatively_correlated_with_time.png",width=7,height=5, dpi=300)
dev.off()


suppressWarnings(rm(prune.dat_top, tax_selected, tabella, corr_top))




##################### \\\\\\ R AND PACKAGES VERSION \\\\\\ #########################


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
package$otherPkgs$DESeq2[1:2]
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


