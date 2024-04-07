##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  # analysis packages
  library("DESeq2")
  library("vegan")
  # utilities
  library("xlsx")  
  library("qiime2R")
}

{dir.create("Data_check")
dir.create("Results")
dir.create("Results/Abundances")
dir.create("Results/Beta_div")
dir.create("Results/DA_DESeq2")
}

options(scipen = 100) # disabling scientific annotation


# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be set as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3")
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3")



####################### IMPORTING DATA #######################

data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample
sample<-substring(sample, first = 10, last= 14) # strtrim to cut until X
sample<-gsub("-","",sample, fixed=T)
sample<-gsub("A","",sample)
sample_names(data)<-sample # update

Metadata <- read.csv2("metadata.csv", header = T)
colnames(Metadata)<-gsub("X","",colnames(Metadata))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,original_names,sample)

sample_data(data)$Time<-factor(sample_data(data)$Time, levels=c("T0","T1")) # decides the order in plots
head(sample_data(data))


# save.image("data.RData")



############################ RAREFACTION ANALYSIS ################################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as gnerated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of he origina samples (and in the same order!!!)
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
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F, ylab = "ASVs")
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)



####################### REMOVING THE OUTLIERS ###########################

if(! "GP22" %in% sample_names(data)){
  cat("\n !!! Looks like the outliers have been already removed...\n\n")
  Sys.sleep(2)}

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}


# TOP 5 Phylum con stacked bar plot
suppressWarnings(rm(top5, others, tabella))
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Phylum)) +
  theme_classic(base_size =12) + 
  facet_grid2(Sample_ID+Time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), 
        strip.text.y = element_text(angle=0, size=9.5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 13 ),
        legend.position="bottom",legend.margin = margin(2,0,3,15)
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Percent abundance", fill="", title="with GP22", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Data_check/TOP5_phyla_withGP22.png",width=6,height=8, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Data_check/TOP5_phyla_withGP22_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top5))[["Phylum"]]))



# TOP 10 Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub("[","",tabella$Genus, fixed=T)
  tabella$Genus<-gsub("]","",tabella$Genus, fixed=T)
  tabella$Genus<-gsub("_group","",tabella$Genus, fixed=T)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Genus)) +
  theme_classic(base_size =12) + 
  facet_grid2(Sample_ID+Time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), 
        strip.text.y = element_text(angle=0, size=9.5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom",legend.margin = margin(-5,0,0,36)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Patients", x="Percent abundance", fill="", title="with GP22", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Data_check/TOP_Genera_withGP22.png",width=6,height=8, dpi=300) 
dev.off()

# means of TOP10 genera (no inoculum)
write.xlsx(file = "Data_check/TOP_Genera_withGP22_Averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))



# PCoA
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  scale_color_manual(values=c("T1"="coral","T0"="chartreuse2")) +
  geom_point(size=3, alpha=0.3) +
  theme_classic(base_size = 11) + 
  labs(color="Time",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples"
  )
ggsave(file="Data_check/PCoA_with_GP22.png", width = 5, height = 4, dpi=300)



# backup
if("GP22" %in% sample_names(data)){data_with_ORT_2_sample<-data}
### --> Removing GP22 (About 60% of its composition are soil/plant bacteria)
data<-subset_samples(data, Sample_ID!="ORT_2")# removing also its pair (this is a paired analysis)

write.csv2(cbind(otu_table(data_with_ORT_2_sample),tax_table(data_with_ORT_2_sample)), file="Data_check/Raw_feature_table_ASV.csv", row.names = T)

rm(data.phy,data.genus,data.genus.prop,data.phy.prop)




#################### FILTERING NOISES FROM DATA SET ####################

unfiltered_data<-data

if("GP22" %in% sample_names(data)){
  stop("\n !!! Looks like the outliers are still in your data set...\n\n")
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genera_under_0005_max_cutoff.csv")

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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Unassigned_domain_checking.csv", row.names = T, quote = F)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


### removing mitochondria and chloroplast
if( "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Mitochondria from the data\n")
  data<-subset_taxa(data, Genus != "Mitochondria")
}
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast from the data\n")
  data<-subset_taxa(data, Genus != "Chloroplast")
}

proof1<-"This is the proof of having performed the filtering steps"




####################### PREPARATION OF THE DATA #######################

if(! "proof1" %in% ls() ){
  stop("\nDid you perform the filtering steps yet?\n")
}


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



#################### % ASSIGNED IN SILVA #########################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Total","Assigned","%","Taxa")
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)


######################## GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(data)==1, data)
length(which(taxa_sums(data)==1)) # if zero there are no singletons
{n<-as.data.frame(otu_table(filter))
  N<-as.data.frame(otu_table(data))
  G<-1-(colSums(n)/colSums(N))
}
con<-file("Data_check/Percentuale_di_singletons_Good's_coverage.txt")
sink(con)
cat("GOOD'S COVERAGE ESTIMATOR \n", fill=TRUE)
cat("1-(n/N) for each sample, where n is number of singletons and N is Total ASV \n \n", fill=TRUE)
G
sink()
close(con)

rm(con, filter, G, n, N)



########################### COUNTS EXPORT ##########################################

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}



###################### ABUNDANCES BAR PLOT ##########################


# TOP 5 Phylum con stacked bar plot
suppressWarnings(rm(top5, others, tabella))
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.prop)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Phylum)) +
  theme_classic(base_size =12) + 
  facet_grid2(Sample_ID+Time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), 
        strip.text.y = element_text(angle=0, size=9.5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 13 ),
        legend.position="bottom",legend.margin = margin(2,0,3,15)
        ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Percent abundance", fill="", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/TOP5_phyla_abundances.png",width=6,height=8, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Abundances/TOP5_phyla_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top5))[["Phylum"]]))



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
tabella$Genus<-gsub("[","",tabella$Genus, fixed=T)
tabella$Genus<-gsub("]","",tabella$Genus, fixed=T)
tabella$Genus<-gsub("_group","",tabella$Genus, fixed=T)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Genus)) +
  theme_classic(base_size =12) + 
  facet_grid2(Sample_ID+Time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), 
        strip.text.y = element_text(angle=0, size=9.5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom",legend.margin = margin(-5,0,0,36)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Patients", x="Percent abundance", fill="", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Abundances/TOP10_genera_abundances.png",width=6,height=8, dpi=300) 
dev.off()

# means of TOP10 genera (no inoculum)
write.xlsx(file = "Results/Abundances/TOP_10_genera_Averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



########################## ALFA DIVERSITY ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated on genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Time", color="Time")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Time, New_data$Sample_ID),]
pAlpha$data<-New_data
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha + theme_classic2(base_size = 9) +
  scale_color_manual(values = c("T1"="chartreuse2", "T0"="coral")) +
  geom_line(aes(group = pAlpha$data$Sample_ID),col="grey",size=0.18) +
  geom_point(alpha=0.5, size=2.8) +
  labs(x="") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=8.5),
        strip.text = element_text(size=9.2)
        ) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=2.8, label.y.npc = "top", vjust=-0.5, hjust=0.4) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Alfa_diversity_GENUS_paired_wilcoxon.png", width = 4.5,height =4, dpi=300)


suppressWarnings( rm(pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor) )



##################### BETA DIVERSITY #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Beta_div/Beta_divers_permanova_Helling.csv",quote=F,row.names = T)


### PLUS: checking it with Bray too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Sample_ID + Time, data=metadata, permutations = 9999, method="bray")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Beta_div/Beta_divers_permanova_BRAY_on_genera.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on genera
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Time)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)



######################## PCoA BETA DIV #########################

# on genera
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# with names
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  scale_color_manual(values=c("T1"="coral","T0"="chartreuse2")) +
  geom_point(size=3, alpha=0.3) +
  theme_classic(base_size = 11) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), color="black", size=2, show.legend = FALSE) +
  labs(color="Time",
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H),
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples"
  )
ggsave(file="Results/Beta_div/PCoA_Beta_div_Hellinger_genera_names.png", width = 5, height = 4, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_ID),col="grey", size=0.15)+
  scale_color_manual(values=c("T1"="coral","T0"="chartreuse2")) +
  geom_point(size=3.2, alpha=0.3) +
  theme_classic(base_size = 11) + 
  labs(color="Time",
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H),
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples"
  )
ggsave(file="Results/Beta_div/PCoA_Beta_div_Hellinger_genera_no_names.png", width = 4.6, height = 4, dpi=300)




############## DIFFERENTIAL ABUNDANCES WITH DESEQ2 #############

if(! "proof1" %in% ls() ){
  stop("\nDid you perform the filtering steps yet?\n")
}

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
system("echo 'In case of pair analysis, the log2foldchange displayes the change from T0 to T1, see https://support.bioconductor.org/p/105981/' > Results/DA_DESeq2/NB.txt ")  

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
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Sample_ID+Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "T0", "T1"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  # res<-res[res$baseMean > 50, ] # arbitrary threshold to avoid the most noisy results
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
    write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_T0_vs_T1.csv"), row.names = F, quote=F, na = "")
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

View(Res_tot)
#write.csv(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Time)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Time)

unique(tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_group","",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("[","",tabella_g$Bacteria, fixed = T)
tabella_g$Bacteria<-gsub("]","",tabella_g$Bacteria, fixed = T)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Time)) +
  theme_classic(base_size = 10.5) +
  scale_color_manual(values = c("T0"="coral","T1"="chartreuse2")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = "Genus")~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time), size=2.2, alpha=0.9) +
  geom_point(aes(color=Time), size=3.8, alpha=0.5) +
  geom_line(aes(group=Sample_ID), size=0.2, color="black", alpha=0.65) +
  theme(strip.text.x=element_text(size=7.9,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        axis.text.y = element_text(size=7.5)
        )+
  scale_x_discrete(labels=rep(unique(levels(tabella_g$Time))), expand=c(0,0.5)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank(),
        panel.grid.major.y = element_line(size=0.2, color="gray"),
  ) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none") +
  scale_y_sqrt(breaks=c(0.01, 0.05, 0.1, 0.25,0.75,seq(0,8,0.5)))
plot_1 +
  labs(y="Percent abundance", fill="Time", x="")
ggsave(filename = "Results/DA_DESeq2/Plot_genera_DeSeq2_Time.png", width = 5.2, height = 4, dpi=300)
dev.off()

plot_2<-ggplot(tabella_2, aes(x= Xaxis, y=Abundance, fill=Time)) +
  theme_classic(base_size = 10.5) +
  scale_color_manual(values = c("T0"="coral","T1"="chartreuse2")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = "Family")~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time), size=2.2, alpha=0.9) +
  geom_point(aes(color=Time), size=3.8, alpha=0.5) +
  geom_line(aes(group=Sample_ID), size=0.2, color="black", alpha=0.65) +
  theme(strip.text.x=element_text(size=7.9,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line"),
        axis.text.y = element_text(size=7.5)
  )+
  scale_x_discrete(labels=rep(unique(levels(tabella_g$Time))), expand=c(0,0.3)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank(),
        panel.grid.major.y = element_line(size=0.2, color="gray"),
  ) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none") +
  scale_y_sqrt(breaks=c(0.01, 0.05, 0.1, 0.25,0.75,seq(0,8,0.5)))
plot_2 +
  labs(y="Percent abundance", fill="Time", x="")

ggsave(filename = "Results/DA_DESeq2/Plot_2_DESeq2_Time_no redundants.png", width = 5, height = 4, dpi=300)
dev.off()



##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results/R_version_and_packages.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"

cat(package$R.version$version.string)
cat("   running on", package$running)
cat("\n", "\n", fill=TRUE)
package$otherPkgs$DESeq2[c(1,4)]
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
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)
#print(package$loadedOnly)


sink() # restore STR OUTPUT to R console
close(con)
suppressWarnings(rm(con))
