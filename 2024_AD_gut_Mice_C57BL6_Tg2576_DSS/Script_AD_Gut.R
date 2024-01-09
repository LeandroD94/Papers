##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
  library("ggplot2")
  library("ggh4x")
  library("vegan")
  library("ggpubr")
  library("dendextend")
  library("DESeq2")
  library("xlsx")
  library("qiime2R")
}

{dir.create("Data_check")
dir.create("Results")
dir.create("Results/Abundances")
dir.create("Results/Hierarchical_clustering")
dir.create("Results/Alpha_diversity")
dir.create("Results/Bray_Curtis_Beta_div")
dir.create("Results/w_Unifrac_Beta_div")
dir.create("Results/DA_DESeq2")
dir.create("Results/PICRUST2_LEFSE_results")
}

options(scipen = 100) # disable scientific annotation


####################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample
sample<-substring(sample, first = 9, last= 20) # strtrim to cut until X
sample_names(data)<-sample

Metadata <- as.data.frame(read.csv(file="Metadata.csv"))

row.names(Metadata)<-Metadata$FASTQ_ID
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
original_length<- original_length - 1 # the sample AD6 has not been sequenced successfully

Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,original_names,sample)

sample_names(data)<-sample_data(data)$Sample
sample_data(data)$Mutation_Treatment<-factor(sample_data(data)$Mutation_Treatment, levels=c("WT_Control","Tg2576_Control","WT_DSS","Tg2576_DSS")) # decides the order in plots
sample_data(data)$Mutation<-factor(sample_data(data)$Mutation, levels=c("WT","Tg2576"))
sample_data(data)$Treatment<-factor(sample_data(data)$Treatment, levels=c("Control","DSS"))

# save.image("data.RData")
write.table(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Data_check/Feature_table_pre_filters.tsv",quote=F, sep="\t")


#################### FILTERING NOISES FROM DATA SET ####################

unfiltered_data<-data

###### cutting under 0.01% to remove noises/contaminants, too conservative but also safe cutoff, see  DOI: 10.1128/mSystems.00290-19
data.genus.temp<-tax_glom(data, taxrank = "Genus")
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.01, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ),
           file="Data_check/Filtered_genus_under_001_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) max(x) <= 0.01, TRUE)))[["Genus"]]
data<-subset_taxa(data, ! Genus %in% filtered)
rm(filtered, data.genus.temp)


######## Checking unassigned in prokaryote kingdom
{Unass<-tax_glom(data, taxrank = "Kingdom") # or domain
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


######################### CHECKING THE OUTLIERS (AD 42) ###############################

dir.create("Data_check/AD42_Outlier")

if("AD42" %in% sample_names(data)){
  data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
} else { cat (" \nLooks like the outlier AD42 is no more in your object already \n ")}

#### Exporting every abundance
dir.create("Data_check/AD42_Outlier/Relative_abundances")
{ write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Data_check/AD42_Outlier/Relative_abundances/counts_genus.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Data_check/AD42_Outlier/Relative_abundances/counts_genus.csv",quote=F)
}

#### Bar plot TOP 8 Genera
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + 
  facet_grid2( Mutation ~ Treatment, scales = "free_x", independent = "x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust = 1, size = 9), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) +
  labs(x="Sample", y="Relative abundance", title = "Eight most abundant genera (with AD42)", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Data_check/AD42_Outlier/TOP_8_Genera.png",width=9,height=7.5,dpi=300)
dev.off()

#### PCoA on ASV
data.prop.labels<-data.prop
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","Tg2576 Control"="coral","WT DSS"="seagreen4","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Data_check/AD42_Outlier/PCoA_Beta_diversity_Bray_Curtis_on_ASV.png", width = 8, height = 6, dpi=300)

#### PCoA ASV in 3D
{Dist<- vegdist(t(otu_table(data.sqrt_prop)), method = "bray")                                                                                                      
  obj<-ecodist::pco(Dist)
  matrix<-obj[["vectors"]]
  row.names(matrix)<-sample_names(data.sqrt_prop)
}
{color_table<-as.data.frame(cbind(as.character(sample_data(data.sqrt_prop)$Mutation_Treatment),as.character(sample_data(data.sqrt_prop)$Mutation_Treatment)))
  colnames(color_table)<-c("Gruppo","Colore")
  unique(color_table$Gruppo)
  colors <- gsub("WT Control","3",color_table$Colore) #green
  colors <- gsub("Tg2576 Control","#FF7F50",colors) # coral
  colors <- gsub("WT DSS","#008B00",colors) # dark green
  colors <- gsub("Tg2576 DSS","#CD2626",colors) # dark red
  shape_table<-as.data.frame(cbind(as.character(sample_data(data.sqrt_prop)$Mutation_Treatment),as.character(sample_data(data.sqrt_prop)$Mutation_Treatment)))
  colnames(shape_table)<-c("Gruppo","shape")
  shapes <- gsub("WT Control","sphere",shape_table$shape)
  shapes <- gsub("Tg2576 Control","sphere",shapes)
  shapes <- gsub("WT DSS","o",shapes)
  shapes <- gsub("Tg2576 DSS","o",shapes)
}
rgl::open3d(windowRect=c(25,25,1200,1200))
pca3d::pca3d(matrix, col=colors, shape = shapes, radius=1.5, 
             show.shadows = T, show.plane = F, show.centroids = F)
rgl::rgl.viewpoint(theta = 25.8, phi = 20.8, fov = 110, zoom = 0.38)
rgl::legend3d("topleft", c("WT","Tg2576", "WT+DSS", "Tg2576+DSS"), 
              col=c(3,"#FF7F50","#008B00","#CD2626"), pch = c(19, 19, 17, 17), magnify=0.7 )
rgl::rgl.snapshot("temp.png")
temp<-png::readPNG("temp.png")
# install.packages("pdftools") # needs also sudo apt install libpoppler-cpp-dev
png::writePNG(temp, "Data_check/AD42_Outlier/PCoA_Bray_Curtis_ASV_3D_Mutation_Treatment.png", dpi = 300)
unlink("temp.png")
rgl::close3d()

#### PCoA on genera
{ data.prop.labels<-data.genus.prop
  sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","Tg2576 Control"="coral","WT DSS"="seagreen4","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional genera", color="Mutation_Treatment", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Data_check/AD42_Outlier/PCoA_Beta_diversity_Bray_Curtis_on_genera.png", width = 8, height = 6, dpi=300)


#### it's better to remove that sample
data<-subset_samples(data, Sample != "AD42")

rm(data.genus, data.prop, tabella_others, top, tabella_others, tabella_top, prune.dat_top, prune.data.others, data.prop.labels, colors, shapes, color_table, shape_table, Dist, obj, matrix, fill_color_5, BC.dist)


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
r<-rarecurve(t(otu_table(data)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)


####################### PREPARATION OF THE DATA #######################

if("AD42" %in% sample_names(data)){
  stop (" \nLooks like the outlier AD42 is still in your object! Remove that sample \n ")
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

rm(con, filter)

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

######################## ABUNDANCES BAR PLOT ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one

# TOP 5 Phyla
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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  facet_grid2( Mutation ~ Treatment, scales = "free_x", space = "free_x", independent = "x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust = 1, size = 9), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Sample", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/TOP_5_phyla.png",width=9,height=7, dpi=300)
dev.off()

rm(top, prune.dat_top,tabella, tabella_top)

# TOP 5 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + 
  facet_grid2( Mutation ~ Treatment, scales = "free_x", independent = "x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust = 1, size = 9), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Sample", y="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Abundances/TOP_5_Genera.png",width=9,height=7,dpi=300)
dev.off()

# TOP 8 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
  facet_grid2( Mutation ~ Treatment, scales = "free_x", independent = "x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust = 1, size = 9), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) +
  labs(x="Sample", y="Relative abundance", title = "Eight most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Abundances/TOP_8_Generi.png",width=9,height=7,dpi=300)
dev.off()

rm(top, prune.dat_top,tabella, tabella_top)

######################## HIERARCHICAL CLUSTERING ###################

#euclidean
{c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
  c<-as.dendrogram(c)
  color_table<-as.data.frame(cbind(as.character(sample_data(data)$Mutation_Treatment),as.character(sample_data(data)$Mutation_Treatment)))
  colnames(color_table)<-c("Gruppo","Colore")
  colors <- gsub("WT_Control","chartreuse",color_table$Colore) #green
  colors <- gsub("Tg2576_Control","coral",colors) #orange
  colors <- gsub("WT_DSS","seagreen4",colors) 
  colors <- gsub("Tg2576_DSS","red3",colors)
  labels_colors(c) <- colors[order.dendrogram(c)] # but sort them based on their order in dendrogram
}
png(file="Results/Hierarchical_clustering/Hierarchical_cluster_Euclidean_sqrt_prop_ASVs.png",width=2500,height=1800, res=300)
par(mar=c(8,4,4,3), cex.lab=0.4, cex.main=1.4, cex.sub=1.1)
plot(c, main="Community structure computed with Euclidean distance \n on sqrt proportional ASVs", 
     sub="WT Control = Green    Tg2576 Control = Orange    WT DSS = Dark Green     Tg2576 DSS = red")
dev.off()

# Bray Curtis
{c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method="bray"))
  c<-as.dendrogram(c)
  color_table<-as.data.frame(cbind(as.character(sample_data(data)$Mutation_Treatment),as.character(sample_data(data)$Mutation_Treatment)))
  colnames(color_table)<-c("Gruppo","Colore")
  colors <- gsub("WT_Control","chartreuse",color_table$Colore) #green
  colors <- gsub("Tg2576_Control","coral",colors) #orange
  colors <- gsub("WT_DSS","seagreen4",colors) 
  colors <- gsub("Tg2576_DSS","red3",colors)
  labels_colors(c) <- colors[order.dendrogram(c)] # but sort them based on their order in dendrogram
}
png(file="Results/Hierarchical_clustering/Hierarchical_cluster_Bray_Curtis_sqrt_prop_ASVs.png",width=2500,height=1800, res=300)
par(mar=c(8,4,3,4), cex.lab=0.4, cex.main=1.4, cex.sub=1.1)
plot(c, cex.lab= 0.2, main="Community structure computed with Bray-Curtis distance \n on sqrt proportional ASVs", 
     sub="WT Control = Green    Tg2576 Control = Orange    WT DSS = Dark Green     Tg2576 DSS = red")
dev.off()


########################## ALFA DIVERSITY ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

{pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Mutation_Treatment")
  pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evanness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
}
pAlpha$data$Mutation_Treatment<-factor(pAlpha$data$Mutation_Treatment, 
                                       levels = c("WT_Control","WT_DSS", "Tg2576_Control", "Tg2576_DSS"))
pAlpha + 
  geom_boxplot(data=pAlpha$data, 
                      aes(x=Mutation_Treatment, y=value, color=Mutation_Treatment), 
                      alpha=0.1, size=0.4) +
  scale_color_manual(values=c("WT_Control"="chartreuse","Tg2576_Control"="coral","WT_DSS"="seagreen4","Tg2576_DSS"="red3")) +
  theme_bw() + 
  labs(x="", title="Alpha diversity between WT, WT+DSS, Tg2576 and Tg2576+DSS groups") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=35, vjust=1, hjust=1, size=11),
        strip.text = element_text(size=11)) +
  stat_compare_means(aes(group = Mutation_Treatment), label="p.format", method = "kruskal.test", label.x= 1.78, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="Results/Alpha_diversity/Alfa_diversity_with_Kruskal.png", width = 9,height =7, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data)
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Mutation_Treatment
kruskal.test(Obser_value$value~factor)

# pairwise comparisons
con <- file("Results/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_Mutation_Treatment.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Mutation_Treatment sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-filter(alphadt, variable=="Observed")
a<-FSA::dunnTest(value ~ Mutation_Treatment, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Shannon")
a<-FSA::dunnTest(value ~ Mutation_Treatment, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Evenness")
a<-FSA::dunnTest(value ~ Mutation_Treatment, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)
rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor)

######################## BETA DIVERSITY BRAY CURTIS #######################

{ASV.prop<-as.data.frame(otu_table(data.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA
metadata<-as(sample_data(data.prop),"data.frame")
library(vegan)

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Mutation_Treatment, data=metadata, permutations = 9999, method="bray")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_Bray<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Mutation_Treatment, data=metadata, permutations = 9999, method="bray")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_Bray<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Mutation_Treatment, data=metadata, permutations = 9999, method="bray")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Mutation_Treatment, data=metadata, permutations = 9999, method="bray")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Mutation_Treatment, data=metadata, permutations = 9999, method="bray")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Mutation_Treatment, data=metadata, permutations = 9999, method="bray")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=metadata$Mutation_Treatment, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV

# exporting beta diversity
rm(con)
con<-file("Results/Bray_Curtis_Beta_div/Beta_diversity_general_and_pairwise_between_Mutation_Treatment_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
  disper<-vegan::betadisper(BC.dist,metadata$Mutation_Treatment)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
rm(con)
con<-file("Results/Bray_Curtis_Beta_div/Beta_dispersion_General_and_Pairwise_between_Mutation_Treatment_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con)

################## GENERAL PCoA BRAY CURTIS #####################

# based on ASV
data.prop.labels<-data.prop
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","Tg2576 Control"="coral","WT DSS"="seagreen4","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Bray_Curtis_Beta_div/PCoA_Beta_diversity_Bray_Curtis_on_ASV.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","Tg2576 Control"="coral","WT DSS"="seagreen4","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Bray_Curtis_Beta_div/PCoA_Beta_diversity_Bray_Curtis_on_ASV_no_ellipses.png", width = 8, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","Tg2576 Control"="coral","WT DSS"="seagreen4","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Bray_Curtis_Beta_div/PCoA_Beta_diversity_Bray_Curtis_on_ASV_points.png", width = 8, height = 6, dpi=300)

# based on ASV (ONLY WT)
data.prop.labels<-data.prop
data.prop.labels<-subset_samples(data.prop.labels, Mutation=="WT")
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","WT DSS"="seagreen4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Treatment",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Bray_Curtis_Beta_div/PCoA_WT_Treatment_ASV.png", width = 8, height = 6, dpi=300)


# based on ASV (ONLY Tg2576)
data.prop.labels<-data.prop
data.prop.labels<-subset_samples(data.prop.labels, Mutation=="Tg2576")
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("Tg2576 Control"="coral","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Treatment",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Bray_Curtis_Beta_div/PCoA_Tg2576_Treatment_ASV.png", width = 8, height = 6, dpi=300)


# ASV in 3D
{Dist<- vegdist(t(otu_table(data.sqrt_prop)), method = "bray")                                                                                                      
  obj<-ecodist::pco(Dist)
  matrix<-obj[["vectors"]]
  row.names(matrix)<-sample_names(data.sqrt_prop)
}
{color_table<-as.data.frame(cbind(as.character(sample_data(data.sqrt_prop)$Mutation_Treatment),as.character(sample_data(data.sqrt_prop)$Mutation_Treatment)))
colnames(color_table)<-c("Gruppo","Colore")
unique(color_table$Gruppo)
colors <- gsub("WT Control","3",color_table$Colore) #green
colors <- gsub("Tg2576 Control","#FF7F50",colors) # coral
colors <- gsub("WT DSS","#008B00",colors) # dark green
colors <- gsub("Tg2576 DSS","#CD2626",colors) # dark red
shape_table<-as.data.frame(cbind(as.character(sample_data(data.sqrt_prop)$Mutation_Treatment),as.character(sample_data(data.sqrt_prop)$Mutation_Treatment)))
colnames(shape_table)<-c("Gruppo","shape")
shapes <- gsub("WT Control","sphere",shape_table$shape)
shapes <- gsub("Tg2576 Control","sphere",shapes)
shapes <- gsub("WT DSS","o",shapes)
shapes <- gsub("Tg2576 DSS","o",shapes)
}
rgl::open3d(windowRect=c(25,25,1200,1200))
pca3d::pca3d(matrix, col=colors, shape = shapes, radius=1.5, 
             show.shadows = T, show.plane = F, show.centroids = F)
rgl::rgl.viewpoint(theta = -17, phi = 20.8, fov = 120, zoom = 0.331)
rgl::legend3d("topleft", c("WT","Tg2576", "WT+DSS", "Tg2576+DSS"), 
              col=c(3,"#FF7F50","#008B00","#CD2626"), pch = c(19, 19, 17, 17), magnify=0.7 )
rgl::rgl.snapshot("temp.png")
temp<-png::readPNG("temp.png")
# install.packages("pdftools") # needs also sudo apt install libpoppler-cpp-dev
png::writePNG(temp, "Results/Bray_Curtis_Beta_div/PCoA_Bray_Curtis_ASV_3D_Mutation_Treatment.png", dpi = 300)
unlink("temp.png")
rgl::close3d()

# again but on genera
{ data.prop.labels<-data.genus.prop
  sample_data(data.prop.labels)$Mutation_Treatment<-gsub("_"," ",sample_data(data.prop.labels)$Mutation_Treatment)
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT Control"="chartreuse","Tg2576 Control"="coral","WT DSS"="seagreen4","Tg2576 DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional genera", color="Mutation_Treatment", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Bray_Curtis_Beta_div/PCoA_Beta_diversity_Bray_Curtis_on_genera.png", width = 8, height = 6, dpi=300)

################ BETA DIVERSITY wUNIFRAC ##############

#### PERMANOVA

metadata<-as(sample_data(data.prop),"data.frame")

rm(perm_ASV, DistBC,data.sqrt_prop_perm)

{data.sqrt_prop_perm<-transform_sample_counts(data.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac") # needing phyloseq to compute unifrac
perm_ASV<- vegan::adonis(DistBC ~Mutation_Treatment, data=metadata, permutations = 9999)
}

{data.sqrt_prop_perm<-transform_sample_counts(data.genus.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_g<- vegan::adonis(DistBC ~Mutation_Treatment, data=metadata, permutations = 9999)
}

{data.sqrt_prop_perm<-transform_sample_counts(data.fam.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_f<- vegan::adonis(DistBC ~Mutation_Treatment, data=metadata, permutations = 9999)
}

{data.sqrt_prop_perm<-transform_sample_counts(data.order.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_o<- vegan::adonis(DistBC ~Mutation_Treatment, data=metadata, permutations = 9999)
}

{data.sqrt_prop_perm<-transform_sample_counts(data.class.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_c<- vegan::adonis(DistBC ~Mutation_Treatment, data=metadata, permutations = 9999)
}

{data.sqrt_prop_perm<-transform_sample_counts(data.phy.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
perm_p<- vegan::adonis(DistBC ~Mutation_Treatment, data=metadata, permutations = 9999)
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

# pairwise beta diversity
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
data.sqrt_prop_perm<-transform_sample_counts(data.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "wunifrac")
pair_ASV<- pairwise.adonis(DistBC, factors=metadata$Mutation_Treatment, perm = 9999)
pair_ASV

# exporting beta diversity
rm(con)
con<-file("Results/w_Unifrac_Beta_div/Beta_diversity_general_and_pairwise_between_Mutation_Treatment_groups")
sink(con, append = TRUE)
cat("General beta diversity wUnifrac \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity wUnifrac (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, con)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on OTU
{data.sqrt_prop_perm<-transform_sample_counts(data.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop_perm, method = "unifrac")
disper<-vegan::betadisper(DistBC,metadata$Mutation_Treatment)
disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
disp_ASV$tab
a<-as.data.frame(disp_ASV$pairwise$permuted)
colnames(a)<-c("permuted_p_value")
a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
rm(con)
con<-file("Results/w_Unifrac_Beta_div/Beta_dispersion_general_and_pairwise_between_Mutation_Treatment_groups")
sink(con, append=TRUE)
cat("General beta dispersion wUnifrac \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con, data.sqrt_prop_perm, DistBC, perm_ASV)


###### PCoA weighted Unifrac

data.prop.labels<-data.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels)) 
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
DistBC = phyloseq::distance(data.sqrt_prop, method = "wunifrac")
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT_Control"="chartreuse","Tg2576_Control"="coral","WT_DSS"="seagreen4","Tg2576_DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with weighted UniFrac distance \n computed on sqrt proportional ASV", color="Mutation_Treatment", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/w_Unifrac_Beta_div/PCoA_Beta_diversity_wUnifrac.png", width = 8, height = 6, dpi=300)

rm(data.sqrt_prop, eigval, DistBC, ordBC,data.prop.labels)

##### PCoA unweighted Unifrac

data.prop.labels<-data.prop # to change labels
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt)
DistBC = phyloseq::distance(data.sqrt_prop, method = "unifrac") # not weighted
ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues 
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Mutation_Treatment") +
  scale_color_manual(values=c("WT_Control"="chartreuse","Tg2576_Control"="coral","WT_DSS"="seagreen4","Tg2576_DSS"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with unweighted UniFrac distance \n computed on sqrt proportional ASV ", color="Mutation_Treatment", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/w_Unifrac_Beta_div/PCoA_Beta_diversity_unweighted_Unifrac.png", width = 8, height = 6, dpi=300)

rm(data.sqrt_prop, eigval, DistBC, ordBC,data.prop.labels)


##################### DA WITH DESEQ2 #######################

# Trimming under sum of 10, see DeSeq2 tutorial
data_pruned<- prune_taxa(taxa_sums(data) > 10, data)

# preparing new taxa vocabularies for plots (some ASV may change glomming after trimming)
{data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
{taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
  tax_table(data_pruned.genus.prop)<-as.matrix(taxa_temp)
  # adding informations to uncultured families
  taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
  for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
    taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  tax_table(data_pruned.fam.prop)<-as.matrix(taxa_temp)
  rm(taxa_temp)
}

# function to automatize the comparison between groups
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Mutation_Treatment)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Mutation_Treatment)
    }
    # DE<-DESeq(d, test="LRT",reduced= ~ 1)
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Mutation_Treatment", nume[i], deno[i]), test="Wald")
      #res<-results(DE, contrast=c("Mutation_Treatment", nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      #res<-lfcShrink(DE, type="normal", contrast=c("Mutation_Treatment",nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        #			cat(paste0(rank,":",nume[i],"/",dena[i],"\n"),file=outfile,append=T)
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("Tg2576_Control","Tg2576_DSS","Tg2576_Control","WT_Control")
deno<-c("WT_Control","WT_DSS","Tg2576_DSS","WT_DSS")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
unlink("DE_results.txt")
DE_results<-DE_results[DE_results$BaseMean>50, ]
system(" echo 'Every result under the arbitrary threshold Basemean=50 has been removed in order to avoid the noisiest results' > Results/DA_DESeq2/NB_results_are_filtered.txt")
write.csv2(DE_results, file="Results/DA_DESeq2/DE_every_results.csv", quote=F, row.names = F)
write.xlsx(DE_results, file="Results/DA_DESeq2/DE_every_results.xlsx",showNA = F, col.names = T, row.names = F)

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
  target<-prune_taxa(target, data_pruned.genus.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Family","ASV"]
  target<-prune_taxa(target, data_pruned.fam.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Order","ASV"]
  target<-prune_taxa(target, data_pruned.order.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Class","ASV"]
  target<-prune_taxa(target, data_pruned.class.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
  target<-prune_taxa(target, data_pruned.phy.prop)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
  rm(tabella)
}

# unique boxplot of DeSeq2 results
{tabella_g$Taxa<-"Genera"
  tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
  colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"
}
{tabella_f$Taxa<-"Families"
  tabella_f[,c("Phylum","Order","Class")]<-NULL
  colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"
}
{tabella_o$Taxa<-"Orders"
  tabella_o[,c("Phylum","Class")]<-NULL
  colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"
}
{tabella_c$Taxa<-"Classes"
  tabella_c[,"Phylum"]<-NULL
  colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"
}
{tabella_p$Taxa<-"Phyla"
  colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"
}

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_c,tabella_o,tabella_p)
tabella_tot$Mutation_Treatment<-gsub("_Control","",tabella_tot$Mutation_Treatment)
tabella_tot$Mutation_Treatment<-factor(tabella_tot$Mutation_Treatment,levels = c("WT","Tg2576","WT_DSS","Tg2576_DSS"))

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Mutation_Treatment)) + 
  facet_grid2(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free", space="free") +
  scale_fill_manual(values=c("WT"="chartreuse","Tg2576"="coral","WT_DSS"="seagreen4","Tg2576_DSS"="red3")) +
  geom_boxplot(width=0.8) + theme_bw( ) + 
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  theme(axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=9),
        legend.margin=margin(-20, 0, 0, 0), legend.position="bottom") +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Mutation_Treatment", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,12,2),seq(12,90,3)))
ggsave(filename = "Results/DA_DESeq2/Every_Results_and_Y_is_zoomed_in.png", width = 12, height = 8, dpi=300)
dev.off()

# to cut off REDUNDANTS
Redund<-c("Lactobacillaceae","Lactobacillales","Cyanobacteria")
tabella_tot2<-subset(tabella_tot, ! Bacteria %in% Redund)
tabella_tot2<-tabella_tot2[!tabella_tot2$Taxa %in% c("Families","Orders","Classes"),]

ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Mutation_Treatment)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("WT"="chartreuse","Tg2576"="coral","WT_DSS"="seagreen4","Tg2576_DSS"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(axis.text.x = element_text(angle = 20, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=9),
        legend.margin=margin(-25, 0, 0, 0), legend.position="bottom") + 
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Mutation_Treatment", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,12,2),seq(12,90,3)))
ggsave(filename = "Results/DA_DESeq2/Every_Results_no_redundants_and_Y_is_zoomed_in.png", width = 12, height = 7.5, dpi=300)
dev.off()

################# now plotting result in pairwise (group vs group)

##### Tg2576_Control vs WT_Control
tabella_tot3<-subset(tabella_tot2, Mutation_Treatment != "WT_DSS")
tabella_tot3<-subset(tabella_tot3, Mutation_Treatment != "Tg2576_DSS")
WT_Control_vs_Tg2576_Control<-DE_results[DE_results$num=="Tg2576_Control" & DE_results$denom=="WT_Control", ]
#NB: the object DE_results has the default group names (with "Control")!

# anything there!

#####  WT_Control vs WT_DSS
tabella_tot3<-subset(tabella_tot2, Mutation_Treatment != "Tg2576")
tabella_tot3<-subset(tabella_tot3, Mutation_Treatment != "Tg2576_DSS")
WT_DSS_vs_WT_Control<-DE_results[DE_results$num=="WT_Control" & DE_results$denom=="WT_DSS", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% WT_DSS_vs_WT_Control [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, WT_DSS_vs_WT_Control[WT_DSS_vs_WT_Control$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)
tabella_tot3[tabella_tot3$Bacteria=="Clostridia_UCG-014","Taxa"]<-"Families"
# the family Clostridia has the same name of the Genus Clostridia then here is corrected as it should be (see Results)

tabella_tot3$Mutation_Treatment<-gsub("_Control","",tabella_tot3$Mutation_Treatment)
ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Mutation_Treatment)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("WT"="chartreuse","WT_DSS"="seagreen4")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=13.5),
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between WT (Control) and WT DSS", y="Proportional Abundance", fill="Mutation_Treatment", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5,0.1,1,seq(2,30,2)))
ggsave(filename = "Results/DA_DESeq2/Results_only_WT_Control_vs_WT_DSS_DeSeq2_no redundants.png", width = 10, height = 8, dpi=300)
dev.off()



#####  Tg2576_Control vs Tg2576_DSS
tabella_tot3<-subset(tabella_tot2, Mutation_Treatment != "WT")
tabella_tot3<-subset(tabella_tot3, Mutation_Treatment != "WT_DSS")
Tg2576_DSS_vs_Tg2576_Control<-DE_results[DE_results$num=="Tg2576_Control" & DE_results$denom=="Tg2576_DSS", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% Tg2576_DSS_vs_Tg2576_Control [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, Tg2576_DSS_vs_Tg2576_Control[Tg2576_DSS_vs_Tg2576_Control$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

tabella_tot3$Mutation_Treatment<-gsub("_Control","",tabella_tot3$Mutation_Treatment)
ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Mutation_Treatment)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("Tg2576"="coral","Tg2576_DSS"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=13.5),
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between Tg2576 (Control) and Tg2576 DSS", y="Proportional Abundance", fill="Mutation_Treatment", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(seq(40,90,2)))
ggsave(filename = "Results/DA_DESeq2/Results_only_Tg2576_Control_vs_Tg2576_DSS_DeSeq2_no redundants.png", width = 10, height = 8, dpi=300)

#####  WT_DSS vs Tg2576_DSS
tabella_tot3<-subset(tabella_tot2, Mutation_Treatment != "WT")
tabella_tot3<-subset(tabella_tot3, Mutation_Treatment != "Tg2576")
DSS_vs_DSS<-DE_results[DE_results$num=="Tg2576_DSS" & DE_results$denom=="WT_DSS", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% DSS_vs_DSS [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, DSS_vs_DSS[DSS_vs_DSS$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)
tabella_tot3[tabella_tot3$Bacteria=="Clostridia_UCG-014","Taxa"]<-"Families"
# the family Clostridia has the same name of the Genus Clostridia then here is corrected as it should be (see Results)

tabella_tot3$Mutation_Treatment<-gsub("_Control","",tabella_tot3$Mutation_Treatment)
ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Mutation_Treatment)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("WT_DSS"="seagreen4","Tg2576_DSS"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=13.5),
        axis.text.y = element_text(size=11.8)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between WT DSS and Tg2576 DSS", y="Proportional Abundance", fill="Mutation_Treatment", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,12,2),seq(12,45,3),seq(45,87,4)))
ggsave(filename = "Results/DA_DESeq2/Results_only_WT_DSS_vs_Tg2576_DSS_DeSeq2_no redundants.png", width = 10, height = 8, dpi=300)
dev.off()

system(" echo 'The following comparison has been done: WT vs Tg, WT_DSS vs Tg_DSS, WT vs WT DSS, Tg vs Tg DSS' >Results/DA_DESeq2/About_the_comparison_performed.txt")
system(" echo '\n\nEvery result under the arbitrary threshold Basemean=50 has been removed in order to avoid the noisiest results' >> Results/DA_DESeq2/About_the_comparison_performed.txt")


############# PLOTTING THE LEFSE RESULT ##################

a <- read.delim("QIIME/PICRUST2_LEFSE/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("QIIME/PICRUST2_LEFSE/Result_LEFSE.res", header=FALSE)
head(Significative_functions_LEFSE, n=2)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=2)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="Results/PICRUST2_LEFSE_results/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)

# to cut names too long
Significative_functions_LEFSE$Pathway<-gsub("","", Significative_functions_LEFSE$Pathway)
# plotting results
{ Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ]
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical sorting
}
Significative_functions_LEFSE$Class_with_highest_mean<-gsub("-Control","",Significative_functions_LEFSE$Class_with_highest_mean)

# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = 0)+
  scale_fill_manual(values=c("WT"="chartreuse",
                             "Tg2576"="coral",
                             "WT-DSS"="seagreen4",
                             "Tg2576-DSS"="red2")) +
  theme_classic(base_size = 15) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_every_result.png", width = 8, height = 12, dpi = 300)
# plotting only results over 3
Significative_functions_LEFSE_3<-Significative_functions_LEFSE[abs(Significative_functions_LEFSE$logLDA_score)>3,]
ggplot(data=Significative_functions_LEFSE_3, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = 0, size=3)+
  scale_fill_manual(values=c("WT"="chartreuse",
                             "Tg2576"="coral",
                             "WT-DSS"="seagreen4",
                             "Tg2576-DSS"="red2")) +
  theme_classic(base_size = 15) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_only_over_3.png", width = 7, height = 5, dpi = 300)

system("echo 'Those are general results of a Kruskal Wallis filtered through a subsequent Wilcoxon test: \neach pathway (of a group) is differently abundant compared to every other group! \n\nI would suggest to focus only on results over 3 logLDAscore.' > Results/PICRUST2_LEFSE_results/NB.txt")

rm(Significative_functions_LEFSE, a)

##################### R and PACKAGES VERSION #########################

package<-sessionInfo()

con <- file("Results/R_version_and_packages.txt") 
sink(con, append = TRUE)

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
print("Dendextend")
packageVersion("dendextend")
cat("\n", "\n", fill=TRUE)
package$otherPkgs$pca3d[c(1,4)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$mixOmics[c(1,4)]

sink()
close(con)
rm(con)