################### PREPARING THE ENVIRONMENT ##############

{library("phyloseq")
library("dendextend")
library("ggplot2")
library("ggh4x")
library("ggpubr")
library("vegan")
library("Hmisc")
library("xlsx")
library("DESeq2")
library("egg")
}

{data<- qiime2R::qza_to_phyloseq(features="QIIME2/table.qza", taxonomy="QIIME2/assigned_taxa.qza")
metadata<- read.delim("metadata_R.tsv", sep="\t")
row.names(metadata)<-metadata$FASTQ
metadata<-metadata[sample_names(data),]
identical(as.numeric(length(which(is.na(rownames(metadata))))),0) # TRUE --> OK
original_length<-length(sample_names(data))
sample_data(data)<-metadata
identical(original_length, length(sample_names(data))) # TRUE, no sample cutted out
sample_names(data)<-sample_data(data)$FASTQ_code
rm(original_length)
}

data_total<-data # original object for subsetting
if(length(unique(sample_data(data_total)$Sampling_Type))!=3){
  cat("\nThere are not three sampling types in the object, please restart the whole script! \n")
  Sys.sleep(2)
} else { rm(data) }

head(sample_data(data_total), n=2)

# <- data saved here
# save.image("data_2022.RData")

dir.create("Results")
setwd("Results")


### function from https://github.com/matteoramazzotti
CircoTax2=function(file,title="CircoTax plot",ramp=c("orange","white","blue"),tax_col=7:11,fc_col=2,sort=c("no","rank","fc","absfc","alpha"),sort_dir="d") {
  data=file
  sort_dir=ifelse(sort_dir == "d",FALSE,TRUE)
  if(length(tax_col) == 5) {
    gplot_labels= unlist(strsplit("PCOFG",""))
  }
  if(length(tax_col) == 6) { #KPCOFG as in RDP
    gplot_labels= unlist(strsplit("KPCOFG",""))
  }
  if(length(tax_col) == 7) { #KPCOFGS as in silva
    gplot_labels= unlist(strsplit("KPCOFGS",""))
  }
  #build the taxa-related variables (tax index and label)
  #ranks are assumed to be in decreasing order from kinkdom(domain) to species
  #y represent the height (from the center of the circle) of the bars => domain=1 (min) species=7 (max)
  y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
  if(sort[1] == "no") {
    #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y
    fc=data[,fc_col]
  }
  if(sort[1] == "rank") {
    #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    o=order(y,decreasing=sort_dir)
    synth=apply(data[o,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "alphalin") {
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    o=order(synth,decreasing=sort_dir)
    synth=synth[o]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "alpha") {
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    o=order(labels)
    labels=labels[o]
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "fc") {
    o=order(data[,fc_col],decreasing=sort_dir)
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    synth=synth[o]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "absfc") {
    #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    o=order(abs(data[,fc_col]),decreasing=sort_dir)
    synth=apply(data[o,tax_col],1,function(x) paste0(x,collapse="-"))
    synth=synth[order(abs(data[,fc_col]),decreasing=sort_dir)]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    fc=data[o,fc_col]
  }
  
  #builds the data.frame for ggplot2 
  df=data.frame("id"=1:dim(data)[1],"name"=labels,"rank"=y,"FC"=fc)
  
  #adds the label angle column
  nbar=dim(df)[1]
  angle=90-360* ((1:nbar)-0.5) / nbar
  angle<-ifelse(angle < -90, angle+180, angle)
  df$angle=angle
  
  #the plot starts here
  ggplot(df, aes(x = id, y = rank)) +
    ggtitle(title) +
    geom_col(
      aes(fill=FC),
      position = "dodge"
    ) + 
    scale_fill_gradient2(
      low="blue",
      mid="white",
      high="orange"
    ) +
    annotate(
      "text",label=gplot_labels, x=rep(0,length(tax_col)), y=1:length(tax_col),size=4.5, vjust=0.58
    ) +
    geom_text(
      aes(
        x= id, 
        y=7, 
        label=name,
      ), 
      color="black", 
      fontface="bold",
      alpha=0.6,
      size=3,
      angle=angle,
    ) +
    geom_hline( 
      yintercept=1:6,
      color="grey"
    ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = NA),
      plot.margin = unit(rep(1,4), "cm"),      # Adjust the margin to make in sort labels are not truncated!
      plot.title = element_text(hjust = 0.5,face="bold", size=18)
    ) +
    coord_polar(start = 0, clip="off") +
    labs(fill="log2FC")
}


##################### EXPORTING GENERAL ABUNDANCES  ######################

{data.genus<- tax_glom(data_total, taxrank = "Genus", NArm = F)
data.fam<- tax_glom(data_total, taxrank = "Family", NArm = F)
data.class<- tax_glom(data_total, taxrank = "Class", NArm = F)
data.order<- tax_glom(data_total, taxrank = "Order", NArm = F)
data.phy<- tax_glom(data_total, taxrank = "Phylum", NArm = F)
}
{data.prop_total <- transform_sample_counts(data_total, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}

dir.create("Every_raw_abundances")
{write.csv2(cbind(as(otu_table(data_total),"matrix"),as(tax_table(data_total),"matrix")),file="Every_raw_abundances/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Every_raw_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Every_raw_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Every_raw_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Every_raw_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Every_raw_abundances/counts_genus.csv",quote=F)
}

options( scipen=100 )
dir.create("Every_percent_abundances")
{write.csv2(cbind(as(otu_table(data.prop_total),"matrix"),as(tax_table(data.prop_total),"matrix")),file="Every_percent_abundances/counts_otu.csv",quote=F)
write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Every_percent_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class.prop),"matrix")),file="Every_percent_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order.prop),"matrix")),file="Every_percent_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam.prop),"matrix")),file="Every_percent_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Every_percent_abundances/counts_genus.csv",quote=F)
}

###################### RAREFACTION ANALYSIS ######################### 

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as generated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of he original samples (and in the same order!!!)
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    check<-"red"
    if(length(x) < 10) {
      points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col="cyan",pch=16,cex=1)
    } else {
      slope<-mean(diff(v[(length(v)-t):length(v)])/dx)
      if(slope<lim) {
        check<-"blue"
        sat = sat+1
      }
      cat(i,slope,check,"\n")
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  legend("bottomright",paste(sat,"saturated samples"),bty="n")}

png(file="Sample_rarefaction.png",width=4500,height=3200, res=300)
r<-rarecurve(t(otu_table(data_total)), step=100,label=F)
evalslopes(r,sample_names(data_total),lim=0.001,cex=1)
dev.off()
rm(r,evalslopes)


######################## GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(data_total)==1, data_total)
length(which(taxa_sums(data_total)==1)) # if zero there are no singletons
{n<-as.data.frame(otu_table(filter))
  N<-as.data.frame(otu_table(data_total))
  G<-1-(colSums(n)/colSums(N))
}
con<-file("Percentuale_di_singletons_Good's_coverage.txt")
sink(con)
cat("GOOD'S COVERAGE ESTIMATOR \n", fill=TRUE)
cat("1-(n/N) for each sample, where n is number of singletons and N is Total ASV \n \n", fill=TRUE)
G
sink()
close(con)

rm(con, filter, G, n, N)


########################### GENERAL PCoA #########################

dir.create("General_PCoA")

#### on ASV

suppressWarnings(rm("data.sqrt"))
data.sqrt<-transform_sample_counts(data.prop_total, sqrt) # sqrt of proportion
{DistBC = phyloseq::distance(data.sqrt, method = "bray")
ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
# to add the sample number
Sample_number<-paste0(sample_data(data_total)$Patient,sample_data(data_total)$Sampling_Type)
Sample_number
{Sample_number<-unlist(strsplit(Sample_number, split = "-"))
Sample_number<-gsub("CRC","",Sample_number)
Sample_number<-gsub("AP","",Sample_number)
Sample_number<-gsub("B","",Sample_number)
Sample_number<-gsub("F","",Sample_number)
Sample_number_BF<-Sample_number
# the object "Sample_number_BF" is needed to connect just B and F --> points connected if they have the same name --> "saliva" names have to be different
Sample_number<-gsub("S","",Sample_number)
}
Sample_number_BF
Sample_number
plot_ordination(data.sqrt, ordBC, color = "Sampling_Type",shape = "Tumor") +
  geom_point(size=2.5) + 
  geom_line(aes(group=Sample_number_BF), color="black", size=0.05) +
  geom_text(aes(label = Sample_number), color="black",show.legend=F, size=2) +
  theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Sampling Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="General_PCoA/PCoA_Bray_Curtis_prop_sqrt_ASV.png", width = 7, height = 5, dpi=300)

#### again but with a zoom on the interest area

# install.packages("ggforce")
library("ggforce")
plot_ordination(data.sqrt, ordBC, color = "Sampling_Type",shape = "Tumor") +
  geom_point(size=2.5) + geom_line(aes(group=Sample_number_BF), color="black", size=0.05) +
  geom_text(aes(label = Sample_number), color="black",show.legend=F, size=2) + theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", color="Sampling Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")) +# si pu? aggiungere shape
  facet_zoom(xlim= c(0.0,0.3)) # i valori sono per zoom sull'asse X
ggsave(file="General_PCoA/PCoA_Bray_with_zoom_and_lines.png", width=8,height=6, dpi=300)

##### again but on biopsies and stools only

data_subset<-subset_samples(data.prop_total, Sampling_Type!="S")

data.sqrt<-transform_sample_counts(data_subset,sqrt)
{DistBC = phyloseq::distance(data.sqrt, method = "bray")
ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
# to add the sample number (new orders)
Sample_number<-paste0(sample_data(data_subset)$Patient,sample_data(data_subset)$Sampling_Type)
Sample_number
{Sample_number<-unlist(strsplit(Sample_number, split = "-"))
  Sample_number<-gsub("CRC","",Sample_number)
  Sample_number<-gsub("AP","",Sample_number)
  Sample_number<-gsub("B","",Sample_number)
  Sample_number<-gsub("F","",Sample_number)
}
plot_ordination(data.sqrt, ordBC, color = "Sampling_Type", shape = "Tumor") +
  geom_point(size=2.5) + geom_line(aes(group=Sample_number), color="black", size=0.05) +
  geom_text(aes(label = Sample_number), color="black",show.legend=F, size=2) +
  theme_classic() + scale_color_manual(values = c("B"="coral","F"="green")) +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance",
       color="Sampling Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="General_PCoA/PCoA_Bray_sqrt_prop_NO_SALIVA.png", width = 6, height = 6, dpi=300)

#### again but only on adenocarcinoma 

data_subset<-subset_samples(data.prop_total, Tumor=="Adenocarcinoma")
{data.sqrt<-transform_sample_counts(data_subset,sqrt)
DistBC = phyloseq::distance(data.sqrt, method = "bray")
ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
# to add the sample number
Sample_number<-paste0(sample_data(data_subset)$Patient,sample_data(data_subset)$Sampling_Type)
{Sample_number<-unlist(strsplit(Sample_number, split = "-"))
  Sample_number<-gsub("CRC","",Sample_number)
  Sample_number<-gsub("AP","",Sample_number)
  Sample_number<-gsub("B","",Sample_number)
  Sample_number<-gsub("F","",Sample_number)
  Sample_number_BF<-Sample_number
  # the object "Sample_number_BF" is needed to connect just B and F --> points connected if they have the same name --> "saliva" names have to be different
  Sample_number<-gsub("S","",Sample_number)
}
plot_ordination(data.sqrt, ordBC, color = "Sampling_Type") +
  geom_point(size=2.5) + theme_classic() + 
  geom_line(aes(group=Sample_number_BF), color="black", size=0.05) +
  geom_text(aes(label = Sample_number), color="black",show.legend=F, size=2) + 
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only subjects with Adenocarcinoma)",
       color="Sampling Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="General_PCoA/PCoA_Bray_sqrt_prop_ONLY_ADENOCARC.png", width = 6, height = 5, dpi=300)


#### again but only on adenoma 

data_subset<-subset_samples(data.prop_total, Tumor=="Adenoma")
{data.sqrt<-transform_sample_counts(data_subset,sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
# to add the sample number
Sample_number<-paste0(sample_data(data_subset)$Patient,sample_data(data_subset)$Sampling_Type)
{Sample_number<-unlist(strsplit(Sample_number, split = "-"))
  Sample_number<-gsub("CRC","",Sample_number)
  Sample_number<-gsub("AP","",Sample_number)
  Sample_number<-gsub("B","",Sample_number)
  Sample_number<-gsub("F","",Sample_number)
  Sample_number_BF<-Sample_number
  # the object "Sample_number_BF" is needed to connect just B and F --> points connected if they have the same name --> "saliva" names have to be different
  Sample_number<-gsub("S","",Sample_number)
}
plot_ordination(data.sqrt, ordBC, color = "Sampling_Type") +
  geom_point(size=2.5) + theme_classic() + geom_line(aes(group=Sample_number_BF), color="black", size=0.05) +
  geom_text(aes(label = Sample_number), color="black",show.legend=F, size=2) + 
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only subjects with Adenoma)",
       color="Sampling Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="General_PCoA/PCoA_Bray_sqrt_prop_ONLY_ADENOMA.png", width = 6, height = 5, dpi=300)

##################### STARTING STOOL ANALYSIS ##############################

if( grepl("_analysis",getwd()) | grepl("_correlations",getwd()) ){setwd("..")}
dir.create("Stool_analysis")
setwd("Stool_analysis")

{to_remove<-ls()
to_remove<-to_remove[ ! to_remove %in% c("data_total","metadata","evalslopes","CircoTax2","Genera.DESEQ2_CRC_AP_stool","Genera.DESEQ2_CRC_AP_biopsy")]
to_remove
rm(list = to_remove) # reset the workspace
}

data<-subset_samples(data_total,Sampling_Type=="F")
head(sample_data(data))
sample_names(data)<-sample_data(data)$Patient
head(sample_data(data))

data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
OTU.genus<-as.data.frame(otu_table(data.genus))
Taxa.genus<-as.data.frame(tax_table(data.genus))

data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
OTU.fam<-as.data.frame(otu_table(data.fam))
Taxa.fam<-as.data.frame(tax_table(data.fam))

data.class<- tax_glom(data, taxrank = "Class", NArm = F)
OTU.class<-as.data.frame(otu_table(data.class))
Taxa.class<-as.data.frame(tax_table(data.class))

data.order<- tax_glom(data, taxrank = "Order", NArm = F)
OTU.order<-as.data.frame(otu_table(data.order))
Taxa.order<-as.data.frame(tax_table(data.order))

data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
OTU.phy<-as.data.frame(otu_table(data.phy))
Taxa.phy<-as.data.frame(tax_table(data.phy))

{data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}

dir.create("CRC_Stages_Analysis")
data.stages<-subset_samples(data, !is.na(Stadiation)) # only samples with adenocarcinoma, in order to confront stadiations
sample_data(data.stages)$Stadiation<-factor(sample_data(data.stages)$Stadiatio, levels=c("I","II","III"))
{data.stages.prop<-transform_sample_counts(data.stages, function(otu) otu/sum(otu)*100)
data.stages.genus<- tax_glom(data.stages, taxrank = "Genus", NArm = F)
data.stages.genus.prop <- transform_sample_counts(data.stages.genus, function(otu) otu/sum(otu)*100)
data.stages.phy<- tax_glom(data.stages, taxrank = "Phylum", NArm = F)
data.stages.phy.prop <- transform_sample_counts(data.stages.phy, function(otu) otu/sum(otu)*100)
}

#################### EXPORTING ABUNDANCES (STOOL) ############################

dir.create("Abundances")

{options( scipen=100 ) 
dir.create("Abundances/Percent_abundances")
write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Percent_abundances/counts_phylum.xlsx",row.names = F, showNA = F)
write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Percent_abundances/counts_genus.xlsx",row.names = F, showNA = F)
write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Percent_abundances/counts_phylum.csv",quote=F)
write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class.prop),"matrix")),file="Abundances/Percent_abundances/counts_class.csv",quote=F)
write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order.prop),"matrix")),file="Abundances/Percent_abundances/counts_order.csv",quote=F)
write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam.prop),"matrix")),file="Abundances/Percent_abundances/counts_family.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Percent_abundances/counts_genus.csv",quote=F)
}

### TOP 5 Orders
{top5 <- names(sort(taxa_sums(data.order.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.order.prop)
others<-taxa_names(data.order.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.order.prop)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Order<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Order)) + 
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every order below rank 5 ") +
  ggtitle("Five most abundant orders and other taxa")
ggsave(file="Abundances/TOP5_Orders.png",width=6.5,height=4.5, dpi=300)
dev.off()

# TOP 5 Phyla 
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
others<-taxa_names(data.phy.prop)
others<-others[!(others %in% top5)]
prune.data.others<-prune_taxa(others,data.phy.prop)
tabella_top<-psmelt(prune.dat_top5)
tabella_others<-psmelt(prune.data.others)
tabella_others$Phylum<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="Abundances/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.genus.prop)
others<-taxa_names(data.genus.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.genus.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Genus<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Genus)) + 
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every genus below rank 5 ") +
  ggtitle("Five most abundant genera and other taxa")
ggsave(file="Abundances/TOP5_Genera.png",width=6.5,height=4.5, dpi=300)
dev.off()

# TOP 5 Families
{top <- names(sort(taxa_sums(data.fam.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.fam.prop)
others<-taxa_names(data.fam.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.fam.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Family<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Family)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every family below rank 5 ") +
  ggtitle("Five most abundant families and other taxa")
ggsave(file="Abundances/TOP5_families.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Classes
{top <- names(sort(taxa_sums(data.class.prop), decreasing=TRUE))[1:5]
prune.dat_top <- prune_taxa(top,data.class.prop)
others<-taxa_names(data.class.prop)
others<-others[!(others %in% top)]
prune.data.others<-prune_taxa(others,data.class.prop)
tabella_top<-psmelt(prune.dat_top)
tabella_others<-psmelt(prune.data.others)
tabella_others$Class<-"Others"
tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Class)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every class below rank 5 ") +
  ggtitle("Five most abundant classes and other taxa")
ggsave(file="Abundances/TOP5_classes.png",width=6.5,height=4.5,dpi=300)
dev.off()

dir.create("CRC_Stages_Analysis/Abundances")

write.xlsx(cbind(as(otu_table(data.stages.phy.prop),"matrix"),as(tax_table(data.stages.phy.prop),"matrix")),file="CRC_Stages_Analysis/Abundances/Percent_counts_phylum.xlsx",row.names = F, showNA = F)
write.xlsx(cbind(as(otu_table(data.stages.genus.prop),"matrix"),as(tax_table(data.stages.genus.prop),"matrix")),file="CRC_Stages_Analysis/Abundances/Percent_counts_genus.xlsx",row.names = F, showNA = F)

# TOP 5 Phyla (Only CRC Stages)
{top5 <- names(sort(taxa_sums(data.stages.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.stages.phy.prop)
  others<-taxa_names(data.stages.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.stages.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Stadiation), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="CRC_Stages_Analysis/Abundances/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Genera (Only CRC Stages)
{top <- names(sort(taxa_sums(data.stages.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.stages.genus.prop)
  others<-taxa_names(data.stages.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.stages.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Genus)) + 
  facet_grid2(cols= vars(Stadiation), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every genus below rank 5 ") +
  ggtitle("Five most abundant genera and other taxa")
ggsave(file="CRC_Stages_Analysis/Abundances/TOP5_Genera.png",width=6.5,height=4.5, dpi=300)
dev.off()

############### ALFA AND BETA DIVERSITY CRC vs AP (STOOL) ###################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1", "Shannon", "Observed"), x="Tumor", color = "Patient") # prendo valori dal grafico
pAlpha
{Chao<-dplyr::filter(pAlpha$data, variable=="Chao1")
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evennes
identical(H$Codice,obs$Codice) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev)
pAlpha$data<-New_data
}
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Tumor, y=value, color=NULL), alpha=0.1) + theme_bw() + labs(x="") +
  guides(fill="none", color="none") +
  theme(axis.text.x = element_text(angle = -20, vjust=1, hjust = 0) ) +
  stat_compare_means(aes(group = Tumor), label="p.format", method = "anova", label.x= 1.4, size=3.2, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alpha_diversity_adenoma_vs_adenoc.png", width = 6,height = 5, dpi=300)

# PLUS:
shapiro.test(pAlpha$data[pAlpha$data$variable=="Observed", "value"]) # not Gaussian
shapiro.test(pAlpha$data[pAlpha$data$variable=="Shannon", "value"]) # not Gaussian

########################## BETA DIVERSITY

{OTU.genus.prop<-as.data.frame(otu_table(data.genus.prop))
OTU.fam.prop<-as.data.frame(otu_table(data.fam.prop))
OTU.class.prop<-as.data.frame(otu_table(data.class.prop))
OTU.order.prop<-as.data.frame(otu_table(data.order.prop))
OTU.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA

sample_OTU<-as.data.frame(t(sqrt(OTU.genus.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permg<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.fam.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permf<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.class.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permc<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.order.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permo<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.phy.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permp<-perm

beta<-rbind(permg$aov.tab[1,],permf$aov.tab[1,],permo$aov.tab[1,],permc$aov.tab[1,],permp$aov.tab[1,])
row.names(beta)<-c("Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta,file="Beta_div_permanova_Bray_ADENOMA_VS_ADENOCARC.csv",quote=F,row.names = T)

##### PCoA

dir.create("PCoA")

# su ASV
{data.sqrt<-transform_sample_counts(data.prop, sqrt)
DistBC = phyloseq::distance(data.sqrt, method = "bray")
ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues
eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Tumor") +
  geom_point(size=2.5) +
  geom_text(aes(label = sample_names(data)), show.legend=F,
            position = (position_nudge(y=0.015)), size=2, color="black") + 
  theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_sqrt_ASV_prop_STOOL.png", width = 6, height = 5, dpi=300)
# with ellipses
plot_ordination(data.sqrt, ordBC, color = "Tumor") +
  geom_point(size=2.5) + stat_ellipse() + theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_sqrt_ASV_prop_STOOL_ELLIPSE.png", width = 6, height = 5, dpi=300)

# Adenocarcinoma only. in order to focus on genders
data_subset<-subset_samples(data.prop, Tumor=="Adenocarcinoma")
{data.sqrt<-transform_sample_counts(data_subset, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +  theme_classic() +
  scale_color_manual(values = c("F"="pink","M"="lightblue"))+
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only patients with adenocarcinoma)", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_ADENOCARCINOMA_GENDER_STOOL.png", width = 6, height = 5, dpi=300)

# Adenoma only, in order to focus on genders
data_subset<-subset_samples(data.prop, Tumor=="Adenoma")
{data.sqrt<-transform_sample_counts(data_subset, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +  theme_classic() +
  scale_color_manual(values = c("F"="pink","M"="lightblue"))+
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only patients with adenoma)", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_ADENOMA_GENDER_STOOL.png", width = 6, height = 5, dpi=300)


###################### DA WITH DESEQ2 CRC vs AP (STOOL) ####################################

dir.create("DESEQ2")

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

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
  
  DEseq_data<-phyloseq_to_deseq2(d, ~Tumor)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Tumor", "Adenocarcinoma", "Adenoma"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("DESEQ2/DA_",t,"_ratio_Adenocarcinoma_vs_Adenoma.csv"), row.names = F, quote=F, na = "")
    r_level<-r
    r_level[, "Taxon"]<- rep(t)
    Res_tot<-rbind.data.frame(Res_tot,r_level)
    # single box plots
    target<-r[[t]]
    colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
    target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
    Table_DE<-psmelt(target)
    colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
    Table_DE$ASV<-NULL
    #Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
    assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
    # appending to unique box plot
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
write.csv(Res_tot, file="DESEQ2/Every_result_DESeq2.csv", row.names = F)

Table_tot$Bacteria<-gsub("Peptostreptococcales-Tissierellales","Peptostreptococcales\nTissierellales", Table_tot$Bacteria)
ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Tumor)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=12.9,colour="black")) + 
  scale_fill_manual(values=c("Adenocarcinoma"="coral","Adenoma"="red4")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  labs(title= "Differently abundant Taxa in Stool samples between patients with Adenocarcinoma and Adenoma", y="Proportional Abundance", 
       fill="Tumor", x="")
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(1,9,1),seq(10,30,2))) # remove the "zoom" if not needed
ggsave(filename = "DESEQ2/DA_Tumor_every_result_asse_normale.png", width = 15, height = 8, dpi=300)
dev.off()

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Tumor)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=12.9,colour="black")) + 
  scale_fill_manual(values=c("Adenocarcinoma"="coral","Adenoma"="red4")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  labs(title= "Differently abundant Taxa in Stool samples between patients with Adenocarcinoma and Adenoma", y="Proportional Abundance", 
       fill="Tumor", x="") +
scale_y_sqrt(breaks=c(0.1,0.5,seq(1,4,1),seq(6,41,2)))
ggsave(filename = "DESEQ2/DA_Tumor_every_result_asse_con_zoom.png", width = 15, height = 8, dpi=300)
dev.off()


Table_tot2<-subset(Table_tot, Taxa == "Genus") # to remove redundant results
Table_tot2$Bacteria<-gsub("_group","",Table_tot2$Bacteria)
Table_tot2$Bacteria<-gsub("[","",Table_tot2$Bacteria, fixed=T)
Table_tot2$Bacteria<-gsub("]","",Table_tot2$Bacteria, fixed = T)
deseq_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Tumor)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5)) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1,0.5,seq(1,4,1),seq(6,41,2))) +
  labs(title= "", y="Proportional Abundance", fill="stool", x="")

deseq_plot +  labs(title= "Differently abundant Taxa in Stool samples between patients with Adenocarcinoma and Adenoma")
ggsave(filename = "DESEQ2/DA_Tumor_no_redundants_zoom.png", width = 13.5, height = 9, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > DESEQ2/NB.txt ")

# for the correlations
Genera.DESEQ2_CRC_AP_stool<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


####### CIRCOTAXA PLOT TO PLOT DA LOG2FOLDCHANGE (STOOL)

Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
Res_tot2<-Res_tot2[Res_tot2$Taxon =="Genus", ] # to remove redundants
Res_tot2$Genus<-gsub("[Eubacterium]_coprostanoligenes_group","Eubacterium_coprostanoligenes",Res_tot2$Genus, fixed=T)
Res_tot2$Genus<-gsub("[Ruminococcus]_torques_group","Ruminococcus_torques",Res_tot2$Genus, fixed=T)

circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")

png(file="DESEQ2/CircoTax_plot.png",width=1800,height=1800, res = 300)
circo + labs(title="Adenocarcinoma vs adenoma \n DA log2foldchange (stool)\n")
dev.off()

circo_no_m<- circo + theme(plot.margin = unit(c(0,0,0,0.8), units = "cm"))

# diff abundances + logfoldchange
png(file="DESEQ2/Final_plot.png",width=4000, height=2150, res = 300)
ggarrange(deseq_plot, circo_no_m, nrow = 1, widths= c(1,0.58), labels = c("A)", "B)"))
dev.off()


################### ALFA AND BETA DIVERSITY ADENONOCARCINOMA STAGES (Stool) ############################

dir.create("CRC_Stages_Analysis/Alpha_diversity")

{pAlpha<-plot_richness(data.stages, measures=c("Shannon", "Observed"), x="Stadiation")
  pAlpha
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
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Stadiation, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Stadiation", title="Alpha diversity between stages I, II and III") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11)) +
  stat_compare_means(aes(group = Stadiation), label="p.format", method = "kruskal.test", label.x= 1.05, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="CRC_Stages_Analysis/Alpha_diversity/Alfa_diversity_with_Kruskal.png", width = 7,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data)
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Stadiation
kruskal.test(Obser_value$value~factor)

# pairwise comparisons
con <- file("CRC_Stages_Analysis/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_Stadiation.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Stadiation sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-filter(alphadt, variable=="Observed")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Shannon")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Evenness")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)
rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor)

###### BETA DIVERSITY

dir.create("CRC_Stages_Analysis")
dir.create("CRC_Stages_Analysis/Beta_diversity")

{ASV.prop<-as.data.frame(otu_table(data.stages.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.stages.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.stages.phy.prop))
}

#### PERMANOVA

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_Bray<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_Bray<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=as(sample_data(data.stages.prop),"data.frame")$Stadiation, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV

# exporting beta diversity
rm(con)
con<-file("CRC_Stages_Analysis/Beta_diversity/Beta_diversity_general_and_pairwise_between_Stadiation_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_p, perm_ASV)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
  disper<-vegan::betadisper(BC.dist,as(sample_data(data.stages.prop),"data.frame")$Stadiation)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
rm(con)
con<-file("CRC_Stages_Analysis/Beta_diversity/Beta_dispersion_General_and_Pairwise_between_Stadiation_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con)

####### PCoA BRAY CURTIS

data.stages.prop.labels<-data.stages.prop
{data.stages.sqrt_prop<-transform_sample_counts(data.stages.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.stages.sqrt_prop, method = "bray")
  ordBC = ordinate(data.stages.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.stages.sqrt_prop, ordBC, color = "Stadiation") +
  scale_color_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Stadiation", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="CRC_Stages_Analysis/Beta_diversity/PCoA_Beta_diversity_Bray_Curtis_on_ASV.png", width = 8, height = 6, dpi=300)

######################## DA WITH DESEQ2 ON CRC STAGES (Stool) #######################

dir.create("CRC_Stages_Analysis/DA_with_DESeq2")
# Trimming under sum of 10, see DeSeq2 tutorial
suppressWarnings(rm(data_pruned))
data_pruned<- prune_taxa(taxa_sums(data.stages) > 10, data.stages)

# preparing new taxa vocabularies for plots (some ASV may change glomming after trimming)
{data.stages.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.stages.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.stages.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.stages.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.stages.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.prop <- transform_sample_counts(data.stages, function(ASV) ASV/sum(ASV)*100)
  data_pruned.phy.prop <- transform_sample_counts(data.stages.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.stages.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.stages.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.stages.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.stages.genus_pruned, function(ASV) ASV/sum(ASV)*100)
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
      d<-phyloseq_to_deseq2(data_pruned, ~Stadiation)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Stadiation)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Stadiation", nume[i], deno[i]), test="Wald")
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("II","I","II")
deno<-c("I","III","III")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
unlink("DE_results.txt")
DE_results<-DE_results[DE_results$BaseMean>50, ]
system(" echo 'Every result under the arbitrary threshold Basemean=50 has been removed in order to avoid the noisiest results' > CRC_Stages_Analysis/DA_with_DESeq2/NB_results_are_filtered.txt")
write.csv2(DE_results, file="CRC_Stages_Analysis/DA_with_DESeq2/DE_every_results.csv", quote=F, row.names = F)
write.xlsx(DE_results, file="CRC_Stages_Analysis/DA_with_DESeq2/DE_every_results.xlsx",showNA = F, col.names = T, row.names = F)

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
  target<-prune_taxa(target, data_pruned.genus.prop)
  tabella<-psmelt(target)
  tabella_g<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Family","ASV"]
  target<-prune_taxa(target, data_pruned.fam.prop)
  tabella<-psmelt(target)
  tabella_f<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Order","ASV"]
  target<-prune_taxa(target, data_pruned.order.prop)
  tabella<-psmelt(target)
  tabella_o<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Class","ASV"]
  target<-prune_taxa(target, data_pruned.class.prop)
  tabella<-psmelt(target)
  tabella_c<-tabella
  rm(tabella)
}

# {target<-DE_results[DE_results$Rank=="Phylum","ASV"]
#   target<-prune_taxa(target, data_pruned.phy.prop)
#   tabella<-psmelt(target)
#   tabella_p<-tabella
#   rm(tabella)
# }

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
# {tabella_p$Taxa<-"Phyla"
#   colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"
# }

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_c,tabella_o)#,tabella_p)

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw( ) + 
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust=1, size=10), 
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))#,seq(22,max(tabella_tot$Abundance),3) ))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/EVERY_RESULT_DeSeq2_Stadiation.png", width = 12, height = 7, dpi=300)
dev.off()

Redund<-c("uncultured_ o Rhodospirillales","Porphyromonadaceae")

tabella_tot2<-subset(tabella_tot, ! Bacteria %in% Redund)
ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=10),
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") + 
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))#,seq(22,max(tabella_tot$Abundance),3)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/EVERY_RESULT_no_redundants_Stadiation.png", width = 10, height = 7, dpi=300)
dev.off()


################# now plotting result in pairwise (group vs group)

##### II vs I
tabella_tot3<-subset(tabella_tot2, Stadiation != "III")
I_vs_II<-DE_results[DE_results$num=="II" & DE_results$denom=="I", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% I_vs_II [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, I_vs_II[I_vs_II$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages I and II", y="Proportional Abundance", fill="Stadiation", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_I_vs_II_no_redundants.png", width = 10, height = 7, dpi=300)

# Results.DESEQ2_I_vs_II<-tabella_tot3[,c("Bacteria","Taxa")]
# Results.DESEQ2_I_vs_II<-Results.DESEQ2_I_vs_II[!duplicated(Results.DESEQ2_I_vs_II$Bacteria),]
# head(Results.DESEQ2_I_vs_II) # just a check

#####  I vs III
tabella_tot3<-subset(tabella_tot2, Stadiation != "II")
III_vs_I<-DE_results[DE_results$num=="I" & DE_results$denom=="III", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% III_vs_I [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, III_vs_I[III_vs_I$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12),
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages I and III", y="Proportional Abundance", fill="Stadiation", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_I_vs_III_DeSeq2_no redundants.png", width = 10, height = 7, dpi=300)

# Results.DESEQ2_I_vs_III<-tabella_tot3[,c("Bacteria","Taxa")]
# Results.DESEQ2_I_vs_III<-Results.DESEQ2_I_vs_III[!duplicated(Results.DESEQ2_I_vs_III$Bacteria),] # one row for each taxa
# head(Results.DESEQ2_I_vs_III) # just a check

##### II vs III
tabella_tot3<-subset(tabella_tot2, Stadiation != "I")
III_vs_II<-DE_results[DE_results$num=="II" & DE_results$denom=="III", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% III_vs_II [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, III_vs_II[III_vs_II$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages II and III", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_II_vs_III_DeSeq2_no redundants.png", width = 10, height = 7, dpi=300)

# Results.DESEQ2_II_vs_III<-tabella_tot3[,c("Bacteria","Taxa")]
# Results.DESEQ2_II_vs_III<-Results.DESEQ2_II_vs_III[!duplicated(Results.DESEQ2_II_vs_III$Bacteria),] # one row for each taxa
# head(Results.DESEQ2_II_vs_III) # just a check


############## ANALYSIS ON FECAL METABOLITES STOOL (CRC vs AP) ###############

dir.create("Analysis_on_metabolites")

Metabolites<-read.csv(file="../../Metabolites.csv", row.names = 1)
Fatty<-read.csv(file="../../Fatty_acid.csv", row.names = 1)
aa<-read.csv(file="../../Aminoacids.csv", row.names = 1)
# the row.names are the "Patient_ID" (identificative number)

head(Metabolites,n=2)
head(Fatty,n=2)
head(aa,n=2)

objs<-list(Metabolites,Fatty,aa)               # NB: it is a list
names(objs)<-c("Metabolites","Fatty","aa")

###### preparing the data sets                 # NB: this script works only with stool sample data! Do not run it with biopsy data (different number of sample)
for(x in 1:length(objs)){
  suppressWarnings(rm(no_microbiota,original_lenght,sample_order))
  original_lenght<-as.numeric(length(!is.na(objs[[x]][1])))
  sample_order<-sample_data(data)$Patient_number[sample_data(data)$Patient_number %in% row.names(objs[[x]])]
  objs[[x]]<-objs[[x]][as.character(sample_order),]
  objs[[x]]<-sapply(objs[[x]], MARGIN = 2, as.numeric)
  row.names(objs[[x]])<-sample_order # they were lost after sapply
  if(identical(as.numeric(length(!is.na(objs[[x]][,1]))),original_lenght)){
    cat("Any sample has been cutted out from the",names(objs)[x],"dataset --> OK\n")
  } else { cat("\n Something has gone wrong ... \n "); Sys.sleep(1) }
  objs[[x]]<-as.data.frame(t(objs[[x]])) # samples are now columns (just as the phyloseq obj)
  #### THIS NORMALIZATION STEP IS OFF! This NRM data is already normalized through PQN
  # objs[[x]]<-apply(objs[[x]], MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
  # if(as.integer(sum(objs[[x]][,1])) != 100){
  #   cat("\nLooks like there has been a problem with the trasformation to relative abundances in",names(objs)[x],"dataset, check it!\n", fill=T)
  #   Sys.sleep(1)} else{ cat("Successfully translated the dataset and trasformed values to relative abundances --> OK\n")}
}

#### Wilcoxon Mann Whithney (To Select the differentially abundant variables)
meta<-sample_data(data)
row.names(meta)<-meta$Patient_number
meta<-meta[as.character(sample_order),"Tumor"] # same order as matrix samples

#dir.create("Analysis_on_metabolites/Processed_datasets_used")

Wilcoxons<-NULL
for(x in 1:length(objs)){
  if(! identical(row.names(meta),colnames(objs[[x]]))){ # just a check
    cat("\nSomething has gone wrong, names are not matching...\n", fill=T)
    Sys.sleep(1)}
  part_Wilcox<-NULL # subset of the "Wilcoxons" object in creation, in order to properly denominates the rows for each dataset
  rows<-rownames(objs[[x]]) # while they are all present
  for(y in 1:length(row.names(objs[[x]]))){
    temp<-wilcox.test(as.numeric(objs[[x]][y,]) ~ meta$Tumor) # samples are on columns --> testing rows (variables)
    part_Wilcox<-rbind.data.frame(part_Wilcox, cbind(unname(temp$statistic),temp$p.value))
  #   if(temp$p.value>0.05){
  #     objs[[x]][y,] <- 0 # the row is "labelled" because it's not significant (it can not be removed during the "y" loop)
  #   }
  }
  # objs[[x]]<-objs[[x]][! rowSums(objs[[x]])==0, ] # removed the not significant variables
  # write.csv(objs[[x]], file=paste("Analysis_on_metabolites/Processed_datasets_used/",names(objs)[x],
  #                                "_significative_according_MANN_WHITNEY.csv") )
  row.names(part_Wilcox)<-rows
  #part_Wilcox$p.adj<-p.adjust(part_Wilcox$V2, "BH")  ### NB: reactivate this row to perform p-value correction done on each subset, not on the whole
  Wilcoxons<-rbind.data.frame(Wilcoxons,part_Wilcox)
}

Wilcoxons$p.adj<-p.adjust(Wilcoxons$V2, "BH") # NB: this correction is done on the whole dataset
colnames(Wilcoxons)<-c("W_statistic","p-value", "p-adjusted_BH")
write.csv(Wilcoxons, file=paste("Analysis_on_metabolites/MANN_WHITNEY_CRC_vs_AP_RESULTS.csv",sep="_"))

suppressWarnings(rm(Significative))
Significative<-row.names(Wilcoxons)[Wilcoxons$`p-adjusted_BH`<0.05]
Significative

aa2 <- cbind.data.frame(aa,Metabolites)
aa2 <- aa2[,Significative]
#aa2<- as.data.frame(t(aa2))


#### now testing correlation between genera and metabolites in CRC vs AP

# vs Five Most abundant Genera
top5 <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
taxa_top5 <- prune_taxa(top5,data.genus.prop)
Names<-as.data.frame(tax_table(taxa_top5))
Names<-Names$Genus
corr_top5 <- t(as.data.frame(otu_table(taxa_top5)))
colnames(corr_top5)<-Names
row.names(corr_top5)<-gsub("CRC","",row.names(corr_top5)) # the "patient_number" variable is inside the sample name itself
row.names(corr_top5)<-gsub("AP","",row.names(corr_top5))

sample_order<-as.character(sample_data(data.genus)$Patient_number[sample_data(data.genus)$Patient_number %in% row.names(aa2)])
corr_top5<-corr_top5[sample_order,]
l<-length(row.names(corr_top5))
l # 28
aa2<-aa2[sample_order, ]
identical(l, length(which(! is.na(corr_top5[,1]))) ) # TRUE --> any sample has been cutted out
identical(row.names(corr_top5),row.names(aa2))

{x<-as.matrix(cbind(corr_top5,aa2))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  identical(data_corr[,1:2],data_pvalue[,1:2])
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Bacteria","Aminoacids","Corr","pvalue")
  data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_top5), ]
  data_corr<-data_corr[data_corr$Aminoacids %in% colnames(aa2), ]
  data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
data_corr$Bacteria<-gsub("[","",data_corr$Bacteria, fixed = T)
data_corr$Bacteria<-gsub("]","",data_corr$Bacteria, fixed = T)
write.csv2(data_corr, file="Analysis_on_metabolites/Spearman_TOP5_abundant_genera_stool_vs_aa.csv", quote=F, row.names = F)

data_corr$simbol<-data_corr$padj
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$Aminoacids, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=8) +
  theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 7)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= simbol), color= "white", size =9, nudge_y = -0.1) +
  labs(title = "Spearman correlation between \n five most abundant bacterial genera in stool samples and \n metabolites in stool differently abundant between CRC and AP samples ",
       y= "", x= "", caption= "* = p-value < 0.05") +
  theme(plot.title = element_text(size=10)) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=8), legend.key.height= unit(1.1, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=10))
ggsave(filename = "Analysis_on_metabolites/Heatmap_TOP5Genera_vs_metab.png", width = 6, height = 4.5, dpi=300)

#### with DESeq2 results
generi_sel <- subset_taxa(data.genus.prop, Genus %in% Genera.DESEQ2_CRC_AP_stool)
identical(row.names(otu_table(generi_sel)),row.names(tax_table(generi_sel)))
{Names<-as.data.frame(tax_table(generi_sel))
  Names<-Names$Genus
  corr_tot<- t(as.data.frame(otu_table(generi_sel)))
  colnames(corr_tot)<-Names
  row.names(corr_tot)<-gsub("CRC","",row.names(corr_tot))
  row.names(corr_tot)<-gsub("AP","",row.names(corr_tot))
}

sample_order<-as.character(sample_data(data.genus)$Patient_number[sample_data(data.genus)$Patient_number %in% rownames(aa2)])
corr_tot<-corr_tot[sample_order,]
l<-length(row.names(corr_tot))
l # 28
aa2<-aa2[sample_order,]
identical(l, length(which(! is.na(row.names(corr_tot)))) ) # TRUE --> any sample has been cutted out
identical(row.names(corr_tot),row.names(aa2))

{x<-as.matrix(cbind(corr_tot,aa2))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  identical(data_corr[,1:2],data_pvalue[,1:2])
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Bacteria","Aminoacids","Corr","pvalue")
  data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_tot), ]
  data_corr<-data_corr[data_corr$Aminoacids %in% colnames(aa2), ]
  data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
data_corr$Bacteria<-gsub("[","",data_corr$Bacteria, fixed = T)
data_corr$Bacteria<-gsub("]","",data_corr$Bacteria, fixed = T)
write.csv2(data_corr, file="Analysis_on_metabolites/Spearman_of_DESeq2_Results_vs_aminoacids.csv", quote=F, row.names = F)

data_corr$simbol<-data_corr$padj
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$Aminoacids, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=8) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 7)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.1) +
  labs(title = "Spearman correlation between \n differently abundant bacterial genera and metabolites in stool sample\n between CRC and AP patients",
       y= "", x= "", caption= "* = p-adjusted < 0.05") +
  theme(plot.title = element_text(size=10)) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=8), legend.key.height= unit(1.1, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=10))
ggsave(filename = "Analysis_on_metabolites/Heatmap_DESeq2_results_vs_metab.png", width = 6, height = 4.5, dpi=300)


############## ANALYSIS ON FECAL METABOLITES STOOL (CRC Stages) ###############

dir.create("CRC_Stages_Analysis/Analysis_on_metabolites")

Metabolites<-read.csv(file="../../Metabolites.csv", row.names = 1)
Fatty<-read.csv(file="../../Fatty_acid.csv", row.names = 1)
aa<-read.csv(file="../../Aminoacids.csv", row.names = 1)
# the row.names are the "Patient_ID" (identificative number)

head(Metabolites,n=2)
head(Fatty,n=2)
head(aa,n=2)

objs<-list(Metabolites,Fatty,aa)
names(objs)<-c("Metabolites","Fatty","aa")

###### preparing the data sets
for(x in 1:length(objs)){
  suppressWarnings(rm(no_microbiota,original_lenght,sample_order))
  sample_order<-sample_data(data.stages)$Patient_number[sample_data(data.stages)$Patient_number %in% row.names(objs[[x]])]
  objs[[x]]<-objs[[x]][as.character(sample_order),]
  objs[[x]]<-sapply(objs[[x]], MARGIN = 2, as.numeric)
  row.names(objs[[x]])<-sample_order # they were lost after sapply
  if(identical(as.numeric(length(!is.na(objs[[x]][,1]))),20)){
    cat("There are 20 samples in the",names(objs)[x],"dataset --> OK\n")
  } else { cat("\n Something has gone wrong ... \n "); Sys.sleep(1) }
  objs[[x]]<-as.data.frame(t(objs[[x]])) # samples are now columns (just as the phyloseq obj)
}

#### Kruskal Wallis (To Select the differentially abundant variables)
meta<-sample_data(data.stages)
row.names(meta)<-meta$Patient_number
meta<-meta[as.character(sample_order),"Stadiation"] # same order as matrix samples

Kruskal<-NULL
for(x in 1:length(objs)){
  if(! identical(row.names(meta),colnames(objs[[x]]))){ # just a check
    cat("\nSomething has gone wrong, names are not matching...\n", fill=T)
    Sys.sleep(1)}
  part_Krusk<-NULL # subset of the "Kruskal" object in creation, in order to properly denominates the rows for each dataset
  rows<-rownames(objs[[x]]) # while they are all present
  for(y in 1:length(row.names(objs[[x]]))){
    temp<-kruskal.test(as.numeric(objs[[x]][y,]) ~ meta$Stadiation) # samples are on columns --> testing rows (variables)
    part_Krusk<-rbind.data.frame(part_Krusk, cbind(unname(temp$statistic),temp$p.value))
    # if(temp$p.value>0.05){
    #   objs[[x]][y,] <- 0 # the row is "labelled" because it's not significant (it can not be removed during the "y" loop)
    # }
  }
  #objs[[x]]<-objs[[x]][! rowSums(objs[[x]])==0, ] # removed the not significant variables
  #write.csv(objs[[x]], file=paste("Analysis_on_metabolites/Processed_datasets_used/",names(objs)[x],
  #                                "_significative_according_MANN_WHITNEY.csv") )
  row.names(part_Krusk)<-rows
  #part_Krusk$p.adj<-p.adjust(part_Krusk$V2, "BH") # NB: p-value correction done on each subset, not on the whole
  Kruskal<-rbind.data.frame(Kruskal,part_Krusk)
}

colnames(Kruskal)<-c("chi-squared","p-value","p-adjusted_BH")
write.csv(Kruskal, file=paste("CRC_Stages_Analysis/Analysis_on_metabolites/KRUSKAL_metabolites_between_stages_RESULTS.csv",sep="_"))

suppressWarnings(rm(Significative))
Significative<-row.names(Kruskal)[Kruskal$`p-value`<0.05]
Significative
# any significative results according to kruskal...

# just checking the post hoc wilcoxon results anyways
pair_W<-NULL
for(x in 1:length(objs)){
  if(! identical(row.names(meta),colnames(objs[[x]]))){ # just a check
    cat("\nSomething has gone wrong, names are not matching...\n", fill=T)
    Sys.sleep(1)}
  part_pair<-NULL # subset of the "pair_W" object in creation, in order to properly denominates the rows for each dataset
  rows<-rownames(objs[[x]]) # while they are all present
  for(y in 1:length(row.names(objs[[x]]))){
    temp<-pairwise.wilcox.test(as.numeric(objs[[x]][y,]),meta$Stadiation) # samples are on columns --> testing rows (variables)
    p_temp<-data.frame("I_vs_II"=temp$p.value[1,1],"I_vs_III"=temp$p.value[2,1],"II_vs_III"=temp$p.value[2,2])
    part_pair<-rbind.data.frame(part_pair, p_temp)
  }
  row.names(part_pair)<-rows
  pair_W<-rbind.data.frame(pair_W,part_pair)
}

write.csv(pair_W, file="CRC_Stages_Analysis/Analysis_on_metabolites/Pairwise_Wilcox_metabolites_between_stages_RESULTS.csv")
Signif_I_vs_II<-pair_W[pair_W$I_vs_II<0.05,"I_vs_II"]
Signif_I_vs_II
Signif_I_vs_III<-pair_W[pair_W$I_vs_III<0.05,"I_vs_III"]
Signif_I_vs_III
Signif_II_vs_III<-pair_W[pair_W$II_vs_III<0.05,"II_vs_III"]
Signif_II_vs_III

system("touch CRC_Stages_Analysis/Analysis_on_metabolites/No_signif_results_for_Kruskal_or_pairwise_Wilcox.txt")

################### ANALYSING LEFSE RESULTS (STOOL) ###########################


suppressWarnings(rm(a))
a <- read.delim("../../QIIME2/PICRUST2_LEFSE_stool/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("../../QIIME2/PICRUST2_LEFSE_stool/Result_LEFSE.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="PICRUST2_LEFSE/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)

# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Adenocarcinoma"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Adenocarcinoma"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Adenoma",1,0), size=3.2) +
  scale_fill_manual(values=c("Adenocarcinoma"="coral", "Adenoma"="deepskyblue")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST2_LEFSE/PICRUST2_LEFSE_plot_diff_METACYC_every_result.png", width = 9, height = 3.5, dpi = 300)

system('touch PICRUST2_LEFSE/no_result_over_3')

rm(Significative_functions_LEFSE, a)


###################### STARTING BIOPSY ANALYSIS ##############################

if( grepl("_analysis",getwd()) | grepl("_correlations",getwd()) ){setwd("..")}
dir.create("Biopsy_analysis")
setwd("Biopsy_analysis")

to_remove<-ls()
to_remove<-to_remove[ ! to_remove %in% c("data_total","metadata","evalslopes", "CircoTax2", "Genera.DESEQ2_CRC_AP_stool","Genera.DESEQ2_CRC_AP_biopsy")]
to_remove
rm(list = to_remove) # reset the workspace

data<-subset_samples(data_total,Sampling_Type=="B")
sample_names(data)<-sample_data(data)$Patient
head(sample_data(data))

### looks like there is a little of eucaryotic DNA contamination here... It's better to check it
if( ! "AP19" %in% sample_names(data) | ! "CRC7" %in% sample_names(data) ){
  cat("Looks like the object is already decontaminated... check it again")
  Sys.sleep(2)
} else{
    data_bio_contam<- data }

############ EXPORTING BIOPSY ABUNDANCES (with contaminants) ################## 

dir.create("Abundances")
dir.create("Abundances/Relat_abundances_with_contaminants")

{ data.genus2<- tax_glom(data_bio_contam, taxrank = "Genus", NArm = F)
data.phy2<- tax_glom(data_bio_contam, taxrank = "Phylum", NArm = F)
data.prop <- transform_sample_counts(data_bio_contam, function(ASV) ASV/sum(ASV)*100)
data.phy.prop <- transform_sample_counts(data.phy2, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus2, function(ASV) ASV/sum(ASV)*100)
}

{options( scipen=100 ) 
  write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Relat_abundances_with_contaminants/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Relat_abundances_with_contaminants/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Relat_abundances_with_contaminants/counts_genus.csv",quote=F)
}

# TOP 5 Phyla 
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="Abundances/Relat_abundances_with_contaminants/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

rm(data.genus2, data.genus.prop, data.phy2, data.phy.prop, data.prop, tabella, top5, others, prune.dat_top5, prune.data.others, tabella_others, tabella_top)

########### CHECKING AND REMOVING PROBABLE EUCARYOTIC CONTAMINANTS (BIOPSY) ##############

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
write.csv2(e[,colnames(e)!="Kingdom"], file="Unassigned_domain_checking.csv", row.names = T, quote = F)

####### exporting ASV names (to check)
human_ASV_biopsy<-as.data.frame(taxa_names(subset_taxa(data_bio_contam,Kingdom=="Unassigned")))
colnames(human_ASV_biopsy)<-"ASV_code"
write.csv(human_ASV_biopsy, file = "human_ASV_biopsy.csv", row.names = F, quote=F)

####### if the domain is Unassigned in SILVA then it's probably a contaminant
data<- subset_taxa(data, Kingdom!="Unassigned")
head(tax_table(data))
write.csv2(tax_table(data), file="Every filtered ASV and taxonomic assignment.csv", row.names = T)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x, human_ASV_biopsy)

##### moreover, the samples CRC54, CRC7 and AP19 have about the 50% of unassigned DNA... better to remove them
data<- subset_samples(data, ! Patient %in% c("CRC7", "AP19", "CRC54"))
if(length(sample_names(data))+3 == length(sample_names(data_bio_contam))){
  cat("Three samples have been removed correctly")
}

### better to check the saturation again
r<-rarecurve(t(otu_table(data)), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
# the samples are still saturated


############# PREPARING THE BIOPSY OBJECTS (AFTER THE DECONTAMINATION) ##################

if( length(sample_names(data)) == length(sample_names(data_bio_contam))){
  cat("\n \n Proceed with the decontamination before!")} else {
    data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
    OTU.genus<-as.data.frame(otu_table(data.genus))
    Taxa.genus<-as.data.frame(tax_table(data.genus))
    data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
    OTU.fam<-as.data.frame(otu_table(data.fam))
    Taxa.fam<-as.data.frame(tax_table(data.fam))
    data.class<- tax_glom(data, taxrank = "Class", NArm = F)
    OTU.class<-as.data.frame(otu_table(data.class))
    Taxa.class<-as.data.frame(tax_table(data.class))
    data.order<- tax_glom(data, taxrank = "Order", NArm = F)
    OTU.order<-as.data.frame(otu_table(data.order))
    Taxa.order<-as.data.frame(tax_table(data.order))
    data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
    OTU.phy<-as.data.frame(otu_table(data.phy))
    Taxa.phy<-as.data.frame(tax_table(data.phy))
    }
  
{data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)
}

dir.create("CRC_Stages_Analysis")
data.stages<-subset_samples(data, !is.na(Stadiation)) # only samples with adenocarcinoma, in order to confront stadiations
sample_data(data.stages)$Stadiation<-factor(sample_data(data.stages)$Stadiatio, levels=c("I","II","III"))
{data.stages.prop<-transform_sample_counts(data.stages, function(otu) otu/sum(otu)*100)
  data.stages.genus<- tax_glom(data.stages, taxrank = "Genus", NArm = F)
  data.stages.genus.prop <- transform_sample_counts(data.stages.genus, function(otu) otu/sum(otu)*100)
  data.stages.phy<- tax_glom(data.stages, taxrank = "Phylum", NArm = F)
  data.stages.phy.prop <- transform_sample_counts(data.stages.phy, function(otu) otu/sum(otu)*100)
}

#################### EXPORTING ABUNDANCES (BIOPSY) ############################

dir.create("Abundances")
dir.create("Abundances/Percent_abundances")

{options( scipen=100 ) 
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Percent_abundances/counts_phylum.xlsx",row.names = F, showNA = F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Percent_abundances/counts_genus.xlsx",row.names = F, showNA = F)
  write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Percent_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class.prop),"matrix")),file="Abundances/Percent_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order.prop),"matrix")),file="Abundances/Percent_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam.prop),"matrix")),file="Abundances/Percent_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Percent_abundances/counts_genus.csv",quote=F)
}

### TOP 5 Orders
{top5 <- names(sort(taxa_sums(data.order.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.order.prop)
  others<-taxa_names(data.order.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.order.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Order<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Order)) + 
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every order below rank 5 ") +
  ggtitle("Five most abundant orders and other taxa")
ggsave(file="Abundances/TOP5_Orders.png",width=6.5,height=4.5, dpi=300)
dev.off()

# TOP 5 Phyla 
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="Abundances/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Genus)) + 
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every genus below rank 5 ") +
  ggtitle("Five most abundant genera and other taxa")
ggsave(file="Abundances/TOP5_Genera.png",width=6.5,height=4.5, dpi=300)
dev.off()

# TOP 5 Families
{top <- names(sort(taxa_sums(data.fam.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.fam.prop)
  others<-taxa_names(data.fam.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.fam.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Family<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Family)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every family below rank 5 ") +
  ggtitle("Five most abundant families and other taxa")
ggsave(file="Abundances/TOP5_families.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Classes
{top <- names(sort(taxa_sums(data.class.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.class.prop)
  others<-taxa_names(data.class.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.class.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Class<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Class)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every class below rank 5 ") +
  ggtitle("Five most abundant classes and other taxa")
ggsave(file="Abundances/TOP5_classes.png",width=6.5,height=4.5,dpi=300)
dev.off()

dir.create("CRC_Stages_Analysis/Abundances")

write.xlsx(cbind(as(otu_table(data.stages.phy.prop),"matrix"),as(tax_table(data.stages.phy.prop),"matrix")),file="CRC_Stages_Analysis/Abundances/Percent_counts_phylum.xlsx",row.names = F, showNA = F)
write.xlsx(cbind(as(otu_table(data.stages.genus.prop),"matrix"),as(tax_table(data.stages.genus.prop),"matrix")),file="CRC_Stages_Analysis/Abundances/Percent_counts_genus.xlsx",row.names = F, showNA = F)

# TOP 5 Phyla (Only CRC Stages)
{top5 <- names(sort(taxa_sums(data.stages.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.stages.phy.prop)
  others<-taxa_names(data.stages.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.stages.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Stadiation), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="CRC_Stages_Analysis/Abundances/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Genera (Only CRC Stages)
{top <- names(sort(taxa_sums(data.stages.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.stages.genus.prop)
  others<-taxa_names(data.stages.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.stages.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Genus)) + 
  facet_grid2(cols= vars(Stadiation), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every genus below rank 5 ") +
  ggtitle("Five most abundant genera and other taxa")
ggsave(file="CRC_Stages_Analysis/Abundances/TOP5_Genera.png",width=6.5,height=4.5, dpi=300)
dev.off()

################## ALFA AND BETA DIVERSITY CRC vs AP (BIOPSY) ###################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1", "Shannon", "Observed"), x="Tumor", color = "Patient") # prendo valori dal grafico
pAlpha
{Chao<-dplyr::filter(pAlpha$data, variable=="Chao1")
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evennes
  identical(H$Codice,obs$Codice) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
}
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Tumor, y=value, color=NULL), alpha=0.1) + theme_bw() + labs(x="") +
  guides(fill="none", color="none") +
  theme(axis.text.x = element_text(angle = -20, vjust=1, hjust = 0) ) +
  stat_compare_means(aes(group = Tumor), label="p.format", method = "wilcox.test", label.x= 1.4, size=3.2, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alpha_diversity_adenoma_vs_adenoc.png", width = 6,height = 5, dpi=300)

########################## BETA DIVERSITY

{OTU.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  OTU.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  OTU.class.prop<-as.data.frame(otu_table(data.class.prop))
  OTU.order.prop<-as.data.frame(otu_table(data.order.prop))
  OTU.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA

sample_OTU<-as.data.frame(t(sqrt(OTU.genus.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permg<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.fam.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permf<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.class.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permc<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.order.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permo<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.phy.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permp<-perm

beta<-rbind(permg$aov.tab[1,],permf$aov.tab[1,],permo$aov.tab[1,],permc$aov.tab[1,],permp$aov.tab[1,])
row.names(beta)<-c("Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta,file="Beta_div_permanova_Bray_ADENOMA_VS_ADENOCARC.csv",quote=F,row.names = T)

##### PCoA

dir.create("PCoA")

# su ASV
{data.sqrt<-transform_sample_counts(data.prop, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Tumor") +
  geom_point(size=2.5) +
  geom_text(aes(label = sample_names(data)), show.legend=F,
            position = (position_nudge(y=0.015)), size=2, color="black") + 
  theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_sqrt_ASV_prop_BIOPSY.png", width = 6, height = 5, dpi=300)
# with ellipses
plot_ordination(data.sqrt, ordBC, color = "Tumor") +
  geom_point(size=2.5) + stat_ellipse() + theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_sqrt_ASV_prop_BIOPSY_ELLIPSES.png", width = 6, height = 5, dpi=300)

# Adenocarcinoma only, in order to focus on genders
data_subset<-subset_samples(data.prop, Tumor=="Adenocarcinoma")
{data.sqrt<-transform_sample_counts(data_subset, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +  theme_classic() +
  scale_color_manual(values = c("F"="pink","M"="lightblue"))+
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only patients with adenocarcinoma)", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_ADENOCARCINOMA_GENDER_BIOPSY.png", width = 6, height = 5, dpi=300)

# Adenoma only, in order to focus on genders
data_subset<-subset_samples(data.prop, Tumor=="Adenoma")
{data.sqrt<-transform_sample_counts(data_subset, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +  theme_classic() +
  scale_color_manual(values = c("F"="pink","M"="lightblue"))+
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only patients with adenoma)", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_ADENOMA_GENDER_BIOPSY.png", width = 6, height = 5, dpi=300)

###################### DA WITH DESEQ2 CRC vs AP (BIOPSY) ####################################

dir.create("DESEQ2")

if(length(sample_names(data))+3 != length(sample_names(data_bio_contam))){
  cat("\nSomething is odd here... Did you perform the decontamination step yet???\n\n")
  Sys.sleep(2)
}

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

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
  
  DEseq_data<-phyloseq_to_deseq2(d, ~Tumor)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Tumor", "Adenocarcinoma", "Adenoma"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("DESEQ2/DA_",t,"_ratio_Adenocarcinoma_vs_Adenoma.csv"), row.names = F, quote=F, na = "")
    r_level<-r
    r_level[, "Taxon"]<- rep(t)
    Res_tot<-rbind.data.frame(Res_tot,r_level)
    # single box plots
    target<-r[[t]]
    colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
    target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
    Table_DE<-psmelt(target)
    colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
    Table_DE$ASV<-NULL
    #Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
    assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
    # appending to unique box plot
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
write.csv(Res_tot, file="DESEQ2/Every_result_DESeq2.csv", row.names = F)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Tumor)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=12.9,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  labs(title= "Differently abundant Taxa in Biopsy samples between patients with Adenocarcinoma and Adenoma", y="Proportional Abundance", 
       fill="Tumor", x="")
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) # remove the "zoom" if not needed
ggsave(filename = "DESEQ2/DA_Tumor_every_result.png", width = 15, height = 8, dpi=300)
dev.off()


Table_tot2<-subset(Table_tot, Taxa %in% c("Order","Family","Genus")) # to remove redundant results
redundants<-c("Yersinaceae","Fusobacteriaceae", "Synergistales","Lactobacillales","Fusobacteriales","Campylobacterales","Gemellaceae")
Table_tot2<-subset(Table_tot2, ! Bacteria %in% redundants)
Table_tot2$Bacteria <-gsub("Peptostreptococcales-Tissierellales","Peptostreptococcales\nTissierellales",Table_tot2$Bacteria) 

ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Tumor)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=15)) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  #scale_y_sqrt(breaks=c(0.1, 0.5,seq(1,4,1),seq(6,20,2), seq(20, max(Table_tot2$Abundance)+2,3))) +
  labs(title= "Differently abundant Taxa in biopsy samples between patients with Adenocarcinoma and Adenoma", y="Proportional Abundance", 
       fill="Tumor", x="")
ggsave(filename = "DESEQ2/DA_Tumor_no_redundants_asse_Y_normale.png", width = 13, height = 8, dpi=300)
dev.off()

deseq_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Tumor)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=15)) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5,seq(1,4,1),seq(6,20,2), seq(20, max(Table_tot2$Abundance)+2,3))) +
  labs(title= "", y="Proportional Abundance", fill="pathological tissue", x="")

deseq_plot +  labs(title= "Differently abundant Taxa in Biopsy samples between patients with Adenocarcinoma and Adenoma")
ggsave(filename = "DESEQ2/DA_Tumor_no_redundants_zoom.png", width = 13.5, height = 8, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > DESEQ2/NB.txt ")


# for the correlations
Genera.DESEQ2_CRC_AP_biopsy<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


####### CIRCOTAXA PLOT TO PLOT DA LOG2FOLDCHANGE (STOOL)

Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
#to remove the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Taxon %in% c("Class","Phylum")), ]
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% redundants & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[! (Res_tot2$Order %in% c(redundants,"Fusobacteriales","Lactobacillales","Campylobacterales") & Res_tot2$Taxon=="Order"), ]
Res_tot2$Order<-gsub("-","\n",Res_tot2$Order, fixed = T)

circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")

png(file="DESEQ2/CircoTax_plot.png",width=1800,height=1800, res = 300)
circo + labs(title="Adenocarcinoma vs adenoma \n DA log2foldchange (biopsy)\n")
dev.off()

circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))
  
# diff abundances + logfoldchange
png(file="DESEQ2/Final_plot.png",width=4560, height=2150, res = 300)
ggarrange(deseq_plot, circo_no_m, nrow = 1, widths= c(1,0.58), labels = c("A)", "B)"))
dev.off()


################### ALFA AND BETA DIVERSITY ADENONOCARCINOMA STAGES (Biopsy) ############################

dir.create("CRC_Stages_Analysis/Alpha_diversity")

{pAlpha<-plot_richness(data.stages, measures=c("Shannon", "Observed"), x="Stadiation")
  pAlpha
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
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Stadiation, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Stadiation", title="Alpha diversity between stages I, II and III") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11)) +
  stat_compare_means(aes(group = Stadiation), label="p.format", method = "kruskal.test", label.x= 1.05, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="CRC_Stages_Analysis/Alpha_diversity/Alfa_diversity_with_Kruskal.png", width = 7,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data)
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Stadiation
kruskal.test(Obser_value$value~factor)

# pairwise comparisons
con <- file("CRC_Stages_Analysis/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_Stadiation.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Stadiation sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-filter(alphadt, variable=="Observed")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Shannon")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Evenness")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)
rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor)

###### BETA DIVERSITY

dir.create("CRC_Stages_Analysis")
dir.create("CRC_Stages_Analysis/Beta_diversity")

{ASV.prop<-as.data.frame(otu_table(data.stages.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.stages.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.stages.phy.prop))
}

#### PERMANOVA

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_Bray<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_Bray<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=as(sample_data(data.stages.prop),"data.frame")$Stadiation, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV

# exporting beta diversity
rm(con)
con<-file("CRC_Stages_Analysis/Beta_diversity/Beta_diversity_general_and_pairwise_between_Stadiation_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_p, perm_ASV)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
  disper<-vegan::betadisper(BC.dist,as(sample_data(data.stages.prop),"data.frame")$Stadiation)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
rm(con)
con<-file("CRC_Stages_Analysis/Beta_diversity/Beta_dispersion_General_and_Pairwise_between_Stadiation_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con)

######## PCoA BRAY CURTIS

# in base ad ASV
data.stages.prop.labels<-data.stages.prop
{data.stages.sqrt_prop<-transform_sample_counts(data.stages.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.stages.sqrt_prop, method = "bray")
  ordBC = ordinate(data.stages.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.stages.sqrt_prop, ordBC, color = "Stadiation") +
  scale_color_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Stadiation", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="CRC_Stages_Analysis/Beta_diversity/PCoA_Beta_diversity_Bray_Curtis_on_ASV.png", width = 8, height = 6, dpi=300)

######################## DA WITH DESEQ2 ON CRC STAGES (Biopsy) #######################

dir.create("CRC_Stages_Analysis/DA_with_DESeq2")
# Trimming under sum of 10, see DeSeq2 tutorial
data_pruned<- prune_taxa(taxa_sums(data.stages) > 10, data.stages)

# preparing new taxa vocabularies for plots (some ASV may change glomming after trimming)
{data.stages.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.stages.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.stages.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.stages.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.stages.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.prop <- transform_sample_counts(data.stages, function(ASV) ASV/sum(ASV)*100)
  data_pruned.phy.prop <- transform_sample_counts(data.stages.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.stages.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.stages.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.stages.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.stages.genus_pruned, function(ASV) ASV/sum(ASV)*100)
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
      d<-phyloseq_to_deseq2(data_pruned, ~Stadiation)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Stadiation)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Stadiation", nume[i], deno[i]), test="Wald")
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("II","I","II")
deno<-c("I","III","III")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
unlink("DE_results.txt")
DE_results<-DE_results[DE_results$BaseMean>50, ]
system(" echo 'Every result under the arbitrary threshold Basemean=50 has been removed in order to avoid the noisiest results' > CRC_Stages_Analysis/DA_with_DESeq2/NB_results_are_filtered.txt")
write.csv2(DE_results, file="CRC_Stages_Analysis/DA_with_DESeq2/DE_every_results.csv", quote=F, row.names = F)
write.xlsx(DE_results, file="CRC_Stages_Analysis/DA_with_DESeq2/DE_every_results.xlsx",showNA = F, col.names = T, row.names = F)

# box plots
{target<-DE_results[DE_results$Rank=="Genus","ASV"]
  target<-prune_taxa(target, data_pruned.genus.prop)
  tabella<-psmelt(target)
  tabella_g<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Family","ASV"]
  target<-prune_taxa(target, data_pruned.fam.prop)
  tabella<-psmelt(target)
  tabella_f<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Order","ASV"]
  target<-prune_taxa(target, data_pruned.order.prop)
  tabella<-psmelt(target)
  tabella_o<-tabella
  rm(tabella)
}

# {target<-DE_results[DE_results$Rank=="Class","ASV"]
#   target<-prune_taxa(target, data_pruned.class.prop)
#   tabella<-psmelt(target)
#   tabella_c<-tabella
#   rm(tabella)
# }

{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
  target<-prune_taxa(target, data_pruned.phy.prop)
  tabella<-psmelt(target)
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
# {tabella_c$Taxa<-"Classes"
#   tabella_c[,"Phylum"]<-NULL
#   colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"
# }
{tabella_p$Taxa<-"Phyla"
  colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"
}

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_p,tabella_o)#,tabella_c)

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw( ) + 
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust=1, size=10), 
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,30,2)))#,seq(22,max(tabella_tot$Abundance),3) ))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/EVERY_RESULT_DeSeq2_Stadiation.png", width = 12, height = 7, dpi=300)
dev.off()

Redund<-c("uncultured_ o Rhodospirillales","Porphyromonadaceae","Paenibacillaceae","Paenibacillales")

tabella_tot2<-subset(tabella_tot, ! Bacteria %in% Redund)
ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=10),
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") + 
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))#,seq(22,max(tabella_tot$Abundance),3)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/EVERY_RESULT_no_redundants_Stadiation.png", width = 10, height = 7, dpi=300)
dev.off()

################# now plotting result in pairwise (group vs group)

##### II vs I
tabella_tot3<-subset(tabella_tot2, Stadiation != "III")
I_vs_II<-DE_results[DE_results$num=="II" & DE_results$denom=="I", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% I_vs_II [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, I_vs_II[I_vs_II$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages I and II", y="Proportional Abundance", fill="Stadiation", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_I_vs_II_no_redundants.png", width = 10, height = 7, dpi=300)

#####  I vs III
tabella_tot3<-subset(tabella_tot2, Stadiation != "II")
III_vs_I<-DE_results[DE_results$num=="I" & DE_results$denom=="III", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% III_vs_I [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, III_vs_I[III_vs_I$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12),
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages I and III", y="Proportional Abundance", fill="Stadiation", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_I_vs_III_DeSeq2_no redundants.png", width = 10, height = 7, dpi=300)

##### II vs III
tabella_tot3<-subset(tabella_tot2, Stadiation != "I")
III_vs_II<-DE_results[DE_results$num=="II" & DE_results$denom=="III", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% III_vs_II [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, III_vs_II[III_vs_II$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages II and III", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_II_vs_III_DeSeq2_no redundants.png", width = 10, height = 7, dpi=300)


################### ANALYSING LEFSE RESULTS (BIOPSY) ###########################

dir.create("PICRUST2_LEFSE")

suppressWarnings(rm(a))
a <- read.delim("../../QIIME2/PICRUST2_LEFSE_biopsy/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("../../QIIME2/PICRUST2_LEFSE_biopsy/Result_LEFSE.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="PICRUST2_LEFSE/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)

# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Adenocarcinoma"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="Adenocarcinoma"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Adenoma",1,0), size=3.2) +
  scale_fill_manual(values=c("Adenocarcinoma"="coral", "Adenoma"="deepskyblue")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.margin = margin(-10,0,0,0)) +
  theme(legend.position = "bottom")
ggsave(filename = "PICRUST2_LEFSE/PICRUST2_LEFSE_plot_diff_METACYC_every_result.png", width = 9, height = 6, dpi = 300)

system('touch PICRUST2_LEFSE/no_result_over_3')

rm(Significative_functions_LEFSE, a)


############## ANALYSIS ON FECAL METABOLITES BIOPSY DATA (CRC vs AP) ###############

dir.create("Analysis_on_metabolites")

Metabolites<-read.csv(file="../../Metabolites.csv", row.names = 1)
Fatty<-read.csv(file="../../Fatty_acid.csv", row.names = 1)
aa<-read.csv(file="../../Aminoacids.csv", row.names = 1)
# the row.names are the "Patient_ID" (identificative number)

head(Metabolites,n=2)
head(Fatty,n=2)
head(aa,n=2)

# the significative metabolites have been already found in "analysis on fecal metabolites STOOL data" paragraph of this script
significative<-c("Leucine","Tyrosine","Alanine","Valine","Lactate")
aa2 <- cbind.data.frame(aa,Metabolites)[,significative]

#### with biopsy DESeq2 results
generi_sel <- subset_taxa(data.genus.prop, Genus %in% Genera.DESEQ2_CRC_AP_biopsy)
identical(row.names(otu_table(generi_sel)),row.names(tax_table(generi_sel)))
{Names<-as.data.frame(tax_table(generi_sel))
  Names<-Names$Genus
  corr_tot<- t(as.data.frame(otu_table(generi_sel)))
  colnames(corr_tot)<-Names
  row.names(corr_tot)<-gsub("CRC","",row.names(corr_tot))
  row.names(corr_tot)<-gsub("AP","",row.names(corr_tot))
}

sample_order<-as.character(sample_data(data.genus)$Patient_number[sample_data(data.genus)$Patient_number %in% rownames(aa2)])
corr_tot<-corr_tot[sample_order,]
l<-length(row.names(corr_tot))
l # 25
if(l!=25){cat("\nWait, did you perform the decontamination step already? Three biopsy microbiota sample should have been removed\n", fill=T)}

aa2<-aa2[sample_order,]
identical(l, length(which(! is.na(row.names(corr_tot)))) ) # TRUE --> any sample has been cutted out
identical(row.names(corr_tot),row.names(aa2))

{x<-as.matrix(cbind(corr_tot,aa2))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  identical(data_corr[,1:2],data_pvalue[,1:2])
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Bacteria","Aminoacids","Corr","pvalue")
  data_corr<-data_corr[data_corr$Bacteria %in% colnames(corr_tot), ]
  data_corr<-data_corr[data_corr$Aminoacids %in% colnames(aa2), ]
  data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
data_corr$Bacteria<-gsub("[","",data_corr$Bacteria, fixed = T)
data_corr$Bacteria<-gsub("]","",data_corr$Bacteria, fixed = T)
write.csv2(data_corr, file="Analysis_on_metabolites/Spearman_of_DESeq2_Results_vs_aminoacids.csv", quote=F, row.names = F)

data_corr$simbol<-data_corr$padj
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Bacteria, y = data_corr$Aminoacids, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=8) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 7)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.1) +
  labs(title = "Spearman correlation between \n differently abundant bacterial genera and metabolites in biopsy sample\n between CRC and AP patients",
       y= "", x= "", caption= "* = p-adjusted < 0.05") +
  theme(plot.title = element_text(size=10)) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=8), legend.key.height= unit(1.1, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=10))
ggsave(filename = "Analysis_on_metabolites/Heatmap_DESeq2_results_vs_metab.png", width = 6, height = 4.5, dpi=300)


################ CORRELATIONS BETWEEN DESEQ2 RESULTS OF BIOPSY AND STOOL DATA #####################

if( grepl("_analysis",getwd()) | grepl("_correlations",getwd()) ){setwd("..")}
dir.create("Biopsy_Stool_DA_correlations")
setwd("Biopsy_Stool_DA_correlations")

to_remove<-ls()
to_remove<-to_remove[ ! to_remove %in% c("data_total","metadata","evalslopes", "CircoTax2", "Genera.DESEQ2_CRC_AP_stool","Genera.DESEQ2_CRC_AP_biopsy")]
to_remove
rm(list = to_remove) # reset the workspace

data<-subset_samples(data_total,Sampling_Type %in% c("B","F"))
both_sampling<-sample_data(data)$Patient[duplicated(sample_data(data)$Patient)]
length(both_sampling)
data<-subset_samples(data, Patient %in% both_sampling)
identical(as.numeric(length(sample_names(data))), length(both_sampling)*2 ) # TRUE
head(sample_data(data), n=2)

### decontaminating the object (just as in biopsy section)
data<- subset_taxa(data, Kingdom!="Unassigned")

### removing from the whole object those samples that are anomalous in the biopsy data set
data<- subset_samples(data, ! Patient %in% c("CRC7", "AP19", "CRC54"))
# NB: only two of those three have also fecal samples --> one (AP19) is already removed
if(length(sample_names(data))+2 == (length(both_sampling)-1)*2 ) {
  cat("Two samples (--> 4 observations) have been removed correctly")
}

data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)


#### extracting rel abundances of Stool DESeq2 results
stool.prop<-subset_samples(data.genus.prop, Sampling_Type == "F")
sample_names(stool.prop)<-sample_data(stool.prop)$Patient
generi_sel_stool <- subset_taxa(stool.prop, Genus %in% c(Genera.DESEQ2_CRC_AP_stool))
identical(row.names(otu_table(generi_sel_stool)),row.names(tax_table(generi_sel_stool)))
{suppressWarnings(rm("Names","corr_tot"))
  Names<-as.data.frame(tax_table(generi_sel_stool))
  Names<-Names$Genus
  corr_tot<- t(as.data.frame(otu_table(generi_sel_stool)))
  colnames(corr_tot)<-Names
  row.names(corr_tot)<-gsub("CRC","",row.names(corr_tot))
  row.names(corr_tot)<-gsub("AP","",row.names(corr_tot))
  corr_stool<-corr_tot
}

#### extracting rel abundances of Biopsy DESeq2 results
biopsy.prop<-subset_samples(data.genus.prop, Sampling_Type == "B")
sample_names(biopsy.prop)<-sample_data(biopsy.prop)$Patient
generi_sel_biopsy <- subset_taxa(biopsy.prop, Genus %in% c(Genera.DESEQ2_CRC_AP_biopsy))
identical(row.names(otu_table(generi_sel_biopsy)),row.names(tax_table(generi_sel_biopsy)))
{suppressWarnings(rm("Names","corr_tot"))
  Names<-as.data.frame(tax_table(generi_sel_biopsy))
  Names<-Names$Genus
  corr_tot<- t(as.data.frame(otu_table(generi_sel_biopsy)))
  colnames(corr_tot)<-Names
  row.names(corr_tot)<-gsub("CRC","",row.names(corr_tot))
  row.names(corr_tot)<-gsub("AP","",row.names(corr_tot))
  corr_biopsy<-corr_tot
}

# ensuring the same order between the two data set before joining
corr_biopsy<-corr_biopsy[row.names(corr_stool),]
identical(rownames(corr_biopsy), row.names(corr_stool)) # TRUE
length(which(is.na(corr_biopsy[,1])))==0 # TRUE

{x<-as.matrix(cbind(corr_stool, corr_biopsy))
  r<-rcorr(x, type = "spearman")
  correlation<-as.data.frame(r$r)
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  identical(data_corr[,1:2],data_pvalue[,1:2])
  data_corr<-cbind(data_corr,data_pvalue[,3])
  colnames(data_corr)<-c("Biopsy_res","Stool_res","Corr","pvalue")
  data_corr<-data_corr[data_corr$Biopsy_res %in% Genera.DESEQ2_CRC_AP_biopsy, ]
  data_corr<-data_corr[data_corr$Stool_res %in% Genera.DESEQ2_CRC_AP_stool, ]
  data_corr$padj<-p.adjust(data_corr$pvalue, method = "BH")
}
data_corr$Stool_res<-gsub("[","",data_corr$Stool_res, fixed = T)
data_corr$Stool_res<-gsub("]","",data_corr$Stool_res, fixed = T)
write.csv2(data_corr, file="Spearman_of_DESeq2_Results_Biopsy_vs_stool.csv", quote=F, row.names = F)

data_corr$simbol<-data_corr$padj
data_corr$simbol[data_corr$simbol<0.05]<-"*"
data_corr$simbol[data_corr$simbol>0.05]<-" "
ggplot(data_corr, aes(x = data_corr$Biopsy_res, y = data_corr$Stool_res, fill = data_corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5,linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=8) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 7)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= simbol), color= "white", size =6, nudge_y = -0.1) +
  labs(y= "Stool DA genera (Stool relative abundances)", 
       x= "Biopsy DA genera (Biopsy relative abundances)", 
       caption= "* = p-adjusted < 0.05") +
  theme(plot.title = element_text(size=10)) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0)) +
  theme(legend.text = element_text(size=8), legend.key.height= unit(1.1, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=10))
ggsave(filename = "Heatmap_DESeq2_results_stool_vs_biopsy.png", width = 6, height = 4.5, dpi=300)



##################### STARTING SALIVA ANALYSIS ############################

if( grepl("_analysis",getwd()) | grepl("_correlations",getwd())  ){setwd("..")}
dir.create("Saliva_analysis")
setwd("Saliva_analysis")

to_remove<-ls()
to_remove<-to_remove[ ! to_remove %in% c("data_total","metadata","CircoTax2")]
to_remove
rm(list = to_remove) # reset the workspace

data<-subset_samples(data_total,Sampling_Type=="S")
head(sample_data(data))
sample_names(data)<-sample_data(data)$Patient
head(sample_data(data))

data.genus<- tax_glom(data, taxrank = "Genus", NArm = F)
OTU.genus<-as.data.frame(otu_table(data.genus))
Taxa.genus<-as.data.frame(tax_table(data.genus))

data.fam<- tax_glom(data, taxrank = "Family", NArm = F)
OTU.fam<-as.data.frame(otu_table(data.fam))
Taxa.fam<-as.data.frame(tax_table(data.fam))

data.class<- tax_glom(data, taxrank = "Class", NArm = F)
OTU.class<-as.data.frame(otu_table(data.class))
Taxa.class<-as.data.frame(tax_table(data.class))

data.order<- tax_glom(data, taxrank = "Order", NArm = F)
OTU.order<-as.data.frame(otu_table(data.order))
Taxa.order<-as.data.frame(tax_table(data.order))

data.phy<- tax_glom(data, taxrank = "Phylum", NArm = F)
OTU.phy<-as.data.frame(otu_table(data.phy))
Taxa.phy<-as.data.frame(tax_table(data.phy))

data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)

dir.create("CRC_Stages_Analysis")
data.stages<-subset_samples(data, !is.na(Stadiation)) # only samples with adenocarcinoma, in order to confront stadiations
sample_data(data.stages)$Stadiation<-factor(sample_data(data.stages)$Stadiatio, levels=c("I","II","III"))
{data.stages.prop<-transform_sample_counts(data.stages, function(otu) otu/sum(otu)*100)
  data.stages.genus<- tax_glom(data.stages, taxrank = "Genus", NArm = F)
  data.stages.genus.prop <- transform_sample_counts(data.stages.genus, function(otu) otu/sum(otu)*100)
  data.stages.phy<- tax_glom(data.stages, taxrank = "Phylum", NArm = F)
  data.stages.phy.prop <- transform_sample_counts(data.stages.phy, function(otu) otu/sum(otu)*100)
}

#################### EXPORTING ABUNDANCES (SALIVA) ############################

dir.create("Abundances")

{options( scipen=100 ) 
  dir.create("Abundances/Percent_abundances")
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Percent_abundances/counts_phylum.xlsx",row.names = F, showNA = F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Percent_abundances/counts_genus.xlsx",row.names = F, showNA = F)
  write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy.prop),"matrix")),file="Abundances/Percent_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class.prop),"matrix")),file="Abundances/Percent_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order.prop),"matrix")),file="Abundances/Percent_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam.prop),"matrix")),file="Abundances/Percent_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus.prop),"matrix")),file="Abundances/Percent_abundances/counts_genus.csv",quote=F)
}

### TOP 5 Orders
{top5 <- names(sort(taxa_sums(data.order.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.order.prop)
  others<-taxa_names(data.order.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.order.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Order<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Order)) + 
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every order below rank 5 ") +
  ggtitle("Five most abundant orders and other taxa")
ggsave(file="Abundances/TOP5_Orders.png",width=6.5,height=4.5, dpi=300)
dev.off()

# TOP 5 Phyla 
{top5 <- names(sort(taxa_sums(data.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.phy.prop)
  others<-taxa_names(data.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="Abundances/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Genera
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  others<-taxa_names(data.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Genus)) + 
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every genus below rank 5 ") +
  ggtitle("Five most abundant genera and other taxa")
ggsave(file="Abundances/TOP5_Genera.png",width=6.5,height=4.5, dpi=300)
dev.off()

# TOP 5 Families
{top <- names(sort(taxa_sums(data.fam.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.fam.prop)
  others<-taxa_names(data.fam.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.fam.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Family<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Family)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every family below rank 5 ") +
  ggtitle("Five most abundant families and other taxa")
ggsave(file="Abundances/TOP5_families.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Classes
{top <- names(sort(taxa_sums(data.class.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.class.prop)
  others<-taxa_names(data.class.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.class.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Class<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Class)) +
  facet_grid2(cols= vars(Tumor), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every class below rank 5 ") +
  ggtitle("Five most abundant classes and other taxa")
ggsave(file="Abundances/TOP5_classes.png",width=6.5,height=4.5,dpi=300)
dev.off()

dir.create("CRC_Stages_Analysis/Abundances")

write.xlsx(cbind(as(otu_table(data.stages.phy.prop),"matrix"),as(tax_table(data.stages.phy.prop),"matrix")),file="CRC_Stages_Analysis/Abundances/Percent_counts_phylum.xlsx",row.names = F, showNA = F)
write.xlsx(cbind(as(otu_table(data.stages.genus.prop),"matrix"),as(tax_table(data.stages.genus.prop),"matrix")),file="CRC_Stages_Analysis/Abundances/Percent_counts_genus.xlsx",row.names = F, showNA = F)

# TOP 5 Phyla (Only CRC Stages)
{top5 <- names(sort(taxa_sums(data.stages.phy.prop), decreasing=TRUE))[1:5]
  prune.dat_top5 <- prune_taxa(top5,data.stages.phy.prop)
  others<-taxa_names(data.stages.phy.prop)
  others<-others[!(others %in% top5)]
  prune.data.others<-prune_taxa(others,data.stages.phy.prop)
  tabella_top<-psmelt(prune.dat_top5)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Phylum)) +
  facet_grid2(cols= vars(Stadiation), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every phylum below rank 5 ") +
  ggtitle("Five most abundant phyla and other taxa")
ggsave(file="CRC_Stages_Analysis/Abundances/TOP5_phyla.png",width=6.5,height=4.5,dpi=300)
dev.off()

# TOP 5 Genera (Only CRC Stages)
{top <- names(sort(taxa_sums(data.stages.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.stages.genus.prop)
  others<-taxa_names(data.stages.genus.prop)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.stages.genus.prop)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
}
ggplot(data=tabella, aes(x=Patient, y=Abundance, fill=Genus)) + 
  facet_grid2(cols= vars(Stadiation), scales = "free_x", space = "free") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size = 11) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 8), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(y="Relative abundance", 
       caption = " 'Others' includes every genus below rank 5 ") +
  ggtitle("Five most abundant genera and other taxa")
ggsave(file="CRC_Stages_Analysis/Abundances/TOP5_Genera.png",width=6.5,height=4.5, dpi=300)
dev.off()

################### ALFA AND BETA DIVERSITY CRC vs AP (SALIVA) ###################

# NO normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Chao1", "Shannon", "Observed"), x="Tumor", color = "Patient") # prendo valori dal grafico
pAlpha
{Chao<-dplyr::filter(pAlpha$data, variable=="Chao1")
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evennes
  identical(H$Codice,obs$Codice) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
}
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Tumor, y=value, color=NULL), alpha=0.1) + theme_bw() + labs(x="") +
  guides(fill="none", color="none") +
  theme(axis.text.x = element_text(angle = -20, vjust=1, hjust = 0) ) +
  stat_compare_means(aes(group = Tumor), label="p.format", method = "wilcox.test", label.x= 1.4, size=3.2, label.y.npc = "top", vjust=-0.5)
ggsave(file="Alpha_diversity_adenoma_vs_adenoc.png", width = 6,height = 5, dpi=300)

########################## BETA DIVERSITY

{OTU.genus.prop<-as.data.frame(otu_table(data.genus.prop))
  OTU.fam.prop<-as.data.frame(otu_table(data.fam.prop))
  OTU.class.prop<-as.data.frame(otu_table(data.class.prop))
  OTU.order.prop<-as.data.frame(otu_table(data.order.prop))
  OTU.phy.prop<-as.data.frame(otu_table(data.phy.prop))
}

#### PERMANOVA

sample_OTU<-as.data.frame(t(sqrt(OTU.genus.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permg<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.fam.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permf<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.class.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permc<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.order.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permo<-perm

sample_OTU<-as.data.frame(t(sqrt(OTU.phy.prop))) # campioni lungo le righe
perm<- adonis(sample_OTU ~Tumor, data=as(sample_data(data),"data.frame"), permutations = 9999, method="bray")
permp<-perm

beta<-rbind(permg$aov.tab[1,],permf$aov.tab[1,],permo$aov.tab[1,],permc$aov.tab[1,],permp$aov.tab[1,])
row.names(beta)<-c("Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta,file="Beta_div_permanova_Bray_ADENOMA_VS_ADENOCARC.csv",quote=F,row.names = T)

##### PCoA

dir.create("PCoA")

# su ASV
{data.sqrt<-transform_sample_counts(data.prop, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Tumor") +
  geom_point(size=2.5) +
  geom_text(aes(label = sample_names(data)), show.legend=F,
            position = (position_nudge(y=0.015)), size=2, color="black") + 
  theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_sqrt_ASV_prop_SALIVA.png", width = 6, height = 5, dpi=300)
# with ellipses
plot_ordination(data.sqrt, ordBC, color = "Tumor") +
  geom_point(size=2.5) + stat_ellipse() + theme_classic() +
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_sqrt_ASV_prop_SALIVA_ELLIPSE.png", width = 6, height = 5, dpi=300)

# Adenocarcinoma only. in order to focus on genders
data_subset<-subset_samples(data.prop, Tumor=="Adenocarcinoma")
{data.sqrt<-transform_sample_counts(data_subset, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +  theme_classic() +
  scale_color_manual(values = c("F"="pink","M"="lightblue"))+
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only patients with adenocarcinoma)", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_ADENOCARCINOMA_GENDER_SALIVA.png", width = 6, height = 5, dpi=300)

# Adenoma only, in order to focus on genders
data_subset<-subset_samples(data.prop, Tumor=="Adenoma")
{data.sqrt<-transform_sample_counts(data_subset, sqrt)
  DistBC = phyloseq::distance(data.sqrt, method = "bray")
  ordBC = ordinate(data.sqrt, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt, ordBC, color = "Gender") +
  geom_point(size=2.5) + stat_ellipse() +  theme_classic() +
  scale_color_manual(values = c("F"="pink","M"="lightblue"))+
  labs(title="PCoA with Bray-Curtis distance \n calculated on squared root ASV percent abundance\n (only patients with adenoma)", 
       color="Tumor", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="PCoA/PCoA_Bray_ADENOMA_GENDER_SALIVA.png", width = 6, height = 5, dpi=300)

###################### DA WITH DESEQ2 CRC vs AP (Saliva) ####################################

dir.create("DESEQ2")

suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

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
  
  DEseq_data<-phyloseq_to_deseq2(d, ~Tumor)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Tumor", "Adenocarcinoma", "Adenoma"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("DESEQ2/DA_",t,"_ratio_Adenocarcinoma_vs_Adenoma.csv"), row.names = F, quote=F, na = "")
    r_level<-r
    r_level[, "Taxon"]<- rep(t)
    Res_tot<-rbind.data.frame(Res_tot,r_level)
    # single box plots
    target<-r[[t]]
    colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
    target<-subset_taxa(d.prop, Aimed_taxa %in% target) # cannot use t %in% target in this function, then it's scripted in this way
    Table_DE<-psmelt(target)
    colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
    Table_DE$ASV<-NULL
    #Table_DE$Abundance<-sqrt(Table_DE$Abundance) # then sqrt of proportion
    assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
    # appending to unique box plot
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

system(" echo 'After removing every result under the arbitrary threshold of basemean=100 there is anything left' > DESEQ2/Nothing.txt ")

################### ALFA AND BETA DIVERSITY ADENONOCARCINOMA STAGES (Saliva) ############################

dir.create("CRC_Stages_Analysis/Alpha_diversity")

{pAlpha<-plot_richness(data.stages, measures=c("Shannon", "Observed"), x="Stadiation")
  pAlpha
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
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Stadiation, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Stadiation", title="Alpha diversity between stages I, II and III") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11)) +
  stat_compare_means(aes(group = Stadiation), label="p.format", method = "kruskal.test", label.x= 1.05, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="CRC_Stages_Analysis/Alpha_diversity/Alfa_diversity_with_Kruskal.png", width = 7,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data)
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
factor<-Obser_value$Stadiation
kruskal.test(Obser_value$value~factor)

# pairwise comparisons
con <- file("CRC_Stages_Analysis/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_Stadiation.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Stadiation sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-filter(alphadt, variable=="Observed")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Shannon")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Evenness")
a<-FSA::dunnTest(value ~ Stadiation, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)
rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor)

###### BETA DIVERSITY

dir.create("CRC_Stages_Analysis")
dir.create("CRC_Stages_Analysis/Beta_diversity")

{ASV.prop<-as.data.frame(otu_table(data.stages.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.stages.genus.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.stages.phy.prop))
}

#### PERMANOVA

{sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
  perm_ASV<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_Bray<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_Bray<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Stadiation, data=as(sample_data(data.stages.prop),"data.frame"), permutations = 9999, method="bray")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta

system(" echo 'The result of PERMANOVA is significant ONLY on ASV level and not after the collapsing at genus level, considering also the low sample number and the environment variability I would define this as a REALLY weak result'   > CRC_Stages_Analysis/Beta_diversity/NB.txt")

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=as(sample_data(data.stages.prop),"data.frame")$Stadiation, p.adjust.m = "BH", sim.method="bray", perm = 9999)
pair_ASV

# exporting beta diversity
rm(con)
con<-file("CRC_Stages_Analysis/Beta_diversity/Beta_diversity_general_and_pairwise_between_Stadiation_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity Bray Curtis \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity Bray-Curtis (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_p, perm_ASV)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
  disper<-vegan::betadisper(BC.dist,as(sample_data(data.stages.prop),"data.frame")$Stadiation)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
rm(con)
con<-file("CRC_Stages_Analysis/Beta_diversity/Beta_dispersion_General_and_Pairwise_between_Stadiation_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion Bray Curtis \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con)

######## PCoA BRAY CURTIS

# in base ad ASV
data.stages.prop.labels<-data.stages.prop
{data.stages.sqrt_prop<-transform_sample_counts(data.stages.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.stages.sqrt_prop, method = "bray")
  ordBC = ordinate(data.stages.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.stages.sqrt_prop, ordBC, color = "Stadiation") +
  scale_color_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="Stadiation", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="CRC_Stages_Analysis/Beta_diversity/PCoA_Beta_diversity_Bray_Curtis_on_ASV.png", width = 8, height = 6, dpi=300)

######################## DA WITH DESEQ2 ON CRC STAGES (Saliva) #########################

dir.create("CRC_Stages_Analysis/DA_with_DESeq2")
# Trimming under sum of 10, see DeSeq2 tutorial
data_pruned<- prune_taxa(taxa_sums(data.stages) > 10, data.stages)

# preparing new taxa vocabularies for plots (some ASV may change glomming after trimming)
{data.stages.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.stages.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.stages.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.stages.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.stages.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.prop <- transform_sample_counts(data.stages, function(ASV) ASV/sum(ASV)*100)
  data_pruned.phy.prop <- transform_sample_counts(data.stages.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.stages.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.stages.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.stages.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.stages.genus_pruned, function(ASV) ASV/sum(ASV)*100)
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
      d<-phyloseq_to_deseq2(data_pruned, ~Stadiation)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Stadiation)
    }
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Stadiation", nume[i], deno[i]), test="Wald")
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
      }
    }
  }
}

nume<-c("II","I","II")
deno<-c("I","III","III")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.05)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
unlink("DE_results.txt")
DE_results<-DE_results[DE_results$BaseMean>50, ]
system(" echo 'Every result under the arbitrary threshold Basemean=50 has been removed in order to avoid the noisiest results' > CRC_Stages_Analysis/DA_with_DESeq2/NB_results_are_filtered.txt")
write.csv2(DE_results, file="CRC_Stages_Analysis/DA_with_DESeq2/DE_every_results.csv", quote=F, row.names = F)
write.xlsx(DE_results, file="CRC_Stages_Analysis/DA_with_DESeq2/DE_every_results.xlsx",showNA = F, col.names = T, row.names = F)

# box plots

# {target<-DE_results[DE_results$Rank=="Genus","ASV"]
#   target<-prune_taxa(target, data_pruned.genus.prop)
#   tabella<-psmelt(target)
#   tabella_g<-tabella
#   rm(tabella)
# }

{target<-DE_results[DE_results$Rank=="Family","ASV"]
  target<-prune_taxa(target, data_pruned.fam.prop)
  tabella<-psmelt(target)
  tabella_f<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Order","ASV"]
  target<-prune_taxa(target, data_pruned.order.prop)
  tabella<-psmelt(target)
  tabella_o<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Class","ASV"]
  target<-prune_taxa(target, data_pruned.class.prop)
  tabella<-psmelt(target)
  tabella_c<-tabella
  rm(tabella)
}

{target<-DE_results[DE_results$Rank=="Phylum","ASV"]
  target<-prune_taxa(target, data_pruned.phy.prop)
  tabella<-psmelt(target)
  tabella_p<-tabella
  rm(tabella)
}

# unique boxplot of DeSeq2 results
# {tabella_g$Taxa<-"Genera"
#   tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
#   colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"
# }
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

tabella_tot<-rbind.data.frame(tabella_f,tabella_p,tabella_o,tabella_c)

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw( ) + 
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust=1, size=10), 
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,30,2)))#,seq(22,max(tabella_tot$Abundance),3) ))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/EVERY_RESULT_DeSeq2_Stadiation.png", width = 12, height = 7, dpi=300)
dev.off()

tabella_tot[tabella_tot$Bacteria=="Gracilibacteria" & tabella_tot$Taxa=="Classes","Bacteria"]<-"redundant"
tabella_tot[tabella_tot$Bacteria=="Absconditabacteriales_(SR1)" & tabella_tot$Taxa=="Orders","Bacteria"]<-"redundant"

tabella_tot2<-subset(tabella_tot, ! Bacteria %in% "redundant")
ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=10),
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") + 
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))#,seq(22,max(tabella_tot$Abundance),3)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/EVERY_RESULT_no_redundants_Stadiation.png", width = 10, height = 7, dpi=300)
dev.off()

################# now plotting result in pairwise (group vs group)

##### II vs I
tabella_tot3<-subset(tabella_tot2, Stadiation != "III")
I_vs_II<-DE_results[DE_results$num=="II" & DE_results$denom=="I", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% I_vs_II [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, I_vs_II[I_vs_II$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("I"="coral","II"="red")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) + 
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages I and II", y="Proportional Abundance", fill="Stadiation", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_I_vs_II_no_redundants.png", width = 10, height = 7, dpi=300)

#####  I vs III
tabella_tot3<-subset(tabella_tot2, Stadiation != "II")
III_vs_I<-DE_results[DE_results$num=="I" & DE_results$denom=="III", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% III_vs_I [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, III_vs_I[III_vs_I$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)
# nothing there

##### II vs III
tabella_tot3<-subset(tabella_tot2, Stadiation != "I")
III_vs_II<-DE_results[DE_results$num=="II" & DE_results$denom=="III", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% III_vs_II [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, III_vs_II[III_vs_II$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Stadiation)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("II"="red","III"="red3")) +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12)) +
  theme(plot.title= element_text(size=16) ,legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=14)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between stages II and III", y="Proportional Abundance", fill="Stadiation", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,1)))
ggsave(filename = "CRC_Stages_Analysis/DA_with_DESeq2/Results_only_II_vs_III_DeSeq2_no redundants.png", width = 10, height = 7, dpi=300)

##################### R AND PACKAGES VERSION #########################

setwd("..")

package<-sessionInfo()

con <- file("R_version_and_packages.txt")
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

sink()
close(con)