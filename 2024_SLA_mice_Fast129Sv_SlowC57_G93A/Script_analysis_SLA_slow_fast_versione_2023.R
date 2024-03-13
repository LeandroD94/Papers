##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggvenn")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  library("dendextend")
  # analysis packages
  library("DESeq2")
  library("vegan")
  # utilities
  library("stringr")
  library("xlsx")  
  library("dplyr")
  library("qiime2R")
}


{dir.create("Data_check")
  dir.create("Data_check/PCoA_test")
  dir.create("Results")
}

Subsets<-c("PS","OS","SY") # to create each folder
for(x in Subsets){
  dir.create(paste0('Results/',x))
  dir.create(paste0('Results/',x,'/Abundances'))
  dir.create(paste0('Results/',x,'/Ratio_Firmi_Bacteroi'))
  dir.create(paste0('Results/',x,'/Hierarchical_clustering'))
  dir.create(paste0('Results/',x,'/Beta_div_Genotype'))
  dir.create(paste0('Results/',x,'/Beta_div_Ntg_vs_G93A'))
  dir.create(paste0('Results/',x,'/DA_DESeq2_Genotype/'))
  dir.create(paste0('Results/',x,'/DA_DESeq2_Ntg_vs_G93A/'))
}
suppressWarnings(rm(Subsets,x))

{ dir.create('Results/Diff_along_time')
  dir.create('Results/Diff_along_time/Abundances')
  dir.create('Results/Diff_along_time/Ratio_Firmi_Bacteroi')
  dir.create('Results/Diff_along_time/Alfa_div')
  dir.create('Results/Diff_along_time/Beta_div')
  dir.create('Results/Diff_along_time/DA_DESeq2')
}

dir.create("Results/General_Beta_div")

options(scipen = 100) # disable scientific annotation



### function from https://github.com/matteoramazzotti
CircoTax2=function(file,title="CircoTax plot",ramp=c("orange","white","blue"),tax_col=7:11,fc_col=2,sort=c("no","rank","fc","absfc","alpha"),sort_dir="d") {
  data=file
  sort_dir=ifelse(sort_dir == "d",FALSE,TRUE)
  if(length(tax_col) == 5) {
    gplot_labels= unlist(strsplit("PCOFG",""))
  }
  if(length(tax_col) == 6) { #KPCOFG as in RDP
    gplot_labels= unlist(strsplit("DPCOFG",""))
  }
  if(length(tax_col) == 7) { #KPCOFGS as in silva
    gplot_labels= unlist(strsplit("DPCOFGS",""))
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
    # geom_text( 
    geom_text_aimed(  # it allows to calculate the ideal angle to have labels perfectly perpendicular to the axis
      aes(
        x= id,
        y=7,
        label=name),
      color="black",
      fontface="bold",
      alpha=0.6,
      size=3
    ) +
    geom_hline( 
      yintercept=1:6,
      color="grey"
    ) +
    coord_polar(start = 0, clip="off") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = NA),
      plot.margin = unit(rep(1,4), "cm"),      # It adjusts the margins of the plot to avoid the trimming of the labels!
      plot.title = element_text(hjust = 0.5,face="bold", size=18)
    ) +
    labs(fill="log2FC")
}



####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample
sample<-substring(sample, first = 10, last= 15) # strtrim to cut until X
sample<-gsub("-A","",sample)
sample<-gsub("-P","",sample)
sample<-gsub("-","",sample, fixed=T)
sample_names(data)<-sample # update

Metadata <- as.data.frame(read.table("Metadata.tsv",sep="\t",header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)

sample_data(data)$Genotype<-factor(sample_data(data)$Genotype, levels=c("C57Ola","129Sv")) # decides the order in plots
head(sample_data(data))


########## IMPORTING ALSO THE CONTAMINATED ONE (FROM HOST DNA)

data_contam<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/eucar_contaminants/contam_taxonomy.qza", tree = "QIIME/eucar_contaminants/rooted-tree_with_contam.qza")
# changing names (already done before, check up there!)
if(identical(original_names,sample_names(data_contam))){
  sample_names(data_contam)<-sample_names(data)
}
if(identical(sample_names(data_contam), row.names(Metadata))) {
  sample_data(data_contam)<-Metadata
  cat("Everything alright there!")} else { cat("\n Sample names are not matching, check them \n\n") }


#################### FILTERING NOISES FROM DATA SET ####################

if(! "filter_proof" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus")
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_mean_0005_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured"] # to avoid the removal of other uncultured genera
data<-subset_taxa(data, ! Genus %in% filtered)
data<-prune_taxa(taxa_sums(data) > 0, data) 
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


# removing chloroplast and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria\n\n")
  data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mitochondria"))
}


filter_proof<- "Marker of the filtering, it is required for the script"

# save.image("data_mice_noise_filtered2023.RData")



################# CHECKING THE COMPOSITION AFTER FILTERING ####################

data.unf.prop<-transform_sample_counts(data_contam, function(x) (x/sum(x))*100)


if(! "filter_proof" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}
data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100)


################# BRAY HELLINGER

##### unfiltered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.unf.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3600, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


################# sqrt BRAY PROPORTIONAL AB

##### unfiltered sqrt BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.unf.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "bray")
  DistBC = sqrt(DistBC)
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.prop.labels, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA sqrt Bray-Curtis (on proportional ASV)\n\n UNfiltered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.prop.labels))
data.prop.labels<-data.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "bray")
  DistBC = sqrt(DistBC)
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_sqrt_BRAY_test.png", width = 3600, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


################# EUCLIDEAN HELLINGER

##### unfiltered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
data.prop.labels<-data.unf.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3600, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


################# wUNIFRAC PROP

##### unfiltered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop, p1))
data.prop.labels<-data.unf.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.prop.labels, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered wUnifrac
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Genotype", shape="Time_point") +
  scale_color_manual(values=c("129Sv"="coral","C57Ola"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$FASTQ_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_test.png", width = 3600, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


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
r<-rarecurve(as(t(otu_table(data)),"matrix"), step=50,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)



############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

if(! "filter_proof" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
Original_number<-sum( c(Original_number$forward.sequence.count, Original_number$reverse.sequence.count))
Remaining_number<-sum(otu_table(data))

con<-file("Data_check/DETAILS_ABOUT_ORIGINAL_AND_PROCESSED_READS_NUMBER.txt")
sink(con, append=TRUE)
cat("Raw reads abundance (sequenced, considering the pairs as one)", fill=TRUE)
cat(Original_number/2)
cat("Reads abundance after every filter", fill=TRUE)
cat(Remaining_number)
cat("Percentual of remaining vs original", fill=TRUE)
cat(paste0(round((Remaining_number/(Original_number/2))*100,digits=2)),"%")
sink()
close(con)
rm(con)

### barplot of reads
Original_read_number<-as.data.frame(read.table(file="QIIME/Original_number_of_reads_for_sample.tsv", header = T, sep="\t", row.names = 1))
DADA2_read_number<-as.data.frame(read.table(file="QIIME/Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
after_filter_number<-as.data.frame(colSums(otu_table(data)))
# updating names from FASTQ codes to sample names
Original_read_number<-Original_read_number[original_names, ] # this vector came from the import section (ordered as "sample", which contains the modified names)
DADA2_read_number<-DADA2_read_number[original_names, ]
row.names(Original_read_number)<-sample # same order, already modified
row.names(DADA2_read_number)<-sample # same order, already modified
# same order of rows
Original_read_number<-Original_read_number[row.names(DADA2_read_number), ]
after_filter_number<-after_filter_number[row.names(DADA2_read_number), ]
# creating the table
table<-cbind.data.frame(Original=Original_read_number$forward.sequence.count, # merged --> just one column, not both of them
                        After_quality_filter=DADA2_read_number$non.chimeric,
                        After_contam_filter=after_filter_number)
# adding factors to re-group the table
table$factor<-as.data.frame(sample_data(data))[sample, ]$Batch
table$factor<-gsub("_o", "o batch", table$factor)
table$Samples<-row.names(Original_read_number)
# plotting
plot<-ggplot(aes(x=Samples, y=Original), data=table) +
  # facet_grid( ~ factor, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Quality read filters"), alpha= 0.8,
            width = 0.6, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Contaminant filters"),
            width= 0.4, stat = "identity") +
  theme_classic( base_size = 10.2) +
  theme(axis.text.x = element_text(size=6, angle = 90, vjust=0.5),
        axis.text.y = element_text(size=5),
        panel.grid.major.y = element_line(size=0.1, color="grey"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size=9.2),
        legend.title = element_text(size=9.8),
        legend.key.height = unit(0.4,"cm")
  ) +
  scale_fill_manual(name='Number of reads:  ',
                    breaks=c('Original number', 'Host read filters', 'Quality read filters', 'Contaminant filters'),
                    values=c('Original number'='green3', 'Host read filters'='blue','Quality read filters'='coral', 'Contaminant filters'='red3')) +
  scale_y_continuous(breaks = c(10000, seq(0, max(table$Original), 25000))) +
  labs(y="Reads abundance", x="FASTQ name")
plot
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 7.2, height = 3.8, dpi=300)

plot +
  facet_grid( ~ factor, space="free_x", scales = "free_x")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters_WITH BATCHES.png", width = 7.2, height = 3.8, dpi=300)

# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)

rm(table, Original_read_number, DADA2_read_number, after_filter_number, plot)



######################## GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(unfiltered_data)==1, unfiltered_data)
length(which(taxa_sums(unfiltered_data)==1)) # if zero there are no singletons
{n<-as.data.frame(otu_table(filter))
  N<-as.data.frame(otu_table(unfiltered_data))
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


####################### % ASSIGNED IN SILVA ############################

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
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database.csv",row.names = F, quote = F)

suppressWarnings(rm(a,b,c,d,e,assigned,data.phy,data.class,data.order,data.fam.data.genus))


########################## GENERAL PCoA ############################

if(! "filter_proof" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

if(length(unique(sample_data(data)[["Time_point"]])) == 1 ){ # It is a control, in case the script has been restarted
  data<-all_data
}


sample_data(data)$Time_point<-factor(sample_data(data)$Time_point, levels=c("PS","OS","SY"))

# New column in which both Time and Genotype
sample_data(data)$Geno_Time<-paste(sample_data(data)$Genotype, sample_data(data)$Time_point, sep="_")

{ data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
  data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
}


####### on ASV
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Geno_Time", shape = "Time_point") +
  scale_color_manual(values=c("C57Ola_PS"="lightgreen","C57Ola_OS"="green3","C57Ola_SY"="chartreuse4",
                              "129Sv_PS"="coral","129Sv_OS"="red","129Sv_SY"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_line(aes(group=Sample_code), color= "grey", size= 0.15) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same sample")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV.png", width = 10, height = 7, dpi=300)
# # without names
# plot_ordination(data.sqrt_prop, ordBC, color = "Geno_Time", shape = "Time_point") +
#   scale_color_manual(values=c("C57Ola_PS"="lightgreen","C57Ola_OS"="green3","C57Ola_SY"="chartreuse4",
#                               "129Sv_PS"="coral","129Sv_OS"="red","129Sv_SY"="red4")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
#   geom_line(aes(group=Sample_code), color= "grey", size= 0.15) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)", 
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        caption="the lines are connecting the same sample")
# ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 10, height = 7, dpi=300)

rm(data.sqrt_prop, eigval, DistBC, ordBC,data.prop.labels)


####### again but on Genera
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Geno_Time", shape = "Time_point") +
  scale_color_manual(values=c("C57Ola_PS"="lightgreen","C57Ola_OS"="green3","C57Ola_SY"="chartreuse4",
                              "129Sv_PS"="coral","129Sv_OS"="red","129Sv_SY"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_line(aes(group=Sample_code), color= "grey", size= 0.15) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same sample")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 10, height = 7, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC,
                color = "Geno_Time",
                shape = "Time_point") +
  scale_color_manual(values=c("C57Ola_PS"="lightgreen","C57Ola_OS"="green3","C57Ola_SY"="chartreuse4",
                              "129Sv_PS"="coral","129Sv_OS"="red","129Sv_SY"="red4")) +
  geom_point(size=3, alpha= 0.2) +
  geom_point(size=2.1, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-10),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.14) +
  geom_line(aes(group=Sample_code), color= "grey", size= 0.08) +
  labs(
    #title="PCoA with Hellinger distance \n(euclidean on Hellinger transformed genera)", 
       color="Genotype", shape="Time point",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same mice at different times")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 7, height = 5, dpi=300)



# #### Genera in 3D
# 
# #setting colors
# {tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Geno_Time),as.character(sample_data(data)$Geno_Time)))
# colnames(tabella_colore)<-c("Gruppo","Colore")
# colors <- gsub("C57Ola_PS","#99FF99",tabella_colore$Colore) # light green
# colors <- gsub("C57Ola_OS","#00EE00FF",colors) # darker
# colors <- gsub("C57Ola_SY","#006600FF",colors) # darker
# colors <- gsub("129Sv_PS", "#FF8080", colors) # light red
# colors <- gsub("129Sv_OS","#CC0000FF",colors) # darker
# colors <- gsub("129Sv_SY","#800000FF",colors) # darker
# }
# {Dist<- vegdist(t(otu_table(data.sqrt_prop)), method = "euclidean") # then Hellinger                                                                                                
#   obj<-ecodist::pco(Dist)
#   matrix<-obj[["vectors"]]
#   row.names(matrix)<-sample_names(data.sqrt_prop)
# }
# rgl::open3d(windowRect=c(25,25,1200,1200))
# pca3d::pca3d(matrix, col=colors, axes.color = "darkgray",
#              radius=1.5, show.shadows = T, show.plane = F, show.centroids = F)
# rgl::rgl.viewpoint(theta = -50.5, phi = 45, fov = 120, zoom = 0.34)
# rgl::legend3d("topleft", c("C57Ola_PS", "C57Ola_OS","C57Ola_SY",
#                            "129Sv_PS","129Sv_OS","129Sv_SY"), 
#               col=c("#99FF99","#00EE00FF","#006600FF",
#                     "#FF8080","#CC0000FF","#800000FF"), pch = 19, magnify=0.7)
# rgl::rgl.snapshot("temp.png")
# temp<-png::readPNG("temp.png")
# # install.packages("pdftools") # it needs also "sudo apt install libpoppler-cpp-dev"
# png::writePNG(temp, "Results/General_Beta_div/PCoA_Hellinger_3D_on_Genera.png", dpi = 300)
# unlink("temp.png")
# rgl::close3d()

rm(data.sqrt_prop, eigval, DistBC, ordBC,data.prop.labels)



################# \\\\\\\\\\ STARTING THE ANALYSIS OF THE PS SUBSET \\\\\\\\\\ ###################

to_reset<-ls() [! ls() %in% c("data","Metadata","unfiltered_data","contam_data","filter_proof","all_data","CircoTax2")]
to_reset<-to_reset [! to_reset %in% ls(pattern="Genera.DESEQ2") ] # to mantain DA results
rm(list = to_reset)


if(! "filter_proof" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
}

if(length(unique(sample_data(data)[["Time_point"]])) > 1 ){   # it works only with the whole dataset
  all_data<-data
  rm(data)
}

data<-subset_samples(all_data, Time_point=="PS")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa
length(sample_names(data))
head(sample_data(data))
tail(sample_data(data))

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

Meta_PS<-Metadata[Metadata$Time_point=="PS",]



########################### COUNTS EXPORT PS ##########################################

dir.create("Results/PS/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/PS/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/PS/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/PS/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/PS/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/PS/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/PS/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/PS/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/PS/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT PS ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one


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
plot<-ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity", position="stack") + 
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(title = element_text(size= 11), 
        axis.text.x=element_text(size= 8.5, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size= 9),
        axis.title.y=element_text(size= 11),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.85, "cm"),
        legend.text = element_text ( size = 12 ),
        legend.margin = margin(3,0,5,0)
  ) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Mice", y="Percentual abundance", title = "Five most abundant phyla at PS time", caption = " 'Others' includes every phylum below rank 5 ")
# 
# plot +
#   facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x")
# ggsave(file="Results/PS/Abundances/TOP_5_phyla_at PS _versione 1.png",width=7,height=5, dpi=300)
# dev.off()

plot +
  facet_nested( ~ Genotype + Strain, space = "free_x", scales = "free_x")
ggsave(file="Results/PS/Abundances/TOP_5_phyla_at PS _versione 2.png",width=7, height=5, dpi=300)
dev.off( )


# means of TOP5 phyla
write.xlsx(file = "Results/PS/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


# # TOP 5 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_5) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
#   labs(x="Mice", y="Relative abundance", title = "Five most abundant genera at PS", caption = " 'Others' includes every genus below rank 5 ")
# ggsave(file="Results/PS/Abundances/TOP_5_genera_PS.png",width=9,height=5,dpi=300)
# dev.off()
# 
# rm(top, tabella, tabella_top)
# 
# 
# # TOP 8 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
#   facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_8) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
#   labs(x="Mices", y="Percentual abundance", title = "Eigth most abundant genera",
#        caption = " 'Others' includes every genus below rank 8 ")
# ggsave(file="Results/PS/Abundances/TOP_8_genera_PS.png",width=9,height=5,dpi=300)
# dev.off()
# 
# suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


### TOP 10 Genera
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
# tabella$XXXX<-gsub("AAAA"," ",tabella$XXXX)    # if it is needed to rename
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
  facet_nested( ~ Genotype + Strain ,scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_10) +
  theme(title = element_text(size= 11), 
        axis.text.x=element_text(size= 8.5, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size= 9),
        axis.title.y=element_text(size= 11),
        legend.key.height = unit(0.18, "cm"),
        legend.key.width = unit(0.75, "cm"),
        legend.text = element_text ( size = 10.5 ),
        legend.position="bottom",
        legend.margin = margin(-8,0,0,-30)
  ) + 
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Mice", y="Percentual abundance",
       title = "Ten most abundant genera at PS time",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/PS/Abundances/TOP_10_genera_at PS_versione 2.png",width=7,height=5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/PS/Abundances/TOP_10_genera_Average_abundances_PS.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



##################### SETTING THE GROUP COLORS ######################

# same function down there will search the colors from here
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Genotype),as.character(sample_data(data)$Genotype)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("C57Ola","chartreuse",tabella_colore$Colore) #green
colors <- gsub("129Sv","coral",colors) #orange

# 

#################### RATIO FIRMICUTES/BACTEROIDES PS (Genotype) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_PS)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_PS[,c("FASTQ_ID","Batch","Strain","Genotype")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Genotype, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Genotype=="129Sv","Ratios_Mean"]<-Ratios_Mean["129Sv"]
ratio_fb[ratio_fb$Genotype=="C57Ola","Ratios_Mean"]<-Ratios_Mean["C57Ola"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Genotype, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Genotype=="129Sv","Ratios_st_err"]<-Ratios_st["129Sv"]/sqrt(length(which(ratio_fb$Genotype=="129Sv")))
ratio_fb[ratio_fb$Genotype=="C57Ola","Ratios_st_err"]<-Ratios_st["C57Ola"]/sqrt(length(which(ratio_fb$Genotype=="C57Ola")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/PS/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Genotype_PS.csv", ratio_fb)

ratio_fb$Genotype<-factor(ratio_fb$Genotype, levels = c("129Sv","C57Ola"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-wilcox.test(Ratio ~ Genotype, paired=F, data=ratio_fb)
p_val<-round(res$p.value, digits = 2)
p_val
# exporting the results
con <- file("Results/PS/Ratio_Firmi_Bacteroi/Genotype_Mann_Whi_Wilcoxon_PS.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Genotype at PS :   V=",res$statistic, "p-value=", p_val, "\n",fill=T)
cat("\n\n Shapiro Wilk test p-value: ",  distr$p.value)
sink()
close(con)

# absolutely not significant
model<-summary(lm(rank(Ratio)~Batch+Strain+Genotype, data=ratio_fb))
model$coefficients[4,4] # Genotype p-value

######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ_ID, fill=Genotype, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("C57Ola"="chartreuse",
                               "129Sv"="coral")) +
  facet_grid2(.~Genotype, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Genotype"))+
  scale_y_sqrt(breaks=c(1,2,3,4,5) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8),
        title = element_text(size=10)) +
  labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Mann Whitney p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/PS/Ratio_Firmi_Bacteroi/Ratio_barplot_Genotype_PS.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Genotype, color=Genotype,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  facet_grid(.~Genotype, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean + Ratios_st_err, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0, # the if else is needed to avoid geom bar below zero
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between genotypes at PS",
       caption = paste("Mann Whitney p-value:",p_val,
                       "\nLinear model (on ranks) p-value:", round(model$coefficients[4,4],2)))
ggsave(file="Results/PS/Ratio_Firmi_Bacteroi/Jitter_with_mean_and_STerr_Genotype.png",width=4,height=3.5, dpi=300) 
dev.off()


####### BOX PLOT
ggplot(data=ratio_fb, aes(x=Genotype, color=Genotype,y=Ratio)) +
  theme_classic(base_size =10) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  facet_grid(.~Genotype, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_boxplot(width = 0.8, size=0.3,
               show.legend = F, aes(color=NULL)) +
  geom_point(size=1.3, aes(color=Genotype)) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between genotypes at PS",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/PS/Ratio_Firmi_Bacteroi/Boxplot_Firm_Bacteroid_Genotype.png",width=3.5,height=3, dpi=300) 
dev.off()



suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY PS (Genotype) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype"
       # , title="Alpha diversity between 129Sv and C57Ola at PS"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9.5)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/PS/Alfa_div_Genotype_ both strain_ at PS.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


### AGAIN, ONLY ON HEALTHY SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Condition=="Healthy"), measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Genotype<-paste(pAlpha$data$Genotype, pAlpha$data$Strain, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv_Ntg"="coral","C57Ola_Ntg"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype" 
       # , title="Alpha diversity at PS (Only Ntg mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/PS/Alfa_div_ Genotype_ ONLY_HEALTHY_PS.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


### AGAIN, ONLY ON SLA SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Condition=="SLA"), measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Genotype<-paste(pAlpha$data$Genotype, pAlpha$data$Strain, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv_G93A"="coral","C57Ola_G93A"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype" 
       # , title="Alpha diversity at PS (Only G93A mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/PS/Alfa_div_Genotype_ ONLY_SLA _PS.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)



##################### BETA DIVERSITY PS (Genotype, with both Ntg and G93A) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[3] # 3 is for the third variable in the design
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[3] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[3]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[3] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[3,],perm_g$aov.tab[3,],perm_f$aov.tab[3,],perm_o$aov.tab[3,],perm_c$aov.tab[3,],perm_p$aov.tab[3,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
# write.csv2(beta, file="Results/PS/Beta_div_Genotype/Beta_diversity_permanova_Hellinger_PS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/PS/Beta_div_Genotype/Beta_diversity_permanova_BOTH STRAINS_at PS.xlsx",col.names = T,row.names = T)


# ### PLUS: checking it with Bray too
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="bray")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/PS/Beta_div_Genotype/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)
# 
# 
# ### PLUS: checking without variable corrections
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Genotype, data=metadata, permutations = 9999, method="euclidean")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/PS/Beta_div_Genotype/Beta_divers_permanova_CHECK_WITHOUT_MODEL_CORRECTIONS.csv",quote=F,row.names = T)
# 

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/PS/Beta_div_Genotype/Beta_disp_permanova_BOTH STRAINS_at PS.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS PS (Genotype, with both Ntg and G93A) #####################

# on Genera
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Strain") +
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset", 
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/PS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_PS.png", width = 9, height = 6, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset",
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/PS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_Genera_no_ellipse.png", width = 9, height = 6, dpi=300)
# # without names
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset", 
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/PS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_points.png", width = 9, height = 6, dpi=300)
# without names and including Ntg/G93A
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape = "Strain") +
  scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/PS/Beta_div_Genotype/Genotypes_PCoA_Hellinger_on_Genera_points.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



##################### BETA DIVERSITY and PCoA at PS (Genotype, Only Ntg) #######################

{
  data_sub<-subset_samples(data, Condition=="Healthy")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/PS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_PS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/PS/Beta_div_Genotype/Beta_div_ C57Ola vs 129sv _ONLY_Ntg_at_PS.xlsx",col.names = T,row.names = T)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/PS/Beta_div_Genotype/Beta_disp_ C57Ola vs 129sv _ONLY_Ntg.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)


########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Genotype<-paste(sample_data(data.prop.labels)$Genotype, sample_data(data.prop.labels)$Strain, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola_Ntg"="chartreuse","129Sv_Ntg"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on Ntg mice at PS", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/PS/Beta_div_Genotype/Genotypes_PCoA_ONLY_Ntg.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



##################### BETA DIVERSITY and PCoA at PS (Genotype, Only G93A) #######################

{
  data_sub<-subset_samples(data, Condition=="SLA")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2] 
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/PS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_G93A_at_PS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/PS/Beta_div_Genotype/Beta_div_ C57Ola vs 129sv _ONLY_G93A_at_PS.xlsx",col.names = T,row.names = T)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/PS/Beta_div_Genotype/Beta_disp_ C57Ola vs 129sv _ONLY_G93A.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)


########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Genotype<-paste(sample_data(data.prop.labels)$Genotype, sample_data(data.prop.labels)$Strain, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on G93A mice at PS", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/PS/Beta_div_Genotype/Genotypes_PCoA_ONLY_G93A.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



################### DA WITH DESEQ2 PS (Genotype) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Strain + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/PS/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_PS.csv"), row.names = F, quote=F, na = "")
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
#write.csv(Res_tot, file="Results/PS/DA_DESeq2_Genotype/Every_result_Genotype_DESeq2_PS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/PS/DA_DESeq2_Genotype/Every_result_Genotype_ both strains_at PS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "Differently abundant Taxa between genotypes at PS",
    y="Proportional Abundance", 
       fill="Genotype", x="")
#ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_every_result_PS.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Actinobacteriota","Coriobacteriia","Coriobacteriales",
          "Lachnospirales","Clostridiales", "Gastranaerophilales",
          "Erysipelotrichales","Marinifilaceae","Rikenellaceae")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "Differently abundant Taxa between genotypes at PS",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
#ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_PS_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola"="chartreuse3","129Sv"="coral2"))
ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_ both strains _at PS.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/PS/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
Genera.DESEQ2_PS<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/PS/DA_DESeq2_Genotype/CircoPlot_ Genotype_ both strains.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins


################### DA WITH DESEQ2 PS (Genotype, ONLY Ntg mice) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition=="Healthy")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/PS/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_PS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/PS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_Ntg_PS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/PS/DA_DESeq2_Genotype/Every_result_Genotype_ ONLY Ntg _at PS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + 
  theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "Differently abundant Taxa between genotypes at PS",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_every_result_PS.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Actinobacteriota","Coriobacteriia","Coriobacteriales",
          "Lachnospirales", "Clostridiales", "Gastranaerophilales",
          "Atopobiaceae", "Erysipelotrichales")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola_Ntg"="chartreuse","129Sv_Ntg"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "DA Taxa between genotypes (only Ntg mice) at PS",
    y="Proportional Abundance %",
    color="Genotype",
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_PS_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype_mutation), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola_Ntg"="chartreuse3","129Sv_Ntg"="coral2"))
ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_ONLY_Ntg_without_redundants.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


# system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/PS/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
# Genera.DESEQ2_PS<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("Coriobacteriaceae_UCG-002","UCG_002", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/PS/DA_DESeq2_Genotype/CircoPlot_Genotpe_ONLY with Ntg_ PS.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins



################### DA WITH DESEQ2 PS (Genotype, ONLY G93A mice) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition=="SLA")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/PS/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_PS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/PS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_G93A_PS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/PS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_G93A_PS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "Differently abundant Taxa between genotypes at PS",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_every_result_PS.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Actinobacteriota","Cyanobacteria",
          "Coriobacteriia","Coriobacteriales", "Vampirivibrionia",
          "Clostridiales", "Marinifilaceae",
          #"Gastranaerophilales",
          "Atopobiaceae")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

unique(Table_tot2$Genotype_mutation)
plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "DA Taxa between genotypes (only G93A mice) at PS",
    color="Genotype",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_PS_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype_mutation), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola_G93A"="chartreuse3","129Sv_G93A"="coral2"))
ggsave(filename = "Results/PS/DA_DESeq2_Genotype/DA_Genotype_ONLY_G93A_without_redundants.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


# system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/PS/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
# Genera.DESEQ2_PS<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("Coriobacteriaceae_UCG-002","UCG_002", Res_tot2$Genus, fixed=T) # more space for the plot


# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/PS/DA_DESeq2_Genotype/CircoPlot_ Genotype _ONLY with G93A_ PS.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins



#################### RATIO FIRMICUTES/BACTEROIDES PS (Ntg vs G93A) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_PS)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_PS[,c("FASTQ_ID","Genotype","Strain","Batch")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Strain, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Strain=="G93A","Ratios_Mean"]<-Ratios_Mean["G93A"]
ratio_fb[ratio_fb$Strain=="Ntg","Ratios_Mean"]<-Ratios_Mean["Ntg"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Strain, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Strain=="G93A","Ratios_st_err"]<-Ratios_st["G93A"]/sqrt(length(which(ratio_fb$Strain=="G93A")))
ratio_fb[ratio_fb$Strain=="Ntg","Ratios_st_err"]<-Ratios_st["Ntg"]/sqrt(length(which(ratio_fb$Strain=="Ntg")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/PS/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Ntg_vs_G93A_PS.csv", ratio_fb)

ratio_fb$Strain<-factor(ratio_fb$Strain, levels = c("Ntg","G93A"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr$p.value
hist(ratio_fb[,"Ratio"])

res_w<-wilcox.test(Ratio ~ Strain, paired=F, data=ratio_fb)
res_w
system("echo 'It is significant according to the Mann Whithney, trying also with a linear model to correct for the batch confounding effect... still significative! ' > Results/PS/Ratio_Firmi_Bacteroi/Note_about_Ratio_H_vs_G93A.txt ")


res<-summary(lm(rank(Ratio)~Batch+Genotype+Strain, data=ratio_fb))
p_val<-as.numeric(round(res$coefficients[,"Pr(>|t|)"][4], digits = 4))
p_val

# PLUS: without Ranks
summary(lm(Ratio~Batch+Genotype+Strain, data=ratio_fb)) # also significant


# exporting the results
con <- file("Results/PS/Ratio_Firmi_Bacteroi/Ntg_vs_G93A_Mann_Whi_Wilcoxon_PS.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Strain at PS :   V=",res_w$statistic, "p-value=", p_val, "\n\n\n\n",fill=T)
cat("Linear Model ~Batch+Genotype+Strain on ranks ... \n")
print(res) # printing linear model too
cat("\n\nPLUS: The model is still significative also without ranks")
cat("\n Shaphiro Wilk test p-value:",distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ_ID, fill=Strain, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Ntg"="deepskyblue",
                               "G93A"="red")) +
  facet_grid2(.~Strain, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Strain"))+
  scale_y_sqrt(breaks=c(1,2,3,4,5,seq(10,80,10)) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("p-value linear model:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/PS/Ratio_Firmi_Bacteroi/Ratio_barplot_Ntg_G93A_PS.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Strain, color=Strain,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(.~Strain, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean + Ratios_st_err, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0, # the if else is needed to avoid geom bar below zero
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio",
       title = "Firmicutes/Bacteroidetes Ratio between Strains at PS",
       caption = paste("Mann Whitney p-value:", round(res_w$p.value,5),
                       "\n Linear Model (on ranks) p-value:", round(res$coefficients[4,4],5) ))
ggsave(file="Results/PS/Ratio_Firmi_Bacteroi/Jitter_with_mean_and_STerr_Ntg_vs_G93A.png",width=4,height=3.5, dpi=300) 
dev.off()


####### BOX PLOT
ggplot(data=ratio_fb, aes(x=Strain, color=Strain,y=Ratio)) +
  theme_classic(base_size =10) +
  scale_color_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(.~Strain, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_boxplot(width = 0.8, size=0.3,
               show.legend = F, aes(color=NULL)) +
  geom_point(size=1.3, aes(color=Strain)) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between strains at PS",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/PS/Ratio_Firmi_Bacteroi/Boxplot_Firm_Bacteroid_STRAIN.png",width=3.5,height=3, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY PS (Ntg vs G93A) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Strain")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Strain, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  theme_classic(base_size = 10) + 
  geom_point(size=1.3, aes(color=Strain)) +
  labs(x="Strain"
       #, title="Alpha diversity between Ntg and G93A at PS"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Strain), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/PS/Alfa_div_Ntg vs G93A_every genotype_PS.png", width = 5,height =4.2, dpi=300)

# # just to test the plotted p value
# alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
# Obser_value<-filter(alphadt, variable=="Observed richness")
# factor<-Obser_value$Strain
# wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)

### AGAIN, ONLY ON C57Ola SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Genotype=="C57Ola"), measures=c("Shannon", "Observed"), x="Strain")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Strain<-paste(pAlpha$data$Strain, pAlpha$data$Genotype, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Strain, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  theme_classic(base_size = 10) + 
  geom_point(size=1.3, aes(color=Strain)) +
  labs(x="Strain"
       #, title="Alpha diversity at PS (Only C57Ola mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Strain), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/PS/Alfa_div_Ntg vs G93A_ONLY C57Ola_PS.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)



### AGAIN, ONLY ON 129Sv SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Genotype=="129Sv"), measures=c("Shannon", "Observed"), x="Strain")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Strain<-paste(pAlpha$data$Strain, pAlpha$data$Genotype, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Strain, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  theme_classic(base_size = 10) + 
  geom_point(size=1.3, aes(color=Strain)) +
  labs(x="Strain"
       # , title="Alpha diversity at PS (Only 129Sv mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Strain), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/PS/Alfa_div_Ntg vs G93A_ONLY 129Sv_PS.png", width = 5,height =4.2, dpi=300)

suppressWarnings(rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor))



##################### BETA DIVERSITY PS (Ntg vs G93A, with both genotypes) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[3] # 3 is for the third variable in the design
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[3]

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[3]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[3] 


sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain,data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[3,],perm_g$aov.tab[3,],perm_f$aov.tab[3,],perm_o$aov.tab[3,],perm_c$aov.tab[3,],perm_p$aov.tab[3,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/PS/Beta_div_Ntg_vs_G93A/Beta_diversity_permanova_Hellinger_PS.csv",quote=F,row.names = T)

# 
# ### PLUS: checking it with Bray too
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="bray")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/PS/Beta_div_Ntg_vs_G93A/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)
# 
# 
# ### PLUS: checking without variable corrections
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Strain, data=metadata, permutations = 9999, method="euclidean")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/PS/Beta_div_Genotype/Beta_divers_permanova_CHECK_WITHOUT_MODEL_CORRECTIONS.csv",quote=F,row.names = T)
# 

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on GENERA (due to the p-value)
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Strain)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/PS/Beta_div_Ntg_vs_G93A/Beta_disp_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# system('echo "THIS p-value is particular, because it is significative ONLY if the confounding factors are included in the model design,\nhowever, the reason becomes clear after checking the extra PCoA plots" > Results/PS/Beta_div_Ntg_vs_G93A/Warning.txt')



########################### PCoA BRAY CURTIS PS (Ntg vs G93A, with both genotypes) #####################

# on Genera
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype_mutation", shape="Strain") + # to color according to both Genotype and Strain
#   scale_color_manual(values=c("C57Ola_G93A"="coral","C57Ola_Ntg"="chartreuse2",
#                               "129Sv_G93A"="red3","129Sv_Ntg"="chartreuse4")) + # 129Sv is the bad genetic background
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset",
#        subtitle = "*hint to read this plot quickly: check the green ones (Ntg) vs the red ones (G93A)",
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        caption = "check the PCoA with only the 129Sv mice to see a difference ALSO in this subset ")
# ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_PS.png", width = 9, height = 6, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset",
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_Genera_no_ellipse.png", width = 9, height = 6, dpi=300)
# # without names
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_points.png", width = 9, height = 6, dpi=300)
# without names but with Genotypes too
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype_mutation", shape="Strain") + # to color according to both Genotype and Strain
  scale_color_manual(values=c("C57Ola_G93A"="coral","C57Ola_Ntg"="deepskyblue",
                              "129Sv_G93A"="red3","129Sv_Ntg"="deepskyblue3")) + # 129Sv is the bad genetic background
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset", 
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_points2.png", width = 7, height = 5, dpi=300)



# # again but on ASV
# {data.prop.labels<-data.prop
#   data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues
#   eigval<- round((eigval/sum(eigval))*100, 1)
# }
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain", shape="Genotype") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 12) + 
#   stat_ellipse(size=0.2) + 
#   #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on PS subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
# ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_ASV.png", width = 9, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



# checking only with 129Sv (genetic background) subset (which is uncertain in the PCoA with the whole dataset)
data.temp<-subset_samples(data.genus, Genotype=="129Sv")
data.temp<-transform_sample_counts(data.temp, function(x) x/sum(x))
{data.sqrt_prop<-transform_sample_counts(data.temp, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
  scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on PS subset", 
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption = "in the PCoA with the *whole* dataset is not that clear, but there is a difference between Ntg and G93A\n ALSO in the 129Sv genetic background subset, just as reported by the PERMANOVA")
ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_ONLY_WITH_129Sv.png", width = 9, height = 6, dpi=300)


##################### BETA DIVERSITY and PCoA at PS (Ntg vs G93A, Only C57Ola) #######################

{
  data_sub<-subset_samples(data, Genotype=="C57Ola")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Strain, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/PS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_PS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/PS/Beta_div_Genotype/Beta_div_ Ntg vs G93A _ ONLY C57Ola_ PS.xlsx",col.names = T,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Strain)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/PS/Beta_div_Genotype/Beta_disp_ Ntg vs G93A _ONLY C57Ola_PS.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)



########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Strain<-paste(sample_data(data.prop.labels)$Strain, sample_data(data.prop.labels)$Genotype, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Strain") + 
  scale_color_manual(values = c("G93A_C57Ola"="coral","Ntg_C57Ola"="deepskyblue")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on C57Ola mice at PS subset", 
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg vs G93A PCoA _but ONLY C57Ola.png", width = 7, height = 5, dpi=300)



##################### BETA DIVERSITY and PCoA at PS (Ntg vs G93A, Only 129Sv) #######################

{
  data_sub<-subset_samples(data, Genotype=="129Sv")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Strain, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/PS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_PS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/PS/Beta_div_Genotype/Beta_div_ Ntg vs G93A _ ONLY 129Sv_ PS.xlsx",col.names = T,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Strain)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/PS/Beta_div_Genotype/Beta_disp_ Ntg vs G93A _ONLY 129Sv_PS.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)



########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Strain<-paste(sample_data(data.prop.labels)$Strain, sample_data(data.prop.labels)$Genotype, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Strain") + 
  scale_color_manual(values = c("G93A_129Sv"="coral","Ntg_129Sv"="deepskyblue")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on 129Sv mice at PS subset", 
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/PS/Beta_div_Ntg_vs_G93A/Ntg vs G93A PCoA _but ONLY 129Sv.png", width = 7, height = 5, dpi=300)



#################### VEEN DIAGRAM (Ntg vs G93A at PS) #########################

### target: very low yet enough abundance (0.1%) in at least 1/4 of total samples (~8 mice)

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

minimum <- round( (length(sample_names(data.venn))/100) * 25 , 0) # at least in more than 25% of PS samples 
minimum # 8
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 (then it is not a noise!) --> "a point"
who<-who[!rowSums(who)> minimum,]
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
who<-who[!is.na(who)] # this is to avoid the NA, because there are also "good" NA
data.venn<-subset_taxa(data.venn, ! Genus %in% who)

G93A<-subset_samples(data.venn, Strain=="G93A")
G93A<-as.character(tax_table(prune_taxa(taxa_sums(G93A)>0, G93A))[,"Genus"])

Ntg<-subset_samples(data.venn, Strain=="Ntg")
Ntg<-as.character(tax_table(prune_taxa(taxa_sums(Ntg)>0, Ntg))[,"Genus"])

ONLY_IN_Ntg<- Ntg[! Ntg %in% G93A]
Ntg_list<-ONLY_IN_Ntg
ONLY_IN_Ntg<- paste(ONLY_IN_Ntg, collapse = ", ")
head(ONLY_IN_Ntg)

ONLY_IN_G93A<- G93A[! G93A %in% Ntg]
ONLY_IN_G93A<- paste(ONLY_IN_G93A, collapse = ", ")
head(ONLY_IN_G93A)

# con<-file("Results/Beta_div/Ntg_vs_G93A_exclusive_bacteria.txt")
# sink(con, append=TRUE)
# cat("\n\nONLY IN G93A", fill=TRUE)
# cat(ONLY_IN_G93A)
# cat("\n\nONLY IN Ntg", fill=TRUE)
# cat(ONLY_IN_Ntg)
# 
# sink()
# close(con)




################### DA WITH DESEQ2 PS (Ntg vs G93A) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype + Strain)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Strain", "Ntg", "G93A"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/PS/DA_DESeq2/DA_",t,"_ratio_Ntg_vs_G93A_PS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/PS/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_PS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/PS/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_PS.xlsx", showNA = F, col.names = T)


plot_table<-ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Strain)) + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance)+2,2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "Differently abundant Taxa between strains at PS",
    y="Proportional Abundance %", 
       fill="Strain", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_PS_v1.png", width = 7.3, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Strain), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red"))
ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_both Genotypes_at PS.png", width = 7.3, height = 7, dpi=300)
dev.off()

# Table_tot2<-subset(Table_tot, Taxa %in% c("Family","Genus")) # to remove redundant results
# ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Strain)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
#   scale_fill_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
#         axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
#         axis.text.y = element_text(size=11.5), 
#         plot.title= element_text(size=18),
#         panel.grid.minor.y= element_blank(),
#         plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
#   #scale_x_discrete(expand=c(-0.2, 1)) +
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
#   #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
#   labs(title= "Differently abundant Taxa between Strains at PS", y="Proportional Abundance", 
#        fill="Strain", x="")
# ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_PS.png", width = 8, height = 9, dpi=300)
# dev.off()


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
Res_tot2$Genus<-gsub("-","_",Res_tot2$Genus)
Res_tot2$Genus<-gsub("001","001             ",Res_tot2$Genus)

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/PS/DA_DESeq2_Ntg_vs_G93A/CircoPlot_Strain_Every genotype_PS.png",width=1650,height=1650, res = 300)
circo
dev.off()


system(" echo 'This model has been corrected by both Batch and Genotype variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/PS/DA_DESeq2_Ntg_vs_G93A/Statistical_Corrections.txt ")


################### DA WITH DESEQ2 PS (Ntg vs G93A, Only C57Ola ) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Genotype=="C57Ola")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~ Batch + Strain)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Strain", "Ntg", "G93A"))
  resultsNames(DE)
  res<-results(DE)
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/PS/DA_DESeq2/DA_",t,"_ratio_Ntg_vs_G93A_PS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/PS/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_PS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/PS/DA_DESeq2_Ntg_vs_G93A/Every result_Ntg vs G93A_ ONLY C57Ola _ PS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Strain)) + 
  theme_classic(base_size = 16) +
  geom_boxplot() +
  scale_fill_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance)+2,2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "Differently abundant Taxa between strains at PS",
    y="Proportional Abundance %", 
       fill="Strain", x="")

Table_tot2<-subset(Table_tot,  Taxa!="Order") # to remove redundant results
Table_tot2$Strain<-paste(Table_tot2$Strain, Table_tot2$Genotype, sep="_")
plot_table<- ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Strain)) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  # geom_boxplot(width=0.8) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("Ntg_C57Ola"="deepskyblue","G93A_C57Ola"="red")) +
  theme(strip.text.x=element_text(size=14,colour="black")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.8), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) 
        ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,12,2) , seq(15, max(Table_tot$Abundance+4),3))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "DA Taxa between Strains (only C57Ola mice) at PS",
    y="Proportional Abundance %",
       fill="Strain", x="")
# ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_PS.png", width = 8, height = 9, dpi=300)
# dev.off()


plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Strain), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("Ntg_C57Ola"="deepskyblue","G93A_C57Ola"="red"))
ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg vs G93A_only C57Ola_no redundants_PS.png", width = 7.3, height = 7, dpi=300)
dev.off()



##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
Res_tot2<-Res_tot2[Res_tot2$Taxon!="Family", ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
Res_tot2$Genus<-gsub("-","_",Res_tot2$Genus)
# Res_tot2$Genus<-gsub("001","001             ",Res_tot2$Genus)

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/PS/DA_DESeq2_Ntg_vs_G93A/CircoPlot_Ntg vs G93A _ONLY C57Ola_PS.png",width=1650,height=1650, res = 300)
circo
dev.off()


################### DA WITH DESEQ2 PS (Ntg vs G93A, Only 129Sv ) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Genotype=="129Sv")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~ Batch + Strain)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Strain", "Ntg", "G93A"))
  resultsNames(DE)
  res<-results(DE)
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/PS/DA_DESeq2/DA_",t,"_ratio_Ntg_vs_G93A_PS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/PS/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_PS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/PS/DA_DESeq2_Ntg_vs_G93A/Every result_Ntg vs G93A_ ONLY 129Sv _ PS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Strain)) + 
  theme_classic(base_size = 16) +
  geom_boxplot() +
  scale_fill_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance)+2,2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between strains at PS",
    y="Proportional Abundance %", 
       fill="Strain", x="")


Table_tot2<-subset(Table_tot,  Taxa!="Family") # to remove redundant results
Table_tot2$Strain<-paste(Table_tot2$Strain, Table_tot2$Genotype, sep="_")
plot_table<- ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Strain)) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  # geom_boxplot(width=0.8) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("Ntg_129Sv"="deepskyblue","G93A_129Sv"="red")) +
  theme(strip.text.x=element_text(size=14,colour="black")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5)
        ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    #title= "DA Taxa between Strains (only 129Sv mice) at PS", 
    y="Proportional Abundance %",
       fill="Strain", x="")
# ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_PS.png", width = 8, height = 9, dpi=300)
dev.off()


plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Strain), size= 1, alpha= 0.65) +
  scale_color_manual(values=c("Ntg_129Sv"="deepskyblue","G93A_129Sv"="red"))
ggsave(filename = "Results/PS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg vs G93A_only 129Sv_no redundants_PS.png", width = 7.3, height = 7, dpi=300)
dev.off()



# ##### CIRCOPLOT TO PLOT THE LOG2FC
# 
# # preparing the object
# Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# Res_tot2<-Res_tot2[Res_tot2$Taxon!="Family", ]
# Res_tot2$Genus<-gsub("-","_",Res_tot2$Genus)
# # Res_tot2$Genus<-gsub("001","001             ",Res_tot2$Genus)
# 
# # plotting
# circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
# png(file="Results/PS/DA_DESeq2_Ntg_vs_G93A/CircoPlot_Ntg vs G93A _ONLY 129Sv_PS.png",width=1650,height=1650, res = 300)
# circo
# dev.off()



################ \\\\\\\\\\ STARTING THE ANALYSIS OF THE OS SUBSET \\\\\\\\\\ ####################

to_reset<-ls() [! ls() %in% c("data","Metadata","unfiltered_data","contam_data","filter_proof","all_data","CircoTax2")]
to_reset<-to_reset [! to_reset %in% ls(pattern="Genera.DESEQ2") ] # to mantain DA results
rm(list = to_reset)


if(! "filter_proof" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

if(length(unique(sample_data(data)[["Time_point"]])) > 1 ){
  all_data<-data
  rm(data)
}

data<-subset_samples(all_data, Time_point=="OS")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa
length(sample_names(data))
head(sample_data(data))


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

Meta_OS<-Metadata[Metadata$Time_point=="OS",]


########################### COUNTS EXPORT OS ##########################################

dir.create("Results/OS/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/OS/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/OS/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/OS/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/OS/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/OS/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/OS/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/OS/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/OS/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT OS ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one


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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Mice", y="Percentual abundance", title = "Five most abundant phyla at OS time", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/OS/Abundances/TOP_5_phyla_OS.png",width=8,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/OS/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


# # TOP 5 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_5) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
#   labs(x="Mice", y="Relative abundance", title = "Five most abundant genera at OS", caption = " 'Others' includes every genus below rank 5 ")
# ggsave(file="Results/OS/Abundances/TOP_5_genera_OS.png",width=9,height=5,dpi=300)
# dev.off()
# 
# rm(top, tabella, tabella_top)
# 
# 
# # TOP 8 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_8) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
#   labs(x="Mice", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
# ggsave(file="Results/OS/Abundances/TOP_8_genera_OS.png",width=9,height=5,dpi=300)
# dev.off()
# 
# suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


### TOP 10 Genera
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
# tabella$XXXX<-gsub("AAAA"," ",tabella$XXXX)    # if it is needed to rename
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.2 ),
        legend.position="bottom",
        legend.margin = margin(-8,0,0,-40)
  ) + 
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Mice", y="Percentual abundance",
       title = "Ten most abundant genera at OS time",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/OS/Abundances/TOP_10_genera_OS.png",width=8,height=5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/OS/Abundances/TOP_10_genera_Average_abundances_OS.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



##################### SETTING THE GROUP COLORS ######################

# same function down there will search the colors from here
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Genotype),as.character(sample_data(data)$Genotype)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("C57Ola","chartreuse",tabella_colore$Colore) #green
colors <- gsub("129Sv","coral",colors) #orange


#################### HIERARCHICAL CLUSTERING OS ###################

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.prop))))) # then Hellinger
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/OS/Hierarchical_clustering/Hierarchical_cluster_Hellinger_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure of OS subset using Hellinger distance",
     sub="C57Ola = green     129Sv = orange")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/OS/Hierarchical_clustering/Hierarchical_cluster_Bray_sqrt_prop_normalized_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure of OS subset using\n Bray-Curtis distance on sqrt proportional ASVs",
     sub="C57Ola = green     129Sv = orange")
dev.off()

suppressWarnings(rm(c,color_table,labels_colors))



#################### RATIO FIRMICUTES/BACTEROIDES OS (Genotype) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_OS)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_OS[,c("FASTQ_ID","Genotype")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Genotype, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Genotype=="129Sv","Ratios_Mean"]<-Ratios_Mean["129Sv"]
ratio_fb[ratio_fb$Genotype=="C57Ola","Ratios_Mean"]<-Ratios_Mean["C57Ola"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Genotype, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Genotype=="129Sv","Ratios_st_err"]<-Ratios_st["129Sv"]/sqrt(length(which(ratio_fb$Genotype=="129Sv")))
ratio_fb[ratio_fb$Genotype=="C57Ola","Ratios_st_err"]<-Ratios_st["C57Ola"]/sqrt(length(which(ratio_fb$Genotype=="C57Ola")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/OS/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Genotype_OS.csv", ratio_fb)

ratio_fb$Genotype<-factor(ratio_fb$Genotype, levels = c("129Sv","C57Ola"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-wilcox.test(Ratio ~ Genotype, paired=F, data=ratio_fb)
p_val<-round(res$p.value, digits = 2)
p_val

# absolutely not significative

# exporting the results
con <- file("Results/OS/Ratio_Firmi_Bacteroi/Genotype_Mann_Whi_Wilcoxon_OS.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Genotype at OS :   V=",res$statistic, "p-value=", p_val, "\n",fill=T)
cat("\n\n Shapiro Wilk test p-value: ",  distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ_ID, fill=Genotype, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("C57Ola"="chartreuse",
                               "129Sv"="coral")) +
  facet_grid2(.~Genotype, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Genotype"))+
  scale_y_sqrt(breaks=c(1,2,3,4,5,seq(10,80,10)) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Mann Whitney p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/OS/Ratio_Firmi_Bacteroi/Ratio_barplot_Genotype_OS.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Genotype, color=Genotype,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  facet_grid(.~Genotype, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean + Ratios_st_err, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0, # the if else is needed to avoid geom bar below zero
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between genotypes at OS",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/OS/Ratio_Firmi_Bacteroi/Jitter_with_mean_and_STerr_Genotype.png",width=4,height=3.5, dpi=300) 
dev.off()


####### BOX PLOT
ggplot(data=ratio_fb, aes(x=Genotype, color=Genotype,y=Ratio)) +
  theme_classic(base_size =10) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  facet_grid(.~Genotype, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_boxplot(width = 0.8, size=0.3,
               show.legend = F, aes(color=NULL)) +
  geom_point(size=1.3, aes(color=Genotype)) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between genotypes at OS",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/OS/Ratio_Firmi_Bacteroi/Boxplot_Firm_Bacteroid_Genotype.png",width=3.5,height=3, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY OS (Genotype) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype"
       #, title="Alpha diversity between 129Sv and C57Ola at OS"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9.5)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/OS/Alfa_div_ Genotype_ both strains _at OS.png", width = 5,height =4.2, dpi=300)

system("echo 'The two bacthes have an almost equilibrated library size, then the OS alpha diversity p-value is quite reliable' > Results/OS/Note_regarding_alpha_div_results.txt ")
# 
# # just to test the plotted p value
# alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
# Obser_value<-filter(alphadt, variable=="Observed richness")
# factor<-Obser_value$Genotype
# wilcox.test(Obser_value$value~factor)
# 
# rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


### AGAIN, ONLY ON HEALTHY SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Condition=="Healthy"), measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Genotype<-paste(pAlpha$data$Genotype, pAlpha$data$Strain, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv_Ntg"="coral","C57Ola_Ntg"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype" 
       #title="Alpha diversity at OS (Only Ntg mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/OS/Alfa_div_Genotype_ONLY_HEALTHY_OS.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


### AGAIN, ONLY ON SLA SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Condition=="SLA"), measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Genotype<-paste(pAlpha$data$Genotype, pAlpha$data$Strain, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv_G93A"="coral","C57Ola_G93A"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype" 
       #title="Alpha diversity at OS (Only G93A mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/OS/Alfa_div_Genotype_ONLY_SLA_OS.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY  OS (Genotype) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[3] # 3 is for the third variable in the design
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[3] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[3]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[3] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[3,],perm_g$aov.tab[3,],perm_f$aov.tab[3,],perm_o$aov.tab[3,],perm_c$aov.tab[3,],perm_p$aov.tab[3,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/OS/Beta_div_Genotype/Beta_diversity_permanova_BOTH STRAINS_OS.csv",quote=F,row.names = T)

# 
# ### PLUS: checking it with Bray too
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Genotype, data=metadata, permutations = 9999, method="bray")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/OS/Beta_div_Genotype/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)
# 
# 
# ### PLUS: checking without variable corrections
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Genotype, data=metadata, permutations = 9999, method="euclidean")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/SY/Beta_div_Genotype/Beta_divers_permanova_CHECK_WITHOUT_MODEL_CORRECTIONS.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/OS/Beta_div_Genotype/Beta_disp_permanova_BOTH STRAINS_OS.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS OS (Genotype) #####################

# on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Genotype_mutation<-gsub("G93A","G93A",sample_data(data.prop.labels)$Genotype_mutation)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Batch") + # to color according to both Genotype and Strain
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset", 
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/OS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_OS_1.png", width = 9, height = 6, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset",
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/OS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_Genera_no_ellipse.png", width = 9, height = 6, dpi=300)
# # without names
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset", 
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/OS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_points.png", width = 9, height = 6, dpi=300)
# without names and including Ntg/G93A
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape = "Strain") +
  scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/OS/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_points2.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


# checking only with 129Sv (genetic background) subset (which is uncertain in the PCoA with the whole dataset)
data.temp<-subset_samples(data.genus, Genotype=="129Sv")
data.temp<-transform_sample_counts(data.temp, function(x) x/sum(x))
{data.sqrt_prop<-transform_sample_counts(data.temp, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
  scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset",
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption = "in the PCoA with the *whole* dataset is not that clear, but there is a SMALL difference between Ntg and G93A\n ALSO in the 129Sv genetic background subset, just as reported by the PERMANOVA")
ggsave(file="Results/OS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA_Hellinger_ONLY_with_129Sv.png", width = 9.5, height = 6, dpi=300)


system("echo 'NB: This difference of beta div is significative but quite borderline, the p-value is too near 0.05!' > Results/OS/Beta_div_Ntg_vs_G93A/Warning.txt")



##################### BETA DIVERSITY and PCoA at OS (Genotype, Only Ntg) #######################

{
  data_sub<-subset_samples(data, Condition=="Healthy")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
#write.csv2(beta, file="Results/OS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_OS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/OS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_OS.xlsx",col.names = T,row.names = T)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/OS/Beta_div_Genotype/Beta_disp_ONLY_Ntg.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)


########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Genotype<-paste(sample_data(data.prop.labels)$Genotype, sample_data(data.prop.labels)$Strain, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola_Ntg"="chartreuse","129Sv_Ntg"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on Ntg mice at OS", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/OS/Beta_div_Genotype/Genotypes_PCoA_ONLY_Ntg.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



##################### BETA DIVERSITY and PCoA at OS (Genotype, Only G93A) #######################

{
  data_sub<-subset_samples(data, Condition=="SLA")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2] 
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/OS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_G93A_at_OS.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/OS/Beta_div_Genotype/Beta_div_Hellinger_ONLY_G93A_at_OS.xlsx",col.names = T,row.names = T)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/OS/Beta_div_Genotype/Beta_disp_ONLY_G93A.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)


########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Genotype<-paste(sample_data(data.prop.labels)$Genotype, sample_data(data.prop.labels)$Strain, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on G93A mice at OS", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/OS/Beta_div_Genotype/Genotypes_PCoA_ONLY_G93A.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### DA WITH DESEQ2 OS (Genotype) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Strain + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/OS/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_OS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/OS/DA_DESeq2_Genotype/Every_result_Genotype_DESeq2_OS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/OS/DA_DESeq2_Genotype/Every_result_Genotype_ both strains _at OS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at OS",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_every_result_OS.png", width = 15, height = 9, dpi=300)
dev.off()

Redund<-c("Actinobacteriota","Coriobacteriia","Coriobacteriales","Vampirivibrionia",
          "Lachnospirales","Clostridiales", "Gastranaerophilales","Cyanobacteria","Campylobacterales",
          "Erysipelotrichales","Marinifilaceae","Helicobacteraceae","Campylobacteria","Campilobacterota")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  theme(strip.text.x=element_text(size=13.5,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at OS",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
#ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_OS_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola"="chartreuse3","129Sv"="coral2"))
ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_ both strains_ no redund _at OS.png", width = 10.5, height = 7, dpi=300)
dev.off()


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("-","_", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/OS/DA_DESeq2_Genotype/CircoPlot_Genotype_both strains.png",width=1650,height=1650, res = 300)
circo
dev.off()


system(" echo 'This model has been corrected by Bacth and Genotype variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/OS/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
Genera.DESEQ2_OS<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


################### DA WITH DESEQ2 OS (Genotype, ONLY Ntg mice) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition=="Healthy")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/OS/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_OS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/OS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_Ntg_OS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/OS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY Ntg_at OS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at OS",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_every_result_OS.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Campilobacterota","Campylobacteria",
          "Coriobacteriia","Campylobacterales",
          "Lachnospirales","Atopobiaceae",
          "Helicobacteraceae", "Marinifilaceae")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola_Ntg"="chartreuse","129Sv_Ntg"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "DA Taxa between genotypes (only Ntg mice) at OS",
    color="Genotype",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_OS_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype_mutation), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola_Ntg"="chartreuse3","129Sv_Ntg"="coral2"))
ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_ONLY_Ntg_without_redundants.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


# system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/OS/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
# Genera.DESEQ2_OS<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/OS/DA_DESeq2_Genotype/CircoPlot_Genotypes_ONLY_Ntg_OS.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins



################### DA WITH DESEQ2 OS (Genotype, ONLY G93A mice) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition=="SLA")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/OS/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_OS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/OS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_G93A_OS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/OS/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_G93A_OS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at OS",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_every_result_OS.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Erysipelotrichales")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "DA Taxa between genotypes (only G93A mice) at OS",
    y="Proportional Abundance %", 
    color="Genotype",
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_OS_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype_mutation), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola_G93A"="chartreuse3","129Sv_G93A"="coral2"))
ggsave(filename = "Results/OS/DA_DESeq2_Genotype/DA_Genotype_ONLY_G93A_without_redundants.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


# system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/OS/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
# Genera.DESEQ2_OS<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/OS/DA_DESeq2_Genotype/CircoPlot_Genotypes_ONLY_G93A_OS.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins



#################### RATIO FIRMICUTES/BACTEROIDES OS (Ntg vs G93A) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_OS)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_OS[,c("FASTQ_ID","Genotype","Strain","Batch")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Strain, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Strain=="G93A","Ratios_Mean"]<-Ratios_Mean["G93A"]
ratio_fb[ratio_fb$Strain=="Ntg","Ratios_Mean"]<-Ratios_Mean["Ntg"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Strain, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Strain=="G93A","Ratios_st_err"]<-Ratios_st["G93A"]/sqrt(length(which(ratio_fb$Strain=="G93A")))
ratio_fb[ratio_fb$Strain=="Ntg","Ratios_st_err"]<-Ratios_st["Ntg"]/sqrt(length(which(ratio_fb$Strain=="Ntg")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/OS/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Ntg_vs_G93A_OS.csv", ratio_fb)

ratio_fb$Strain<-factor(ratio_fb$Strain, levels = c("Ntg","G93A"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res_w<-wilcox.test(Ratio ~ Strain, paired=F, data=ratio_fb)
res_w$p.value

res<-summary(lm(rank(Ratio)~Batch+Genotype+Strain, data=ratio_fb))
p_val<-as.numeric(round(res$coefficients[,"Pr(>|t|)"][4], digits = 2))
p_val

system("echo 'It is NOT significant according to Wilcoxon, trying also with a linear model to correct for the batch confounding effect... even less significant! ' > Results/OS/Ratio_Firmi_Bacteroi/Note_about_Ratio_H_vs_G93A.txt ")

# exporting the results
con <- file("Results/OS/Ratio_Firmi_Bacteroi/Ntg_vs_G93A_Mann_Whi_Wilcoxon_OS.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Strain at OS :   V=",res_w$statistic, "p-value=", p_val, "\n\n\n\n",fill=T)
cat("Linear Model ~Batch+Genotype+Strain on ranks ... \n")
print(res) # printing linear model too
cat("\n\nPLUS: The model is still NOT significative also without ranks")
cat("\n\n Shapiro Wilk test p-value: ",  distr$p.value)
sink()
close(con)



######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ_ID, fill=Strain, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Ntg"="deepskyblue",
                               "G93A"="red")) +
  facet_grid2(.~Strain, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Strain"))+
  scale_y_sqrt(breaks=c(1,2,3,4,5,seq(10,80,10)) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("p-value linear model:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/OS/Ratio_Firmi_Bacteroi/Ratio_barplot_Ntg_G93A_OS.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Strain, color=Strain,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(.~Strain, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean + Ratios_st_err, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0, # the if else is needed to avoid geom bar below zero
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between Strains at OS",
       caption = paste("p-value linear model:",p_val) )
ggsave(file="Results/OS/Ratio_Firmi_Bacteroi/Jitter_with_mean_and_STerr_Ntg_vs_G93A.png",width=4,height=3.5, dpi=300) 
dev.off()


####### BOX PLOT
ggplot(data=ratio_fb, aes(x=Strain, color=Strain,y=Ratio)) +
  theme_classic(base_size =10) +
  scale_color_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(.~Strain, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_boxplot(width = 0.8, size=0.3,
               show.legend = F, aes(color=NULL)) +
  geom_point(size=1.3, aes(color=Strain)) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between strains at OS",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/OS/Ratio_Firmi_Bacteroi/Boxplot_Firm_Bacteroid_STRAIN.png",width=3.5,height=3, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY OS (Ntg vs G93A) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Strain")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Strain, y=value, color=NULL),
               size= 0.25, alpha=0.1) + 
  theme_classic(base_size = 10) + 
  geom_point(size=1.3, aes(color=Strain)) +
  scale_color_manual(values = c("G93A"="coral","Ntg"="deepskyblue")) +
  labs(x="Strain" 
       # title="Alpha diversity between Ntg and G93A at OS"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9.5)) +
  stat_compare_means(aes(group = Strain), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/OS/Alfa_diversity_Ntg vs G93A_ both Genotypes _OS.png", width = 5,height =4.2, dpi=300)


# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Strain
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY OS (Ntg vs G93A) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[3] # 3 is for the third variable in the design
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[3]

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[3]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[3] 


sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[3,],perm_g$aov.tab[3,],perm_f$aov.tab[3,],perm_o$aov.tab[3,],perm_c$aov.tab[3,],perm_p$aov.tab[3,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/OS/Beta_div_Ntg_vs_G93A/Beta_diversity_permanova_Hellinger_OS.csv",quote=F,row.names = T)

# 
# ### PLUS: checking it with Bray too
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="bray")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/OS/Beta_div_Ntg_vs_G93A/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)
# 

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on GENERA (due to the p-value)
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist, metadata$Strain)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/OS/Beta_div_Ntg_vs_G93A/Beta_disp_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS OS (Ntg vs G93A) #####################

# on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Genotype_mutation<-gsub("G93A","G93A",sample_data(data.prop.labels)$Genotype_mutation)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype_mutation", shape="Strain") + # to color according to both Genotype and Strain
#   scale_color_manual(values=c("C57Ola_G93A"="coral","C57Ola_Ntg"="chartreuse2",
#                               "129Sv_G93A"="red3","129Sv_Ntg"="chartreuse4")) + # 129Sv is the bad genetic background
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        subtitle = "*hint: bright red vs bright green, then dark red vs dark green")
# ggsave(file="Results/OS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_OS.png", width = 9, height = 6, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset",
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/OS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_Genera_no_ellipse.png", width = 9, height = 6, dpi=300)
# # without names
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/OS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_points.png", width = 9, height = 6, dpi=300)
# without names but with Genotypes too
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype_mutation", shape="Strain") + # to color according to both Genotype and Strain
  scale_color_manual(values=c("C57Ola_G93A"="coral","C57Ola_Ntg"="deepskyblue",
                              "129Sv_G93A"="red3","129Sv_Ntg"="deepskyblue3")) + # 129Sv is the bad genetic background
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on OS subset", 
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/OS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_points.png", width = 7, height = 5, dpi=300)


# # again but on ASV
# {data.prop.labels<-data.prop
#   data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues
#   eigval<- round((eigval/sum(eigval))*100, 1)
# }
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on OS subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
# ggsave(file="Results/OS/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_ASV.png", width = 9, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



################### DA WITH DESEQ2 OS (Ntg vs G93A) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

##### STARTING THE DIFFERENTIAL ANALYSIS
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype + Strain)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Strain", "Ntg", "G93A"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/OS/DA_DESeq2/DA_",t,"_ratio_Ntg_vs_G93A_OS.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/OS/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_OS.csv", row.names = F)
write.xlsx(Res_tot, file="Results/OS/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_OS.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Strain)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  labs(title= "Differently abundant Taxa between Strains at OS", y="Proportional Abundance", 
       fill="Strain", x="")
# ggsave(filename = "Results/OS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_Strain_every_result_OS.png", width = 18, height = 9, dpi=300)
dev.off()

Table_tot2<-subset(Table_tot, Taxa %in% "Order") # to remove the redundant result
plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Strain)) +
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  labs(
    # title= "Differently abundant Taxa between strains at OS",
    y="Proportional Abundance %", 
       fill="Strain", x="")


plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/OS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_OS_v1.png", width = 7.3, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Strain), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red"))
ggsave(filename = "Results/OS/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_OS_v2.png", width = 7.3, height = 7, dpi=300)
dev.off()


system(" echo 'This model has been corrected by both Batch and Genotype variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/OS/DA_DESeq2_Ntg_vs_G93A/Statistical_Corrections.txt ")


################# \\\\\\\\\\ STARTING THE ANALYSIS OF THE SY SUBSET \\\\\\\\\\ #######################

to_reset<-ls() [! ls() %in% c("data","Metadata","unfiltered_data","contam_data","filter_proof","all_data","CircoTax2")]
to_reset<-to_reset [! to_reset %in% ls(pattern="Genera.DESEQ2") ] # to mantain DA results
rm(list = to_reset)


if(! "filter_proof" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

if(length(unique(sample_data(data)[["Time_point"]])) > 1 ){
  all_data<-data
  rm(data)
}

data<-subset_samples(all_data, Time_point=="SY")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa
length(sample_names(data))
head(sample_data(data))


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

Meta_SY<-Metadata[Metadata$Time_point=="SY",]


########################### COUNTS EXPORT SY ##########################################

dir.create("Results/SY/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/SY/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/SY/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/SY/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/SY/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/SY/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/SY/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/SY/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/SY/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT SY ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one


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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Mice", y="Percentual abundance", title = "Five most abundant phyla at SY time", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/SY/Abundances/TOP_5_phyla_SY.png",width=8,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/SY/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


# # TOP 5 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_5) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
#   labs(x="Mice", y="Relative abundance", title = "Five most abundant genera at SY", caption = " 'Others' includes every genus below rank 5 ")
# ggsave(file="Results/SY/Abundances/TOP_5_genera_SY.png",width=9,height=5,dpi=300)
# dev.off()
# 
# rm(top, tabella, tabella_top)
# 
# # TOP 8 Genera
# {top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:8]
#   prune.dat_top <- prune_taxa(top,data.genus.prop)
#   tax_selected<-as.data.frame(tax_table(prune.dat_top))
#   tax_selected<-Taxa.genus.update[row.names(tax_selected),]
#   tax_table(prune.dat_top)<-as.matrix(tax_selected)
#   others<-taxa_names(data.genus.prop)
#   others<-others[!(others %in% top)]
#   prune.data.others<-prune_taxa(others,data.genus.prop)
#   tabella_top<-psmelt(prune.dat_top)
#   tabella_others<-psmelt(prune.data.others)
#   tabella_others$Genus<-"Others"
#   tabella<-rbind.data.frame(tabella_top,tabella_others)
#   tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
#   tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
# }
# ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
#   geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
#   scale_fill_manual(values=fill_color_8) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5), 
#         legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
#   theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
#   labs(x="Mice", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
# ggsave(file="Results/SY/Abundances/TOP_8_genera_SY.png",width=9,height=5,dpi=300)
# dev.off()
# 
# suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


### TOP 10 Genera
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
# tabella$XXXX<-gsub("AAAA"," ",tabella$XXXX)    # if it is needed to rename
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Genotype),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.1 ),
        legend.position="bottom",
        legend.margin = margin(-8,0,0,-40)
  ) + 
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Mice", y="Percentual abundance",
       title = "Ten most abundant genera at SY time",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/SY/Abundances/TOP_10_genera_SY.png",width=8,height=5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/SY/Abundances/TOP_10_genera_Average_abundances_SY.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



##################### SETTING THE GROUP COLORS ######################

# same function down there will search the colors from here !
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Genotype),as.character(sample_data(data)$Genotype)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("C57Ola","chartreuse",tabella_colore$Colore) #green
colors <- gsub("129Sv","coral",colors) #orange


#################### HIERARCHICAL CLUSTERING SY ###################

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.prop))))) # then Hellinger
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/SY/Hierarchical_clustering/Hierarchical_cluster_Hellinger_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure of SY subset using Hellinger distance",
     sub="C57Ola = green     129Sv = orange")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/SY/Hierarchical_clustering/Hierarchical_cluster_Bray_sqrt_prop_normalized_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure of SY subset using\n Bray-Curtis distance on sqrt proportional ASVs",
     sub="C57Ola = green     129Sv = orange")
dev.off()

suppressWarnings(rm(c,color_table,labels_colors))



#################### RATIO FIRMICUTES/BACTEROIDES SY (Genotype) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_SY)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_SY[,c("FASTQ_ID","Batch","Strain","Genotype")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Genotype, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Genotype=="129Sv","Ratios_Mean"]<-Ratios_Mean["129Sv"]
ratio_fb[ratio_fb$Genotype=="C57Ola","Ratios_Mean"]<-Ratios_Mean["C57Ola"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Genotype, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Genotype=="129Sv","Ratios_st_err"]<-Ratios_st["129Sv"]/sqrt(length(which(ratio_fb$Genotype=="129Sv")))
ratio_fb[ratio_fb$Genotype=="C57Ola","Ratios_st_err"]<-Ratios_st["C57Ola"]/sqrt(length(which(ratio_fb$Genotype=="C57Ola")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/SY/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Genotype_SY.csv", ratio_fb)

ratio_fb$Genotype<-factor(ratio_fb$Genotype, levels = c("129Sv","C57Ola"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-wilcox.test(Ratio ~ Genotype, paired=F, data=ratio_fb)
p_val<-round(res$p.value, digits = 5)
p_val

model<-summary(lm( rank(Ratio) ~ Batch + Strain + Genotype, data=ratio_fb))
model$coefficients[4,4] # genotype p-value

# exporting the results
con <- file("Results/SY/Ratio_Firmi_Bacteroi/Genotype_Mann_Whi_Wilcoxon_SY.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Genotype at SY :   V=",res$statistic, "p-value=", p_val, "\n\n\n",fill=T)
cat("Linear Model ~Batch+Strain+Genotype on ranks ... \n")
print(model) # printing linear model too
cat("\n\nPLUS: The model is still significative also without ranks")
cat("\n Shaphiro Wilk test p-value:",distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ_ID, fill=Genotype, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("C57Ola"="chartreuse",
                               "129Sv"="coral")) +
  facet_grid2(.~Genotype, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Genotype"))+
  scale_y_sqrt(breaks=c(1,2,3,4,5,seq(10,80,10)) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Mann Whitney p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/SY/Ratio_Firmi_Bacteroi/Ratio_barplot_Genotype_SY.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Genotype, color=Genotype,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  facet_grid(.~Genotype, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean + Ratios_st_err, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0, # the if else is needed to avoid geom bar below zero
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between genotypes at SY",
       caption = paste("Mann Whitney p-value:", round(res$p.value,5),
                       "\n Linear Model (on ranks) p-value:", round(model$coefficients[4,4],5) ))
ggsave(file="Results/SY/Ratio_Firmi_Bacteroi/Jitter_with_mean_and_STerr_Genotype.png",width=4,height=3.5, dpi=300) 
dev.off()


####### BOX PLOT
ggplot(data=ratio_fb, aes(x=Genotype, color=Genotype,y=Ratio)) +
  theme_classic(base_size =10) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  facet_grid(.~Genotype, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_boxplot(width = 0.8, size=0.3,
               show.legend = F, aes(color=NULL)) +
  geom_point(size=1.3, aes(color=Genotype)) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between genotypes at SY",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/SY/Ratio_Firmi_Bacteroi/Boxplot_Firm_Bacteroid_Genotype.png",width=3.5,height=3, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY SY (Genotype) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv"="coral","C57Ola"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype"
       # title="Alpha diversity between 129Sv and C57Ola at SY"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9.5)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/SY/Alfa_div_ Genotype_ both strains _at SY.png", width = 5,height =4.2, dpi=300)

system("echo 'The two bacthes have equilibrated library size, then the SY alpha diversity p-value is quite reliable' > Results/SY/Note_regarding_alpha_div_results.txt ")


# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Genotype
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


### AGAIN, ONLY ON HEALTHY SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Condition=="Healthy"), measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Genotype<-paste(pAlpha$data$Genotype, pAlpha$data$Strain, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv_Ntg"="coral","C57Ola_Ntg"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype"
       # title="Alpha diversity at SY (Only Ntg mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/SY/Alfa_div_Genotype_ONLY HEALTHY_SY.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


### AGAIN, ONLY ON SLA SUBJECTS
pAlpha<-plot_richness(subset_samples(data, Condition=="SLA"), measures=c("Shannon", "Observed"), x="Genotype")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Genotype<-paste(pAlpha$data$Genotype, pAlpha$data$Strain, sep="_")
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Genotype, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  geom_point(size=1.3, aes(color=Genotype)) +
  scale_color_manual(values = c("129Sv_G93A"="coral","C57Ola_G93A"="chartreuse")) +
  theme_classic(base_size = 10) + 
  labs(x="Genotype" 
       #title="Alpha diversity at SY (Only G93A mice)"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Genotype), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/SY/Alfa_div_Genotype_ONLY SLA_SY.png", width = 5,height =4.2, dpi=300)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY  SY (Genotype) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[3] # 3 is for the third variable in the design
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[3] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[3]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[3] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Strain+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[3,],perm_g$aov.tab[3,],perm_f$aov.tab[3,],perm_o$aov.tab[3,],perm_c$aov.tab[3,],perm_p$aov.tab[3,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/SY/Beta_div_Genotype/Beta_diversity_permanova_BOTH STRAINS_at SY.csv",quote=F,row.names = T)

# 
# ### PLUS: checking it with Bray too
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Genotype, data=metadata, permutations = 9999, method="bray")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/SY/Beta_div_Genotype/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)
# 

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/SY/Beta_div_Genotype/Beta_disp_permanova_BOTH STRAINS_SY.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS SY (Genotype) #####################

# on Genera
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape="Strain") +
  scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_SY.png", width = 9, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 13) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(
    #title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset",
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA__Hellinger_Genera_no_ellipse.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_points.png", width = 9, height = 6, dpi=300)
# without names and including Ntg/G93A
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype", shape = "Strain") +
  scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_Genera_points2.png", width = 7, height = 5, dpi=300)

# # Genera in 3D
# {Dist<- vegdist(t(otu_table(data.sqrt_prop)), method = "euclidean")                                                                                                      
#   obj<-ecodist::pco(Dist)
#   matrix<-obj[["vectors"]]
#   row.names(matrix)<-sample_names(data.sqrt_prop)
# }
# rgl::open3d(windowRect=c(25,25,1200,1200))
# pca3d::pca3d(matrix, col=colors, axes.color = "darkgray",
#              radius=1.5, show.shadows = T, show.plane = F, show.centroids = F)
# rgl::rgl.viewpoint(theta = -40.8, phi = 30.8, fov = 120, zoom = 0.35)
# rgl::legend3d("topleft", c("129Sv", "C57Ola"), 
#               col=c(2,3), pch = c(19, 19), magnify=0.7)
# rgl::rgl.snapshot("temp.png")
# temp<-png::readPNG("temp.png")
# # install.packages("pdftools") # it needs also "sudo apt install libpoppler-cpp-dev"
# png::writePNG(temp, "Results/SY/Beta_div_Genotype/PCoA_Hellinger_3D_on_Genera_Genotype_SY.png", dpi = 300)
# unlink("temp.png")
# rgl::close3d()
# 
# 
# # again but on ASV
# {data.prop.labels<-data.prop
#   data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues
#   eigval<- round((eigval/sum(eigval))*100, 1)
# }
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
#   scale_color_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on SY subset", 
#        color="Genotype", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
# ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA__Hellinger_on_ASV.png", width = 9, height = 6, dpi=300)


suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


##################### BETA DIVERSITY and PCoA at SY (Genotype, Only Ntg) #######################

{
  data_sub<-subset_samples(data, Condition=="Healthy")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/SY/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_SY.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/SY/Beta_div_Genotype/Beta_div_Hellinger_ONLY_Ntg_at_SY.xlsx",col.names = T,row.names = T)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/SY/Beta_div_Genotype/Beta_disp_ONLY_Ntg.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)


########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Genotype<-paste(sample_data(data.prop.labels)$Genotype, sample_data(data.prop.labels)$Strain, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola_Ntg"="chartreuse","129Sv_Ntg"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on Ntg mice at SY", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA_ONLY_Ntg.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



##################### BETA DIVERSITY and PCoA at SY (Genotype, Only G93A) #######################

{
  data_sub<-subset_samples(data, Condition=="SLA")
  data_sub.prop<-transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
  data_sub.genus<-tax_glom(data_sub.prop, taxrank = "Genus", NArm = F)
  data_sub.phylum<-tax_glom(data_sub.prop, taxrank = "Phylum", NArm = F)
}
suppressWarnings(rm(ASV.prop))
{ASV.prop<-as.data.frame(otu_table(data_sub.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data_sub.genus))
  ASV.phy.prop<-as.data.frame(otu_table(data_sub.phylum))
}


#### PERMANOVA
metadata<-as(sample_data(data_sub.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2] 
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Phyla")
beta
# write.csv2(beta, file="Results/SY/Beta_div_Genotype/Beta_div_Hellinger_ONLY_G93A_at_SY.csv",quote=F,row.names = T)
write.xlsx(beta, file="Results/SY/Beta_div_Genotype/Beta_div_Hellinger_ONLY_G93A_at_SY.xlsx",col.names = T,row.names = T)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Genotype)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/SY/Beta_div_Genotype/Beta_disp_ONLY_G93A.csv",quote=F,row.names = T)

rm(beta, perm_g, BC.dist, perm_p, perm_ASV)


########################### PCoA

# on Genera
data.prop.labels<-data_sub.genus # prop
sample_data(data.prop.labels)$Genotype<-paste(sample_data(data.prop.labels)$Genotype, sample_data(data.prop.labels)$Strain, sep="_")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype") +
  scale_color_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on G93A mice at SY", 
       color="Genotype", x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/SY/Beta_div_Genotype/Genotypes_PCoA_ONLY_G93A.png", width = 7, height = 5, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



################### DA WITH DESEQ2 SY (Genotype) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

##### STARTING THE DIFFERENTIAL ANALYSIS
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Strain + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/SY/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_SY.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/SY/DA_DESeq2_Genotype/Every_result_Genotype_DESeq2_SY.csv", row.names = F)
write.xlsx(Res_tot, file="Results/SY/DA_DESeq2_Genotype/Every_result_Genotype_DESeq2_SY.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at SY",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_every_result_SY.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Bifidobacteriaceae","Bacteroidaceae", "Atopobiaceae",
          "Rikenellaceae","Rickettsiales","Bifidobacteriales","Bifidobacteriaceae",
          "Actinobacteria","Campylobacteria","Campylobacterales","Campilobacterota",
          "Coriobacteriia","Clostridiales", "Vampirivibrionia","Cyanobacteria","Erysipelotrichaceae")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund)         # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=2,b=5, l=19) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at SY",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_SY_v1.png", width = 11, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola"="chartreuse3","129Sv"="coral2"))
ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_ both strains_no redund _at SY.png", width = 11, height = 7, dpi=300)
dev.off()

rm(plot_table)



##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("-","_", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/SY/DA_DESeq2_Genotype/CircoPlot_Genotype_ both strains_SY.png",width=1650,height=1650, res = 300)
circo
dev.off()


system(" echo 'This model has been corrected by Bacth and Genotype variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/SY/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
Genera.DESEQ2_SY<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


################### DA WITH DESEQ2 SY (Genotype, ONLY Ntg mice) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition=="Healthy")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/SY/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_SY.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/SY/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_Ntg_SY.csv", row.names = F)
write.xlsx(Res_tot, file="Results/SY/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_Ntg_SY.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola"="chartreuse","129Sv"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at SY",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_every_result_SY.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Deferribacterota", "Deferribacteres",
          "Deferribacterales", "Deferribacteraceae"
          )
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola_Ntg"="chartreuse","129Sv_Ntg"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=8) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "DA Taxa between genotypes (only Ntg mice) at SY",
    color="Genotype",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_SY_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype_mutation), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola_Ntg"="chartreuse3","129Sv_Ntg"="coral2"))
ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_ONLY_Ntg_without_redundants.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


# system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/SY/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
# Genera.DESEQ2_SY<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/SY/DA_DESeq2_Genotype/CircoPlot_Genotypes_ONLY Ntg_SY.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins



################### DA WITH DESEQ2 SY (Genotype, ONLY G93A mice) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data_sub, data.genus_pruned))
data_sub<-subset_samples(data, Condition=="SLA")
data_pruned<- prune_taxa(taxa_sums(data_sub) > 10, data_sub) 
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Genotype", "C57Ola", "129Sv"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/SY/DA_DESeq2/DA_",t,"_ratio_C57Ola_vs_129Sv_SY.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/SY/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_G93A_SY.csv", row.names = F)
write.xlsx(Res_tot, file="Results/SY/DA_DESeq2_Genotype/Every_result_Genotype_ONLY_G93A_SY.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  scale_fill_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "Differently abundant Taxa between genotypes at SY",
    y="Proportional Abundance", 
       fill="Genotype", x="")
# ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_every_result_SY.png", width = 18, height = 9, dpi=300)
dev.off()

Redund<-c("Actinobacteriota","Campilobacterota",
          "Campylobacteria","Actinobacteria",
          "Coriobacteriales",
          "Coriobacteriia","Bifidobacteriales",
          "Campylobacterales","Erysipelotrichales",
          "Atopobiaceae", "Bacteroidaceae",
          "Helicobacteraceae", "Bifidobacteriaceae",
          "Atopobiaceae")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results

plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Genotype_mutation)) + 
  theme_classic(base_size = 16) +
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("C57Ola_G93A"="chartreuse","129Sv_G93A"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        panel.grid.major.y = element_line(size=0.15, color="grey"),
        legend.key.size=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=9.5), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=0,b=5, l=37) ) +  
  #scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(
    # title= "DA Taxa between genotypes (only G93A mice) at SY",
    color="Genotype",
    y="Proportional Abundance %", 
       fill="Genotype", x="")

plot_table + # version 1
  geom_boxplot(width=0.8)
# ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_WITHOUT_redundants_SY_v1.png", width = 10.5, height = 7, dpi=300)
dev.off()

plot_table + # version 2
  geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
             aes(color=Genotype_mutation), size= 0.7, alpha= 0.65) +
  scale_color_manual(values=c("C57Ola_G93A"="chartreuse3","129Sv_G93A"="coral2"))
ggsave(filename = "Results/SY/DA_DESeq2_Genotype/DA_Genotype_ONLY_G93A_without_redundants.png", width = 10.5, height = 7, dpi=300)
dev.off()

rm(plot_table)


# system(" echo 'This model has been corrected by Batch and Strain variables,\nmoreover, every result under the arbitrary threshold of basemean=100 has been removed in order to avoid the most noisy results' > Results/SY/DA_DESeq2_Genotype/Statistical_Corrections.txt ")

# for further comparisons ...
# Genera.DESEQ2_SY<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])


##### CIRCOPLOT TO PLOT THE LOG2FC

# preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
# removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
# label settings
Res_tot2$Genus <- gsub("[","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("]_","", Res_tot2$Genus, fixed=T) # more space for the plot
Res_tot2$Genus <- gsub("_group","", Res_tot2$Genus, fixed=T) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/SY/DA_DESeq2_Genotype/CircoPlot_plot_ONLY_G93A_SY.png",width=1650,height=1650, res = 300)
circo
dev.off()

# circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))  # no margins



#################### RATIO FIRMICUTES/BACTEROIDES SY (Ntg vs G93A) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_SY)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_SY[,c("FASTQ_ID","Genotype","Strain","Batch")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Strain, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Strain=="G93A","Ratios_Mean"]<-Ratios_Mean["G93A"]
ratio_fb[ratio_fb$Strain=="Ntg","Ratios_Mean"]<-Ratios_Mean["Ntg"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Strain, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Strain=="G93A","Ratios_st_err"]<-Ratios_st["G93A"]/sqrt(length(which(ratio_fb$Strain=="G93A")))
ratio_fb[ratio_fb$Strain=="Ntg","Ratios_st_err"]<-Ratios_st["Ntg"]/sqrt(length(which(ratio_fb$Strain=="Ntg")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/SY/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Ntg_vs_G93A_SY.csv", ratio_fb)

ratio_fb$Strain<-factor(ratio_fb$Strain, levels = c("Ntg","G93A"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res_w<-wilcox.test(Ratio ~ Strain, paired=F, data=ratio_fb)
res_w

res<-summary(lm(rank(Ratio)~Batch+Genotype+Strain, data=ratio_fb))
p_val<-as.numeric(round(res$coefficients[,"Pr(>|t|)"][4], digits = 2))
p_val

system("echo 'It is NOT significant according to Wilcoxon, trying also with a linear model to correct for the batch confounding effect... still NOT significant! ' > Results/SY/Ratio_Firmi_Bacteroi/Note_about_Ratio_H_vs_G93A.txt ")

# exporting the results
con <- file("Results/SY/Ratio_Firmi_Bacteroi/Ntg_vs_G93A_Mann_Whi_Wilcoxon_SY.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Strain at SY :   V=",res_w$statistic, "p-value=", p_val, "\n\n\n\n",fill=T)
cat("Linear Model ~Batch+Genotype+Strain on ranks ... \n")
print(res) # printing linear model too
cat("\n\nPLUS: The model is still significative also without ranks")
cat("\n Shaphiro Wilk test p-value:",distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ_ID, fill=Strain, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Ntg"="deepskyblue",
                               "G93A"="red")) +
  facet_grid2(.~Strain, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Strain"))+
  scale_y_sqrt(breaks=c(1,2,3,4,5,seq(10,80,10)) ) +
  #scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("p-value linear model:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/SY/Ratio_Firmi_Bacteroi/Ratio_barplot_Ntg_G93A_SY.png",width=8,height=4, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Strain, color=Strain,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(.~Strain, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean, # to add the mean line
                    ymax = Ratios_Mean),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean + Ratios_st_err, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean - Ratios_st_err < 0, # the if else is needed to avoid geom bar below zero
                                 0, Ratios_Mean - Ratios_st_err)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between Strains at SY",
       caption = paste("Mann Whitney p-value:", round(res_w$p.value,2),
                       "\n Linear Model (on ranks) p-value:", round(res$coefficients[4,4],2) ))
ggsave(file="Results/SY/Ratio_Firmi_Bacteroi/Jitter_with_mean_and_STerr_Ntg_vs_G93A.png",width=4,height=3.5, dpi=300) 
dev.off()


####### BOX PLOT
ggplot(data=ratio_fb, aes(x=Strain, color=Strain,y=Ratio)) +
  theme_classic(base_size =10) +
  scale_color_manual(values = c("G93A"="red","Ntg"="deepskyblue")) +
  facet_grid(.~Strain, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,2,3,4,5,seq(10,80,10)) ) +
  geom_boxplot(width = 0.8, size=0.3,
               show.legend = F, aes(color=NULL)) +
  geom_point(size=1.3, aes(color=Strain)) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
        title = element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio between strains at SY",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/SY/Ratio_Firmi_Bacteroi/Boxplot_Firm_Bacteroid_STRAIN.png",width=3.5,height=3, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY SY (Ntg vs G93A) ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Strain")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Strain, y=value, color=NULL),
               size= 0.25, alpha=0.1) +
  theme_classic(base_size = 10) + 
  geom_point(size=1.3, aes(color=Strain)) +
  scale_color_manual(values = c("G93A"="coral","Ntg"="deepskyblue")) +
  labs(x="Strain"
       # , title="Alpha diversity between Ntg and G93A at SY"
       ) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=25, vjust=1, hjust=1, size=9.5)) +
  stat_compare_means(aes(group = Strain), label="p.format",
                     method = "wilcox.test", label.x= 1.35,
                     size=3.1, label.y.npc = "top", vjust=-0.38, hjust=0.32)
ggsave(file="Results/SY/Alfa_diversity_Ntg_vs_G93A_Mann_Withn_Wilcox_SY.png", width = 5,height =4.2, dpi=300)


# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Strain
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY SY (Ntg vs G93A) #######################

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

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[3] # 3 is for the third variable in the design
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[3]

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[3]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[3] 


sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain,data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[3,],perm_g$aov.tab[3,],perm_f$aov.tab[3,],perm_o$aov.tab[3,],perm_c$aov.tab[3,],perm_p$aov.tab[3,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/SY/Beta_div_Ntg_vs_G93A/Beta_diversity_permanova_Hellinger_SY.csv",quote=F,row.names = T)

# 
# ### PLUS: checking it with Bray too
# suppressWarnings(rm(sample_OTU,perm_ASV))
# sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
# perm_ASV<- vegan::adonis(sample_OTU ~Batch+Genotype+Strain, data=metadata, permutations = 9999, method="bray")
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/SY/Beta_div_Ntg_vs_G93A/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)
# 

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS SY (Ntg vs G93A) #####################

# on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Genotype_mutation<-gsub("G93A","G93A",sample_data(data.prop.labels)$Genotype_mutation)
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Genotype_mutation", shape="Strain") + # to color according to both Genotype and Strain
#   scale_color_manual(values=c("C57Ola_G93A"="coral","C57Ola_Ntg"="chartreuse2",
#                               "129Sv_G93A"="red3","129Sv_Ntg"="chartreuse4")) + # 129Sv is the bad genetic background
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        subtitle = "*hint: dark green vs dark red, then bright green vs bright red")
# ggsave(file="Results/SY/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_SY.png", width = 9, height = 6, dpi=300)
# # without ellipses
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + 
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset",
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/SY/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_Genera_no_ellipse.png", width = 9, height = 6, dpi=300)
# # without names
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
# ggsave(file="Results/SY/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_points.png", width = 9, height = 6, dpi=300)
# without names but with Genotypes too
plot_ordination(data.sqrt_prop, ordBC, color = "Genotype_mutation", shape="Strain") + # to color according to both Genotype and Strain
  scale_color_manual(values=c("C57Ola_G93A"="coral","C57Ola_Ntg"="deepskyblue",
                              "129Sv_G93A"="red3","129Sv_Ntg"="deepskyblue3")) + # 129Sv is the bad genetic background
  geom_point(size=4.2, alpha= 0.3) +
  geom_point(size=2.4, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on SY subset", 
       color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/SY/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_Genera_points2.png", width = 7, height = 5, dpi=300)


# # again but on ASV
# {data.prop.labels<-data.prop
#   data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues
#   eigval<- round((eigval/sum(eigval))*100, 1)
# }
# plot_ordination(data.sqrt_prop, ordBC, color = "Strain") +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse(size=0.2) + 
#   #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on SY subset", 
#        color="Strain", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
# ggsave(file="Results/SY/Beta_div_Ntg_vs_G93A/Ntg_vs_G93A_PCoA__Hellinger_on_ASV.png", width = 9, height = 6, dpi=300)
# 
# suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### DA WITH DESEQ2 SY (Ntg vs G93A) #####################

if(! "filter_proof" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


##### STARTING THE DIFFERENTIAL ANALYSIS
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
    tax_table(d.prop)<-as.matrix(taxa_temp)
    rm(taxa_temp) }
  
  ### starting the analysis
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch + Genotype + Strain)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Strain", "Ntg", "G93A"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    #write.csv2(r, file=paste0("Results/SY/DA_DESeq2/DA_",t,"_ratio_Ntg_vs_G93A_SY.csv"), row.names = F, quote=F, na = "")
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
# write.csv(Res_tot, file="Results/SY/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_SY.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/SY/DA_DESeq2_Ntg_vs_G93A/Every_result_Ntg_vs_G93A_DESeq2_SY.xlsx", showNA = F, col.names = T)
# 
# ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Strain)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 14) +
#   theme(strip.text.x=element_text(size=12.5,colour="black")) + 
#   scale_fill_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=16),
#         axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
#         axis.text.y = element_text(size=12), 
#         plot.title= element_text(size=18),
#         panel.grid.minor.y= element_blank() ) +   
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
#   labs(
#     # title= "Differently abundant Taxa between Strains at SY",
#     y="Proportional Abundance", 
#        fill="Strain", x="")
# # ggsave(filename = "Results/SY/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_Strain_every_result_SY.png", width = 15, height = 9, dpi=300)
# dev.off()
# 
# Table_tot2<-subset(Table_tot, Taxa %in% "Genus") # to remove redundant results
# # moreover, removing also the Cyanobacteria... a PHYLUM result with just 0.5 relative abundance is too weak (and anomalous... "Cyanobacteria!)...
# 
# plot_table<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Strain)) +
#   theme_classic(base_size = 16) +
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   scale_fill_manual(values=c("Ntg"="deepskyblue","G93A"="red")) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
#         panel.grid.major.y = element_line(size=0.15, color="grey"),
#         legend.key.size=unit(0.8,"cm"),
#         legend.text=element_text(size=16),
#         axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=15), 
#         axis.text.y = element_text(size=9.5), 
#         plot.title= element_text(size=18),
#         panel.grid.minor.y= element_blank(),
#         plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
#   labs(
#     # title= "Differently abundant Taxa between strains at SY",
#     y="Proportional Abundance %", 
#        fill="Strain", x="")
# 
# plot_table + # version 1
#   geom_boxplot(width=0.8)
# # ggsave(filename = "Results/SY/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_SY_v1.png", width = 7.3, height = 7, dpi=300)
# dev.off()
# 
# plot_table + # version 2
#   geom_boxplot(width=0.85, size= 0.5, alpha= 0.05, outlier.alpha = 0) +
#   geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.4),
#              aes(color=Strain), size= 0.7, alpha= 0.65) +
#   scale_color_manual(values=c("Ntg"="deepskyblue","G93A"="red"))
# ggsave(filename = "Results/SY/DA_DESeq2_Ntg_vs_G93A/DA_Ntg_vs_G93A_WITHOUT_redundants_SY.png", width = 7.3, height = 7, dpi=300)
# dev.off()


system(" echo 'If both the variables Batch and Genotype are taken in account then there are no significant results between Ntg and G93A at SY' > Results/SY/DA_DESeq2_Ntg_vs_G93A/Statistical_Corrections.txt ")



########## \\\\\\\\\\ STARTING THE ANALYSIS OF Time_point ALONG THE TIME \\\\\\\\\\ #############

to_reset<-ls() [! ls() %in% c("data","Metadata","unfiltered_data","contam_data","filter_proof","all_data","CircoTax2")]
to_reset<-to_reset [! to_reset %in% ls(pattern="Genera.DESEQ2") ] # to mantain DA results
rm(list = to_reset)


if(! "filter_proof" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
}

if(length(unique(sample_data(data)[["Time_point"]])) > 1 ){
  all_data<-data
  rm(data)
}

if("all_data" %in% ls()){
  data<-all_data
}
tail(sample_data(data))

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


########################### COUNTS EXPORT (TIMES) ##########################################

dir.create("Results/Diff_along_time/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Diff_along_time/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Diff_along_time/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Diff_along_time/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Diff_along_time/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Diff_along_time/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Diff_along_time/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/Diff_along_time/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Diff_along_time/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


######################## ABUNDANCES BAR PLOT (TIMES) ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one


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
tabella$Time_point<-factor(tabella$Time_point, levels=c("PS","OS","SY"))

# too long, splitting in two parts
tabella1<-tabella[tabella$Genotype=="129Sv",]
t1<-ggplot(data=tabella1, aes(x=Abundance, y=Sample, fill=Phylum)) + theme_classic(base_size =11.5) + 
  facet_grid2(Strain + Sample_code+Time_point~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"),
        strip.text.y = element_text(size=8.5, angle = 0)
        )+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_x_continuous (expand = c(0,1) ) + # to space the bars from the strips
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme( axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         #axis.text.y=element_text(angle=0, vjust=0.5, size=7),
        legend.key.size = unit(0.3, "cm"),
        title = element_text(size=12.5)) + 
  guides(fill="none") +
  labs(y="", x="Percentual abundance", 
       title = "Five most abundant phyla along times", subtitle = "129Sv mice")#caption = " 'Others' includes every phylum below rank 5 ")

tabella2<-tabella[tabella$Genotype=="C57Ola",]
t2<-ggplot(data=tabella2, aes(x=Abundance, y=Sample, fill=Phylum)) + theme_classic(base_size =11.5) + 
  facet_grid2(Strain + Sample_code+Time_point~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"),
        strip.text.y = element_text(size=8.5, angle = 0))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_x_continuous (expand = c(0,1) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        #axis.text.y=element_text(angle=0, vjust=0.5, size=7),
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.65, "cm"),
        legend.margin = margin(0,35,0,10),
        legend.title = element_text(size = 13.2),
        legend.text = element_text ( size = 12.5 )) + 
  theme(legend.position="right") +
  guides(fill=guide_legend(nrow=6)) +
  labs(y="", x="Percentual abundance",
       title = "", subtitle = "C57Ola mice",caption = " 'Others' includes every phylum below rank 5 ")

png(filename = "Results/Diff_along_time/Abundances/TOP_5_phyla.png", width = 5000, height = 2900, res=300)
ggarrange(t1, t2, nrow = 1)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Diff_along_time/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, prune.dat_top,tabella, tabella_top, tabella1, t1, tabella2, t2)



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
tabella$Time_point<-factor(tabella$Time_point, levels=c("PS","OS","SY"))

# too long, splitting in two parts
tabella1<-tabella[tabella$Genotype=="129Sv",]
t1<-ggplot(data=tabella1, aes(x=Abundance, y=Sample, fill=Genus)) + theme_classic(base_size =11.5) + 
  facet_grid2(Strain + Sample_code+Time_point~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"),
        strip.text.y = element_text(size=8.5, angle = 0))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_x_continuous (expand = c(0,1) ) +
  scale_fill_manual(values=fill_color_10) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        #axis.text.y=element_text(angle=0, vjust=0.5, size=7), 
        legend.key.size = unit(0.3, "cm"),
        title = element_text(size=12.5)) + 
  guides(fill="none") +
  labs(y="", x="Percentual abundance", 
       title = "Ten most abundant genera along times",
       subtitle = "129Sv mice")#caption = " 'Others' includes every phylum below rank 5 ")

tabella2<-tabella[tabella$Genotype=="C57Ola",]
tabella2$Genus<-gsub("_group","",tabella2$Genus)
tabella2$Genus<-factor(tabella2$Genus, levels = c(unique(tabella2$Genus)[unique(tabella2$Genus)!="Others"],"Others"))
t2<-ggplot(data=tabella2, aes(x=Abundance, y=Sample, fill=Genus)) +
  theme_classic(base_size =11.5) + 
  facet_grid2(Strain + Sample_code+Time_point~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"),
        strip.text.y = element_text(size=8.5, angle = 0))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_x_continuous (expand = c(0,1) ) +
  scale_fill_manual(values=fill_color_10) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        #axis.text.y=element_text(angle=0, vjust=0.5, size=7),
        legend.title = element_text(size = 13.2),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.65, "cm"),
        legend.margin = margin(0,-2,0,0),
        legend.text = element_text ( size = 10.8 )) + 
  theme(legend.position="right") +
  guides(fill=guide_legend(nrow=11)) +
  labs(y="", x="Percentual abundance", 
       title = "",
       subtitle = "C57Ola mice",
       caption = " 'Others' includes every genus below rank 10 ")

png(filename = "Results/Diff_along_time/Abundances/TOP_10_genera.png", width = 5000, height = 2900, res=300)
ggarrange(t1, t2, nrow = 1)
dev.off()


# means of TOP10 genera
write.xlsx(file = "Results/Diff_along_time/Abundances/TOP10_genera_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genera"= as.data.frame(tax_table(prune.dat_top))[["Genus"]] ))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top, t1,t2, tabella1,tabella2))



#################### RATIO FIRMICUTES/BACTEROIDES  (TIMES) ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-cbind.data.frame(otu_table(data_fb),tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) )
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<-Metadata$FASTQ_ID
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb,Metadata[,c("Sample_code","Time_point")])

# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Time_point, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Time_point=="PS","Ratios_Mean"]<-Ratios_Mean["PS"]
ratio_fb[ratio_fb$Time_point=="OS","Ratios_Mean"]<-Ratios_Mean["OS"]
ratio_fb[ratio_fb$Time_point=="SY","Ratios_Mean"]<-Ratios_Mean["SY"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Time_point, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Time_point=="PS","Ratios_st_err"]<-Ratios_st["PS"]/sqrt(length(which(ratio_fb$Time_point=="PS")))
ratio_fb[ratio_fb$Time_point=="OS","Ratios_st_err"]<-Ratios_st["OS"]/sqrt(length(which(ratio_fb$Time_point=="OS")))
ratio_fb[ratio_fb$Time_point=="SY","Ratios_st_err"]<-Ratios_st["SY"]/sqrt(length(which(ratio_fb$Time_point=="SY")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/Diff_along_time/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio.csv", ratio_fb)

ratio_fb$Time_point<-factor(ratio_fb$Time_point, levels = c("PS","OS","SY"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

# ratio_fb<-ratio_fb[order(ratio_fb$Sample_number,ratio_fb$Time),] # if paired
res<-kruskal.test(Ratio ~ Time_point, data=ratio_fb)
p_val<-round(res$p.value, digits = 2)
p_val

model<-summary(lm(rank(Ratio)~Sample_code+Time_point, data=ratio_fb))
model

# exporting the results
con <- file("Results/Diff_along_time/Ratio_Firmi_Bacteroi/Kruskal_Wallis AND linear_model results.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Kruskal Wallis test \n", fill=T)
cat("Ratio~Time_point :   V=",res$statistic, "p-value=", p_val, "\n",fill=T) # da adattare al kruskal
cat("\n Linear Model (on ranks):")
print(model)
cat("\n\n Shapiro Wilk test p-value: ",distr$p.value)
sink()
close(con)


# ######## BAR PLOT
# ggplot(data=ratio_fb, aes(x=Sample_code, fill=Time_point, y=Ratio)) +
#   theme_bw(base_size =12) + 
#   scale_fill_manual(values = c("PS"="orange","OS"="red","SY"="red4")) +
#   facet_grid2(.~Time_point, scales = "free_x", 
#               space="free", strip = strip_nested(size="constant"))+
#   theme(panel.spacing.x = unit(2,"pt"))+
#   scale_x_discrete (expand = c(0.01,0) ) +
#   geom_line(aes(y= ratio_fb$Ratios_Mean, group="Time_point"))+
#   scale_y_sqrt(breaks=c(1,5,seq(10,80,10)) ) +
#   # scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)) +
#   geom_bar(stat="identity", position="stack", width = 0.8) +
#   guides(fill="none") +
#   theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
#   labs(x="Mice", y="Firmicutes/Bacteroidetes Ratio",
#        subtitle = paste("Kruskal Wallis p-value:",p_val),
#        caption = "the line is the average ratio of the group")
# ggsave(file="Results/Diff_along_time/Ratio_Firmi_Bacteroi/Ratio_barplot.png",width=8,height=4, dpi=300) 
# dev.off()


####### LINE PLOT
ggplot(data=ratio_fb, aes(x=Time_point, color=Time_point,y=Ratio)) +
  theme_bw(base_size =10) +
  geom_point(size=3, alpha=0.5) +
  geom_point(size=1.8) +
  geom_line(aes(group=Sample_code), size= 0.02, color="black") +
  scale_color_manual(values = c("PS"="orange","OS"="red","SY"="red4")) +
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,1.5,3,5,seq(10,80,10)) ) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=10),
        axis.text.y=element_text(size=7)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio",
       caption = paste("Kruskal Wallis p-value:",p_val) )
ggsave(file="Results/Diff_along_time/Ratio_Firmi_Bacteroi/Ratios_line_plot.png",width=5,height=4, dpi=300) 
dev.off()

suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


########################## ALFA DIVERSITY  (TIMES) - WITH EVERY SAMPLE ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

{pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Time_point", color="Time_point")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
identical(H$Sample_code, obs$Sample_code) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
New_data<-rbind.data.frame(obs,H,ev)
pAlpha$data<-New_data
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
pAlpha$data$Time_point<-factor(pAlpha$data$Time_point, levels = c("PS","OS","SY"))
}

pAlpha +theme_classic() + 
  labs(x="Time_point"
       # title="Alpha diversity at PS, OS and SY time"
       ) +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  guides(fill="none", color="none") + 
  geom_point(size=1.3, alpha=0.5) +
  geom_line(aes(group=Sample_code), size=0.05, color= "grey")+
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11.5)) +
  stat_compare_means(aes(group = Time_point), label="p.format",
                     method = "kruskal.test", label.x= 1.3, size=3.1,
                     label.y.npc = "top", vjust=-0.2, hjust=0)
ggsave(file="Results/Diff_along_time/Alfa_div/EVERY_SAMPLE_Alfa_diversity_with_Kruskal.png",
       width = 5,height =4.2, dpi=300)


##### Trying again but as a liner model with more indipendent variables
alpha_index<-alphadt<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Observed richness",])
o<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
o
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Shannon",])
s<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
s
alpha_index<-alphadt<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Evenness",])
e<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
e
# exporting the results
con <- file("Results/Diff_along_time/Alfa_div/EVERY_SAMPLE_Linear_model_results.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("USING EVERY SAMPLE IN THE DATASET: \n\n", fill=T)
cat("In addiction to the Kruskal Wallis also more complex models have been computed (also not significant)\n\n", fill=T)
cat("rank(value) ~ Sample_code + Time_points     # Observed Richness\n", fill=T)
print(o$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Shannon\n", fill=T)
print(s$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Evenness\n", fill=T)
print(s$coefficients)
sink()
close(con)


########################## ALFA DIVERSITY  (TIMES) - SUBSETS ############################

########## 129Sv
data.temp<-subset_samples(data.genus, Genotype=="129Sv")
{pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time_point", color="Time_point")
  pAlpha
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
  pAlpha$data$Time_point<-factor(pAlpha$data$Time_point, levels = c("PS","OS","SY"))
}
pAlpha +theme_classic() + 
  labs(x="Time_point"
       # , title="Alpha diversity at PS, OS and SY (only 129Sv)"
       ) +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  guides(fill="none", color="none") + 
  geom_point(size=1.3, alpha=0.5) +
  geom_line(aes(group=Sample_code), size=0.05, color= "grey")+
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11.5)) +
  stat_compare_means(aes(group = Time_point), label="p.format",
                     method = "kruskal.test", label.x= 1.3, size=3.1,
                     label.y.npc = "top", vjust=-0.2, hjust=0)
ggsave(file="Results/Diff_along_time/Alfa_div/ONLY_129Sv_Alfa_diversity_with_Kruskal.png",
       width = 5,height =4.2, dpi=300)

# Trying again but as a liner model with more indipendent variables
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Observed richness",])
o<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
o
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Shannon",])
s<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
s
alpha_index<-alphadt<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Evenness",])
e<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
e
# exporting the results
con <- file("Results/Diff_along_time/Alfa_div/ONLY_129Sv_Linear_model_results.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("USING ONLY 129Sv SAMPLES: \n\n", fill=T)
cat("In addiction to the Kruskal Wallis also more complex models have been computed (also not significant)\n\n", fill=T)
cat("rank(value) ~ Sample_code + Time_points     # Observed Richness\n", fill=T)
print(o$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Shannon\n", fill=T)
print(s$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Evenness \n", fill=T)
print(s$coefficients)
sink()
close(con)

rm(pAlpha, con)


########## C57Ola
data.temp<-subset_samples(data.genus, Genotype=="C57Ola")
{pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time_point", color="Time_point")
  pAlpha
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
  pAlpha$data$Time_point<-factor(pAlpha$data$Time_point, levels = c("PS","OS","SY"))
}
pAlpha +theme_classic() + 
  labs(x="Time_point"
       # , title="Alpha diversity at PS, OS and SY (only C57Ola)"
       ) +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  guides(fill="none", color="none") + 
  geom_point(size=1.3, alpha=0.5) +
  geom_line(aes(group=Sample_code), size=0.05, color= "grey")+
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11.5)) +
  stat_compare_means(aes(group = Time_point), label="p.format",
                     method = "kruskal.test", label.x= 1.3, size=3.1,
                     label.y.npc = "top", vjust=-0.2, hjust=0)
ggsave(file="Results/Diff_along_time/Alfa_div/ONLY_C57Ola_Alfa_diversity_with_Kruskal.png",
       width = 5,height =4.2, dpi=300)

# Trying again but as a liner model with more indipendent variables
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Observed richness",])
o<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
o
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Shannon",])
s<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
s
alpha_index<-alphadt<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Evenness",])
e<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
e
# exporting the results
con <- file("Results/Diff_along_time/Alfa_div/ONLY_C57Ola_Linear_model_results.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("USING ONLY C57Ola SAMPLES: \n\n", fill=T)
cat("In addiction to the Kruskal Wallis also more complex models have been computed (also not significant)\n\n", fill=T)
cat("rank(value) ~ Sample_code + Time_points     # Observed Richness\n", fill=T)
print(o$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Shannon\n", fill=T)
print(s$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Evenness\n", fill=T)
print(s$coefficients)
sink()
close(con)

rm(pAlpha, con)


########## G93A
data.temp<-subset_samples(data.genus, Strain=="G93A")
{pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time_point", color="Time_point")
  pAlpha
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
  pAlpha$data$Time_point<-factor(pAlpha$data$Time_point, levels = c("PS","OS","SY"))
}
pAlpha +theme_classic() + 
  labs(x="Time_point"
       # , title="Alpha diversity at PS, OS and SY (only G93A)"
       ) +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  guides(fill="none", color="none") + 
  geom_point(size=1.3, alpha=0.5) +
  geom_line(aes(group=Sample_code), size=0.05, color= "grey")+
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11.5)) +
  stat_compare_means(aes(group = Time_point), label="p.format",
                     method = "kruskal.test", label.x= 1.3, size=3.1,
                     label.y.npc = "top", vjust=-0.2, hjust=0)
ggsave(file="Results/Diff_along_time/Alfa_div/ONLY_G93A_Alfa_diversity_with_Kruskal.png",
       width = 5,height =4.2, dpi=300)

# Trying again but as a liner model with more indipendent variables
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Observed richness",])
o<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
o
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Shannon",])
s<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
s
alpha_index<-alphadt<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Evenness",])
e<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
e
# exporting the results
con <- file("Results/Diff_along_time/Alfa_div/ONLY_G93A_Linear_model_results.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("USING ONLY G93A SAMPLES: \n\n", fill=T)
cat("In addiction to the Kruskal Wallis also more complex models have been computed (also not significant)\n\n", fill=T)
cat("rank(value) ~ Sample_code + Time_points     # Observed Richness\n", fill=T)
print(o$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Shannon\n", fill=T)
print(s$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Evenness\n", fill=T)
print(s$coefficients)
sink()
close(con)

rm(pAlpha, con)


########## Ntg (Healthy)
data.temp<-subset_samples(data.genus, Strain=="Ntg")
{pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time_point", color="Time_point")
  pAlpha
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  identical(H$Sample_code, obs$Sample_code) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
  pAlpha$data$Time_point<-factor(pAlpha$data$Time_point, levels = c("PS","OS","SY"))
}
pAlpha +theme_classic() + 
  labs(x="Time_point"
       # , title="Alpha diversity at PS, OS and SY (only Ntg)"
       ) +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  guides(fill="none", color="none") + 
  geom_point(size=1.3, alpha=0.5) +
  geom_line(aes(group=Sample_code), size=0.05, color= "grey")+
  theme(axis.text.x= element_text(angle=0, vjust=1, hjust=0.5, size=11.5)) +
  stat_compare_means(aes(group = Time_point), label="p.format",
                     method = "kruskal.test", label.x= 1.3, size=3.1,
                     label.y.npc = "top", vjust=-0.2, hjust=0)
ggsave(file="Results/Diff_along_time/Alfa_div/ONLY_Ntg_Alfa_diversity_with_Kruskal.png",
       width = 5,height =4.2, dpi=300)

# Trying again but as a liner model with more indipendent variables
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Observed richness",])
o<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
o
alpha_index<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Shannon",])
s<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
s
alpha_index<-alphadt<- as.data.frame(pAlpha$data[pAlpha$data$variable=="Evenness",])
e<-summary(lm(rank(value)~Sample_code + Time_point, data= alpha_index))
e
# exporting the results
con <- file("Results/Diff_along_time/Alfa_div/ONLY_Ntg_Linear_model_results.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("USING ONLY Ntg SAMPLES: \n\n", fill=T)
cat("In addiction to the Kruskal Wallis also more complex models have been computed (also not significant)\n\n", fill=T)
cat("rank(value) ~ Sample_code + Time_points     # Observed Richness\n", fill=T)
print(o$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Shannon\n", fill=T)
print(s$coefficients)
cat("rank(value) ~ Sample_code + Time_points     # Evenness\n", fill=T)
print(s$coefficients)
sink()
close(con)



######################## BETA DIVERSITY  (TIMES) #######################

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
  perm_ASV<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="euclidean")
  perm_ASV$aov.tab$`Pr(>F)`[2]
  perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="euclidean")
  perm_g$aov.tab$`Pr(>F)`[2]
  perm_g_H<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="euclidean")
}

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta


# PLUS: checking it with Bray too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Sample_code+Time_point, data=metadata, permutations = 9999, method="bray")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Diff_along_time/Beta_div/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)

suppressWarnings(rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV))


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
{BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,metadata$Time_point)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/Diff_along_time/Beta_div/Beta_disp_General_and_Pairwise_between_Time_point_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con)


####################### PCoA BETA DIV  (TIMES) #########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Time_point<-factor(sample_data(data.prop.labels)$Time_point, levels=c("PS","OS","SY"))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
#   geom_line(aes(group=Sample_code), color="black", size=0.05)+
#   scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
#   geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2)  +
#   geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Time_point",
#        x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
#        caption="Significant differences just between samples, NOT along the time")
# ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 9, height = 6, dpi=300)
# without names, ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  geom_point(size=3.9, alpha= 0.3) +
  geom_point(size=2.5, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  labs(
    # title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
    color="Time point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Beta_diversity_Hellinger_on_GENERA_points.png", width = 5, height = 4.2, dpi=300)
# # without names, ellipses and lines
# plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
#   scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
#   geom_point(size=3.9, alpha= 0.3) +
#   geom_point(size=2.5, alpha= 0.8) +
#   theme_classic(base_size = 12) +
#   theme(legend.margin = margin(0,2,0,-12),
#         plot.caption = element_text(size=8),
#         plot.title = element_text(size=12)
#   ) +
#   stat_ellipse(size=0.15) +
#   geom_line(aes(group="Time_point"), size=0.1, color="grey") +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Time point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Beta_diversity_Hellinger_on_GENERA_no_ellipses.png", width = 5, height = 4.2, dpi=300)

# # again but on ASV
# { data.prop.labels<-data.prop
#   data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
#   DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
#   ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
#   eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
#   eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
# }
# plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
#   scale_color_manual(values=c("PS"="chartreuse","OS"="deepskyblue","SY"="coral")) +
#   geom_line(aes(group=Sample_code), color="black", size=0.05) +
#   geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2)  +
#   labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", color="Time_point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
# ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV.png", width = 8, height = 6, dpi=300)


##### TRYING AGAIN ON SUB-SUBSETS
# on 129Sv Ntg
data.prop.labels<-subset_samples(data.genus, Genotype=="129Sv")
data.prop.labels<-subset_samples(data.prop.labels, Strain=="Ntg")
sample_data(data.prop.labels)$Time_point<-factor(sample_data(data.prop.labels)$Time_point, levels=c("PS","OS","SY"))
data.prop.labels<-transform_sample_counts(data.prop.labels, function(x) (x/sum(x)*100))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  geom_point(size=3.9, alpha= 0.3) +
  geom_point(size=2.5, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  geom_line(aes(group=Sample_code), color="black", size=0.05) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n on 129Sv subset only",
    color="Time point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Hellinger genus_ 129Sv_Ntg _only.png", width = 5, height = 4.2, dpi=300)

# on 129Sv G93A
data.prop.labels<-subset_samples(data.genus, Genotype=="129Sv")
data.prop.labels<-subset_samples(data.prop.labels, Strain=="G93A")
sample_data(data.prop.labels)$Time_point<-factor(sample_data(data.prop.labels)$Time_point, levels=c("PS","OS","SY"))
data.prop.labels<-transform_sample_counts(data.prop.labels, function(x) (x/sum(x)*100))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  geom_point(size=3.9, alpha= 0.3) +
  geom_point(size=2.5, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  geom_line(aes(group=Sample_code), color="black", size=0.05) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n on 129Sv subset only",
    color="Time point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Hellinger genus_ 129Sv_G93A _only.png", width = 5, height = 4.2, dpi=300)


# on C57Ola G93A
data.prop.labels<-subset_samples(data.genus, Genotype=="C57Ola")
data.prop.labels<-subset_samples(data.prop.labels, Strain=="G93A")
sample_data(data.prop.labels)$Time_point<-factor(sample_data(data.prop.labels)$Time_point, levels=c("PS","OS","SY"))
data.prop.labels<-transform_sample_counts(data.prop.labels, function(x) (x/sum(x)*100))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  geom_point(size=3.9, alpha= 0.3) +
  geom_point(size=2.5, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  geom_line(aes(group=Sample_code), color="black", size=0.05) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n on C57Ola subset only",
    color="Time point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Hellinger genus_ C57Ola_G93A _only.png", width = 5, height = 4.2, dpi=300)


# on C57Ola G93A
data.prop.labels<-subset_samples(data.genus, Genotype=="C57Ola")
data.prop.labels<-subset_samples(data.prop.labels, Strain=="Ntg")
sample_data(data.prop.labels)$Time_point<-factor(sample_data(data.prop.labels)$Time_point, levels=c("PS","OS","SY"))
data.prop.labels<-transform_sample_counts(data.prop.labels, function(x) (x/sum(x)*100))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time_point") +
  scale_color_manual(values=c("PS"="orange","OS"="red","SY"="red4")) +
  geom_point(size=3.9, alpha= 0.3) +
  geom_point(size=2.5, alpha= 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.margin = margin(0,2,0,-12),
        plot.caption = element_text(size=8),
        plot.title = element_text(size=12)
  ) +
  stat_ellipse(size=0.15) +
  geom_line(aes(group=Sample_code), color="black", size=0.05) +
  labs(
    # title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n on C57Ola subset only",
    color="Time point", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Diff_along_time/Beta_div/PCoA_Hellinger genus_ C57Ola_Ntg _only.png", width = 5, height = 4.2, dpi=300)


##################### DA WITH DESEQ2 (TIMES)  #######################

if(! "unfiltered_data" %in% ls() ){
  cat("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}


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
test1<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.01) {  # pt is the confidence level (-->p-value)
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Sample_code+Time_point)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Time_point)
    }
    # DE<-DESeq(d, test="LRT",reduced= ~ 1)
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Time_point", nume[i], deno[i]), test="Wald")
      #res<-results(DE, contrast=c("Time_point", nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      #res<-lfcShrink(DE, type="normal", contrast=c("Time_point",nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
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

nume<-c("PS","OS","PS")
deno<-c("OS","SY","SY")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
test1(data_pruned,ranks,nume,deno,"DE_results.txt",1,0.01)   #(automatically computes all DESeq2 analyses)
DE_results <- read.delim("DE_results.txt", header=FALSE)
DE_results<-DE_results[,-c(1,17)]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
unlink("DE_results.txt")
DE_results<-DE_results[DE_results$BaseMean>100, ]
system(" echo 'Every result under the arbitrary threshold Basemean=100 has been removed in order to avoid the noisiest results' > Results/Diff_along_time/DA_DESeq2/NB_results_are_filtered.txt")
write.csv2(DE_results, file="Results/Diff_along_time/DA_DESeq2/DE_every_results.csv", quote=F, row.names = F)
write.xlsx(DE_results, file="Results/Diff_along_time/DA_DESeq2/DE_every_results.xlsx",showNA = F, col.names = T, row.names = F)

{target<-DE_results[DE_results$Rank=="Class","ASV"] # only one result and at this taxonomic rank
  target<-prune_taxa(target, data_pruned.class.prop)
  tabella<-psmelt(target)
  tabella$Taxa<-"Classes"
  tabella[,"Phylum"]<-NULL
  colnames(tabella)[colnames(tabella)=="Class"]<-"Bacteria"
  tabella_tot<-tabella
}

# building segment plot basics
tabella_tot<-tabella_tot[order(match(tabella_tot$Time_point, levels(tabella_tot$Time_point))), ]
tabella_tot$Xaxis<-paste0(tabella_tot$Bacteria, tabella_tot$Time_point)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_tot$Xaxis<-factor(tabella_tot$Xaxis, levels = unique(tabella_tot$Xaxis))
# to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
levels(tabella_tot$Xaxis)

plot_1<-ggplot(tabella_tot, aes(x= Xaxis, y=Abundance, fill=Time_point)) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("PS"="orange","OS"="red","SY"="red4")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = "Classes")~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time_point), size=2.2) +
  geom_point(aes(color=Time_point), size=4.8, alpha=0.4) +
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,2))) +
  geom_line(aes(group=Sample_code), size=0.08, color="grey") +
  theme(strip.text.x=element_text(size=14,colour="black"), 
        strip.switch.pad.wrap = unit(10,"line")  ) + 
  theme(axis.text.y = element_text(size=12))+
  scale_x_discrete(labels=unique(levels(tabella_tot$Time_point)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=18)) +
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
plot_1 +
  labs(
    # title= "Differently abundant taxa along times",
    y="Proportional Abundance", fill="Time_point", 
       x="", 
    # caption="Significant difference between PS (higher) and SY (lower)\np-adj:0.006 , log2FC: 1.5"
    )
ggsave(filename = "Results/Diff_along_time/DA_DESeq2/DA_taxa_along_times.png", width = 8, height = 7, dpi=300)
dev.off()


##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

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
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))

