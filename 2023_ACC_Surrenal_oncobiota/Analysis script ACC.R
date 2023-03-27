##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
library("ggplot2")
library("vegan")
library("ggpubr")
library("ggh4x")
library("egg")
library("DESeq2")
library("stringr") #
library("xlsx")  
library("dplyr")
library("qiime2R")
}

{dir.create("Data_check")
dir.create("Data_check/PCoA_test")
dir.create("Results")
}

Steps<-c("ACC_stages","Healthy_vs_ACC","Secretion") # to create each folder
for(x in Steps){
  dir.create(paste0('Results/',x))
  dir.create(paste0('Results/',x,'/Abundances'))
  dir.create(paste0('Results/',x,'/Beta_div'))
  dir.create(paste0('Results/',x,'/DA_DESeq2/'))
}
rm(Steps)


options(scipen = 100) # disable scientific annotation


####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample_names(data)<-gsub("-88-ML18-","-88-ML38-",sample_names(data)) # duplicated names between batches
sample_names(data)<-gsub("-87-ML17-","-87-ML37-",sample_names(data))
sample<-sample_names(data)
original_names<-sample
sample
{sample<-substring(sample, first = 11, last= 20)
sample<-gsub("-[A-Z][0-9]-[A-Z][0-9][0-9]","",sample)
sample<-gsub("-[A-Z][0-9]-[A-Z][0-9]","",sample)
sample<-gsub("-PL2-[A-Z]","",sample)
sample<-gsub("-[0-9][0-9]-","",sample)
sample<-gsub("-[0-9]-","",sample)
sample<-gsub("-A2-[A-B]","",sample)
sample<-gsub("-","",sample)
}
sample_names(data)<-sample # update

Metadata <- as.data.frame(readxl::read_excel("Metadata.xlsx"))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,sample)

sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels=c("Healthy","ACC")) # decides the order in plots
head(sample_data(data))


########## IMPORTING ALSO THE CONTAMINATED ONE (FROM HOST DNA)

data_contam<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/eucar_contaminants/contam_taxonomy.qza", tree = "QIIME/eucar_contaminants/rooted-tree_with_contam.qza")
sample_names(data_contam)<-gsub("-88-ML18-","-88-ML38-",sample_names(data_contam)) # duplicated names between batches
sample_names(data_contam)<-gsub("-87-ML17-","-87-ML37-",sample_names(data_contam))
# changing names (already done before, check up there!)
if(identical(original_names,sample_names(data_contam))){
  sample_names(data_contam)<-sample_names(data)
  }
if(identical(sample_names(data_contam), row.names(Metadata))) {
  sample_data(data_contam)<-Metadata
  cat("Everything alright there!")} else { cat("\n Sample names are not matching, check them \n\n") }


# save.image("data_pre_noise_filtering.RData")


############### CHECKING THE HOST CONTAMINATION IN THE RAW FASTQ (NO BOWTIE) ###################

### if simil Bayes classification (unclassified are NA phyla)
data_contam_k<- tax_glom(data_contam, taxrank = "Phylum", NArm = F) # without Bowtie
data_contam_k<- transform_sample_counts(data_contam_k, function(x) (x/sum(x))*100)
data_contam_k<- subset_taxa(data_contam_k, is.na(tax_table(data_contam_k)[,"Phylum"]) )
write.csv2(file="Data_check/Host_Contamin_proportion_in_raw_FASTQ.csv", cbind.data.frame(tax_table(data_contam_k)[,c("Kingdom","Phylum")],otu_table(data_contam_k) ))

rm(data_contam_k)


#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus")
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see  PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Kingdom","Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_mean_0005_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured"] # to avoid the removal of other uncultured genera
data<-subset_taxa(data, ! Genus %in% filtered)
data<-prune_taxa(taxa_sums(data) > 0, data) 
rm(filtered, data.genus.temp)


######## Checking unassigned in prokaryote kingdom
{Unass<-tax_glom(unfiltered_data, taxrank = "Kingdom") # or domain
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
####### now filtering out every Unassigned ASV
data<- subset_taxa(data, ! Kingdom %in% c("Unassigned","d__Eukaryota") ) # just some ASV, but filtering them out anyway
head(tax_table(data))
write.csv2(tax_table(data), file="Data_check/Every_filtered_ASV_and_taxonomic_assignment.csv", row.names = T)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


proof1<- "Marker of the filtering, it is required for the script"



################# CHECKING THE COMPOSITION AFTER FILTERING ####################

#data.unf.prop<-transform_sample_counts(unfiltered_data, function(x) (x/sum(x))*100) # if BOWTIE STEP IS NOT PERFORMED
data.unf.prop<-transform_sample_counts(data_contam, function(x) (x/sum(x))*100)


if(! "proof1" %in% ls()){
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
p1<-plot_ordination(data.sqrt_prop, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  guides(color="none",shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3200, height = 1800, res=300)
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
p1<-plot_ordination(data.prop.labels, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  guides(color="none",shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA sqrt Bray-Curtis (on proportional ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.prop.labels))
data.prop.labels<-data.prop
{DistBC = phyloseq::distance(data.prop.labels, method = "bray")
  DistBC = sqrt(DistBC)
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_sqrt_BRAY_test.png", width = 3200, height = 1800, res=300)
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
p1<-plot_ordination(data.sqrt_prop, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  guides(color="none",shape="none") +
  geom_point(size=2.5) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  geom_point(size=2.5) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3200, height = 2000, res=300)
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
p1<-plot_ordination(data.prop.labels, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


############ COMPARISON BETWEEN THE TWO BATCHES

#### BATCH 1
suppressWarnings(rm(data.prop.labels1, data.sqrt_prop, p1))
data.prop.labels1<-subset_samples(data, Batch=="1_o")
data.sqrt_prop<-transform_sample_counts(data.prop.labels1, function(x) (x/sum(x))*100) # square root of proportion
data.sqrt_prop<-transform_sample_counts(data.sqrt_prop, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p1<-plot_ordination(data.sqrt_prop, ordBC, shape="Secretion", color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels1)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Hellinger \n\n 1o Batch", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### BATCH 2
data.prop.labels2<-subset_samples(data, Batch=="2_o")
data.sqrt_prop<-transform_sample_counts(data.prop.labels2, function(x) (x/sum(x))*100) # square root of proportion
data.sqrt_prop<-transform_sample_counts(data.sqrt_prop, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC,  color = "Condition") +
  scale_color_manual(values=c("ACC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.prop.labels2)$Sample_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n\n 2o Batch", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_between_batches.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels1.data.prop.labels2))


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
    check<-"red4"
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


################## ///// STARTING HEALTHY VS ACC ANALYSIS ///// #######################

# sample_names(data)<-sample_data(data)$Sample_ID

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


########################### COUNTS EXPORT HEALTHY_VS_ACC ##########################################

dir.create("Results/Healthy_vs_ACC/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Raw_counts/counts_genus.csv",quote=F)
  }

options(scipen = 100)
dir.create("Results/Healthy_vs_ACC/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  }
{write.xlsx2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  write.xlsx2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_family.xlsx",showNA = F, col.names = T )
  write.xlsx2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_class.xlsx",showNA = F, col.names = T)
  write.xlsx2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Healthy_vs_ACC/Abundances/Relative_abundances/counts_order.xlsx",showNA = F, col.names = T)
  }


###################### ABUNDANCES WITHOUT CONTAMINATION ############################

suppressWarnings(rm(top, others, tabella))
data.temp<-tax_glom(data_contam, taxrank = "Phylum", NArm=F)
data.temp<-transform_sample_counts(data.temp, function(x) (x/sum(x))*100 )
{top <- names(sort(taxa_sums(data.temp), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.temp)
  others<-taxa_names(data.temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}
### plot phylum Healthy vs ACC
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Phylum)) + 
  facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Data_check/TOP_5_phyla_without_decontam.png",width=8,height=5.5, dpi=300)
dev.off()


###################### ABUNDANCES BAR PLOT HEALTHY_vs_ACC ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one

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
### plot phylum Healthy vs ACC
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Phylum)) + 
  facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Healthy_vs_ACC/Abundances/TOP_5_phyla.png",width=9,height=5, dpi=300)
dev.off()
### plot phylum BETWEEN BATCHES (TO CHECK)
tabella$Batch<-gsub("1_o","First batch",tabella$Batch)
tabella$Batch<-gsub("2_o","Second batch",tabella$Batch)
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Phylum)) + 
  facet_grid(cols= vars(Batch),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla between batches", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Data_check/TOP_5_phyla_between_Batches.png",width=9,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Healthy_vs_ACC/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
### plot genera Condition
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Healthy_vs_ACC/Abundances/TOP_5_genera.png",width=9,height=5,dpi=300)
dev.off()

rm(top, tabella, tabella_top)


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
### plot Genera 8 Condition
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Healthy_vs_ACC/Abundances/TOP_8_genera.png",width=9,height=5,dpi=300)
### plot Genera 8 BETWEEN BATCH (TO CHECK)
tabella$Batch<-gsub("1_o","First batch",tabella$Batch)
tabella$Batch<-gsub("2_o","Second batch",tabella$Batch)
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Batch),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera (Between Batches)", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Data_check/TOP_8_genera_between_batches.png",width=9,height=5,dpi=300)


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


########################## ALFA DIVERSITY Healthy_vs_ACC ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Condition", title="Alpha diversity between ACC and Healthy patients") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 0.7, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="Results/Healthy_vs_ACC/Alfa_diversity_with_Mann_Withn_Wilcox.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


######################## BETA DIVERSITY Healthy_vs_ACC #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Condition, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Condition, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Condition, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Healthy_vs_ACC/Beta_div/Beta_diversity_permanova_Helling.csv",quote=F,row.names = T)


### PLUS: checking it without BATCH correction
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Healthy_vs_ACC/Beta_div/Beta_divers_WITHOUT_BATCH_CORRECTION.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Condition)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Healthy_vs_ACC/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


######################## PCoA BRAY HEALTHY_vs_ACC #####################

# on ASV
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC,shape = "Batch" , color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Healthy_vs_ACC/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_BATCH.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, shape="Batch",color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
ggsave(file="Results/Healthy_vs_ACC/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 8, height = 6, dpi=300)
# again but without Batches
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       caption="NB: the batch variable has not been considered here",
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Healthy_vs_ACC/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_points2.png", width = 9, height = 6, dpi=300)

rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels)


# on Genera (main target, due to the different batches)
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
# without Batches
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       caption="NB: the batch variable has not been considered here",
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Healthy_vs_ACC/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 8, height = 6, dpi=300)
# with Batches
plot_ordination(data.sqrt_prop, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Healthy_vs_ACC/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_BATCH.png", width = 8, height = 6, dpi=300)
# with Batches, no names
plot_ordination(data.sqrt_prop, ordBC, shape="Batch", color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Healthy_vs_ACC/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_Points.png", width = 8, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### DA WITH DESEQ2 Healthy_vs_ACC #####################

if(! "unfiltered_data" %in% ls() ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch+Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Healthy", "ACC"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
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
    # write.csv2(r, file=paste0("Results/Healthy_vs_ACC/DA_DESeq2/DA_",t,"_ratio_Healthy_vs_ACC.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Healthy_vs_ACC/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/Healthy_vs_ACC/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.5, 1,2,3,5,8,seq(10,max(Table_tot$Abundance+5),5))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance", 
       fill="Condition", x="")
ggsave(filename = "Results/Healthy_vs_ACC/DA_DESeq2/DA_Condition_every_result.png", width = 12, height = 8, dpi=300)
dev.off()

Redund<- c("Erwiniaceae","Pseudomonadaceae","Yersiniaceae","Rhodanobacteraceae","Enterobacterales")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  scale_fill_manual(values=c("Healthy"="chartreuse","ACC"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 30, vjust=1, hjust=1, size=11), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.5, 1,2,3,5,8,seq(10,max(Table_tot$Abundance+5),5))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa between Healthy and ACC", y="Proportional Abundance", 
       fill="Condition", x="")
ggsave(filename = "Results/Healthy_vs_ACC/DA_DESeq2/DA_Condition_no_redundants.png", width = 10, height = 7, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Healthy_vs_ACC/DA_DESeq2/NB.txt ")


###################### ///// STARTING ACC_Stages Analysis ///// #############################

to_reset<-ls()[! ls() %in% c("data","original_names","proof1") ]
rm(list = to_reset)

if(! "proof1" %in% ls() ){
  cat("\n\n The noises filtering proof is missing (?) \n\n")
  Sys.sleep(2)
}

# preparing the base object
data_sub<-subset_samples(data, Condition!="Healthy")
data_sub<-prune_taxa(taxa_sums(data_sub)>0, data_sub)
sample_data(data_sub)$ACC_Stage<-factor(sample_data(data_sub)$ACC_Stage, levels=c("I-II","III-IV"))

{data.phy = tax_glom(data_sub, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data_sub, taxrank = "Class", NArm = F)
  data.order = tax_glom(data_sub, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data_sub, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data_sub, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
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


################### COUNTS EXPORT ACC_stages (ONLY PROPORTIONALS) ##########################################

options(scipen = 100)
dir.create("Results/ACC_stages/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/ACC_stages/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT ACC_stages ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one


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
### plot phylum
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Phylum)) + 
  facet_grid(cols= vars(ACC_Stage),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/ACC_stages/Abundances/TOP_5_phyla.png",width=9,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/ACC_stages/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
### plot genera
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(ACC_Stage),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/ACC_stages/Abundances/TOP_5_genera.png",width=9,height=5,dpi=300)
dev.off()

rm(top, tabella, tabella_top)


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
### plot Genera 8 ACC_Stage
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(ACC_Stage),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/ACC_stages/Abundances/TOP_8_genera.png",width=9,height=5,dpi=300)
### plot Genera 8 ACC_Stages BETWEEN BATCHES
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + 
  facet_grid(cols= vars(Batch),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera between batches", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/ACC_stages/Abundances/TOP_8_genera_BATCHES_STAGES.png",width=9,height=5,dpi=300)


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



########################## ALFA DIVERSITY ACC_stages ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data_sub, measures=c("Shannon", "Observed"), x="ACC_Stage")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=ACC_Stage, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="ACC_stage", title="Alpha diversity between I-II and III-IV stages") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = ACC_Stage), label="p.format", method = "wilcox.test", label.x= 0.85, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="Results/ACC_stages/Alfa_diversity_with_Mann_Withn_Wilcox.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$ACC_Stage
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


######################## BETA DIVERSITY ACC_stages #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Batch+ACC_Stage, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+ACC_Stage, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+ACC_Stage, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+ACC_Stage, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+ACC_Stage, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+ACC_Stage, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/ACC_stages/Beta_div/Beta_diversity_permanova_Helling.csv",quote=F,row.names = T)


### PLUS: checking it without BATCH correction
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~ACC_Stage, data=metadata, permutations = 9999, method="euclidean")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/ACC_stages/Beta_div/Beta_divers_WITHOUT_BATCH_CORRECTION.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist, metadata$ACC_Stage)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/ACC_stages/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


######################## PCoA BRAY ACC_stages #####################

# on ASV
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC,shape = "Batch" , color = "ACC_Stage") +
  scale_color_manual(values=c("I-II"="coral","III-IV"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       color="ACC_stage", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/ACC_stages/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_BATCH.png", width = 9, height = 6, dpi=300)
# again but with Gender
plot_ordination(data.sqrt_prop, ordBC, shape="Gender", color = "ACC_Stage") +
  scale_color_manual(values=c("I-II"="coral","III-IV"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       caption="NB: no batch variable included",
       color="ACC_stage", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/ACC_stages/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_GENDER.png", width = 9, height = 6, dpi=300)
# again but with Secretion
plot_ordination(data.sqrt_prop, ordBC, shape="Secretion",color = "ACC_Stage") +
  scale_color_manual(values=c("I-II"="coral","III-IV"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       color="ACC_stage", caption="NB: no batch variable included",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
ggsave(file="Results/ACC_stages/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_Secretion.png", width = 8, height = 6, dpi=300)

rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels)

# on Genera (main target, due to the different batches)
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, shape="Batch",color = "ACC_Stage") +
  scale_color_manual(values=c("I-II"="coral","III-IV"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       caption="NB: no batch variable included",
       color="ACC_stage", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/ACC_stages/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 8, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### DA WITH DESEQ2 ACC_stages #####################

if(! "proof1" %in% ls() ){
  cat("\n\n The noises filtering proof is missing (?) \n\n")
  Sys.sleep(2)
}

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch+ACC_Stage)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("ACC_Stage", "I-II", "III-IV"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("Results/ACC_stages/DA_DESeq2/DA_",t,"_ratio_I-II_vs_III-IV.csv"), row.names = F, quote=F, na = "")
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
# 
# View(Res_tot)
# write.csv(Res_tot, file="Results/ACC_stages/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/ACC_stages/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)
# 
# ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=ACC_stage)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   scale_fill_manual(values=c("I-II"="coral","III-IV"="red4")) +
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
#         axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
#         axis.text.y = element_text(size=12), plot.title= element_text(size=18),
#         panel.grid.minor.y= element_blank() ) +   
#   scale_x_discrete(expand=c(-0.2, 1)) +
#   #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
#   #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
#   labs(title= "Differently abundant Taxa", y="Proportional Abundance", 
#        fill="ACC_stage", x="")
# ggsave(filename = "Results/ACC_stages/DA_DESeq2/DA_ACC_stage_every_result.png", width = 15, height = 8, dpi=300)
# dev.off()
# 
# # Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
# ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=ACC_stage)) + 
#   facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
#   geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
#   scale_fill_manual(values=c("I-II"="coral","III-IV"="red4")) +
#   theme(strip.text.x=element_text(size=14,colour="black")) + 
#   guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
#         legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
#         axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
#         axis.text.y = element_text(size=12), 
#         plot.title= element_text(size=18),
#         panel.grid.minor.y= element_blank(),
#         plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
#   scale_x_discrete(expand=c(-0.2, 1)) +
#   #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
#   #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
#   labs(title= "Differently abundant Taxa", y="Proportional Abundance", 
#        fill="ACC_stage", x="")
# ggsave(filename = "Results/ACC_stages/DA_DESeq2/DA_ACC_stage_no_redundants.png", width = 12, height = 8, dpi=300)
# dev.off()

system(" echo 'There are no DA results here!' > Results/ACC_stages/DA_DESeq2/Anything_there.txt ")


###################### ///// STARTING Ormone Secretion (ACC) Analysis ///// #############################

to_reset<-ls()[! ls() %in% c("data","original_names","proof1") ]
rm(list = to_reset)

if(! "proof1" %in% ls() ){
  cat("\n\n The noises filtering proof is missing (?) \n\n")
  Sys.sleep(2)
}

# preparing the base object
data_sub<-subset_samples(data, Condition!="Healthy")
data_sub<-subset_samples(data_sub, Secretion %in% c("None","Cortisol"))
data_sub<-prune_taxa(taxa_sums(data_sub)>0, data_sub)
sample_data(data_sub)$Secretion<-factor(sample_data(data_sub)$Secretion, levels=c("None","Cortisol"))

{data.phy = tax_glom(data_sub, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data_sub, taxrank = "Class", NArm = F)
  data.order = tax_glom(data_sub, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data_sub, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data_sub, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data_sub, function(ASV) ASV/sum(ASV)*100)
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


################### COUNTS EXPORT Secretion (ONLY PROPORTIONALS) ##########################################

options(scipen = 100)
dir.create("Results/Secretion/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Secretion/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT Secretion ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one


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
### plot phylum
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Phylum)) + 
  facet_grid(cols= vars(Secretion),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Secretion/Abundances/TOP_5_phyla.png",width=9,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Secretion/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
### plot genera
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Secretion),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Secretion/Abundances/TOP_5_genera.png",width=9,height=5,dpi=300)
dev.off()

rm(top, tabella, tabella_top)


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
### plot Genera 8 Secretion
ggplot(data=tabella, aes(x=Sample_ID, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Secretion),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Secretion/Abundances/TOP_8_genera.png",width=9,height=5,dpi=300)


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



########################## ALFA DIVERSITY Secretion ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data_sub, measures=c("Shannon", "Observed"), x="Secretion")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Secretion, y=value, color=NULL), alpha=0.1) + theme_bw() + 
  labs(x="Secretion", title="Alpha diversity between None versus Cortisol Secretion") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=40, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = Secretion), label="p.format", method = "wilcox.test", label.x= 1.1, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=-0.4)
ggsave(file="Results/Secretion/Alfa_diversity_with_Mann_Withn_Wilcox.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Secretion
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


######################## BETA DIVERSITY Secretion #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Batch+Secretion, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Batch+Secretion, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Batch+Secretion, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Batch+Secretion, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Batch+Secretion, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Batch+Secretion, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Secretion/Beta_div/Beta_diversity_permanova_Helling.csv",quote=F,row.names = T)


### PLUS: checking it without BATCH correction
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Secretion, data=metadata, permutations = 9999, method="euclidean")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Secretion/Beta_div/Beta_divers_WITHOUT_BATCH_CORRECTION.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist, metadata$Secretion)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Secretion/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


######################## PCoA BRAY Secretion #####################

# on ASV
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC,shape = "Batch" , color = "Secretion") +
  scale_color_manual(values=c("None"="coral","Cortisol"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       color="Secretion", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Secretion/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_BATCH.png", width = 9, height = 6, dpi=300)
# again but without_Batches
plot_ordination(data.sqrt_prop, ordBC, color = "Secretion") +
  scale_color_manual(values=c("None"="coral","Cortisol"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed ASV)", 
       caption="NB: no batch variable included",
       color="Secretion", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Secretion/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_2.png", width = 9, height = 6, dpi=300)

rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels)

# on Genera (main target, due to the different batches)
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, shape="Batch",color = "Secretion") +
  scale_color_manual(values=c("None"="coral","Cortisol"="red4")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  # geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", 
       caption="NB: no batch variable included",
       color="Secretion", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Secretion/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 8, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


################### DA WITH DESEQ2 Secretion #####################

if(! "proof1" %in% ls() ){
  cat("\n\n The noises filtering proof is missing (?) \n\n")
  Sys.sleep(2)
}

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Batch+Secretion)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Secretion", "None", "Cortisol"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("Results/Secretion/DA_DESeq2/DA_",t,"_ratio_None_vs_Cortisol.csv"), row.names = F, quote=F, na = "")
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
# 
# View(Res_tot)
# write.csv(Res_tot, file="Results/Secretion/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
# write.xlsx(Res_tot, file="Results/Secretion/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

system(" echo 'There are no DA results here!' > Results/Secretion/DA_DESeq2/Anything_there.txt ")


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
cat("\n \n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
