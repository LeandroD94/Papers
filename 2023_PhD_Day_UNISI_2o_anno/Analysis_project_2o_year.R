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
  library("dplyr")
  library("qiime2R")
}

{dir.create("Data_check")
dir.create("Data_check/PCoA_test")
dir.create("Results")
dir.create("Results/General_data_analysis")
}
for(f in c("Time_series_in_R1","Time_series_in_R2","R1_versus_R2")){
  dir.create(paste0("Results/",f))
  dir.create(paste0("Results/",f,"/Abundances"))
  dir.create(paste0("Results/",f,"/Beta_div"))
  if(f %in% c("Time_series_in_R1","Time_series_in_R2")){
    dir.create(paste0("Results/",f,"/PROPORTIONAL_counts_correlated_with_time"))
    dir.create(paste0("Results/",f,"/LOG_RATIO_counts_correlated_with_time"))
    dir.create(paste0("Results/",f,"/LOG_RATIO_counts_correlated_with_parameters"))
    
  }
}
dir.create("Results/R1_versus_R2/DA_DESeq2")
dir.create("Results/Adherent_vs_Suspended_R1_differential_abundances/")
dir.create("Results/Adherent_vs_Suspended_R2_differential_abundances/")


options(scipen = 100) # disable scientific annotation


####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("15.*F","",sample)
sample_names(data)<-sample # update

Metadata <- as.data.frame(read.csv("metadata.csv", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)]) - 1 # there is a sample (F3) collected but NOT sequenced
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,sample)

sample_names(data)<-sample_data(data)$PCR_ID

sample_data(data)$Type<-factor(sample_data(data)$Type, levels=c("R1","R2")) # decides the order in plots
sample_data(data)$Experiment_state<-factor(sample_data(data)$Experiment_state, levels=c("Still_preparing","Started","Almost_stabilized","Adherent")) # decides the order in plots
head(sample_data(data))


############ EXCLUDING THE SAMPLE COMING FROM ANOTHER PROJECT #################

if( "GRANULI_LOTTI" %in% Metadata$PCR_ID){
  data_complete_with_lotti<-data
  metadata_complete_with_lotti<-Metadata
  data<-subset_samples(data, PCR_ID!="GRANULI_LOTTI")
  Metadata<-Metadata[Metadata$PCR_ID!="GRANULI_LOTTI", ]
}

# save.image("data_PhD_including_LOTTI_sample.RData")


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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Unassigned_domain_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


### removing mitochondria
if( "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Mitochondria from the data\n")
  data<-subset_taxa(data, Genus != "Mitochondria")
}


if( ! "GRANULI_LOTTI" %in% Metadata$PCR_ID ){
  proof1<- "Marker of the filtering, it is required for the script"
  } else { stop("Wait! The other project is still among the data!")
}


############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

if(! "proof1" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

Original_number<-read.table("QIIME/Original_number_of_reads_for_sample.tsv", sep= "\t", header=T)
Original_number<-sum( c(Original_number$forward.sequence.count, Original_number$reverse.sequence.count))
Remaining_number<-sum(otu_table(data_complete_with_lotti))

con<-file("Data_check/DETAILS_ABOUT_ORIGINAL_AND_PROCESSED_READS_NUMBER.txt")
sink(con, append=TRUE)
cat("Raw reads abundance (sequenced)", fill=TRUE)
cat(Original_number)
cat("Reads abundance after every filter", fill=TRUE)
cat(Remaining_number)
cat("Percentual of remaining vs original", fill=TRUE)
cat(paste0(round(Remaining_number/Original_number*100,digits=2)),"%")
sink()
close(con)
rm(con)


################# CHECKING THE COMPOSITION AFTER FILTERING ####################

data.unf.prop<-transform_sample_counts(unfiltered_data, function(x) (x/sum(x))*100)

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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Bray-Curtis (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered BRAY
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3200, height = 1800, res=300)
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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  guides(color="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA Euclidean (on Hellinger transformed ASV)\n\n UNfiltered data", 
       color="Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered EUCLIDEAN
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_hellinger_test.png", width = 3200, height = 1800, res=300)
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
p1<-plot_ordination(data.prop.labels, ordBC, color = "Type") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data", 
       color="Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered wUNIFRAC
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Type") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Type", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_wUnifrac_test.png", width = 3200, height = 1800, res=300)
ggarrange(p1,p2, nrow = 1)
dev.off()

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


#################### FURTHER CHECKING THE HOST CONTAMINATION PROPORTION ###################

# ### if simil Bayes classification (unclassified are NA phyla)
data_contam_k<- tax_glom(data, taxrank = "Phylum", NArm = F) # without Bowtie
data_contam_k<- transform_sample_counts(data_contam_k, function(x) (x/sum(x))*100)
data_contam_k<- subset_taxa(data_contam_k, is.na(tax_table(data_contam_k)[,"Phylum"]) )
write.csv2(file="Data_check/Host_Contamin_proportion_in_raw_FASTQ.csv", cbind.data.frame(tax_table(data_contam_k)[,c("Kingdom","Phylum")],otu_table(data_contam_k) ))

rm(data_contam_k)



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
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)


################## REMOVING THE UNSATURATED SAMPLE ####################

if("S13" %in% Metadata$PCR_ID){
data<-subset_samples(data, PCR_ID!="S13")
Metadata<-Metadata[Metadata$PCR_ID!="S13", ]
proof2<-"The unsaturated sample (S13) has been removed"
}


############### BACKUP OF THE ORIGINAL (COMPLETE) DATA  #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering steps yet?\n", fill = T)
}

if( length(Metadata$Experiment_state)>3)  {
  data<-subset_samples(data, PCR_ID!="GRANULI_LOTTI")
  Metadata<-Metadata[Metadata$PCR_ID!="GRANULI_LOTTI", ]
  data_complete<-prune_taxa(taxa_sums(data)>0,data)
  metadata_complete<-Metadata
  rm(data)
}


#################### % ASSIGNED IN SILVA #########################

data<-data_complete
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


####################### GENERAL PCoA BEFORE SUBSETTING (EVERY SAMPLE) ##########################

suppressWarnings(rm(data.prop_temp, data.sqrt_prop))

if( ! "proof1" %in% ls() | ! "proof2" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}
data.prop_temp<-transform_sample_counts(data_complete, function(x) (x/sum(x))*100)

data.prop.labels<-data.prop_temp
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
# type
plot_ordination(data.sqrt_prop, ordBC, color = "Type", shape = "Experiment_state") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=2.5) +
  theme_bw(base_size = 8.5) +
  geom_line(aes(group=Sampling_date), color="gray", size= 0.12) +
  stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="General PCoA",
       color="Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_TYPE1.png", width = 6, height = 4.5, dpi=300)
# type _ only 2 elipses
plot_ordination(data.sqrt_prop, ordBC, color = "Type", shape = "Experiment_state") +
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=2.5) +
  theme_bw(base_size = 8.5) +
  geom_line(aes(group=Sampling_date), color="gray", size= 0.12) +
  stat_ellipse(size=0.2, aes(group=Type)) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="General PCoA",
       color="Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_TYPE2.png", width = 6, height = 4.5, dpi=300)
# Experiment state
plot_ordination(data.sqrt_prop, ordBC, shape = "Type", color = "Experiment_state") +
  geom_point(size=2.5) +
  theme_bw(base_size = 8.5) +
  geom_line(aes(group=Sampling_date), color="gray", size= 0.12) +
  stat_ellipse(size=0.2) + 
  geom_text(aes(label= PCR_ID), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="General PCoA",
       color="Experiment_state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_SAMPLING.png", width = 6, height = 4.5, dpi=300)
# no ellipses, R1 vs R2
plot_ordination(data.sqrt_prop, ordBC, shape = "Type", color = "Experiment_state") +
  scale_color_manual(values = c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red","Adherent"="red2")) +
  geom_point(size=2.5) +
  theme_bw(base_size = 8.5) +
  geom_line(aes(group=Sampling_date), color="black", size= 0.12) +
  geom_text(aes(label= PCR_ID), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="General PCoA",
       color="Experiment_state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption = "Lines connect the samples from R1 and R2 collected the same day")
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_no_ellipses1.png", width = 6, height = 4.5, dpi=300)
# no ellipses, same reactor ("Along time")
Line_order1<- as.character(sample_data(data.sqrt_prop)$Type)
Line_order1[sample_data(data.sqrt_prop)$Experiment_state=="Adherent" & sample_data(data.sqrt_prop)$Type=="R1"] <- "Adherent_R1"
Line_order1[sample_data(data.sqrt_prop)$Experiment_state=="Adherent" & sample_data(data.sqrt_prop)$Type=="R2"] <- "Adherent_R2"
plot_ordination(data.sqrt_prop, ordBC, shape = "Type", color = "Experiment_state") +
  scale_color_manual(values = c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red","Adherent"="red2")) +
  theme_bw(base_size = 8.4) +
  theme(legend.text = element_text(size=8.2)) +
  geom_point(size=2.7) +
  geom_line(aes(group= Line_order1), color="black", size= 0.12) +
  geom_line(aes(group= paste0(metadata_complete$Type,metadata_complete$Sampling_date)), color="red", size= 0.12) +
  geom_text(aes(label= PCR_ID), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="General PCoA (with adherent biomass sampling)",
       caption = "Gray lines connect the samples from the same reactor along the time;\nRed lines connect the adherent biomass with the suspended one of the same day\n\nNB: there are not sample F13 and S3 --> the related sudden shift is actually a gap.",
       color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_ALONG_TIME1.png", width = 6, height = 4.6, dpi=300)


#### focusing on suspended biomass only, time shift
data.prop.labels<-subset_samples(data.prop_temp, Experiment_state!="Adherent")
{sample_data(data.prop.labels)$Numeration<-as.factor(gsub("[A-Z]","",sample_names(data.prop.labels)))
sample_data(data.prop.labels)$State_for_plot<-gsub("_"," ",sample_data(data.prop.labels)$Experiment_state)
sample_data(data.prop.labels)$State_for_plot<-gsub("Almost","~",sample_data(data.prop.labels)$State_for_plot)
sample_data(data.prop.labels)$State_for_plot<-gsub("S","s",sample_data(data.prop.labels)$State_for_plot)
sample_data(data.prop.labels)$State_for_plot<-factor(sample_data(data.prop.labels)$State_for_plot,
                                                     levels=c("still preparing","started","~ stabilized" ))
}
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, shape = "Type", color = "Experiment_state") +
  scale_color_manual(values = c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  theme_bw(base_size = 8.4) +
  theme(legend.text = element_text(size=8.2)) +
  geom_point(size=2.7) +
  geom_text(aes(label= Numeration), 
            color="black", size=2,
            show.legend = FALSE) +
  geom_line(aes(group=as.factor(Type)), col= "darkgray", size = 0.12) +
           # arrow=arrow(length =unit(0.23,"cm"), type = "closed") ) +
  labs(title="General PCoA",
       caption = "Gray lines connect the samples from the same reactor along the time; \nNB: there are not sample F13 and S3 --> the related sudden shift is actually a gap.",
       color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_ALONG_TIME2.png", width = 6, height = 4.6, dpi=300)
# again, but smaller (poster)
plot_ordination(data.sqrt_prop, ordBC, shape = "Type", color = "State_for_plot") +
  scale_color_manual(values = c("still preparing"="orange","started"="coral2","~ stabilized"="red")) +
  theme_bw(base_size = 7) +
  theme(legend.text = element_text(size=6.5),
        legend.margin = margin(2,-1,2,-4)
        )+
  geom_point(size=1.8) +
  geom_point(size=3.3, alpha= 0.4) +
  geom_text(aes(label= Numeration), 
            color="black", size=2,
            show.legend = FALSE) +
  geom_line(aes(group=as.factor(Type)), col= "darkgray", size = 0.12) +
  # arrow=arrow(length =unit(0.23,"cm"), type = "closed") ) +
  labs(color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/GENERAL_PCoA_Beta_diversity_Hellinger_on_genera_ALONG_TIME3.png", width = 4, height = 2.8, dpi=300)

suppressWarnings(rm(p1,p2,data.sqrt_prop,eigval,ordBC,DistBC,data.prop.labels))


############################# GENERAL COUNTS EXPORT ##############################

data<-data_complete
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

dir.create("Results/General_data_analysis/General_Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/General_data_analysis/General_Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}



################# STUDING THE ABUNDANCES OF THE PHA ACCUMULATORS ######################

#fill_color_13<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3","blue", "gray", "yellow4","red") # "others" will be setted as the last one
fill_color_17<-c("wheat3","deeppink","darkmagenta","coral","cyan","yellow2","brown","firebrick3","springgreen2","violet","bisque","deepskyblue2","darkslategray3","blue", "gray", "yellow4","red") # "others" will be setted as the last one

table_PHA<-read.table(file="Lista_accumulatori_PHA_da_articoli_5aprile2023.tsv", fill = T, header = F, sep="\t")
colnames(table_PHA)<-"Bacteria"
lista_PHA<-table_PHA[,1]     # they derive from papers


data.genus.temp<-tax_glom(data_complete, taxrank ="Genus", NArm=F)
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))


tabella$Experiment_state<-gsub("_"," ",tabella$Experiment_state)
tabella$Experiment_state<-factor(tabella$Experiment_state, levels = c("Still preparing","Started","Almost stabilized","Adherent"))
# factoring also the sample names to ensure a certain order in the plot (=time) ; in the sample_names they are already ordered
tabella$Sample<-factor(tabella$Sample, levels=sample_names(data_complete))
tabella$Time<-gsub("[A-Z][A-Z]","",tabella$Sample)
tabella$Time<-gsub("[A-Z]","",tabella$Time)
tabella$Time<-factor(tabella$Time, levels = seq(1,length(tabella$Time),1))
tabella$Genus<-gsub("Candidatus_","",tabella$Genus) # synonyms

tabella2<-tabella[tabella$Experiment_state!="Adherent", ]
ggplot(data=tabella2, aes( x=Time, y=Abundance, fill=Genus)) + 
  facet_grid2( Type~Experiment_state, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_17)) +
  theme(axis.text.x=element_text(angle=25, hjust=1,vjust=1, size=10),
        axis.title.y = element_text(size=10), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 9.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom", legend.margin =  margin(3,14,0,0)) +
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="Sample code (~ Time)", y="Relative % abundance (on every bacteria) in the sample", 
       fill="",
       title = paste("PHA producing genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"))
ggsave(file="Results/General_data_analysis/PHA_accumulators_along_the_time.png",width=8,height=5.8,dpi=300)
dev.off()

# again but focusing on adherent biomass
tabella2<-tabella[tabella$Experiment_state %in% c("Started","Adherent","Almost stabilized"), ]
tabella2<-tabella2[tabella2$PCR_ID!="F13", ] # to pair the number of samples between two groups (necessary for the color trick after the facet grid)
tabella2$Sampling_date<-gsub("_23","",tabella2$Sampling_date) # the following lines allow to properly"order" the time column as factor
unique(tabella2$Sampling_date)
tabella2$Sampling_date[tabella2$Experiment_state=="Adherent"]<-paste0(tabella2$Sampling_date[tabella2$Experiment_state=="Adherent"],"_Adherent")
tabella2$Sampling_date[tabella2$Experiment_state!="Adherent"]<-paste0(tabella2$Sampling_date[tabella2$Experiment_state!="Adherent"],"") # ... add nothing! Row left to easily modify the plot
tabella2<-tabella2[order(tabella2$Sampling_date) , ] # days
order_time<-tabella2$Sampling_date[c(grep("_04",tabella2$Sampling_date),grep("_05",tabella2$Sampling_date),grep("_06",tabella2$Sampling_date))] # months
tabella2$Sampling_date<-factor(tabella2$Sampling_date, levels = unique(order_time))
color_set<-ifelse(grepl("Adherent",levels(tabella2$Sampling_date)),"red","black")

ggplot(data=tabella2, aes( x=Sampling_date, y=Abundance, fill=Genus)) + 
  facet_grid( .~Type, scales = "free", space="free")+
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_13)) +
  theme(axis.text.x=element_text(angle=38, hjust=1,vjust=1, size=8.2,
                                 color=color_set
                                 ),
        axis.title.y = element_text(size=10), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 9.5 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="bottom", 
        legend.margin = margin(3,33,0,0)) +
  guides(fill=guide_legend(nrow=3)
         ) + 
  labs(x="Sampling date", y="Relative % abundance (on every bacteria) in the sample", 
       fill="",
       title = paste("PHA producing genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"),
       caption= "the red labels mark the adherent biomass")
ggsave(file="Results/General_data_analysis/PHA_accumulators_Adherents_focused.png",width=8,height=5.5,dpi=300)
dev.off()


# only R1 (poster version)
tabella3<-tabella[tabella$Type=="R1",  ]
tabella3<-tabella3[tabella3$Experiment_state!="Adherent",  ]
tabella3$TimeS <- paste0("S",tabella3$Time)
tabella3$TimeS<- factor(tabella3$TimeS, levels = c(paste0("S",1:12),paste0("S",14:23)) )
ggplot(data=tabella3, aes( x=TimeS, y=Abundance, fill=Genus)) + 
  facet_grid2( .~Experiment_state, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =11.5) +
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=c(fill_color_17)) +
  theme(axis.text.x=element_text(angle=35, hjust=1,vjust=1, size=8),
        axis.title.x = element_text(size=11), 
        axis.title.y = element_text(size=10.5), 
        axis.text.y = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.x = unit(0.05, "cm"),
        legend.key.height = unit(0.42, "cm"),
        legend.text = element_text ( size = 8.1 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="right",
        # plot.margin = margin(0,0,0,0),
        legend.margin =  margin(0,-3,0,-10)
        ) +
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample code (~ Time)", y="Percentual abundance", 
       fill="",
       #title = paste("PHA producing genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )")
       )
ggsave(file="Results/Time_series_in_R1/Abundances/PHA_accumulators_along_the_time_R1.png",width=5,height=3,dpi=300)
dev.off()



#######################  !!!!!  STARTING THE ANALYSIS OF R1 vs R2 !!!!!   #############################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
to_remove<-ls()
to_remove<-to_remove[! to_remove %in% c("proof1","proof2","data_complete","metadata_complete",
                                        "unfiltered_data","data_complete_with_lotti","metadata_complete_with_lotti",
                                        "Genera.DESEQ2_R1_vs_R2","significative_abundances_R1","significative_abundances_R2")]
rm(list=to_remove)
}

data<-subset_samples(data_complete, Experiment_state %in% c("Started","Almost_stabilized"))
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


###################### ABUNDANCES BAR PLOT _ R1 versus R2 ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
# fill_color_8<-c("wheat3","darkmagenta","coral","yellow2","firebrick3","springgreen2","violet","deepskyblue2","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one


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
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
unique(tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Phylum)) + theme_classic(base_size =14) + 
  facet_grid2(Sampling_date+Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 14 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y="Sample", x="Percentual abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/R1_versus_R2/Abundances/TOP5_phyla_abundances.png",width=7,height=10.3, dpi=300) 
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/R1_versus_R2/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
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
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
unique(tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(Sampling_date+Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_10) +
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(y="Sample", x="Percentual abundance", title = "Ten most abundant genera", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/R1_versus_R2/Abundances/TOP10_genera_abundances.png",width=7,height=10.3, dpi=300) 
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/R1_versus_R2/Abundances/TOP_10_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))


########################## ALFA DIVERSITY _ R1 versus R2 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

data.genus.temp<-subset_samples(data.genus, PCR_ID!="F13") # the S13 (pair) is unsaturated
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Type", color="Type")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evenness
{identical(H$Sampling_date, obs$Sampling_date) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Type, New_data$Sampling_date),]
pAlpha$data<-New_data
pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
# with points only
pAlpha + theme_classic2() +
  scale_color_manual(values = c("R1"="chartreuse", "R2"="coral")) +
  geom_point(size=3.1, alpha= 0.4) +
  geom_line(aes(group = pAlpha$data$Sampling_date),col="grey",size=0.15) +
  labs(x="Type", title="Alpha diversity between R1 and R2 reactors") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Type), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.47) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/R1_versus_R2/Alfa_diversity_GENUS_paired_wilcoxon.png", width = 6,height =5.2, dpi=300)
# with ID
pAlpha + theme_classic2() +
  scale_color_manual(values = c("R1"="chartreuse", "R2"="coral")) +
  geom_point(size=3.1, alpha= 0.4) +
  geom_line(aes(group = pAlpha$data$Sampling_date),col="grey",size=0.15) +
  geom_text(aes(label=PCR_ID), size= 2.5, color= "black") +
  labs(x="Type", title="Alpha diversity between R1 and R2 reactors") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Type), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.47) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/R1_versus_R2/Alfa_diversity_GENUS_paired_wilcoxon_with_ID.png", width = 6,height =5.2, dpi=300)
# poster version
pAlpha + theme_classic2(base_size = 7) +
  scale_color_manual(values = c("R1"="chartreuse", "R2"="coral")) +
  geom_point(size=2.8, alpha= 0.3) +
  geom_line(aes(group = pAlpha$data$Sampling_date),col="grey",size=0.07) +
  labs(x="Type", title="Alpha diversity between R1 and R2 reactors") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=7),
        axis.text.y= element_text(angle=0, vjust=0, hjust=0.5, size=4),
        ) +
  stat_compare_means(aes(group = Type), label="p.format", method = "wilcox.test", paired = T,
                     label.x= 1.5, size=2.2, label.y.npc = "top", vjust=-0.1, hjust=0.47) +
  theme(panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/R1_versus_R2/Alfa_diversity_GENUS_paired_wilcoxon_poster_version.png", width = 3,height =2.5, dpi=300)


##################### BETA DIVERSITY _ R1 versus R2 #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Sampling_date + Type, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[2]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Sampling_date + Type, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[2]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[2] # needed later for the plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Sampling_date + Type, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Sampling_date + Type, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Sampling_date + Type, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Sampling_date + Type, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[2,],perm_g$aov.tab[2,],perm_f$aov.tab[2,],perm_o$aov.tab[2,],perm_c$aov.tab[2,],perm_p$aov.tab[2,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/R1_versus_R2/Beta_div/Beta_divers_permanova_Helling.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Type)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/R1_versus_R2/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)


######################## PCoA BETA DIV _ R1 versus R2 #########################

# on genera
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  geom_line(aes(group=Sampling_date),col="grey", size=0.15)+
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=3, alpha=0.3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sampling_date), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Type", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), 
       subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H),
       caption="lines connect paired samples")
ggsave(file="Results/R1_versus_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 7, height = 5, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  geom_line(aes(group=Sampling_date),col="grey", size=0.15)+
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=3, alpha=0.3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sampling_date), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Type", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples")
ggsave(file="Results/R1_versus_R2/Beta_div/PCoA_Beta_div_Hellinger_genera_no_ellipse.png", width = 7, height = 5, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Type") +
  geom_line(aes(group=Sampling_date),col="grey", size=0.15)+
  scale_color_manual(values=c("R2"="coral","R1"="chartreuse")) +
  geom_point(size=3, alpha=0.3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Type", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="lines connect paired samples")
ggsave(file="Results/R1_versus_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_points.png", width = 7, height = 5, dpi=300)


############## DA WITH DESEQ2 _ R1 versus R2 #############

if(! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\nDid you perform the filtering steps yet?\n", fill = T)
}


suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
data_pruned<- subset_samples(data_pruned, PCR_ID!="F13")
# removing the unpaired sample
system(" echo 'In case of pair analysis, the log2foldchange displayes the change from R1 to R2, see https://support.bioconductor.org/p/105981/ \n' > Results/R1_versus_R2/DA_DESeq2/NB.txt ")  
system(" echo 'Every result with a baseMean value lower than 50 has been arbitrarly filtered.' >> Results/R1_versus_R2/DA_DESeq2/NB.txt ")  

Table_tot<-NULL
Res_tot<-NULL

# correcting also the time (with more time, the difference increases)
sample_data(data_pruned)$time<-as.factor(gsub("^[0-9][0-9]_","",sample_data(data_pruned)$Sampling_date))
### --> error if included in the model: too much auto-correlations among the variables (the time itself is the paired factor already!)
# then, the april samples (low difference) will be discarded to avoid this bias
data_pruned2<- subset_samples(data_pruned, time!="04_23")


for( t in c("Genus","Family","Class","Order","Phylum") ){
  cat("\nWorking on",t,"level...\n")
  suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level")))
  d <- tax_glom(data_pruned2, taxrank = t, NArm = F)
  d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
  
  # if the paired samples are seen as numbers (it depends from on levels name)
  sample_data(d)[["Sampling_date"]]<-as.character(sample_data(d)[["Sampling_date"]])
  
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sampling_date +Type) # interaction: more time passed = even more differences!
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Type", "R1", "R2"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    write.csv2(r, file=paste0("Results/R1_versus_R2/DA_DESeq2/DA_",t,"_ratio_R1_vs_R2.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/R1_versus_R2/DA_DESeq2/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/R1_versus_R2/DA_DESeq2/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

Table_tot$Type<-factor(Table_tot$Type, levels = unique(Table_tot$Type))
Table_tot<-Table_tot[order(match(Table_tot$Type, levels(Table_tot$Type))), ]
# View(Table_tot)

tabella_g<-Table_tot[Table_tot$Taxa=="Genus",]
tabella_2<-Table_tot[Table_tot$Taxa%in% c("Family","Order","Class","Phylum"),]
# building segment plot basics
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Type)
tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Type)
# to prevent the reorder from ggplot2 (NB: to do after the re-order of table_tot)
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis))
tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis))

# to further check the same order of the factors (NB: otherwise the plot can invert the group labels)
levels(tabella_g$Xaxis)
levels(tabella_g$Type)


# marking the PHA accumulators
lista_PHA<-read.table(file="Lista_accumulatori_PHA_da_articoli_5aprile2023.tsv", fill = T, header = F, sep="\t")
colnames(lista_PHA)<-"Bacteria"
lista_PHA<-lista_PHA[,1]     # they derive from papers

# plotting
unique(tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_ f ","\nf_",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_ o ","\no_",tabella_g$Bacteria)
tabella_g$Bacteria<-gsub("_marine_group","\nmarine_group",tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Type)) +
  theme_bw(base_size = 9.1) +
  scale_color_manual(values = c("R2"="coral","R1"="chartreuse")) +
  facet_wrap2(nrow=4,~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Type), size=0.5, alpha=1) +
  geom_point(aes(color=Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=8.4,colour="black"),
        strip.switch.pad.wrap = unit(3,"line")  ) + 
  theme(axis.text.y = element_text(size=8))+
  scale_x_discrete(labels=unique(levels(tabella_g$Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
plot_1 +
  geom_line(aes(group=Sampling_date), size=0.2) +
  labs(title= "Differently abundant genera", y="Percentual Abundance", fill="Type", x="")
ggsave(filename = "Results/R1_versus_R2/DA_DESeq2/Plot_genera_DeSeq2_Type.png", width = 10.8, height = 7, dpi=300)
# again but with name
plot_1 +
  geom_text(aes(label=PCR_ID), size=2.3, color="black") +
  geom_line(aes(group=Sampling_date), size=0.2, col="darkgray") +
  labs(title= "Differently abundant genera", y="Percentual Abundance", fill="Type", x="")
ggsave(filename = "Results/R1_versus_R2/DA_DESeq2/Plot_genera_DeSeq2_Type_with_name.png", width = 10.8, height = 7, dpi=300)
dev.off()
# again but poster version
tabella_g2<-tabella_g
tabella_g2$Bacteria<-gsub("[a-z]_","",tabella_g$Bacteria)
ggplot(tabella_g2, aes(x= Xaxis, y=Abundance, fill=Type)) +
  theme_classic(base_size = 6) +
  scale_color_manual(values = c("R2"="coral","R1"="chartreuse")) +
  facet_wrap2(nrow=5,~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 20)),
              scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Type), size=0.05, alpha=1) +
  geom_point(aes(color=Type), size=1.2, alpha=0.2) +
  theme(strip.text.x=element_text(size=4.01,colour="black"),
        strip.switch.pad.wrap = unit(0.2,"line")  ) + 
  theme(axis.text.y = element_text(size=4))+
  scale_x_discrete(labels=unique(levels(tabella_g$Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=8)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  # scale_y_sqrt(breaks=c(0.1, 0.25,0.75,1,seq(2,max(tabella_g$Abundance+2),2))) +
  theme(legend.margin=margin(0, 0, 0, 0),
        legend.position = "none") +
  geom_line(aes(group=Sampling_date), size=0.045, color="darkgrey") +
  labs(title= "Differently abundant genera between reactor R1 and R2", y="Percentual Abundance", fill="Type", x="")
ggsave(filename = "Results/R1_versus_R2/DA_DESeq2/Plot_genera_DeSeq2_poster_version.png", width = 4.52, height = 3.2, dpi=300)


##### removing redundant results
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

Redund<-auto_redund
#Redund<-c("") # to select redundants (same ASV, same results at different levels)
tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)

tabella_3<-subset(tabella_2, Taxa %in% c("Phylum","Class"))
plot_3<-ggplot(tabella_3, aes(x= Xaxis, y=Abundance, fill=Type)) +
  theme_bw(base_size = 9) +
  scale_color_manual(values = c("R2"="coral","R1"="chartreuse")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Type), size=0.7, alpha=1) +
  geom_point(aes(color=Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=8.4,colour="black"), 
        strip.switch.pad.wrap = unit(2,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_3$Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  labs(title= "Differently abundant families and orders", y="Proportional Abundance", fill="Type", x="") +
  geom_line(aes(group=Sampling_date), size=0.2) 

tabella_4<-subset(tabella_2, ! Taxa %in% c("Phylum","Class"))
tabella_4$Bacteria<-gsub("-","\n",tabella_4$Bacteria, fixed = F)
# tabella_4$Bacteria<-gsub("Peptostreptococcales","Peptostreptococc.",tabella_4$Bacteria)
plot_4<-ggplot(tabella_4, aes(x= Xaxis, y=Abundance, fill=Type)) +
  theme_bw(base_size = 9) +
  scale_color_manual(values = c("R2"="coral","R1"="chartreuse")) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phylum","Class","Order","Family"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
  geom_point(aes(color=Type), size=0.7, alpha=1) +
  geom_point(aes(color=Type), size=2.5, alpha=0.3) +
  theme(strip.text.x=element_text(size=7.8,colour="black"), 
        strip.switch.pad.wrap = unit(2,"line")  ) + 
  theme(axis.text.y = element_text(size=8.5))+
  scale_x_discrete(labels=unique(levels(tabella_4$Type)), expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14)) +
  #theme(panel.spacing.x = unit(1, "pt"))+
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") +
  labs(title= "Differently abundant families and orders", y="Proportional Abundance", fill="Type", x="") +
  geom_line(aes(group=Sampling_date), size=0.2) 

# now an unique plot (no dates)
head_plot<-plot_4 + # refining first plot with the title
  labs(title= "Differently abundant families, orders, class and phyla", y="Percentual abundance", x="") + 
  theme(plot.title = element_text(size=12.5))
png(filename = "Results/R1_versus_R2/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants.png", width = 2750, height = 1500, res=300)
grid.arrange(head_plot,
             plot_3 + labs(title= "", y="Percentual abundance", x=""))
dev.off()

# again, an unique plot (with dates)
head_plot<-plot_4 + # refining first plot with the title
  labs(title= "Differently abundant families, orders, class and phyla", y="Percentual abundance", x="") + 
  theme(plot.title = element_text(size=12.5)) + 
  geom_text(aes(label=PCR_ID), size=2, color="black")
plot3_bis<- plot_3 +   geom_text(aes(label=PCR_ID), size=2, color="black") +
png(filename = "Results/R1_versus_R2/DA_DESeq2/Plot_TOTAL_DESeq2_Type_no redundants_WITH_DATE.png", width = 2750, height = 1500, res=300)
grid.arrange(head_plot,
             plot3_bis + labs(title= "", y="Percentual abundance", x=""))
dev.off()


# for comparisons...
Genera.DESEQ2_R1_vs_R2<-unique(tabella_g[tabella_g$Taxa=="Genus","Bacteria"])
rm(tabella_2, plot_2, plot_1, head_plot)


#######################   !!!!!  STARTING THE ANALYSIS OF TIME SERIES in R1 !!!!!   #############################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","data_complete","metadata_complete",
                                          "unfiltered_data","data_complete_with_lotti","metadata_complete_with_lotti",
                                          "Genera.DESEQ2_R1_vs_R2","significative_abundances_R1","significative_abundances_R2")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Type == "R1" )
data<-subset_samples(data, Experiment_state != "Adherent" )
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

sample_data(data)$Experiment_state<-factor(sample_data(data)$Experiment_state, levels=c("Still_preparing","Started","Almost_stabilized")) # decides the order in plots

days_ordered<-sample_data(data)$Sampling_date[order(sample_data(data)$Sampling_date)]
order_time<-days_ordered[c( grep("_03",days_ordered), grep("_04",days_ordered),grep("_05",days_ordered),grep("_06",days_ordered))] # months
sample_data(data)$Sampling_date<-factor(sample_data(data)$Sampling_date, levels = unique(order_time))

metadata<-metadata_complete[metadata_complete$Type=="R1" & metadata_complete$Experiment_state!="Adherent", ]


######################## ABUNDANCES BAR PLOT _ TIME SERIE R1 ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
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
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percentual abundance", title = "Five most abundant phyla (R1)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Time_series_in_R1/Abundances/TOP_5_phyla.png",width=7,height=4.5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Time_series_in_R1/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, prune.dat_top,tabella, tabella_top)


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
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percentual abundance", title = "Ten most abundant genera (R1)", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Time_series_in_R1/Abundances/TOP_10_Generi.png",width=7,height=4.5,dpi=300)
dev.off()

# top10 genera again, but poster version
tabella$Sample<-factor(tabella$Sample,
                       levels = c(paste0("S",1:12),paste0("S",14:23))
                       )
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + 
    theme_classic(base_size=11.5) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.title.x=element_text(size=11),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=8),
        axis.title.y = element_text(size=10.5), 
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.x = unit(0.05, "cm"),
        legend.key.height = unit(0.42, "cm"),
        legend.text = element_text ( size = 8.1 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="right",
        legend.margin =  margin(0,-3,0,-10)) +
    guides(fill=guide_legend(ncol=1)) + 
    labs(x="Sample code (~ Time)", 
         fill="",
         y="Percentual abundance")
ggsave(file="Results/Time_series_in_R1/Abundances/TOP_10_Generi_poster_version.png",width=5,height=3,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/Time_series_in_R1/Abundances/TOP_10_genera_of_all_dataset_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))


# TOP 10 Genera of only "ALMOST STABILIZED" group
data.genus.subset<-subset_samples(data.genus.prop, Experiment_state=="Almost_stabilized")
{top <- names(sort(taxa_sums(data.genus.subset), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.genus.subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=8),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percentual abundance", title = "Ten most abundant genera (R1)\n considering only the 'almost stabilized' samples", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Time_series_in_R1/Abundances/TOP_10_Generi_ONLY_STABILIZED_SAMPLES.png",width=7,height=4.5,dpi=300)
dev.off()



########### LINE PLOTS OF MOST ABUNDANT GENERA ##############


fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# TOP 5 Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =14) + 
  geom_point(aes(color=Genus), size =1.8) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_continuous(breaks = seq(0, max(tabella$Abundance), 10)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "5 most abundant genera through time (R1)")
ggsave(file="Results/Time_series_in_R1/Abundances/TOP5_Genera_along_time.png",width=7,height=5, dpi=300) 
dev.off()


########################## ALFA DIVERSITY _ TIME SERIE R1 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

dir.create("Results/Time_series_in_R1/Alpha_diversity/")
pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Experiment_state", color="Experiment_state")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
{ H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
  obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
  # adding evenness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Experiment_state<-factor(pAlpha$data$Experiment_state, levels = levels(sample_data(data)$Experiment_state))
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Experiment_state, y=value, color=NULL), alpha=0) + theme_classic2() + 
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.2, alpha= 0.4) + 
  labs(x="Experiment state", title="Alpha diversity between time groups") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=28, vjust=1, hjust=1, size=9.8)) +
  stat_compare_means(aes(group = Experiment_state), label="p.format", method = "kruskal.test",
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.4, hjust=0.18)
ggsave(file="Results/Time_series_in_R1/Alpha_diversity/Alfa diversity_on genera_Kruskal test.png", width = 6,height =5.2, dpi=300)
# again, but smaller (poster version)
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Experiment_state, y=value, color=NULL),
               alpha=0, size = 0.25) +
  theme_classic2(base_size = 8.6) + 
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=2.9, alpha= 0.3) + 
  labs(x="", title="") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5),
        axis.title.y = element_text(size=7)
        ) +
  stat_compare_means(aes(group = Experiment_state), label="p.format", method = "kruskal.test",
                     label.x= 1.5, size=2.3, label.y.npc = "top", vjust=-0.1, hjust=0.22)
ggsave(file="Results/Time_series_in_R1/Alpha_diversity/Alfa diversity_on genera_Kruskal_poster.png", width = 3.7,height =2.4, dpi=300)


# Statistical analyses
alphadt<- as.data.frame(pAlpha$data)
Obser_value<-alphadt[alphadt$variable=="Observed richness", ]
factor<-Obser_value$Experiment_state
kruskal.test(Obser_value$value~factor) # just to test the plotted p value

# pairwise comparisons
con <- file("Results/Time_series_in_R1/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_of_experiment_state.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Experiment_state sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Shannon", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Evenness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)


# correlation through time (Observed r)
Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=PCR_ID), size= 3.7) +
  geom_line(aes(group=Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Observed richness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 2))
       )
ggsave("Results/Time_series_in_R1/Alpha_diversity/Dot_plot_Observed_richness_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Shannon)
Filter_value<-alphadt[alphadt$variable=="Shannon", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=PCR_ID), size= 3.7) +
  geom_line(aes(group=Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Shannon",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/Time_series_in_R1/Alpha_diversity/Dot_plot_Shannon_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Evenness)
Filter_value<-alphadt[alphadt$variable=="Evenness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=PCR_ID), size= 3.7) +
  geom_line(aes(group=Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Evenness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/Time_series_in_R1/Alpha_diversity/Dot_plot_Evenness_through_time.png", width = 5, height = 3, dpi=300)

rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor, Resulting_corr)



######################## BETA DIVERSITY _ TIME SERIE R1 #######################

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
  perm_ASV<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_H<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_genus<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
pair_genus<- pairwise.adonis(sample_genus, factors=metadata$Experiment_state, p.adjust.m = "BH", sim.method="euclidean", perm = 9999)
pair_genus

# exporting beta diversity
suppressWarnings(rm(con))
con<-file("Results/Time_series_in_R1/Beta_div/Beta_divers_general_and_pairwise_between_Experiment_state_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity on Hellinger \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity on Hellinger (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_genus
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
{BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,metadata$Experiment_state)
  disp_genus<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_genus$tab
  a<-as.data.frame(disp_genus$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/Time_series_in_R1/Beta_div/Beta_dispersion_General_and_Pairwise_between_Experiment_state_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger \n")
disp_genus$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_genus,a, con)


####################### PCoA BETA DIV _ TIME SERIE R1 #########################

### on Genera
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) + stat_ellipse(size=0.2)  +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment_state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 6.5, height = 4.5, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment state", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_no_ellipses.png", width = 6.5, height = 4.5, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_line(aes(group=Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.2,"cm"), type = "closed") 
            ) +
  theme_classic(base_size = 8.5) +
  theme(title=element_text(size=8),
        axis.title.x = element_text(size=7.8),
        legend.margin = margin(2,-2,2,-7)
  ) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera) on reactor R1", color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation","\n            (the X axis corresponds with the passing of time, but it is NOT computed to be so)"),
       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
       )
ggsave(file="Results/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 5, height = 3.8, dpi=300)
# again, but smaller (poster)
sample_data(data.sqrt_prop)$Experiment_state<-gsub("_"," ",sample_data(data.sqrt_prop)$Experiment_state)
sample_data(data.sqrt_prop)$Experiment_state<-gsub("Almost s","~ S",sample_data(data.sqrt_prop)$Experiment_state, fixed = T)
sample_data(data.sqrt_prop)$Experiment_state<-factor(sample_data(data.sqrt_prop)$Experiment_state,
                                                     levels = c("Still preparing","Started","~ Stabilized"))
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values = c("Still preparing"="orange","Started"="coral2","~ Stabilized"="red2")) +
  theme_classic(base_size = 8.5) +
  theme(legend.text = element_text(size=6.5),
        legend.title = element_text(size=7.5),
        legend.margin = margin(2,-1,2,-8)
  )+
  geom_point(size=2.5) +
  geom_point(size=4.1, alpha= 0.4) +
  geom_line(aes(group=as.factor(Type)), col= "lightgray",
            size = 0.3, alpha=1,
            arrow=arrow(length =unit(0.28,"cm"), type = "closed") ) +
  labs(color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation"),
       y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Time_series_in_R1/Beta_div/PCoA_Beta_diversity_along_time_POSTER_VERSION.png", width = 4.3, height = 2.8, dpi=300)


#################### VEEN DIAGRAM _ TIME SERIE R1 ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.1% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Still_preparing<-subset_samples(data.venn, Experiment_state=="Still_preparing")
Still_preparing<-as.character(tax_table(prune_taxa(taxa_sums(Still_preparing)>0, Still_preparing))[,"Genus"])

Started<-subset_samples(data.venn, Experiment_state=="Started")
Started<-as.character(tax_table(prune_taxa(taxa_sums(Started)>0, Started))[,"Genus"])

Almost_stabilized<-subset_samples(data.venn, Experiment_state=="Almost_stabilized")
Almost_stabilized<-as.character(tax_table(prune_taxa(taxa_sums(Almost_stabilized)>0, Almost_stabilized))[,"Genus"])


ONLY_IN_Started<- Started[! Started %in% Almost_stabilized & ! Started %in% Still_preparing]
ONLY_IN_Started<- paste(ONLY_IN_Started, collapse = ", ")
head(ONLY_IN_Started)

ONLY_IN_Still_preparing<- Still_preparing[! Still_preparing %in% Almost_stabilized & ! Still_preparing %in% Started]
ONLY_IN_Still_preparing<- paste(ONLY_IN_Still_preparing, collapse = ", ")
head(ONLY_IN_Still_preparing)

ONLY_IN_Almost_stabilized<- Almost_stabilized[! Almost_stabilized %in% Still_preparing & ! Almost_stabilized %in% Started]
ONLY_IN_Almost_stabilized<- paste(ONLY_IN_Almost_stabilized, collapse = ", ")
head(ONLY_IN_Almost_stabilized)

IN_Still_preparing_AND_Almost_stabilized_BUT_NOT_Started <- Still_preparing[Still_preparing %in% Almost_stabilized & ! Still_preparing %in% Started ]
head(IN_Still_preparing_AND_Almost_stabilized_BUT_NOT_Started)


IN_Still_preparing_AND_Started_BUT_NOT_Almost_stabilized <- Still_preparing[! Still_preparing %in% Almost_stabilized & Still_preparing %in% Started ]
head(IN_Still_preparing_AND_Started_BUT_NOT_Almost_stabilized)

IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing <- Almost_stabilized[Almost_stabilized %in% Started & ! Almost_stabilized %in% Still_preparing ]
head(IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing)


x<-list(Still_preparing=Still_preparing,Started=Started,Almost_stabilized=Almost_stabilized)
ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram \n(only genera with minimal abundance > 0.1% \n at least in three samples) of R1 reactor",
       caption = paste("Not in 'Almost_stabilized': \n",
                       paste(IN_Still_preparing_AND_Started_BUT_NOT_Almost_stabilized,collapse = "   "),
                       "\n\n Not in 'Still preparing': \n",IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing)
  )
ggsave(filename = "Results/Time_series_in_R1/Beta_div/Venn_Diagramm.png", width = 4.85, height = 4.5, dpi=300, bg = "white")
dev.off()

suppressWarnings(rm(ONLY_IN_Started, ONLY_IN_Still_preparing, ONLY_IN_Almost_stabilized, 
                    x, con, Started, Still_preparing, Almost_stabilized, data.venn, who))


################ CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R1 (PROPORTION) #################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "PCR_ID"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances))) # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")

write.csv2(Corr_results, "Results/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1<-row.names(temp_for_save)


con <- file("Results/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R1 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


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
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most positely correlated with time (R1) \naccording to spearman correlation on % counts")
ggsave(file="Results/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
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
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most negatively correlated with time (R1)\n according to spearman correlation on % counts")
ggsave(file="Results/Time_series_in_R1/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_5, prune.dat_top, tax_selected, tabella, corr_6_top))


########## CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R1 (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "PCR_ID"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances))) # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")

write.csv2(Corr_results, "Results/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1<-row.names(temp_for_save)


con <- file("Results/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R1 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

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
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most positely correlated with time (R1)\n according to spearman correlation on CLR counts")
ggsave(file="Results/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
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
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most negatively correlated with time (R1)\n according to spearman correlation on CLR counts")
ggsave(file="Results/Time_series_in_R1/LOG_RATIO_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_5, prune.dat_top, tax_selected, tabella, corr_6_top))


####################### CORRELATIONS WITH REACTORS PARAMETERS _ R1 (CLR ABUNDANCES)  ####################

# loading the csv with reactor parameters
parameters<-as.data.frame(read.csv("Parametri_reattori_da_Serena.csv", header = T, na.strings = ""))
parameters_R1<-parameters[grepl("S", parameters$PCR_ID) & !grepl("A",parameters$PCR_ID), ]
parameters_R1<-parameters_R1[parameters_R1$PCR_ID %in% sample_names(data) , ] # removing missing sequencing
row.names(parameters_R1)<-parameters_R1$PCR_ID
parameters_R1<-parameters_R1[ , ! colnames(parameters_R1) %in% c("PCR_ID","Sampling_date")]
parameters_R1<-apply(parameters_R1, MARGIN=2, round, digits= 0) # avoiding digits to not inflate spearman ranks
parameters_R1[, "COD_in"]<- round(parameters_R1[, "COD_in"], digits = -1) # those column numbers are even more bigger --> they needs a further binning
parameters_R1<-t(parameters_R1)

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)


# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# avoiding uncultured and NA
abundances <- abundances[ ! grepl("uncultured", row.names(abundances)) & ! grepl("NA", row.names(abundances)) , ]

# same order between object and phyloseq data
abundances<-abundances[ , colnames(parameters_R1)]


# Computing the correlations...
{
  x<-as.matrix(t(rbind.data.frame(abundances,parameters_R1)))
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Abundant_bacteria","Parameter","Corr","pvalue")
}
corr<-correlation [ correlation$Abundant_bacteria %in% row.names(abundances), ]
corr<-corr [ corr$Parameter %in% row.names(parameters_R1), ]
colnames(corr)<-c("Abundant_bacteria","Parameter","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(y = corr$Abundant_bacteria, x = corr$Parameter, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=9.8) +
  theme(axis.text.x=element_text(angle = -18, hjust = 0, size= 8.5)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =7.5, vjust= 0.72) +
  labs(title = "", x= "R1 parameters", y= "R1 identified Genera (over 0.1% and in >50% samples)", 
       caption= "\n adjusted p-value (BH) lower than 0.01 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=8),
        legend.margin = margin(0,-2,0,-5),
        legend.key.height= unit(1.8, 'cm'),
        legend.key.width= unit(0.7, 'cm'), legend.title = element_text(size=10))
ggsave(file="Results/Time_series_in_R1/LOG_RATIO_counts_correlated_with_parameters/R1_genera_versus_parameters.png", dpi=300, width = 5, height = 5)
write.csv2(corr, file="Results/Time_series_in_R1/LOG_RATIO_counts_correlated_with_parameters/R1_genera_versus_parameters.csv", quote = F, na="", row.names = F)

# parameters_R1["OUR", ] # only three points ...



#######################   !!!!!  STARTING THE ANALYSIS OF TIME SERIES in R2 !!!!!   #############################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","data_complete","metadata_complete",
                                          "unfiltered_data","data_complete_with_lotti","metadata_complete_with_lotti",
                                          "Genera.DESEQ2_R1_vs_R2","significative_abundances_R1","significative_abundances_R2")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Type == "R2" )
data<-subset_samples(data, Experiment_state != "Adherent" )
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

sample_data(data)$Experiment_state<-factor(sample_data(data)$Experiment_state, levels=c("Still_preparing","Started","Almost_stabilized")) # decides the order in plots

days_ordered<-sample_data(data)$Sampling_date[order(sample_data(data)$Sampling_date)]
order_time<-days_ordered[c( grep("_03",days_ordered), grep("_04",days_ordered),grep("_05",days_ordered),grep("_06",days_ordered))] # months
sample_data(data)$Sampling_date<-factor(sample_data(data)$Sampling_date, levels = unique(order_time))

metadata<-metadata_complete[metadata_complete$Type=="R2" & metadata_complete$Experiment_state!="Adherent", ]


######################## ABUNDANCES BAR PLOT _ TIME SERIE R2 ##########################

# choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
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
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 11 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percentual abundance", title = "Five most abundant phyla (R2)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Time_series_in_R2/Abundances/TOP_5_phyla.png",width=7,height=4.5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Time_series_in_R2/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(top, prune.dat_top,tabella, tabella_top)


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
  tabella_top$Genus<-gsub("uncultured_ f Paludibacteraceae","uncultured_ f Paludibact.",tabella_top$Genus)
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 10.3 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percentual abundance", title = "Ten most abundant genera (R2)", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Time_series_in_R2/Abundances/TOP_10_Generi.png",width=7,height=4.5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/Time_series_in_R2/Abundances/TOP_10_genera_of_all_dataset_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))


# TOP 10 Genera of only "ALMOST STABILIZED" group
data.genus.subset<-subset_samples(data.genus.prop, Experiment_state=="Almost_stabilized")
{top <- names(sort(taxa_sums(data.genus.subset), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.genus.subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.genus.subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.genus.subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella_top$Genus<-gsub("uncultured_ f Paludibacteraceae","uncultured_ f Paludibact.",tabella_top$Genus)
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sampling_date<-factor(tabella$Sampling_date, levels(sample_data(data)$Sampling_date))
ggplot(data=tabella, aes(x=Sampling_date, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Experiment_state),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=13) + 
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=8.5),
        axis.text.x=element_text(angle=35, vjust=1, hjust=1, size=9),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 10.5 )) + 
  theme(legend.position="bottom",legend.margin = margin(5,0,5,0)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Samples", y="Percentual abundance", title = "Ten most abundant genera (R2)\n considering only the 'almost stabilized' samples", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Time_series_in_R2/Abundances/TOP_10_Generi_ONLY_STABILIZED_SAMPLES.png",width=7,height=4.8,dpi=300)
dev.off()



########### LINE PLOTS OF MOST ABUNDANT GENERA ##############


fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# TOP 5 Genera
suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.genus.prop), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.genus.prop)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =14) + 
  geom_point(aes(color=Genus), size =1.8) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13.5 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_continuous(breaks = seq(0, max(tabella$Abundance), 10)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "5 most abundant genera through time (R2)")
ggsave(file="Results/Time_series_in_R2/Abundances/TOP5_Genera_along_time.png",width=7,height=5, dpi=300) 
dev.off()


########################## ALFA DIVERSITY _ TIME SERIE R2 ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

dir.create("Results/Time_series_in_R2/Alpha_diversity/")
pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Experiment_state", color="Experiment_state")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
{ H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
  obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
  # adding evenness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Experiment_state<-factor(pAlpha$data$Experiment_state, levels = levels(sample_data(data)$Experiment_state))
pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Experiment_state, y=value, color=NULL), alpha=0) + theme_classic2() + 
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.2, alpha= 0.4) + 
  labs(x="Experiment state", title="Alpha diversity between time groups") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=28, vjust=1, hjust=1, size=9.8)) +
  stat_compare_means(aes(group = Experiment_state), label="p.format", method = "kruskal.test",
                     label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.4, hjust=0.18)
ggsave(file="Results/Time_series_in_R2/Alpha_diversity/Alfa diversity_on genera_Kruskal test.png", width = 6,height =5.2, dpi=300)


# Statistical analyses
alphadt<- as.data.frame(pAlpha$data)
Obser_value<-alphadt[alphadt$variable=="Observed richness", ]
factor<-Obser_value$Experiment_state
kruskal.test(Obser_value$value~factor) # just to test the plotted p value

# pairwise comparisons
con <- file("Results/Time_series_in_R2/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons_of_experiment_state.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Experiment_state sub groups \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Shannon", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-alphadt[alphadt$variable=="Evenness", ]
a<-FSA::dunnTest(value ~ Experiment_state, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)


# correlation through time (Observed r)
Filter_value<-alphadt[alphadt$variable=="Observed richness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=PCR_ID), size= 3.7) +
  geom_line(aes(group=Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Observed richness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 2))
  )
ggsave("Results/Time_series_in_R2/Alpha_diversity/Dot_plot_Observed_richness_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Shannon)
Filter_value<-alphadt[alphadt$variable=="Shannon", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=PCR_ID), size= 3.7) +
  geom_line(aes(group=Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Shannon",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/Time_series_in_R2/Alpha_diversity/Dot_plot_Shannon_through_time.png", width = 5, height = 3, dpi=300)

# correlation through time (Evenness)
Filter_value<-alphadt[alphadt$variable=="Evenness", ]
row.names(Filter_value)<-Filter_value$Sampling_date
Filter_value<- Filter_value[levels(sample_data(data)$Sampling_date), ] # ordered by sampling date ... the from 1o to last
Resulting_corr<-cor.test(Filter_value$value, 1:length(Filter_value$Sampling_date), method = "spearman")
ggplot( data=Filter_value, aes(y=Filter_value$value, x= 1:length(Filter_value$Sampling_date)) ) +
  geom_point(col="blue", size= 0.5,shape=15) +
  geom_text(aes(label=PCR_ID), size= 3.7) +
  geom_line(aes(group=Type), size=0.5, col="red", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x= "Time passing", y= "Evenness",
       caption = paste("spearman correlation through time:   rho", round(Resulting_corr$estimate,digits = 2), "  p-value", round(Resulting_corr$p.value,digits = 3))
  )
ggsave("Results/Time_series_in_R2/Alpha_diversity/Dot_plot_Evenness_through_time.png", width = 5, height = 3, dpi=300)

rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor, Resulting_corr)



######################## BETA DIVERSITY _ TIME SERIE R2 #######################

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
  perm_ASV<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_H<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Experiment_state, data=metadata, permutations = 9999, method="euclidean")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta

### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_genus<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
pair_genus<- pairwise.adonis(sample_genus, factors=metadata$Experiment_state, p.adjust.m = "BH", sim.method="euclidean", perm = 9999)
pair_genus

# exporting beta diversity
suppressWarnings(rm(con))
con<-file("Results/Time_series_in_R2/Beta_div/Beta_divers_general_and_pairwise_between_Experiment_state_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity on Hellinger \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity on Hellinger (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_genus
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)

# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on Genera
{BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,metadata$Experiment_state)
  disp_genus<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_genus$tab
  a<-as.data.frame(disp_genus$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/Time_series_in_R2/Beta_div/Beta_dispersion_General_and_Pairwise_between_Experiment_state_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger \n")
disp_genus$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_genus,a, con)


####################### PCoA BETA DIV _ TIME SERIE R2 #########################

### on Genera
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) + stat_ellipse(size=0.2)  +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment_state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Time_series_in_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 6.5, height = 4.5, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.8, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)", color="Experiment state", 
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Time_series_in_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_no_ellipses.png", width = 6.5, height = 4.5, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_line(aes(group=Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.2,"cm"), type = "closed") 
  ) +
  theme_classic(base_size = 8.5) +
  theme(title=element_text(size=8),
        axis.title.x = element_text(size=7.8),
        legend.margin = margin(2,-2,2,-7)
  ) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera) on reactor R2", color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation","\n            (the X axis corresponds with the passing of time, but it is NOT computed to be so)"),
       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/Time_series_in_R2/Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 5, height = 3.8, dpi=300)


#################### VEEN DIAGRAM _ TIME SERIE R2 ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.1% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Still_preparing<-subset_samples(data.venn, Experiment_state=="Still_preparing")
Still_preparing<-as.character(tax_table(prune_taxa(taxa_sums(Still_preparing)>0, Still_preparing))[,"Genus"])

Started<-subset_samples(data.venn, Experiment_state=="Started")
Started<-as.character(tax_table(prune_taxa(taxa_sums(Started)>0, Started))[,"Genus"])

Almost_stabilized<-subset_samples(data.venn, Experiment_state=="Almost_stabilized")
Almost_stabilized<-as.character(tax_table(prune_taxa(taxa_sums(Almost_stabilized)>0, Almost_stabilized))[,"Genus"])


ONLY_IN_Started<- Started[! Started %in% Almost_stabilized & ! Started %in% Still_preparing]
ONLY_IN_Started<- paste(ONLY_IN_Started, collapse = ", ")
head(ONLY_IN_Started)

ONLY_IN_Still_preparing<- Still_preparing[! Still_preparing %in% Almost_stabilized & ! Still_preparing %in% Started]
ONLY_IN_Still_preparing<- paste(ONLY_IN_Still_preparing, collapse = ", ")
head(ONLY_IN_Still_preparing)

ONLY_IN_Almost_stabilized<- Almost_stabilized[! Almost_stabilized %in% Still_preparing & ! Almost_stabilized %in% Started]
ONLY_IN_Almost_stabilized<- paste(ONLY_IN_Almost_stabilized, collapse = ", ")
head(ONLY_IN_Almost_stabilized)

IN_Still_preparing_AND_Almost_stabilized_BUT_NOT_Started <- Still_preparing[Still_preparing %in% Almost_stabilized & ! Still_preparing %in% Started ]
head(IN_Still_preparing_AND_Almost_stabilized_BUT_NOT_Started)


IN_Still_preparing_AND_Started_BUT_NOT_Almost_stabilized <- Still_preparing[! Still_preparing %in% Almost_stabilized & Still_preparing %in% Started ]
head(IN_Still_preparing_AND_Started_BUT_NOT_Almost_stabilized)

IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing <- Almost_stabilized[Almost_stabilized %in% Started & ! Almost_stabilized %in% Still_preparing ]
head(IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing)


x<-list(Still_preparing=Still_preparing,Started=Started,Almost_stabilized=Almost_stabilized)
ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram \n(only genera with minimal abundance > 0.1% \n at least in three samples) of R2 reactor",
       caption = paste("Not in 'Almost_stabilized': \n",
                       paste(IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing[1:3],collapse = "     "),"\n",
                       paste(IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing[4:6],collapse = "     "),"\n",
                       paste(IN_Almost_stabilized_AND_Started_BUT_NOT_Still_preparing[7:8],collapse = "     "),"\n",
                       "\n\n Not in 'Still preparing': \n",
                       paste(IN_Still_preparing_AND_Started_BUT_NOT_Almost_stabilized,collapse = "     ")
                      )
  )
ggsave(filename = "Results/Time_series_in_R2/Beta_div/Venn_Diagramm.png", width = 4.85, height = 4.5, dpi=300, bg = "white")
dev.off()

suppressWarnings(rm(ONLY_IN_Started, ONLY_IN_Still_preparing, ONLY_IN_Almost_stabilized, 
                    x, con, Started, Still_preparing, Almost_stabilized, data.venn, who))


################ CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R2 (PROPORTION) #################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "PCR_ID"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances))) # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")

write.csv2(Corr_results, "Results/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R2<-row.names(temp_for_save)


con <- file("Results/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R2 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!


# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[1:6]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most positely correlated with time (R2) \naccording to spearman correlation on % counts")
ggsave(file="Results/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[(length(significative_abundances_R2)-5):length(significative_abundances_R2)]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
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
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most negatively correlated with time (R2)\n according to spearman correlation on % counts")
ggsave(file="Results/Time_series_in_R2/PROPORTIONAL_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_5, prune.dat_top, tax_selected, tabella, corr_6_top))


########## CORRELATIONS ABUNDANCES vs TIME _ TIME SERIE R2 (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sampling_date
order_of_samples<-metadata[levels(sample_data(data)$Sampling_date), "PCR_ID"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]

Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances))) # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_holm<-p.adjust(Corr_results$pvalue, method = "holm")
Corr_results$sign<-ifelse(Corr_results$padj_holm<0.05,"*","")

write.csv2(Corr_results, "Results/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/Correlations_of_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R2<-row.names(temp_for_save)


con <- file("Results/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/Only_certain_genera_have_been_tested.txt")
sink(con, append = TRUE)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the R2 samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))


########### LINE PLOTS OF STRONGEST CORRELATIONS

fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[1:6]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella$Genus<-gsub("uncultured_","uncult. ", tabella$Genus)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12.6 ),
        legend.position="bottom", legend.margin = margin(1,0,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most positely correlated with time (R2)\n according to spearman correlation on CLR counts")
ggsave(file="Results/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/6_Genera_most_positively_correlated.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R2[(length(significative_abundances_R2)-5):length(significative_abundances_R2)]
suppressWarnings(rm(top, others, tabella))
{ tax_table(data.genus.prop)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(data.genus.prop, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella$Genus<-gsub("uncultured_","uncult. ", tabella$Genus)
  tabella$Genus<-gsub("Rhodobacteraceae","Rhodobacterac. ", tabella$Genus)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Sampling_date<-gsub("_23","",tabella$Sampling_date)
tabella<-tabella[order(tabella$Sampling_date) , ] # days
order_time<-tabella$Sampling_date[c(grep("_03",tabella$Sampling_date),grep("_04",tabella$Sampling_date),grep("_05",tabella$Sampling_date),grep("_06",tabella$Sampling_date))] # months
tabella$Sampling_date<-factor(tabella$Sampling_date, levels = unique(order_time))
ggplot(data=tabella, aes(y=Abundance, x=Sampling_date, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.6 ),
        legend.position="bottom",
        legend.margin = margin(1,25,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 5)) +
  labs(x="Time of sampling", y="Percentual abundance", title = "the 6 genera most negatively correlated with time (R2)\n according to spearman correlation on CLR counts")
ggsave(file="Results/Time_series_in_R2/LOG_RATIO_counts_correlated_with_time/6_Genera_most_negatively_correlated.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_5, prune.dat_top, tax_selected, tabella, corr_6_top))


####################### CORRELATIONS WITH REACTORS PARAMETERS _ R2 (CLR ABUNDANCES)  ####################

# loading the csv with reactor parameters
parameters<-as.data.frame(read.csv("Parametri_reattori_da_Serena.csv", header = T, na.strings = ""))
parameters_R2<-parameters[grepl("F", parameters$PCR_ID) & !grepl("A",parameters$PCR_ID), ]
parameters_R2<-parameters_R2[parameters_R2$PCR_ID %in% sample_names(data) , ] # removing missing sequencing
row.names(parameters_R2)<-parameters_R2$PCR_ID
parameters_R2<-parameters_R2[ , ! colnames(parameters_R2) %in% c("PCR_ID","Sampling_date")]
parameters_R2<-apply(parameters_R2, MARGIN=2, round, digits= 0) # avoiding digits to not inflate spearman ranks
parameters_R2[, "COD_in"]<- round(parameters_R2[, "COD_in"], digits = -1) # those column numbers are even more bigger --> they needs a further binning
parameters_R2<-t(parameters_R2)

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)


# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# avoiding uncultured and NA
abundances <- abundances[ ! grepl("uncultured", row.names(abundances)) & ! grepl("NA", row.names(abundances)) , ]

# same order between object and phyloseq data
abundances<-abundances[ , colnames(parameters_R2)]


# Computing the correlations...
{
  x<-as.matrix(t(rbind.data.frame(abundances,parameters_R2)))
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Abundant_bacteria","Parameter","Corr","pvalue")
}
corr<-correlation [ correlation$Abundant_bacteria %in% row.names(abundances), ]
corr<-corr [ corr$Parameter %in% row.names(parameters_R2), ]
colnames(corr)<-c("Abundant_bacteria","Parameter","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(y = corr$Abundant_bacteria, x = corr$Parameter, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=9.8) +
  theme(axis.text.x=element_text(angle = -18, hjust = 0, size= 8.5)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =7.5, vjust= 0.72) +
  labs(title = "", x= "R2 parameters", y= "R2 identified Genera (over 0.1% and in >50% samples)", 
       caption= "\n adjusted p-value (BH) lower than 0.01 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=8),
        legend.margin = margin(0,-2,0,-5),
        legend.key.height= unit(1.8, 'cm'),
        legend.key.width= unit(0.7, 'cm'), legend.title = element_text(size=10))
ggsave(file="Results/Time_series_in_R2/LOG_RATIO_counts_correlated_with_parameters/R2_genera_versus_parameters.png", dpi=300, width = 5, height = 5)
write.csv2(corr, file="Results/Time_series_in_R2/LOG_RATIO_counts_correlated_with_parameters/R2_genera_versus_parameters.csv", quote = F, na="", row.names = F)


################### DIFFERENTIAL ABUNDANCES WITH DESEQ2 _ Adherent vs Suspension (R1) #####################

if(! "unfiltered_data" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","data_complete","metadata_complete",
                                          "unfiltered_data","data_complete_with_lotti","metadata_complete_with_lotti",
                                          "Genera.DESEQ2_R2_vs_R2","significative_abundances_R2","significative_abundances_R2")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Type == "R1" )
data<-subset_samples(data, Sampling_date %in% c("06_06_23","08_06_23","09_06_23",
                                                "22_05_23","22_05_23",
                                                "15_06_23","13_06_23","19_06_23")
) # only the data around the sampling of adherent biomass (only one in pair) have been selected

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
sample_data(data_pruned)$Sampling_area<-ifelse(sample_data(data_pruned)$Experiment_state=="Adherent","Adherent","Suspension")

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sampling_area)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sampling_area", "Suspension", "Adherent"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Suspension_vs_Adherent.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Adherent_vs_Suspended_R1_differential_abundances/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/Adherent_vs_Suspended_R1_differential_abundances/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Sampling_area)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Suspension"="chartreuse","Adherent"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", 
       fill="Sampling area", x="")
#ggsave(filename = "Results/DA_DESeq2/DA_Sampling_area_every_result.png", width = 13, height = 7, dpi=300)
#dev.off()

Redund<- c("Fusibacteraceae" ,"Desulfovibrionales","Verrucomicrobiota","Acholeplasmatales")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
Table_tot2<-subset(Table_tot2, ! Taxa %in% c("Class","Phylum")) # to remove redundant results
Table_tot2$Bacteria<-gsub("Peptostreptococcales-Tissierellales","Peptostreptococ.-T", Table_tot2$Bacteria)
ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Sampling_area) ) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8, alpha=0.25, size=0.22) +
  geom_point(aes(color=Sampling_area), position = position_dodge(width=0.8), size=0.7, alpha=0.7) +
  theme_bw(base_size = 8.5) +
  scale_fill_manual(values=c("Suspension"="chartreuse","Adherent"="coral")) +
  scale_color_manual(values=c("Suspension"="chartreuse4","Adherent"="coral")) +
  theme(strip.text.x=element_text(size=8,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=8),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=7.5), 
        axis.text.y = element_text(size=5.8), 
        plot.title= element_text(size=10),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+3),2))) +
  labs(title= "Differently abundant Taxa (R1)", y="Proportional Abundance %", 
       fill="Sampling area", x="") +
  guides(color="none")
ggsave(filename = "Results/Adherent_vs_Suspended_R1_differential_abundances/DA_Sampling_area_no_redundants.png", width = 5.5, height = 3.5, dpi=300)
dev.off()

system(" echo 'Only the sample with the sampling data around the one of the adherent biomass have been selected (only one adherent sample is in pair!). Moreover, every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Adherent_vs_Suspended_R1_differential_abundances/NB.txt ")


################### DIFFERENTIAL ABUNDANCES WITH DESEQ2 _ Adherent vs Suspension (R2) #####################

if(! "unfiltered_data" %in% ls() ){
  stop("\nDid you perform the filtering step yet?\n", fill = T)
  Sys.sleep(2)
}

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","data_complete","metadata_complete",
                                          "unfiltered_data","data_complete_with_lotti","metadata_complete_with_lotti",
                                          "Genera.DESEQ2_R2_vs_R2","significative_abundances_R2","significative_abundances_R2")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Type == "R2" )
data<-subset_samples(data, Sampling_date %in% c("06_06_23","08_06_23","09_06_23",
                                                "22_05_23","22_05_23",
                                                "15_06_23","13_06_23","19_06_23")
) # only the data surrounding the sampling of adherent biomass (only one in pair) have been selected

##### STARTING THE DIFFERENTIAL ANALYSIS
suppressWarnings(rm(data_pruned, data.genus_pruned))
data_pruned<- prune_taxa(taxa_sums(data) > 10, data) 
# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)
sample_data(data_pruned)$Sampling_area<-ifelse(sample_data(data_pruned)$Experiment_state=="Adherent","Adherent","Suspension")

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
  DEseq_data<-phyloseq_to_deseq2(d, ~Sampling_area)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Sampling_area", "Suspension", "Adherent"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.01) & abs(res$log2FoldChange)>1,]
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
    # write.csv2(r, file=paste0("Results/DA_DESeq2/DA_",t,"_ratio_Suspension_vs_Adherent.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Adherent_vs_Suspended_R2_differential_abundances/Every_result_DESeq2.csv", row.names = F)
write.xlsx(Res_tot, file="Results/Adherent_vs_Suspended_R2_differential_abundances/Every_result_DESeq2.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Sampling_area)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Suspension"="chartreuse","Adherent"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
  #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5))) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", 
       fill="Sampling_area", x="")
#ggsave(filename = "Results/DA_DESeq2/DA_Sampling_area_every_result.png", width = 13, height = 7, dpi=300)
#dev.off()

Redund<- c("Gammaproteobacteria","Verrucomicrobiota")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
Table_tot2<-subset(Table_tot2, ! Taxa %in% c("Family","Order")) # to remove redundant results
ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Sampling_area) ) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8, alpha=0.25, size=0.22) +
  geom_point(aes(color=Sampling_area), position = position_dodge(width=0.8), size=0.7, alpha=0.7) +
  theme_bw(base_size = 8.5) +
  scale_fill_manual(values=c("Suspension"="chartreuse","Adherent"="coral")) +
  scale_color_manual(values=c("Suspension"="chartreuse4","Adherent"="coral")) +
  theme(strip.text.x=element_text(size=8,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-15, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=8),
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=7.5), 
        axis.text.y = element_text(size=5.8), 
        plot.title= element_text(size=10),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,18,4),3.5,24,30,seq(36,max(Table_tot$Abundance+3),10))) +
  labs(title= "Differently abundant Taxa (R2)", y="Proportional Abundance %", 
       fill="Sampling area", x="") +
  guides(color="none")
ggsave(filename = "Results/Adherent_vs_Suspended_R2_differential_abundances/DA_Sampling_area_no_redundants.png", width = 4.8, height = 3.5, dpi=300)
dev.off()

system(" echo 'Only the sample with the sampling data around the one of the adherent biomass have been selected (only one adherent sample is in pair!). Moreover, every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Adherent_vs_Suspended_R2_differential_abundances/NB.txt ")


############### !!!!! EXTRA : OLD REACTOR VS NEW REACTOR PCoA !!!!! #################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","data_complete","metadata_complete",
                                          "unfiltered_data","data_complete_with_lotti","metadata_complete_with_lotti",
                                          "Genera.DESEQ2_R1_vs_R2","significative_abundances_R1","significative_abundances_R2")]
  rm(list=to_remove)
}

{
old_data<-load("../../Primo sequenziamento DNA_aero VS anaero_e_conservaz/Analisi_dati_1/data_backup_after_filters.RData")
# sample_data(data) # this is now the old object!
data_old<-subset_samples(data, Owner=="Serena")
sample_names(data_old)<-gsub("S","N_S",sample_names(data_old))
sample_data(data_old)$Experiment_state<-rep("N_limitation (old R)",length(sample_names(data_old)))
# recreating the same object but without tree
data_old<- phyloseq(otu_table(otu_table(data_old)), tax_table(tax_table(data_old)) , sample_data(sample_data(data_old)) )
rm(data)
}

data_merged<-merge_phyloseq(data_complete, data_old)
data_merged<-tax_glom(data_merged, taxrank = "Genus")
data_merged<-transform_sample_counts(data_merged, function(x) (x/sum(x)) * 100)
# sample_data(data_merged)$Experiment_state   # OK

### on Genera
data.prop.labels<-data_merged
data.prop.labels<-subset_samples(data.prop.labels, Experiment_state!="Adherent")
sample_data(data.prop.labels)$Experiment_state<-factor(sample_data(data.prop.labels)$Experiment_state,
                                                       levels = c("N_limitation (old R)","Still_preparing","Started","Almost_stabilized"))
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_state") +
  scale_color_manual(values=c("N_limitation (old R)"="lightblue","Still_preparing"="orange","Started"="coral2","Almost_stabilized"="red")) +
  geom_point(size=3.4, alpha=0.4) +
  geom_point(size=1, alpha=1) +
  theme_classic(base_size = 11) + stat_ellipse(size=0.2)  +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera) \n with old reactor (N lim) and the new one (P lim)", color="Experiment state",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/EXTRA_OLD_reactor_vs_NEW.png", width = 6.5, height = 4.5, dpi=300)


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
print("Dendextend")
packageVersion("dendextend")
cat("\n", "\n", fill=TRUE)
print("egg")
packageVersion("egg")
cat("\n", "\n", fill=TRUE)
print("ggh4x")
packageVersion("ggh4x")
cat("\n", "\n", fill=TRUE)
package$otherPkgs$pca3d[c(1,4)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$mixOmics[c(1,4)]
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink() # restore STR OUTPUT to R console
close(con)
suppressWarnings(rm(con))
