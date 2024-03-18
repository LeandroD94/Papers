##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggh4x")
  library("egg")
  # analysis packages
  library("vegan")
  library("ggpubr")
  # utilities
  library("xlsx")  
  library("qiime2R")
}

{dir.create("Data_check")
dir.create("Results")
dir.create("Results/General_data_analysis")
dir.create("Results/General_data_analysis/Alpha_divers_EVERY_SAMPLE")
dir.create("Results/General_data_analysis/Alpha_divers_with_Nlim_and_Plim")
dir.create("Results/General_data_analysis/PCoA_EVERY_SAMPLE")
dir.create("Results/General_data_analysis/PCoA_with_Nlim_and_Plim")
}
{
dir.create("Results/N_lim_reactor")
dir.create("Results/N_lim_reactor/Abundances")
dir.create("Results/N_lim_reactor/Inoculum_comparison_with_CuoioDep")
dir.create("Results/N_lim_reactor/PHA_Accumulators")
dir.create("Results/N_lim_reactor/Time_series_in_N_lim")
}
dir.create("Results/P_lim_reactor")
for(f in c("Time_series_in_R1","Time_series_in_R2","R1_versus_R2")){
  dir.create(paste0("Results/P_lim_reactor/",f))
  dir.create(paste0("Results/P_lim_reactor/",f,"/Abundances"))
  dir.create(paste0("Results/P_lim_reactor/",f,"/Beta_div"))
  if(f %in% c("Time_series_in_R1","Time_series_in_R2")){
    dir.create(paste0("Results/P_lim_reactor/",f,"/PROPORTIONAL_counts_correlated_with_time"))
    dir.create(paste0("Results/P_lim_reactor/",f,"/LOG_RATIO_counts_correlated_with_time"))
  }
}
dir.create("Results/P_lim_reactor/R1_versus_R2/DA_DESeq2")


options(scipen = 100) # disable scientific annotation


colors<-c("CuoioDepur Aero"="deepskyblue",   # groups color in PCoA plots
          "N limiting"="lightblue3",
          "P limiting"="coral"
)
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one



####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^14.*F","",sample)
sample<-gsub("^15.*F","",sample)
sample_names(data)<-sample # update

Metadata <- as.data.frame(read.csv("metadata.csv", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)]) - 1 # there is a sample (F3) collected but NOT sequenced
original_order<-Metadata$Sample_name # to maintein this sample order in plots
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

# adding a space to some day (technical replicate) to differentiate the labels in the plot
Metadata$Experiment_day[Metadata$Sample_name=="P2"]<-paste0(Metadata$Experiment_day[Metadata$Sample_name=="P2"]," ")
Metadata$Experiment_day[Metadata$Sample_name=="P7"]<-paste0(Metadata$Experiment_day[Metadata$Sample_name=="P7"]," ")

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length)

sample_names(data)<-sample_data(data)$Sample_name
sample_data(data)$Sample_name<-factor(sample_data(data)$Sample_name,
                                      levels = original_order)
sample_data(data)$Experiment_day<-factor(sample_data(data)$Experiment_day,
                                         levels = unique(sample_data(data)$Experiment_day))



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

proof1<-"This is the proof of having performed the filtering steps"



############ ANNOTATING THE READS NUMBER BEFORE AND AFTER THE PROCESSING ##################

if(! "proof1" %in% ls()){
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
DADA2_read_number<-as.data.frame(read.table(file="QIIME/DADA2_Initial_and_processed_reads.tsv", header = T, sep="\t", row.names = 1))
after_filter_number<-as.data.frame(colSums(otu_table(data)))
# updating names from FASTQ codes to sample names
Original_read_number<-Original_read_number[original_names, ] # this vector came from the import section (ordered as "sample", which contains the modified names)
DADA2_read_number<-DADA2_read_number[original_names, ]
row.names(Original_read_number)<-sample # same order, already modified
row.names(DADA2_read_number)<-sample # same order, already modified
row.names(after_filter_number)<-sample # same order, already modified
# creating the table
table<-cbind.data.frame(Original=Original_read_number$forward.sequence.count, # merged --> just one column, not both of them
                        After_quality_filter=DADA2_read_number$non.chimeric,
                        After_contam_filter=after_filter_number[,1])
table$FASTQ<-row.names(Original_read_number)
# re-naming and adding groups
if(identical( table$FASTQ , sample_data(data)$FASTQ_ID) ){ # same order
  table$Samples<-sample_data(data)$Sample_name
  table$Experiment_day<-sample_data(data)$Experiment_day
  table$factor <- sample_data(data)$Reactor_Type
  table$factor <- gsub("CuoioD_Aero","CuoioDepur\nAero tank",table$factor)
  table$factor <- gsub("N_limitation","N limiting",table$factor)
  table$factor <- gsub("P_limitation","P limiting",table$factor)
  table$factor <- as.factor(table$factor)
  cat("\n\nOK!\n\n")
}
# plotting
ggplot(aes(x=Samples, y=Original), data=table) +
  facet_grid( ~ factor, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Read quality filters"), alpha= 0.8,
            width = 0.6, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
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
                    breaks=c('Original number', 'Read quality filters', 'Relative abundance filters'),
                    values=c('Original number'='green3', 'Read quality filters'='coral', 'Relative abundance filters'='red3')) +
  scale_y_continuous(breaks = c(10000, seq(0, max(table$Original), 25000))) +
  labs(y="Reads abundance", x="FASTQ")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 10, height = 3.8, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)


### again, without the other 4 inoculum samples and the adherent ones
table2<-table[table$factor %in% c("N limiting","P limiting") & !(table$Samples%in%c("AS1","AS2","AS3","AF1","AF2","AF3")), ]
ggplot(aes(x=Samples, y=Original), data=table2) +
  facet_grid( ~ factor, space="free_x", scales = "free_x") +
  geom_bar( aes(fill="Original number") ,  width=0.9,  stat = "identity", alpha= 0.5) +
  geom_bar( aes(y= After_quality_filter, fill="Read quality filters"), alpha= 0.8,
            width = 0.6, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
            width= 0.2, stat = "identity") +
  theme_classic( base_size = 10.2) +
  theme(axis.text.x = element_text(size=7, angle = 50, vjust=1, hjust=1),
        axis.text.y = element_text(size=5),
        panel.grid.major.y = element_line(size=0.1, color="grey"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size=9.2),
        legend.title = element_text(size=9.8),
        legend.key.height = unit(0.4,"cm")
  ) +
  scale_fill_manual(name='Number of reads:  ',
                    breaks=c('Original number', 'Read quality filters', 'Relative abundance filters'),
                    values=c('Original number'='green3', 'Read quality filters'='coral', 'Relative abundance filters'='red3')) +
  scale_y_continuous(breaks = c(10000, seq(0, max(table$Original)+25000, 25000))) +
  labs(y="Reads abundance", x="FASTQ")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters_without_CuoioD_or_adherent.png", width = 7.2, height = 3.8, dpi=300)

write.csv2(table2, file="Data_check/Number_of_reads_pre_and_post_filters_without_CuoioD_or_adherent.csv", quote=F, row.names = F)


rm(table, table2, Original_read_number, DADA2_read_number, after_filter_number)



#################### FURTHER CHECKING THE HOST CONTAMINATION PROPORTION ###################

# ### if simil Bayes classification (unclassified are NA phyla)
data_contam_k<- tax_glom(data, taxrank = "Phylum", NArm = F) # without Bowtie
data_contam_k<- transform_sample_counts(data_contam_k, function(x) (x/sum(x))*100)
data_contam_k<- prune_taxa(taxa_sums(data_contam_k)>0, data_contam_k)
data_contam_k<- subset_taxa(data_contam_k, is.na(tax_table(data_contam_k)[,"Phylum"]) )
write.csv2(file="Data_check/Host_Contamin_proportion_in_raw_FASTQ.csv", cbind.data.frame(tax_table(data_contam_k)[,c("Kingdom","Phylum")],otu_table(data_contam_k) ) )

rm(data_contam_k)

{
  data<-subset_taxa(data, ! is.na(tax_table(data)[,"Phylum"]) )
  proof2<-"This is the proof of removing NA phyla (probable contaminants and NOT bacteria, considering the classificator nature)"
}


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


# genus number without CuoioD or adherent samples (published ones only)
no_CuoioD <- subset_samples(data, ! Sample_name %in% c("A1","A2","A3","A4","AS1","AS2","AS3","AF1","AF2","AF3") )
data.genus_nc = tax_glom(no_CuoioD, taxrank = "Genus", NArm = F)
Taxa.genus_nc<-as.data.frame(tax_table(data.genus_nc))
assigned<-cbind.data.frame(length(Taxa.genus_nc$Genus),length(which(!is.na(Taxa.genus_nc$Genus))),length(which(!is.na(Taxa.genus_nc$Genus)))/length(Taxa.genus_nc$Genus),"Genus")
colnames(assigned)<-c("Total","Assigned","%","Taxa")
assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters_noCuoioD_noAdherents.csv",row.names = F, quote = F)
rm(a,e,assigned,data.genus_nc, data.phy_nc, Taxa.genus_nc, Taxa.phy_nc)



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


# removing the other inoculum clones and the "adherent biomass"
no_CuoioD <- subset_samples(data, ! Sample_name %in% c("A1","A2","A3","A4","AS1","AS2","AS3","AF1","AF2","AF3") )
png(file="Data_check/Rarefaction_curve_noCuoioD.png",width=3000,height=2100, res=300)
#r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
r<-rarecurve(t(as(otu_table(no_CuoioD),"matrix")), step=100,label=F)
evalslopes(r,sample_names(no_CuoioD),lim=0.001,cex=1)
dev.off()
rm(r)



### REMOVING THE UNSATURATED SAMPLES 

if("S13" %in% Metadata$Sample_name){
data<-subset_samples(data, Sample_name!="S13")
Metadata<-Metadata[Metadata$Sample_name!="S13", ]
proof3<-"The unsaturated sample (S13) has been removed"
}



############### \\\\ BACKUP OF THE ORIGINAL (COMPLETE) DATA FOR SUBSETES \\\\ #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls() | ! "proof3" %in% ls()){
  stop("\nDid you perform the filtering steps yet?\n")
}

if( length(Metadata$Experiment_state)>3)  { # then still complete
  data_complete<-prune_taxa(taxa_sums(data)>0,data)
  metadata_complete<-Metadata
  rm(data)
}


####################### GENERAL PCoA BEFORE SUBSETTING (EVERY SAMPLE) ##########################

suppressWarnings(rm(data.prop_temp, data.sqrt_prop))

if( ! "proof1" %in% ls() | ! "proof2" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}
data.temp<-subset_samples(data_complete, Experiment_state!="Adherent")
data.prop_temp<-transform_sample_counts(data.temp, function(x) (x/sum(x))*100)
data.prop.labels<-data.prop_temp

sample_data(data.prop.labels)$Reactor_Type<-gsub("CuoioD_Aero","CuoioDepur Aero",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Reactor_Type<-gsub("N_limitation","N limiting",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Reactor_Type<-gsub("P_limitation","P limiting",sample_data(data.prop.labels)$Reactor_Type)

{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}

sample_names(data.sqrt_prop)<-as.factor(sample_names(data.sqrt_prop))

# Type 1 (everything)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.2) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.12) +
  stat_ellipse(size=0.15) + 
  geom_text(aes(label= Sample_name), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_EVERY_SAMPLE/GENERAL_PCoA_Hellinger_on_genera_Ellypses_Names_Lines.png", width = 6, height = 4.5, dpi=300)
# Type 2 (no ellypses)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.8) +
  geom_point(size=4.8, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.12) +
  geom_text(aes(label= Sample_name), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_EVERY_SAMPLE/GENERAL_PCoA_Hellinger_on_genera_Names_Lines.png", width = 6, height = 4.5, dpi=300)
# Type 3 (no ellypses, no names)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=1.8) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.12) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_EVERY_SAMPLE/GENERAL_PCoA_Hellinger_on_genera_NO_Names.png", width = 6, height = 4.5, dpi=300)
# Type 4 (no ellypses, no lines)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=3) +
  geom_point(size=4.8, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_text(aes(label= Sample_name), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_EVERY_SAMPLE/GENERAL_PCoA_Hellinger_on_genera_NO_Lines_NO_Ellypses.png", width = 6, height = 4.5, dpi=300)
# Type 5 (no lines)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.2) +
  geom_point(size=4, alpha= 0.5) +
  theme_classic(base_size = 9) +
  stat_ellipse(size=0.15) + 
  geom_text(aes(label= Sample_name), 
            color="black", size=2,
            show.legend = FALSE) +
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_EVERY_SAMPLE/GENERAL_PCoA_Hellinger_on_genera_NO_Lines.png", width = 6, height = 4.5, dpi=300)
# Type 6 (ONLY Ellypses)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.5) +
  geom_point(size=4.5, alpha= 0.5) +
  theme_classic(base_size = 9) +
  stat_ellipse(size=0.15) + 
  labs(title="",
       color="Sample Type",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/General_data_analysis/PCoA_EVERY_SAMPLE/GENERAL_PCoA_Hellinger_on_genera_ONLY_Ellypses.png", width = 6, height = 4.5, dpi=300)




############################# GENERAL COUNTS EXPORT ##############################

data<-data_complete
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}
{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

dir.create("Results/General_data_analysis/General_Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/General_data_analysis/General_Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/General_data_analysis/General_Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}



### again but only with the samples used in the publication
no_CuoioD <- subset_samples(data_complete, ! Sample_name %in% c("A1","A2","A3","A4","AS1","AS2","AS3","AF1","AF2","AF3") )
data.genus_nc = tax_glom(no_CuoioD, taxrank = "Genus", NArm = F)
write.csv2(cbind(as(otu_table(no_CuoioD),"matrix"),as(tax_table(no_CuoioD),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_otu_noCuoio_noAdherents.csv",quote=F)
write.csv2(cbind(as(otu_table(data.genus_nc),"matrix"),as(tax_table(data.genus_nc),"matrix")),file="Results/General_data_analysis/General_Raw_counts/counts_genus_noCuoio_noAdherents.csv",quote=F)
rm(data.genus_nc)
  
  

################## ABUNDANCES BAR PLOTS (EVERY SAMPLE) ##################


### TOP 5 Phyla
data.target<-tax_glom(data_complete, taxrank = "Phylum")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.target)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("CuoioD_Aero","CuoioD",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("P_limitation","P limiting",tabella$Reactor_Type)
tabella$Reactor_Type[tabella$Sample_Type=="R2"]<-"P limiting R2"
tabella$Reactor_Type[tabella$Sample_Type=="R1"]<-"P limiting R1"
levels(tabella$Sample_name)

ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.2), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12 ),
        legend.position="bottom",
        #legend.margin = margin(0,0,0,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Sample", y="Percent abundance", 
       title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/General_data_analysis/TOP_5_phyla_EVERY_SAMPLE_CuoioDep_Nlim_Plim.png",width=8.6,height=5, dpi=300)
dev.off()

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)


### TOP 10 Genera
data.target<-tax_glom(data_complete, taxrank = "Genus")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("CuoioD_Aero","CuoioD",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("P_limitation","P limiting",tabella$Reactor_Type)
tabella$Reactor_Type[tabella$Sample_Type=="R2"]<-"P limiting R2"
tabella$Reactor_Type[tabella$Sample_Type=="R1"]<-"P limiting R1"

ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.2), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.8 ),
        legend.position="bottom", 
        legend.margin = margin(0,0,0,0)
  ) + 
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="Sample", y="Percent abundance",
       title = "Ten most abundant genera",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/General_data_analysis/TOP_10_Genera_EVERY_SAMPLE_CuoioDep_Nlim_Plim.png",width=8.6,height=5,dpi=300)
dev.off()



################## ABUNDANCES BAR PLOTS (ONLY INOCULUM) ##################


### TOP 5 Phyla Inoculum
data.target<-subset_samples(data_complete, Experiment_day %in% c("1", "1 ") )
data.target<-tax_glom(data.target, taxrank = "Phylum")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

suppressWarnings(rm(top, others, tabella))
{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:5]
  prune.dat_top <- prune_taxa(top,data.target)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Phylum<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Phylum<-factor(tabella$Phylum, levels = c(unique(tabella$Phylum)[unique(tabella$Phylum)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("CuoioD_Aero","CuoioD",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)
#levels(tabella$Sample_name)

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.2), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12 ),
        legend.position="bottom",
        #legend.margin = margin(0,0,0,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment_day", y="Percent abundance", 
       title = "Five most abundant phyla (only inoculum)", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/General_data_analysis/TOP_5_phyla_INOCULUM.png",width=8.6,height=5, dpi=300)
dev.off()

write.xlsx(file = "Results/General_data_analysis/TOP_5_phyla_INOCULUM_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)



### TOP 10 Genera Inoculum
data.target<-subset_samples(data_complete, Experiment_day %in% c("1", "1 ") )
data.target<-tax_glom(data.target, taxrank = "Genus")
data.target<-transform_sample_counts(data.target, function(x) (x/sum(x))*100 )

{top <- names(sort(taxa_sums(data.target), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.target)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.target)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.2), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.8 ),
        legend.position="bottom", 
        legend.margin = margin(0,0,0,0)
  ) + 
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="Experiment_day", y="Percent abundance",
       title = "Ten most abundant genera (only inoculum)",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/General_data_analysis/TOP_10_Genera_INOCULUM.png",width=8.6,height=5,dpi=300)
dev.off()

write.xlsx(file = "Results/General_data_analysis/TOP_10_Genera_INOCULUM_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))



################### ALPHA DIVERSITIES VALUES (OF EVERY SAMPLE) #####################

# selecting only the last 4 sample for each reactor (the one that SHOULD be at least near the steady state)
data.genus.temp<-subset_samples(data_complete, Sample_name %in% c("A2","A3","A4","P1", # not using P2 being a technical replicate
                                                                  "P9","P10","P11","P12",
                                                                  "S19","S20","S21","S22",
                                                                  "F19","F20","F21","F22"))
data.genus.temp<-tax_glom(data.genus.temp, taxrank ="Genus", NArm=F)
pAlpha<-plot_richness(data.genus.temp, measures=c("Shannon", "Observed"), x="Reactor_Type", color="Reactor_Type")
# plot_richness( ) compute diversity like estimate_diversity( )
H <- pAlpha$data[pAlpha$data$variable=="Shannon", ]
obs <- pAlpha$data[pAlpha$data$variable=="Observed", ]
# adding evenness
{identical(H$Sampling_date, obs$Sampling_date) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Sample_Type, New_data$Sampling_date),]
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha$data$Reactor_Type<-gsub("CuoioD_Aero","CuoioDepur",pAlpha$data$Reactor_Type)
pAlpha$data$Reactor_Type<-gsub("N_limitation","N limiting",pAlpha$data$Reactor_Type)
pAlpha$data$Reactor_Type<-gsub("P_limitation","P limiting",pAlpha$data$Reactor_Type)
pAlpha$data$Reactor_Type[pAlpha$data$Sample_Type=="R2"]<-"P limiting R2"
pAlpha$data$Reactor_Type[pAlpha$data$Sample_Type=="R1"]<-"P limiting R1"

pAlpha$data$Reactor_Type[pAlpha$data$Sample_name=="P1"]<-"CuoioDepur" # this is the starting inoculum, then it is treated as one of CuoioD samples

colors4<-c("CuoioDepur"="deepskyblue",
           "N limiting"="lightblue3",
           "P limiting R1"="coral",
           "P limiting R2"="red"
)

# with points only
pAlpha +
  theme_classic2() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL), size= 0.4, alpha= 0) +
  geom_point(size=3.1, alpha= 0.4) +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "kruskal.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top",
                     vjust=-0.4, hjust=-0.05) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_EVERY_SAMPLE/Alfa_diversity_EVERY_SAMPLE_with_Kruskal.png", width = 7,height =5.2, dpi=300)
# with ID
pAlpha + theme_classic2() +
  scale_color_manual(values = colors4 ) +
  geom_boxplot(data=pAlpha$data, aes(x=Reactor_Type, y=value, color=NULL), size= 0.4, alpha= 0) +
  geom_text(aes(label=Sample_name), size= 2.5, color= "black") +
  labs(x="", title="",
       caption = "Only the last four samples for each reactor have been included") +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=20, vjust=1, hjust=1, size=9)) +
  stat_compare_means(aes(group = Reactor_Type), label="p.format", method = "kruskal.test", paired = T,
                     label.x= 1.5, size=3.5, label.y.npc = "top",
                     vjust=-0.4, hjust=-0.05) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/General_data_analysis/Alpha_divers_EVERY_SAMPLE/Alfa_diversity_EVERY_SAMPLE_with_Kruskal_and_sampleID.png", width = 7,height =5.2, dpi=300)


######## pairwise comparisons
alphadt<-pAlpha$data
Filter_value_obs<-alphadt[alphadt$variable=="Observed richness", ]
Filter_value_sha<-alphadt[alphadt$variable=="Shannon", ]
Filter_value_eve<-alphadt [alphadt$variable=="Evenness", ]

# with Dunnett (one factor only)
con <- file("Results/General_data_analysis/Alpha_divers_EVERY_SAMPLE/Alfa_diversity_groups_pairwise_comparisons.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on sub groups Alpha diversity  \n P-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_obs, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_sha, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
a<-FSA::dunnTest(value ~ Reactor_Type, data=Filter_value_eve, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a$res$Sign<-a$res$P.adj
a$res$Sign[a$res$P.adj<=0.05]<-"*"
a$res$Sign[a$res$P.adj>0.05]<-""
a$res
sink()
close(con)



############## SAVING FOR THE 2nd PART OF THE ANALYSIS #################

save(list=c("data_complete","metadata_complete","unfiltered_data","proof1","proof2","proof3"),
     file="Data_prepared_by_1st_part_of_the_script.RData"
)


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
print("egg")
packageVersion("egg")
cat("\n", "\n", fill=TRUE)
print("ggh4x")
packageVersion("ggh4x")
cat("\n", "\n", fill=TRUE)
print("ggvenn")
packageVersion("ggvenn")
cat("\n", "\n", fill=TRUE)
print("vegan")
packageVersion("vegan")
cat("\n", "\n", fill=TRUE)
print("microbiome")
packageVersion("microbiome")
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
