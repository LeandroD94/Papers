##################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  library("microbiome")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  # analysis packages
  library("DESeq2")
  library("vegan")
  # utilities
  library("xlsx")  
  library("Hmisc")
  library("qiime2R")
}

{dir.create("Data_check")
dir.create("Results")
dir.create("Results/N_and_P")
dir.create("Results/N_lim_reactor")
dir.create("Results/N_lim_reactor/Abundances")
dir.create("Results/N_lim_reactor/Spearman_CLR_in_N_lim")
dir.create("Results/P_lim_reactor")
dir.create("Results/P_lim_reactor/Abundances")
dir.create("Results/P_lim_reactor/Spearman_CLR_in_P_lim")
}
options(scipen = 100) # disable scientific annotation


colors<-c("CuoioDepur Aero"="deepskyblue",   # groups color in PCoA plots
          "N limiting"="lightblue3",
          "P limiting"="coral"
)
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one


table_PHA<-read.table(file="PHA_Accumulators_19marzo2024.tsv", fill = T, header = T, sep="\t")



####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")

# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^14.*F","",sample)
sample<-gsub("^15.*F","",sample)
sample_names(data)<-sample # update

Metadata <- as.data.frame(read.csv2("Metadata_N_P_R1.csv", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
head(Metadata)
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
original_order<-Metadata$Sample_name # to maintein this sample order in plots
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

# adding a space to some day (technical replicate) to differentiate the labels in the plot (otherwise they will be aggregated by ggplot2!)
Metadata$Experiment_day[Metadata$Sample_name=="N2"]<-paste0(Metadata$Experiment_day[Metadata$Sample_name=="N2"]," ")
Metadata$Experiment_day[Metadata$Sample_name=="N7"]<-paste0(Metadata$Experiment_day[Metadata$Sample_name=="N7"]," ")

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
if(identical( table$FASTQ , as.character(sample_data(data)$FASTQ_ID) ) ){ # same order
  table$Samples<-sample_data(data)$Sample_name
  table$Experiment_day<-sample_data(data)$Experiment_day
  table$factor <- sample_data(data)$Reactor_Type
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
            width = 0.65, stat = "identity") +
  geom_bar( aes(y= After_contam_filter, fill="Relative abundance filters"),
            width= 0.35, stat = "identity") +
  theme_classic( base_size = 10.2) +
  theme(axis.text.x = element_text(size=6, angle = 90, vjust=0.5),
        axis.text.y = element_text(size=5),
        panel.grid.major.y = element_line(size=0.1, color="grey"),
        legend.position = "bottom",
        legend.margin = margin(0, 30, 0, 0),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9.5),
        legend.key.height = unit(0.4,"cm")
  ) +
  scale_fill_manual(name='Number of reads:  ',
                    breaks=c('Original number', 'Read quality filters', 'Relative abundance filters'),
                    values=c('Original number'='green3', 'Read quality filters'='coral', 'Relative abundance filters'='red3')) +
  scale_y_continuous(breaks = c(10000, seq(0, max(table$Original+10000), 25000))) +
  labs(y="Reads abundance", x="FASTQ")
ggsave(file="Data_check/Number_of_reads_pre_and_post_filters.png", width = 6.5, height = 3.8, dpi=300)
# saving also the table itself
write.csv2(table, file="Data_check/Number_of_reads_pre_and_post_filters.csv", quote=F, row.names = F)


suppressWarnings( rm(table, table2, Original_read_number, DADA2_read_number, after_filter_number) )


 
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
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F, ylab = "ASVs", xlab= "Reads amount")
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)



### REMOVING THE UNSATURATED SAMPLES 

if("S13" %in% Metadata$Sample_name){
data<-subset_samples(data, Sample_name!="S13")
Metadata<-Metadata[Metadata$Sample_name!="S13", ]
proof3<-"The unsaturated sample (S13) has been removed"
}



############### \\\\\\ PREPARING THE ORIGINAL (COMPLETE) OBJECTS WITH P AND N LIM \\\\\\ #######################

if(! "proof1" %in% ls() | ! "proof2" %in% ls() | ! "proof3" %in% ls()){
  stop("\nDid you perform the filtering steps yet?\n")
}

if( length(Metadata$Reactor_Type) >1 )  { # then the objects are still complete
  data_complete<-prune_taxa(taxa_sums(data)>0,data)
  metadata_complete<-Metadata
  rm(data)
  
  data_complete.prop <- transform_sample_counts(data_complete, function(ASV) ASV/sum(ASV)*100)
  
  data.genus_complete<-tax_glom(data_complete, taxrank ="Genus", NArm=F)
  data.genus_complete.prop<-transform_sample_counts(data.genus_complete, function(x) (x/sum(x))*100)
  
  data.phy_complete = tax_glom(data_complete, taxrank = "Phylum", NArm = F)
  data.phy_complete.prop <- transform_sample_counts(data.phy_complete, function(ASV) ASV/sum(ASV)*100)
    
  }



####################### GENERAL PCoA BEFORE SUBSETTING (EVERY SAMPLE) ##########################

suppressWarnings(rm(data.prop_temp, data.sqrt_prop))

if( ! "proof1" %in% ls() | ! "proof2" %in% ls() | ! "proof3" %in% ls()  ){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

data.prop.labels<-data.genus_complete.prop
sample_data(data.prop.labels)$Reactor_Type<-gsub("N_limitation","N limiting",sample_data(data.prop.labels)$Reactor_Type)
sample_data(data.prop.labels)$Reactor_Type<-gsub("P_limitation","P limiting",sample_data(data.prop.labels)$Reactor_Type)

{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}

sample_names(data.sqrt_prop)<-as.factor(sample_names(data.sqrt_prop))

# Type 1 (sample names)
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
ggsave(file="Results/N_and_P/PCoA_Helling_on_genera_Names.png", width = 6, height = 4.5, dpi=300)
# Type 2 (experiment day)
plot_ordination(data.sqrt_prop, ordBC, color = "Reactor_Type") +
  scale_color_manual(values=colors) +
  geom_point(size=2.6) +
  geom_point(size=5, alpha= 0.5) +
  theme_classic(base_size = 9) +
  geom_path(aes(group=Sample_Type), color="gray", size= 0.14) +
  geom_text(aes(label= Experiment_day), 
            color="black", size=1.6,
            alpha=0.9,
            show.legend = FALSE) +
  labs(
    color="",
    shape="",
    x=paste("PC1: ",eigval[1],"% variation"),
    y=paste("PC2: ",eigval[2],"% variation")
  )
ggsave(file="Results/N_and_P/PCoA_Helling_on_genera_Exp_Day.png", width = 5.5, height = 4.2, dpi=300)




############################# GENERAL COUNTS EXPORT ##############################


dir.create("Results/N_and_P/General_Raw_counts")
{write.csv2(cbind(as(otu_table(data_complete),"matrix"),as(tax_table(data_complete),"matrix")),file="Results/N_and_P/General_Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy_complete),"matrix"),as(tax_table(data.phy_complete),"matrix")),file="Results/N_and_P/General_Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus_complete),"matrix"),as(tax_table(data.genus_complete),"matrix")),file="Results/N_and_P/General_Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/N_and_P/General_Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy_complete.prop),"matrix"),as(tax_table(data.phy_complete),"matrix")),file="Results/N_and_P/General_Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus_complete.prop),"matrix"),as(tax_table(data.genus_complete),"matrix")),file="Results/N_and_P/General_Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus_complete.prop),"matrix"),as(tax_table(data.genus_complete),"matrix")),file="Results/N_and_P/General_Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy_complete.prop),"matrix"),as(tax_table(data.phy_complete),"matrix")),file="Results/N_and_P/General_Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}



################## ABUNDANCES BAR PLOTS (EVERY SAMPLE) ##################


### TOP 5 Phyla
data.target<-data.phy_complete.prop

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

tabella$Reactor_Type<-gsub("N_limitation","N limiting",tabella$Reactor_Type)
tabella$Reactor_Type<-gsub("P_limitation","P limiting",tabella$Reactor_Type)
levels(tabella$Sample_name)

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
  labs(x="Experiment day", y="Percent abundance", 
       title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/N_and_P/TOP_5_phyla_EVERY_SAMPLE_CuoioDep_Nlim_Plim.png",width=7.5,height=5, dpi=300)
dev.off()

write.xlsx(file = "Results/N_and_P/TOP_5_phyla_EVERY_SAMPLE.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)



### TOP 10 Genera
data.target<-data.genus_complete.prop

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
tabella$Reactor_Type<-gsub("P_limitation","P limiting",tabella$Reactor_Type)

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
  labs(x="Experiment day", y="Percent abundance",
       title = "Ten most abundant genera",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/N_and_P/TOP_10_Genera_EVERY_SAMPLE_CuoioDep_Nlim_Plim.png",width=7.5,height=5,dpi=300)
dev.off()

write.xlsx(file = "Results/N_and_P/TOP_10_genera_EVERY_SAMPLE.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)



################## ABUNDANCES BAR PLOTS (ONLY INOCULUM) ##################


### TOP 5 Phyla Inoculum
data.target<-subset_samples(data.phy_complete.prop, Experiment_day %in% c("1", "1 ") )

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

tabella$Reactor_Type<-gsub("N_limitation","Inoculum ",tabella$Reactor_Type)
#levels(tabella$Sample_name)

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Phylum)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.2), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.5 ),
        legend.position="bottom",
        #legend.margin = margin(0,0,0,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Experiment day", y="Percent abundance", 
       title = "Five most abundant phyla (only inoculum)", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/N_and_P/TOP_5_phyla_INOCULUM.png",width=7,height=5, dpi=300)
dev.off()

write.xlsx(file = "Results/N_and_P/TOP_5_phyla_INOCULUM_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others, data.target)



### TOP 10 Genera Inoculum
data.target<-subset_samples(data.genus_complete.prop, Experiment_day %in% c("1", "1 ") )

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
  tabella$Genus[tabella$Genus=="uncultured"]<-paste("uncultured_f ", tabella$Family[tabella$Genus=="uncultured"] )
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}

tabella$Reactor_Type<-gsub("N_limitation","Inoculum ",tabella$Reactor_Type)

ggplot(data=tabella, aes(x=Experiment_day, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Reactor_Type),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =12) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size= 8.2), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11 ),
        legend.position="bottom", 
        legend.margin = margin(-5,0,2,0)
  ) + 
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Experiment day", y="Percent abundance",
       title = "Ten most abundant genera (only inoculum)",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/N_and_P/TOP_10_Genera_INOCULUM.png",width=7,height=5,dpi=300)
dev.off()

write.xlsx(file = "Results/N_and_P/TOP_10_Genera_INOCULUM_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))




################# STUDING THE ABUNDANCES OF THE PHA ACCUMULATORS (both N lim AND P lim) ######################

lista_PHA<-table_PHA[,1]     # they derive from papers

data.genus.temp<-data.genus_complete.prop
prune.dat_top <- subset_taxa(data.genus.temp, Genus %in% lista_PHA)
prune.dat_top <- filter_taxa(prune.dat_top, function(x) max(x)>0, TRUE) # per scartare quelli del tutto assenti
tax_selected<-as.data.frame(tax_table(prune.dat_top))
tax_table(prune.dat_top)<-as.matrix(tax_selected)
tabella<-psmelt(prune.dat_top)
tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
tabella$Genus<-gsub("Candidatus_","",tabella$Genus)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))

# fill_color_20<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","firebrick3","springgreen4","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","darkgreen","lightgrey","coral","deepskyblue2") 
# fill_color_19<-c("wheat3","deeppink","darkmagenta","bisque2","cyan","yellow2","brown","darkgreen","firebrick3","violet","darkslategray3","blue", "gray", "pink3","yellow4","red","springgreen3","coral","deepskyblue2") 
fill_color_19<-c("grey92","lightblue4","wheat","darkgreen","lightgreen","grey20","deepskyblue2","deeppink","cyan","coral","wheat4","darkmagenta","brown4","springgreen","yellow2","bisque3","firebrick3","darkslategray3","blue", "yellow4","red","black","green3","gold","gray","pink3","violet") 

tabella$Sample_Type<-gsub("Inoculum","N limiting", tabella$Sample_Type)
tabella$Sample_Type<-gsub("R1","P limiting", tabella$Sample_Type)
# tabella$Sample_Type<-gsub("R2","P limiting (R2)", tabella$Sample_Type)

ggplot(data=tabella, aes( x=Experiment_day, y=Abundance, fill=Genus)) + 
  facet_grid2( ~Sample_Type, scales = "free", space="free",
               strip = strip_nested(size="constant"))+
  geom_bar(stat="identity", position="stack", width = 0.95, alpha= 0.85) +
  theme_classic(base_size =12) +
  scale_fill_manual(values=c(fill_color_19)) +
  theme(axis.text.x=element_text(angle=40, hjust=1,vjust=1, size=8.8),
        axis.title.y = element_text(size=11), 
        axis.text.y = element_text(size=9.5), 
        strip.text.x = element_text(size=12.5), 
        legend.key.height = unit(0.42, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text ( size = 9.1 ),
        panel.grid.major.y =  element_line(linewidth=0.25, color = "gray"),
        panel.grid.minor.y =  element_line(linewidth=0.05, color = "gray"),
        legend.position="right",
        plot.margin = margin(2,2,2,2),
        legend.margin =  margin(-36,2,0,-5)
  ) +
  guides(fill=guide_legend(nrow=27)) + 
  scale_y_continuous(lim=c(0,100)) +
  labs(x="Experiment day", y="Percent abundance", 
       fill="",
       title = paste("PHA accumulating genera found \n (",length(unique(taxa_names(prune.dat_top))),"founded on",length(unique(table_PHA[,2])),"known names searched )"))
ggsave(file="Results/N_and_P/PHA_Accumulators_in_both_N_and_Plim.png",width=7.25,height=5.5,dpi=300)
dev.off()




#################### \\\\\\ STARTING THE ANALYSIS OF N lim's SUBSET \\\\\ #######################

data<-subset_samples(data_complete, Reactor_Type=="N_limitation")
data<-prune_taxa(taxa_sums(data)>0, data) # to clean the other taxa

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

metadata<-metadata_complete[metadata_complete$Reactor_Type=="N_limitation",  ]



###################### ABUNDANCES BAR PLOT _ N_lim Only ##########################

second_half <- median(as.numeric(metadata$Experiment_day))
# plot on every sample, computation only on the second half (mature community)


# TOP 5 Phylum stacked bar plot
suppressWarnings(rm(top, others, tabella))
# data.target<-subset_samples(data.phy.prop, ! Experiment_day %in% c("1","1 ","3") ) # no inoculum
data.target<-data.phy.prop
logical <- as.numeric(levels(sample_data(data.target)$Experiment_day)) >= second_half
data.target_second_half <- prune_samples(logical , data.target )
{top <- names(sort(taxa_sums(data.target_second_half), decreasing=TRUE))[1:5]
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
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Phylum)) +
  theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 14.1 ),
        legend.position="bottom",legend.margin = margin(2,0,5,0)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Experiment day", y="Percent abundance", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/N_lim_reactor/Abundances/TOP5_phyla_abundances_second_half.png",width=7,height=6.5, dpi=300) 
dev.off()

# means of TOP5 phyla (no inoculum)
write.xlsx(file = "Results/N_lim_reactor/Abundances/TOP5_phyla_second_half_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))



# TOP 10 genera stacked bar plot
suppressWarnings(rm(top, others, tabella))
# data.target<-subset_samples(data.genus.prop, ! Experiment_day %in% c("1","1 ","3") ) # no inoculum
data.target<-data.genus.prop
logical <- as.numeric(levels(sample_data(data.target)$Experiment_day)) >= second_half
data.target_second_half <- prune_samples(logical , data.target )
{top <- names(sort(taxa_sums(data.target_second_half), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
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
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Genus)) + theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_10) +
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12.5 )) + 
  theme(legend.position="bottom",
        legend.margin = margin(0,0,0,0)) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Experiment day", y="Percent abundance", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/N_lim_reactor/Abundances/TOP10_genera_abundances_second_half.png",width=7,height=6.5, dpi=300) 
dev.off()

# means of TOP10 genera (no inoculum)
write.xlsx(file = "Results/N_lim_reactor/Abundances/TOP_10_genera_second_half_Averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



######################## PCoA BETA DIV _ N limiting #########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day))

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC , color = "Sample_Type") +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  scale_color_manual(values = c("lightblue3")) +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none") +
  theme_classic(base_size = 12.5) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.8, show.legend = FALSE) +
  labs(title="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/N_lim_reactor/Beta_divers_Hellinger_on_genera_WITH_NAMES.png", width = 6.8, height = 5, dpi=300)
# with shade of color AND also exp day
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2") +
  scale_color_gradient2(low="lightblue", mid = "deepskyblue", high="blue", lim=c(NA,NA),
                        breaks=c(1,30,90,150,210,262),
                        midpoint = 90) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(title=" ",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/N_lim_reactor/Beta_divers_Hellinger_on_genera_shade_and_ExpDayText.png", width = 4.2, height = 3.8, dpi=300)



########## CORRELATIONS ABUNDANCES vs TIME _ N limiting (CENTERED LOG RATIO) ################

# selecting only genera with at least 0.1% abundance and at least in 50% of the sample)
data.genus.temp<-data.genus.prop
data.genus.temp <- subset_samples(data.genus.temp, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
min_prevalence<-round(50*length(sample_names(data.genus.temp))/100,digit=0)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
### abundance AND prevalence filter (0.1% abundance at least in 50% of the samples)
who<-as.data.frame(otu_table(data.genus.temp))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 --> "a point"
who<-who[!rowSums(who)>=min_prevalence,]
who<-as.vector(tax_table(data.genus.temp)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)

# switching to log ratio transformation
data.genus.temp <- subset_samples(data.genus, ! Sample_name %in% c("N2","N7") ) # removing the technical replicates
data.genus.temp<-microbiome::transform(data.genus.temp, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sample_name
order_of_samples<-metadata[levels(sample_data(data)$Sample_name), "Sample_name"]   # the sample data has been ordered already during the preparation
order_of_samples <- order_of_samples[! order_of_samples %in% c("N2","N7") ]
abundances<-abundances[ , order_of_samples]


Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "spearman") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_bh<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_bh<0.05,"*","")

write.csv2(Corr_results, "Results/N_lim_reactor/Spearman_CLR_in_N_lim/Correlations_of_CLR_bacteria_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_R1_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )

con <- file("Results/N_lim_reactor/Spearman_CLR_in_N_lim/README_Only_certain_genera_have_been_tested.txt")
sink(con, append = T)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the N_lim samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))



########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R1_pos
suppressWarnings(rm(top, others, tabella))
no_replicates<- subset_samples(data.genus.prop, ! Sample_name %in% c("N2","N7") )
{ tax_table(no_replicates)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(no_replicates, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella<-tabella[order(tabella$Sample_name) , ] # days
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom", legend.margin = margin(1,20,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 0.5)) +
  labs(x="Experiment day", y="Percent abundance", title = "genera positely correlated with time (N limiting)\n according to spearman correlation on CLR counts")
ggsave(file="Results/N_lim_reactor/Spearman_CLR_in_N_lim/Genera_CLR_positively_correlated_with_time.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
#corr_6_top<-significative_abundances_R1[(length(significative_abundances_R1)-5):length(significative_abundances_R1)]
corr_6_top<-significative_abundances_R1_neg
suppressWarnings(rm(top, others, tabella))
no_replicates<- subset_samples(data.genus.prop, ! Sample_name %in% c("N2","N7") )
{ tax_table(no_replicates)<-as.matrix(Taxa.genus.update)
  prune.dat_top <- subset_taxa(no_replicates, Genus %in% corr_6_top)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella<-tabella[order(tabella$Sample_name) , ] # days
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom", legend.margin = margin(1,20,0,0)) +
  guides(color=guide_legend(nrow=2)) + 
  scale_y_sqrt(breaks = seq(0, max(tabella$Abundance), 0.5)) +
  labs(x="Experiment day", y="Percent abundance", title = "genera negatively correlated with time (N limiting)\n according to spearman correlation on CLR counts")
ggsave(file="Results/N_lim_reactor/Spearman_CLR_in_N_lim/Genera_CLR_negatively_correlated_with_time.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_6, prune.dat_top, tax_selected, tabella, corr_6_top))




#################### \\\\\\ STARTING THE ANALYSIS OF P lim's SUBSET \\\\\ #######################

if("data_complete" %in% ls() & "proof2" %in% ls() ){
  to_remove<-ls()
  to_remove<-to_remove[! to_remove %in% c("proof1","proof2","proof3","data_complete","metadata_complete",
                                          "colors","fill_color_5","fill_color_10",
                                          "table_PHA",
                                          "unfiltered_data")]
  rm(list=to_remove)
}

data<-subset_samples(data_complete, Reactor_Type == "P_limitation" )
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
  data.class = tax_glom(data, taxrank = "Class", NArm = F)
  data.order = tax_glom(data, taxrank = "Order", NArm = F)
  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
  data.class.prop <- transform_sample_counts(data.class, function(x) x/sum(x)*100)
  data.order.prop <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
  data.fam.prop <- transform_sample_counts(data.fam, function(x) x/sum(x)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
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

metadata<-metadata_complete[metadata_complete$Reactor_Type=="P_limitation" , ]



###################### ABUNDANCES BAR PLOT _ P_lim Only #########################

second_half <- median(as.numeric(metadata$Experiment_day))
# plot on every sample, computation only on the second half (mature community)


# TOP 5 Phylum stacked bar plot
suppressWarnings(rm(top, others, tabella))
# data.target<-subset_samples(data.phy.prop, ! Experiment_day %in% c("1","1 ","3") ) # no inoculum
data.target<-data.phy.prop
logical <- as.numeric(levels(sample_data(data.target)$Experiment_day)) >= second_half
data.target_second_half <- prune_samples(logical , data.target )
{top <- names(sort(taxa_sums(data.target_second_half), decreasing=TRUE))[1:5]
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
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Phylum)) +
  theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_fill_manual(values=fill_color_5) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 14.1 ),
        legend.position="bottom",legend.margin = margin(2,0,5,0)
  ) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="Experiment day", y="Percent abundance", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/P_lim_reactor/Abundances/TOP5_phyla_abundances.png",width=7,height=6.5, dpi=300) 
dev.off()

# means of TOP5 phyla (no inoculum)
write.xlsx(file = "Results/P_lim_reactor/Abundances/TOP5_phyla_second_half_averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))



# TOP 10 Genera stacked bar plot
suppressWarnings(rm(top, others, tabella))
# data.target<-subset_samples(data.phy.prop, ! Experiment_day %in% c("1","1 ","3") ) # no inoculum
data.target<-data.genus.prop
logical <- as.numeric(levels(sample_data(data.target)$Experiment_day)) >= second_half
data.target_second_half <- prune_samples(logical , data.target )
{top <- names(sort(taxa_sums(data.target_second_half), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.target)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
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
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, fill=Genus)) + theme_classic(base_size =14) + 
  # facet_grid2(Sampling_date+Sample_Type~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_fill_manual(values=fill_color_10) +
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 0.95) +
  theme(axis.text.x=element_text(angle=50, vjust=1, hjust = 1, size=11),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12.5 )) + 
  theme(legend.position="bottom",
        legend.margin = margin(0,0,0,0)) +
  guides(fill=guide_legend(nrow=4)) + 
  labs(x="Experiment day", y="Percent abundance", caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/P_lim_reactor/Abundances/TOP10_genera_abundances_second_half_Averages.png",width=7,height=6.5, dpi=300) 
dev.off()

# means of TOP10 genera (no inoculum)
write.xlsx(file = "Results/P_lim_reactor/Abundances/TOP_10_genera_second_half_Averages.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))


suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))




######################## PCoA BETA DIV _ P limiting #########################

# on genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Experiment_day2<-as.numeric(as.character(sample_data(data.prop.labels)$Experiment_day))

# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC , color = "Sample_Type") +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  scale_color_manual(values = c("coral")) +
  geom_point(size=3.5, alpha=0.5,  show.legend = F) +
  guides(color="none") +
  theme_classic(base_size = 12.5) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.8, show.legend = FALSE) +
  labs(title="",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/P_lim_reactor/Beta_divers_Hellinger_on_genera_WITH_NAMES.png", width = 6.8, height = 5, dpi=300)
# with shade of color AND experiment day
plot_ordination(data.sqrt_prop, ordBC, color = "Experiment_day2") +
  scale_color_gradient2(low="orange", mid = "coral2", high="red3", lim=c(315,409), midpoint = 350,
                        breaks=c(315, 335, 365, 390, 409)
  ) +
  geom_point(size=4.3, alpha=0.3) +
  geom_point(size=2.3, alpha=1) +
  geom_path(aes(group=Sample_Type), col= "darkgray", size = 0.3,
            arrow=arrow(length =unit(0.25,"cm"), type = "closed") 
  ) + 
  theme_classic(base_size = 8.8) +
  theme(title=element_text(size=7.5),
        axis.title.x = element_text(size=8.8),
        axis.title.y = element_text(size=8.8),
        legend.margin = margin(0,0,0,-10),
        legend.title = element_text(size=9),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.height = unit(0.3,"cm"),
        legend.key.width  = unit(0.8,"cm")
  ) +
  geom_text(aes(label=sample_data(data.sqrt_prop)$Experiment_day2), color="black", vjust= -0.2, size=2.2, show.legend = FALSE) +
  labs(title=" ",
       color="Experiment day    ",
       x=paste("PC1: ",eigval[1],"% variation"),       y=paste("PC2: ",eigval[2],"% variation"),
       #caption="NB: the points are NOT MANUALLY placed along the plot"
  )
ggsave(file="Results/P_lim_reactor/Beta_divers_Hellinger_on_genera_shade_and_ExpDayText.png", width = 4.2, height = 3.8, dpi=300)




########## CORRELATIONS ABUNDANCES vs TIME _ P limiting (CENTERED LOG RATIO) ################

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
data.genus.temp<-microbiome::transform(data.genus, "clr") # clr using natural logarithm and adding minimal pseudocounts to rows (half of lowest relative percentace among the samples)
# head(otu_table(data.genus.temp))
# head(otu_table(data.genus))
data.genus.temp<-subset_taxa(data.genus.temp, ! Genus %in% who)
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)

# extracting the abundances
abundances<-otu_table(data.genus.temp)
row.names(abundances)<-as.character(tax_table(data.genus.temp)[,"Genus"])

# ordering the sample names according to the time flow 
row.names(metadata)<-metadata$Sample_name
order_of_samples<-metadata[levels(sample_data(data)$Sample_name), "Sample_name"]   # the sample data has been ordered already during the preparation
abundances<-abundances[ , order_of_samples]


Corr_results<-NULL
for(i in 1:length(row.names(abundances))){
  save<-cor.test(as.numeric(abundances[i,]), 1:length(colnames(abundances)), method = "spearman") # correlated with the flow of time (the samples have been ordered accordingly)
  new_row<-cbind.data.frame( row.names(abundances)[i] , save$estimate , save$p.value )
  Corr_results<-rbind.data.frame(Corr_results, new_row)
}
Corr_results<- Corr_results[! Corr_results$`row.names(abundances)[i]` %in% c("uncultured_ o uncultured","NA_ o NA") , ]
{
  row.names(Corr_results)<-Corr_results$`row.names(abundances)[i]`
  Corr_results<-Corr_results[ , -1]
  colnames(Corr_results)<-c("rho","pvalue")
}
Corr_results$padj_BH<-p.adjust(Corr_results$pvalue, method = "BH")
Corr_results$sign<-ifelse(Corr_results$padj_BH<0.05,"*","")

write.csv2(Corr_results, "Results/P_lim_reactor/Spearman_CLR_in_P_lim/Correlations_of_CLR_bacteria_abundances_with_time.csv", row.names = T, quote = F)

temp_for_save<-Corr_results[Corr_results$sign=="*", ]
temp_for_save<-temp_for_save[order(temp_for_save$rho, decreasing = T), ] # NB: the first are positively correlated, while the last are negatively
significative_abundances_R1_pos<-row.names(temp_for_save[temp_for_save$rho>0 , ] )
significative_abundances_R1_neg<-row.names(temp_for_save[temp_for_save$rho<0 , ] )


con <- file("Results/P_lim_reactor/Spearman_CLR_in_P_lim/README_Only_certain_genera_have_been_tested_in_Spearman.txt")
sink(con, append = T)
cat("Only the genera with at least 0.1% minimal abundance at least in ~50% of the P_lim samples (-->",min_prevalence,"samples) have been correlated with the passing of time.")
sink()
close(con)

suppressWarnings(rm(con, temp_for_save, Corr_results, abundances, order_of_samples, data.genus.temp, min_prevalence, save, new_row))



########### LINE PLOTS OF STRONGEST CORRELATIONS

# fill_color_6<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # 5 ... but 6 actually!

# THE 6 GENERA MOST POSIVELY CORRELATED (rho +) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R1_pos
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

fill_color_corr<-c("grey92","lightblue4","wheat","darkgreen","lightgreen","grey20","deepskyblue2","deeppink","cyan","coral","wheat4","darkmagenta","brown4","springgreen","yellow2","bisque3","firebrick3","darkslategray3","blue", "yellow4","red","black","green3","gold","gray","pink3","violet") 
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10.5),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 9 ),
        legend.position="bottom", legend.margin = margin(1,20,0,0)) +
  guides(color=guide_legend(nrow=6)) +
  scale_color_manual(values=fill_color_corr) +
  scale_y_sqrt(breaks = c( 1, seq(0, max(tabella$Abundance), 5)) ) +
  labs(x="Experiment day", y="Percent abundance",
       color="",
       title = "genera positely correlated with time (P limiting)\n according to spearman correlation on CLR counts")
ggsave(file="Results/P_lim_reactor/Spearman_CLR_in_P_lim/Genera_CLR_positively_correlated_with_time.png",width=7,height=5, dpi=300) 
dev.off()


# THE 6 GENERA MOST NEGATIVELY CORRELATED (rho -) WITH THE PASSING TIME
corr_6_top<-significative_abundances_R1_neg
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
ggplot(data=tabella, aes(y=Abundance, x=Experiment_day, color=Genus)) +
  theme_classic(base_size =13.8) + 
  geom_point(aes(color=Genus), size =2) +
  geom_point(aes(color=Genus), size =3.8, alpha=0.4) +
  geom_line(aes(group=Genus), size= 0.4, alpha= 0.4) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10.5),
        axis.text.x=element_text(angle=30, vjust=1, hjust=1, size=10),
        panel.grid.major.y = element_line(linewidth=0.05, color="gray"),
        panel.grid.minor.y = element_line(linewidth=0.02, color="gray"),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 9 ),
        legend.position="bottom", legend.margin = margin(1,20,0,0)) +
  guides(color=guide_legend(nrow=6)) + 
  scale_y_sqrt(breaks = c( 1, seq(0, max(tabella$Abundance), 5)) ) +
  labs(x="Experiment day",
       y="Percent abundance",
       color="",
       title = "genera negatively correlated with time (P limiting)\n according to spearman correlation on CLR counts")
ggsave(file="Results/P_lim_reactor/Spearman_CLR_in_P_lim//Genera_CLR_negatively_correlated_with_time.png",width=7,height=5, dpi=300) 
dev.off()

suppressWarnings(rm(fill_color_6, prune.dat_top, tax_selected, tabella, corr_6_top))



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
