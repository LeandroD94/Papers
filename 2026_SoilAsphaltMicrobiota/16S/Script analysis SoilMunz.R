#################### PREPARING THE ENVIRONMENT #################

{ 
  library("phyloseq")
  # graphical packages
  library("ggplot2")
  library("ggpubr")
  library("ggh4x")
  # analysis packages
  library("vegan")
  library("ecodist")
  library("ggvenn")
  # utilities
  library("reshape")
  library("xlsx")  
  library("qiime2R")
}

{dir.create("Data_check")
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Beta_Diversity")
  dir.create("Results/UniqueBacteria")
}
options(scipen = 100) # disable scientific annotation




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy_SILVA.qza")

# changing names
sample<-sample_names(data)
original_names<-sample
sample<-gsub("^.*F","",sample)
sample_names(data)<-sample # update

# NB: the order of the row in the Metadata is employed also to set the order of the samples in the plots
metadata <- as.data.frame(read.table("metadata_EliaSuolo.txt", header = T, sep="\t"))
row.names(metadata)<-metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)
original_length<-length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])
original_order<-metadata$Sample_name # to maintein this sample order in plots
metadata<-metadata[sample, ]
identical(as.numeric(length(metadata$FASTQ_ID[!is.na(metadata$FASTQ_ID)])),as.numeric(original_length))

# metadata$Sample_name<- factor(metadata$Sample_name, levels = unique(metadata$Sample_name) )
metadata$Sample_name <- factor(metadata$Sample_name, levels = metadata$Sample_name )

sample_data(data)<-metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

# restoring the original row order of the metadata (it sets also the order in the plots!)
row.names(metadata)<-metadata$Sample_name
metadata<-metadata[original_order, ]

rm(original_length)




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
target_filter<- filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE)
filtered<-taxa_names( target_filter )
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_0005_max_cutoff.csv")


# if the ASVs are filtered at ASV level then only the representative of that genera will be filtered ...
# otherwise, if the genera name are filtered then also every eventual NA and uncultured will be filtered also in other families
# then, a more specific taxonomy table (with a partial path included in the unidentified genus) will be now build

taxa_temp<-as.data.frame(tax_table(target_filter))
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
  genera_to_filter<-taxa_temp[["Genus"]]  # every NA or uncultured will now have a more detailed path in the genus name
}

unf_data2<-unfiltered_data
taxa_temp<-as.data.frame(tax_table(unf_data2))
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
  tax_table(unf_data2)<-as.matrix(taxa_temp)
}

unf_data2<-subset_taxa(unf_data2, ! Genus %in% genera_to_filter)
ASVs_code_to_maintain<- taxa_names(unf_data2)  # NB: taxa_names here means ASVs codes
data<-prune_taxa(ASVs_code_to_maintain, data)  # working on ASVs searched with specific names between different levels of taxonomy (raw vs genus) 
data<-prune_taxa(taxa_sums(data) > 0, data) 
rm(filtered, data.genus.temp, unf_data2, genera_to_filter, target_filter)



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
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}


png(file="Data_check/Rarefaction_curve.png",width=2200,height=1800, res=300)
r<-rarecurve(t(as(otu_table(data.genus),"matrix")), step=100,label=F, ylab = "Genera", xlab= "Reads amount")
evalslopes(r,sample_data(data.genus)$Sample_name,lim=0.001,cex=1)
dev.off()
rm(r)



proof2<-"No unsaturated sample to remove"




#################### \\\\\\ STARTING THE ANALYSIS \\\\\ #######################


if(! "proof1" %in% ls() ){
  stop("\nDid you perform the filtering steps yet?\n")
}



### everything ready, let's start with the analysis
{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
#  data.class = tax_glom(data, taxrank = "Class", NArm = F)
#  data.order = tax_glom(data, taxrank = "Order", NArm = F)
#  data.fam = tax_glom(data, taxrank = "Family", NArm = F)
  data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{ data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
  data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
#  data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
#  data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
#  data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{ Taxa.genus<-as.data.frame(tax_table(data.genus))
#  Taxa.fam<-as.data.frame(tax_table(data.fam))
  Taxa.phy<-as.data.frame(tax_table(data.phy))
#  Taxa.class<-as.data.frame(tax_table(data.class))
#  Taxa.order<-as.data.frame(tax_table(data.order))
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


# save.image("data_prepared_for_analysis.RData")




########################### COUNTS EXPORT ##########################################

dir.create("Results/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Abundances/Raw_counts/ASV_counts.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Raw_counts/counts_phylum.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Raw_counts/counts_class.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Raw_counts/counts_order.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

dir.create("Results/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Abundances/Relative_abundances/counts_class.csv",quote=F)
  #write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Abundances/Relative_abundances/counts_order.csv",quote=F)
#  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  #write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
  #write.xlsx(cbind(as(otu_table(data.prop),"matrix"),as(tax_table(data.prop),"matrix")),file="Results/Abundances/Relative_abundances/counts_species.xlsx",showNA = F, col.names = T)
}


write.csv2(cbind(as(otu_table(unfiltered_data),"matrix"),as(tax_table(unfiltered_data),"matrix")),file="Data_check/ASVs_TAXA_counts_before_filters.csv",quote=F)


# This one is designed to be quick to use...
temp <- data.genus.prop
sample_names(temp)<- sample_data(temp)$Sample_name
write.xlsx(cbind(round(as(otu_table(temp),"matrix"), digits = 4),
                 as(tax_table(temp),"matrix")[ , c("Kingdom","Phylum","Class","Order","Family","Genus")]
                 ),
           file="Results/Abundances/Relative_abundances/Counts_genus_as_percentages_WITH_taxonomy.xlsx",showNA = F, col.names = T)
rm(temp)




###################### ABUNDANCES BAR PLOTS ##########################

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
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bacter.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatimon", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Pedosphaer.","Pedosphaeraceae", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_name<- factor(tabella$Sample_name, levels = c("Compost","R1","R2","R3","R4",
                                                             "S1","S2","S3","S4",
                                                             "T1","T2","T3","T4")
)

custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue3",
                 "slateblue3","cyan","black","brown4",
                 "greenyellow", "red",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                 "yellow", "chartreuse3","yellow3", "chocolate4","coral2","springgreen3",
                 "slateblue1", "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.17, "cm"),
  legend.text = element_text ( size = 7.5 ),
  legend.position="bottom",
  legend.margin = margin(-13 ,35, 1 ,5),
  plot.margin = margin(2,1,1,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 30 ")
ggsave(file="Results/Abundances/Most_Abundant_Genera___AveragesEverySample.png",width=5.5,height=4, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



### Again, to be plot side by side with other graphics...
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
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bacter.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatimon", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Saprospir.","Saprospiraceae", tabella$Genus)
  tabella$Genus<-gsub ("Pedosphaer.","Pedosphaeraceae", tabella$Genus)
  tabella$Genus<-gsub ("67-14","Solirubrobacterac. 67-14", tabella$Genus, fixed=T)
  tabella$Genus<-gsub ("KCM-B-112","Acidithiobacil. KCM-B-112", tabella$Genus, fixed=T)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_name<- factor(tabella$Sample_name, levels = c("Compost","R1","R2","R3","R4",
                                                             "S1","S2","S3","S4",
                                                             "T1","T2","T3","T4")
                             )
custom_colors<-c( rev(c("navyblue","gray80","brown4",
                         "darkgreen", "slateblue4",
                        "gold3","deepskyblue2",
                        "indianred",  "coral2",
                        "slateblue1", "grey10", "pink", 
                        "yellow","lightgreen", "darkmagenta","chocolate4","cyan",
                        "red", "green2")),
                  "darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=20,
                                 vjust=1,
                                 hjust=1,
                                 size= 7
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=9),
  axis.title.x = element_text(vjust=2.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=10),
  legend.key.height = unit(0.12, "cm"),
  legend.key.width = unit(0.18, "cm"),
  legend.text = element_text ( size = 7.85 ),
  legend.position="bottom",
  legend.margin = margin(-13 ,35, 1 ,5),
  plot.margin = margin(2,1,1,1)) +
  guides(fill=guide_legend(nrow=5)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="")
ggsave(file="Results/Abundances/Most_Abundant_Genera___AveragesEverySample2.png",width=4.8,height=4, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



### AGAIN, BUT EXCLUDING THE COMPOST SAMPLE
suppressWarnings(rm(top, others, tabella, unass_data))
data_subset <- subset_samples(data.genus.prop, Sample_name!="Compost")
{top <- names(sort(taxa_sums(data_subset), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bact.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Pedosphaer.","Pedosphaeraceae", tabella$Genus)
  tabella$Genus<-gsub ("_terrestrial_group"," terrestrial", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue3",
                 "slateblue3","cyan","black","brown4",
                 "greenyellow", "red",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                 "yellow", "chartreuse3","yellow3", "chocolate4","coral2","springgreen3",
                 "slateblue1", "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.17, "cm"),
  legend.text = element_text ( size = 7.15 ),
  legend.position="bottom",
  legend.margin = margin(-14 ,35, 1 ,5),
  plot.margin = margin(2,1,1,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 30 ")
ggsave(file="Results/Abundances/Most_Abundant_Genera___AveragesExcludingCompost.png",width=5.5,height=4, dpi=300)
dev.off()

# means of TOP Genera
to_save<- cbind.data.frame("Overall Mean"=as.numeric( apply(otu_table(prune.dat_top),1,mean) ),
                           "Mean in R" = as.numeric( apply(otu_table(subset_samples(prune.dat_top,Group=="R")),1,mean) ) ,
                           "Mean in T" = as.numeric( apply(otu_table(subset_samples(prune.dat_top,Group=="T")),1,mean) ),
                           "Mean in S" = as.numeric( apply(otu_table(subset_samples(prune.dat_top,Group=="S")),1,mean) ),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]])
to_save<-to_save[order(to_save$`Overall Mean`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/Most_Abundant_Genera___AveragesExcludingCompost_VALUES.xlsx", row.names = F, 
           x = to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



### AGAIN, BUT ONLY GROUP R
suppressWarnings(rm(top, others, tabella, unass_data))
data_subset <- subset_samples(data.genus.prop, Group=="R" & Sample_name!="Compost")
{top <- names(sort(taxa_sums(data_subset), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bact.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Pedosphaer.","Pedosphaeraceae", tabella$Genus)
  tabella$Genus<-gsub ("_terrestrial_group"," terrestrial", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue3",
                 "slateblue3","cyan","black","brown4",
                 "greenyellow", "red",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                 "yellow", "chartreuse3","yellow3", "chocolate4","coral2","springgreen3",
                 "slateblue1", "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.17, "cm"),
  legend.text = element_text ( size = 7.35 ),
  legend.position="bottom",
  legend.margin = margin(-12 ,35, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 30 ")
ggsave(file="Results/Abundances/Most_Abundant_Genera___OnlyGroupR.png",width=5.2,height=4, dpi=300)
dev.off()

to_save<- cbind.data.frame("Overall Mean in R"=as.numeric( apply(otu_table(prune.dat_top),1,mean) ),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
)
to_save<-to_save[order(to_save$`Overall Mean in R`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/Most_Abundant_Genera___OnlyGroupR.xlsx", row.names = F, 
           x = to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



### AGAIN, BUT ONLY GROUP S
suppressWarnings(rm(top, others, tabella, unass_data))
data_subset <- subset_samples(data.genus.prop, Group=="S")
{top <- names(sort(taxa_sums(data_subset), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bact.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Vicinamibacter.","Vicinamibacter", tabella$Genus, fixed = T)
  tabella$Genus<-gsub ("Rhodanobacter.","Rhodanobact.", tabella$Genus, fixed = T)
  tabella$Genus<-gsub ("Saprospir.","Saprospiraceae", tabella$Genus, fixed = T)
  tabella$Genus<-gsub ("_terrestrial_group"," terrestrial", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue3",
                 "slateblue3","cyan","black","brown4",
                 "greenyellow", "red",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                 "yellow", "chartreuse3","yellow3", "chocolate4","coral2","springgreen3",
                 "slateblue1", "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.145, "cm"),
  legend.text = element_text ( size = 6.8 ),
  legend.position="bottom",
  legend.margin = margin(-12 ,35, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 30 ")
ggsave(file="Results/Abundances/Most_Abundant_Genera___OnlyGroupS.png",width=5.2,height=4, dpi=300)
dev.off()

to_save<- cbind.data.frame("Overall Mean in S"=as.numeric( apply(otu_table(prune.dat_top),1,mean) ),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
)
to_save<-to_save[order(to_save$`Overall Mean in S`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/Most_Abundant_Genera___OnlyGroupS.xlsx", row.names = F, 
           x = to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



### AGAIN, BUT ONLY GROUP T
suppressWarnings(rm(top, others, tabella, unass_data))
data_subset <- subset_samples(data.genus.prop, Group=="T")
{top <- names(sort(taxa_sums(data_subset), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bact.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Vicinamibacter.","Vicinamibacter", tabella$Genus, fixed = T)
  tabella$Genus<-gsub ("Rhodanobacter.","Rhodanobact.", tabella$Genus, fixed = T)
  tabella$Genus<-gsub ("_terrestrial_group"," terrestrial", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue3",
                 "slateblue3","cyan","black","brown4",
                 "greenyellow", "red",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                 "yellow", "chartreuse3","yellow3", "chocolate4","coral2","springgreen3",
                 "slateblue1", "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.145, "cm"),
  legend.text = element_text ( size = 6.85 ),
  legend.position="bottom",
  legend.margin = margin(-12 ,35, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 30 ")
ggsave(file="Results/Abundances/Most_Abundant_Genera___OnlyGroupT.png",width=5.2,height=4, dpi=300)
dev.off()

to_save<- cbind.data.frame("Overall Mean in T"=as.numeric( apply(otu_table(prune.dat_top),1,mean) ),
                           "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]
)
to_save<-to_save[order(to_save$`Overall Mean in T`, decreasing=T), ]
write.xlsx(file = "Results/Abundances/Most_Abundant_Genera___OnlyGroupT.xlsx", row.names = F, 
           x = to_save )

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, top, others, tabella_top, tabella_top))



###################### ABUNDANCES OF ARCHAEA ##########################

### TOP Genera Archaea
suppressWarnings(rm(top, others, tabella, unass_data))
data.temp <- subset_taxa(data.genus.prop, Kingdom=="d__Archaea")
{top <- names(sort(taxa_sums(data.temp), decreasing=TRUE))[1:10]
  prune.dat_top <- prune_taxa(top,data.temp)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data.temp)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data.temp)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_name<- factor(tabella$Sample_name, levels = c("Compost","R1","R2","R3","R4",
                                                             "S1","S2","S3","S4",
                                                             "T1","T2","T3","T4")
)
custom_colors<-c("chartreuse","black", "firebrick4",
                 "orange","red", 
                 "darkblue","grey65", "magenta",
                 "yellow","darkgreen","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  # scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.17, "cm"),
  legend.text = element_text ( size = 7.5 ),
  legend.position="bottom",
  legend.margin = margin(-13 ,35, 1 ,5),
  plot.margin = margin(2,1,1,1)) +
  guides(fill=guide_legend(nrow=4)) +
  labs(x="", y="Percentual abundance of clades (Archaea)",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " 'Others' includes every genus below rank 30 ")
ggsave(file="Results/Abundances/Most_Abundant_Archaea.png",width=4.8,height=4, dpi=300)

dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))




########################## ALFA DIVERSITY ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Group", color="Group")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
{ H<-pAlpha$data[pAlpha$data$variable=="Shannon", ]
  obs<-pAlpha$data[pAlpha$data$variable=="Observed", ]
  # adding evenness
  identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  New_data$Sample_label <- gsub("[A-Z]","", New_data$Sample_name)
  New_data$Sample_label[New_data$Sample_name=="Compost"] <- "Comp"
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha + 
  geom_boxplot(data=pAlpha$data, aes(x=Group, y=value, color=NULL), alpha=0, linewidth=0.2, color="gray30") +
  theme_classic2(base_size = 8.2) + 
  scale_color_manual(values = c("R"="coral", 
                                "S"="royalblue1",
                                "T"="chartreuse2")) +
  labs(x=""
       # title="Alpha diversity in Healthy and CCHS",
  ) +
  geom_text(aes(label=pAlpha$data$Sample_label), color="white", size=2.1, show.legend = FALSE) +
  geom_text(aes(label=pAlpha$data$Sample_label), color="black", size=1.95, vjust=0.54, show.legend = FALSE) +
  guides(fill="none", color="none") +
  theme(axis.text.x= element_text(angle=0, vjust=0.5, hjust=0.5, size=7.5))
ggsave(file="Results/Alfa_diversity_of_Genera.png", width = 3.35,height =2.5, dpi=300)

suppressWarnings ( rm(con, pAlpha, alphadt,H, ev, obs, a, New_data) )




####################### PCoA BETA DIV #########################

### on Genera
data.prop.labels<-data.genus.prop
sample_data(data.prop.labels)$Labels <- gsub("[A-Z]","",sample_data(data.prop.labels)$Sample_name) 
sample_data(data.prop.labels)$Labels[sample_data(data.prop.labels)$Sample_name=="Compost"] <- "Comp  "
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Group") +
  scale_color_manual(values = c("R"="coral", 
                                "S"="royalblue1",
                                "T"="chartreuse")) +
  geom_point(size=3.5, alpha= 0.5) +
  theme_classic(base_size = 9) +
  theme( legend.text = element_text(size=9)) +
  geom_text(aes(label=Labels), color="black", size=2.5, show.legend = FALSE) +
  labs( # title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)",
    # subtitle= paste0("PERMANOVA Pr(>F) Healthy - CCHS: ", perm_g_H["Condition2", "Pr(>F)"] ),
    color="",
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_Diversity/PCoA_EVERY_SAMPLE_computed_on_genera.png", width = 3.5, height = 3.2, dpi=300)



### AGAIN, EXCLUDING THE COMPOST
data.prop.labels<-subset_samples( data.genus.prop , Sample_name!="Compost" )
sample_data(data.prop.labels)$Labels <- gsub("[A-Z]","",sample_data(data.prop.labels)$Sample_name) 
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Group") +
  scale_color_manual(values = c("R"="coral", 
                                "S"="royalblue1",
                                "T"="chartreuse")) +
  geom_point(size=2.3, alpha= 0.4) +
  theme_classic(base_size = 8.25) +
  theme( legend.text = element_text(size=9)) +
  stat_ellipse( aes(group=Group) , linewidth=0.12) +
  geom_text(aes(label=Labels), color="gray80", size=2.8, show.legend = FALSE) +
  geom_text(aes(label=Labels), color="black", size=2.6, show.legend = FALSE) +
  labs( # title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed Genera)",
    # subtitle= paste0("PERMANOVA Pr(>F) Healthy - CCHS: ", perm_g_H["Condition2", "Pr(>F)"] ),
    color="",
    x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_Diversity/PCoA_EXCLUDING_COMPOST_computed_on_genera.png", width = 3, height = 2.65, dpi=300)




#################### VEEN DIAGRAM (COMMON AND UNIQUE BACTERIA) ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-subset_samples( data.genus.temp , Sample_name!="Compost" )

### abundance AND prevalence filter (0.05% abundance at least in 2 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.05, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>1,] # more than 2 "points" --> at least in 2 samples 
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


R<-subset_samples(data.venn, Group=="R")
R<-as.character(tax_table(prune_taxa(taxa_sums(R)>0, R))[,"Genus"])

GroupT<-subset_samples(data.venn, Group=="T")
GroupT<-as.character(tax_table(prune_taxa(taxa_sums(GroupT)>0, GroupT))[,"Genus"])

S<-subset_samples(data.venn, Group=="S")
S<-as.character(tax_table(prune_taxa(taxa_sums(S)>0, S))[,"Genus"])


ONLY_IN_T<- GroupT[! GroupT %in% S & ! GroupT %in% R]
ONLY_IN_T_characters<- ONLY_IN_T
ONLY_IN_T<- paste(ONLY_IN_T, collapse = ", ")
head(ONLY_IN_T)

ONLY_IN_R<- R[! R %in% S & ! R %in% GroupT]
ONLY_IN_R_characters<- ONLY_IN_R
ONLY_IN_R<- paste(ONLY_IN_R, collapse = ", ")
head(ONLY_IN_R)

ONLY_IN_S<- S[! S %in% R & ! S %in% GroupT]
ONLY_IN_S_characters<- ONLY_IN_S
ONLY_IN_S<- paste(ONLY_IN_S, collapse = ", ")
head(ONLY_IN_S)

con<-file("Results/UniqueBacteria/Unique_Bacteria_see_Venn_Diagramm.txt")
sink(con, append=TRUE)
cat("ONLY IN S", fill=TRUE)
cat(ONLY_IN_S)
cat("\n\nONLY IN R", fill=TRUE)
cat(ONLY_IN_R)
cat("\n\nONLY IN T", fill=TRUE)
cat(ONLY_IN_T)
sink()
close(con)

# plot
x<-list(FlaskR=R,FlaskT=GroupT,FlaskS=S)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("chartreuse","coral","deepskyblue")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram \n(only genera with minimal abundance > 0.05% \n at least in 2 samples for group)")
ggsave(filename = "Results/UniqueBacteria/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()



### PLOTTING THESE UNIQUE BACTERIA
suppressWarnings(rm(top, others, tabella, unass_data))
data_subset <- subset_taxa(data.venn, Genus %in% c(ONLY_IN_R_characters,ONLY_IN_S_characters,ONLY_IN_T_characters) )
{prune.dat_top <-data_subset
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  tabella<-psmelt(prune.dat_top)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
  tabella$Genus<-gsub ("bacterales","bact.", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
  tabella$Genus<-gsub ("aceae",".", tabella$Genus)
  tabella$Genus<-gsub ("Pirellul","Pirellulaceae", tabella$Genus, fixed=T)
  tabella$Genus<-gsub ("f Unknown_Family","Gammaproteobacteria", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = unique(tabella$Genus))
}
custom_colors<-c("green","deeppink",
                 "orange","darkmagenta", 
                 "navyblue","grey25",
                 "wheat2","cyan","chocolate4",
                 "red",
                 "deepskyblue2","darkgreen",
                 "slateblue4", "darkred" , "violet", "blue1",
                 "yellow2", "chartreuse2","yellow3", "coral2","springgreen4",
                 "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) , breaks = c(0,1,2,3,4,5) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  panel.grid.major.y = element_line( colour="gray", linewidth=0.15 ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.17, "cm"),
  legend.text = element_text ( size = 7.15 ),
  legend.position="bottom",
  legend.margin = margin(-10 ,28, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " Only unique genera with % abund over 0.05% in at least 2 sample in each group ")
ggsave(file="Results/UniqueBacteria/Unique_Bacteria_In_Each_Group.png",width=5,height=4, dpi=300)



### PLOTTING MOST ABUNDANT COMMON BACTERIA
suppressWarnings(rm(top, others, tabella, unass_data))
data_subset <- subset_taxa(data.venn, Genus %in% GroupT[ GroupT %in% S & GroupT %in% R] )
{ top <- names(sort(taxa_sums(data_subset), decreasing=TRUE))[1:29]
  prune.dat_top <- prune_taxa(top,data_subset)
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  others<-taxa_names(data_subset)
  others<-others[!(others %in% top)]
  prune.data.others<-prune_taxa(others,data_subset)
  tabella_top<-psmelt(prune.dat_top)
  tabella_others<-psmelt(prune.data.others)
  tabella_others$Genus<-"Others"
  tabella<-rbind.data.frame(tabella_top,tabella_others)
  tabella<-tabella[order(sort(tabella$Abundance, decreasing = T)), ]
}
tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
tabella$Genus<-gsub ("uncultured","uncult", tabella$Genus)
tabella$Genus<-gsub ("bacterales","bact.", tabella$Genus)
tabella$Genus<-gsub ("Gemmatimonad","Gemmatim", tabella$Genus)
tabella$Genus<-gsub ("aceae",".", tabella$Genus)
tabella$Genus<-gsub ("o Vicinamibact.","o Vicinamib", tabella$Genus)
tabella$Genus<-gsub ("Pedosphaer.","Pedosphaeraceae", tabella$Genus)
tabella$Genus<-gsub ("_terrestrial_group"," terrestrial", tabella$Genus)
tabella$Genus<-gsub ("Anaeroline."," Anaerolinaceae", tabella$Genus, fixed = T)
tabella$Genus<-gsub (" Caldiline."," Caldilin.", tabella$Genus, fixed = T)
tabella$Genus<-gsub ("Sandaracin.","Sandaracinac.", tabella$Genus, fixed = T)
tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue3",
                 "slateblue3","cyan","black","brown4",
                 "greenyellow", "red",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                 "yellow", "chartreuse3","yellow3", "chocolate4","coral2","springgreen3",
                 "slateblue1", "orangered","darkslategray3")
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0,1) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  panel.grid.major.y = element_line( colour="gray", linewidth=0.15 ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.15, "cm"),
  legend.text = element_text ( size = 6.94 ),
  legend.position="bottom",
  legend.margin = margin(-10 ,32, 1 ,5),
  plot.margin = margin(4,1,4,1)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " Here ONLY common genera across groups have been included ('others' here are also common bacteria) ")
ggsave(file="Results/UniqueBacteria/Common_Bacteria_Across_Groups.png",width=5.18,height=4, dpi=300)


suppressWarnings(rm(ONLY_IN_T, ONLY_IN_R, ONLY_IN_S, 
                    x, con, GroupT, R, S, who,
                    data.venn)
)




################# SEARCHING FOR AOB AND NOB BACTERIA ####################

# AOBs derive from from MIDAS, also in other articles I can't find other names...
AOB_list <- c("Nitrosomonas", "Nitrosococcus", "Nitrosomonadales",   # synonymous (see MIDAS) or higher levels
              "Nitrosospira",  "Nitrosolobus", "Nitrosovibrio",   # synonymous, see MIDAS
              "966-1",     # https://www.sciencedirect.com/science/article/pii/S0045653523014649
              "Ellin6067" )   # https://www.sciencedirect.com/science/article/pii/S0960852423008970

# adding also AOA, from PMID: 24559743
AOA_list <- c("Nitrososphaera", "Nitrososphaerota", 
              "Nitrosocaldus", 
              "Nitrosotalea",
              "Nitrocosmicus")  
AOA_list <- c(AOA_list,
              "Nitrosarchaeum", # DOI: 10.1099/ijsem.0.002926
              "Nitrosacidococcus", # https://doi.org/10.1016/j.wroa.2022.100157
              "Nitrosocosmicus",  # https://www.sciencedirect.com/science/article/pii/S0038071720300225
              "Nitrosocaldus", "Nitrosothermus" , "Nitrosocaldaceae", # https://doi.org/10.3389/fmicb.2020.608832
              "Nitrosoabyssus", "Nitrosymbiomonas", # doi: 10.1186/s40168-025-02146-2
              "Nitrosopelagicus", # doi.org/10.1002/9781118960608.gbm01969
              "Nitrosopolaris", # https://doi.org/10.1093/femsmc/xtac019
              "Nitrosopumilales", "Nitrosopumilaceae",  # https://doi.org/10.1002/9781118960608.obm00122
              "Nitrososphaeraceae"   # yes, also as "Genus", ref from https://www.sciencedirect.com/science/article/pii/S0038071720300225
)

# All of Nitrosomonadaceae are AOB ... REF: https://link.springer.com/rwe/10.1007/978-3-642-30197-1_372
Nitrosomonadaceae <- Taxa.genus[ Taxa.genus$Family=="Nitrosomonadaceae" , "Genus"]
Nitrosomonadaceae <- Nitrosomonadaceae[!is.na(Nitrosomonadaceae)]
Nitrosomonadaceae <- Nitrosomonadaceae[Nitrosomonadaceae!="uncultured"]
                                       
# NOBs
NOB_list <- c("Nitrospira","Nitrospirae", "Nitrospirae", "Nitrospirales","Nitrospiraceae", "Nitrospirota",  # alcuni da MIDAS, ma tutti elencati in DOI: 10.1016/S0076-6879(11)86005-2
              # few Nitrospira names are quite similar, these are synonims or, at least, they still are of Nitrospirales order which is still known for AO activity, see https://www.nature.com/articles/s41396-023-01518-6
              "Nitromaritima",
              "Nitrospina","Nitrospinae","Nitrospinaceae","Nitrospinota",
              "Nitrobacter",
              "Nitrococcus",
              "Nitronereus" , # DOI: 10.1038/s41396-023-01518-6
              "Nitrohelix", "Nitronauta", "Nitrocaldera", "Nitrotheca", # doi: 10.1021/acs.est.3c00636
              "Nitrolancea", #doi: 10.1099/ijs.0.062232-0.
              "Nitrotoga")

{ top_data <- data.genus.prop
  vector_without_Ca <- gsub("Candidatus_","", tax_table( top_data )[ ,"Genus"])  # No "Candidatus" but same order, to a easier use of the logical below
  prune.dat_top <- subset_taxa(top_data , vector_without_Ca %in% c(AOB_list,AOA_list,NOB_list,Nitrosomonadaceae) )
  tax_selected<-as.data.frame(tax_table(prune.dat_top))
  if( identical(taxa_names(data.genus) , taxa_names(top_data)) ){
    tax_selected<-Taxa.genus.update[row.names(tax_selected),]
  } else {
    stop("Something went wrong when updating the taxonomy")
  }
  tax_table(prune.dat_top)<-as.matrix(tax_selected)
  # prune.dat_top <- tax_glom(prune.dat_top, taxrank = "Genus", NArm=T ) 
  tabella_top<-psmelt(prune.dat_top)
  tabella<-tabella_top[order(sort(tabella_top$Abundance, decreasing = T)), ]
  tabella$Genus2 <- gsub("Candidatus_","", tabella$Genus, fixed = T)
  
  tabella$Genus[tabella$Genus2 %in% c(AOB_list,Nitrosomonadaceae)] <- paste(tabella$Genus[tabella$Genus2 %in% c(AOB_list,Nitrosomonadaceae)], "(AOB)" )
  tabella$Genus[tabella$Genus2 %in% c(AOA_list)] <- paste(tabella$Genus[tabella$Genus2 %in% c(AOA_list)], "(AOA)" )
  tabella$Genus[tabella$Genus2 %in% c(NOB_list)] <- paste(tabella$Genus[tabella$Genus2 %in% c(NOB_list)], "(NOB)  " )
  # Adding AOB/NOB group column
  tabella$AOBNOB <- rep("none")
  tabella$AOBNOB[tabella$Genus2 %in% c(AOB_list,Nitrosomonadaceae)] <- "AOB"
  tabella$AOBNOB[tabella$Genus2 %in% c(AOA_list)] <- "AOA"
  tabella$AOBNOB[tabella$Genus2 %in% c(NOB_list)] <- "NOB"
  
  tabella$Genus<-gsub ("Candidatus_","Ca. ", tabella$Genus)
  tabella$Genus<-gsub ("uncultured","uncult. ", tabella$Genus)
  tabella$Genus<-factor(tabella$Genus, levels = c(unique(tabella$Genus)[unique(tabella$Genus)!="Others"],"Others"))
}
tabella$Sample_name<- factor(tabella$Sample_name, levels = c("Compost","R1","R2","R3","R4",
                                                             "S1","S2","S3","S4",
                                                             "T1","T2","T3","T4")
)
custom_colors<-c("gray30","deepskyblue","navyblue",
                 "royalblue1","red2","cadetblue3",
                 "tomato","black","skyblue1",
                 "slateblue2","cyan","deepskyblue3" )
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=Genus)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values=custom_colors) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0.01,0.01) , breaks = seq(0.5,10,0.5), limits = c(0, 7 ) ) +
  theme(axis.text.x=element_text(angle=20,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  panel.grid.major.y = element_line( colour="gray", linewidth=0.15 ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.15, "cm"),
  legend.text = element_text ( size = 6.85 ),
  legend.position="bottom",
  legend.margin = margin(-13 ,33, 3 ,5),
  plot.margin = margin(3,1,1,1)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="",
       caption = " AOA and AOB (Archaea and Bacteria oxidizing NH4+), NOB (bacteria oxidizing NO2-) ")
ggsave(file="Results/AOA_AOB_NOB_genera.png",width=4.5,height=3.5, dpi=300)

to_save <- round(otu_table(prune.dat_top), 2)
temp_meta <- metadata
row.names(temp_meta) <- temp_meta$FASTQ_ID
colnames(to_save) <- temp_meta[colnames(to_save) ,  "Sample_name"]
write.xlsx(cbind.data.frame( to_save, tax_table(prune.dat_top)[,c("Kingdom","Family","Genus")] ),
           file="Results/AOA_AOB_NOB_PERCENTAGES_AND_TAXONOMY.xlsx"
           )


# Again, only AOA, AOB and NOB groups
ggplot(data=tabella, aes(x=Sample_name, y=Abundance, fill=AOBNOB)) +
  geom_bar( stat="identity", position="stack", na.rm = F) + 
  scale_fill_manual(values= c("darkblue","royalblue2","red" ) ) +
  facet_nested(~ Group, scales = "free_x", space="free") +
  theme_classic(base_size =9) +
  scale_y_continuous(expand=c(0.01,0.01) , breaks = seq(0.5,10,0.5), limits = c(0, 7 ) ) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 8.2
  ),
  panel.grid.major.y = element_line( colour="gray", linewidth=0.15 ),
  axis.text.y=element_text(size=7),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(vjust=6.5),
  axis.title =element_text(size=10),
  strip.text = element_text(size=11),
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.5, "cm"),
  legend.text = element_text ( size = 9.5 ),
  legend.position="bottom",
  legend.margin = margin(-15 ,5, 3 ,5),
  plot.margin = margin(3,1,1,1)) +
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Percentual abundance of clades",
       #title = "Most abundant identified phyla (every sample)",
       fill="")
ggsave(file="Results/AOA_AOB_NOB_genera2.png",width=4.5,height=3.8, dpi=300)




################### SEARCHING FOR COMMON TRENDS ######################

dir.create("Results/Gradual_increases_along_samples")

this_table <- as.data.frame( as( otu_table(data.genus.prop) , "matrix" ) )
head(this_table, n=2)
this_table <- round(this_table, digits = 2)
# NB: yes, rounding to two digits infos regarding extremely rare taxa are lose (e.g. 0.004 --> 0) , but this is actually good (if they are THAT rare, then let's consider them as zero)

meta2 <- metadata
row.names(meta2) <- metadata$FASTQ_ID
actual_names <- meta2[ colnames(this_table) , "Sample_name" ]
colnames(this_table) <- actual_names

Table_R <- this_table[ , c("R1","R2","R3","R4")]
Table_S <- this_table[ , c("S1","S2","S3","S4")]
Table_T <- this_table[ , c("T1","T2","T3","T4")]

# fill_color_thrends<-c("navyblue","firebrick3","springgreen2","blue1","cadetblue2","deeppink2",
#                       "wheat","coral2","cyan","slateblue","black", "yellow3","darkmagenta","indianred3","magenta",
#                       "darkred","pink3", "gray90","yellow","darkgreen","violet","gray30","chocolate4",
#                       "chartreuse3","deepskyblue1","green","red1","royalblue2","darkslategray3","slategray"
#                       )
fill_color_thrends<-c("firebrick3","springgreen2","blue1","deeppink2",
                      "wheat","coral","cyan","slateblue1","black","darkmagenta",
                      "darkred","pink3", "gold","darkgreen","violet","gray30","chocolate4",
                      "chartreuse2","deepskyblue1","green","red1","royalblue2","darkslategray3","slategray2"
                      )

to_loop <- list( Table_R, Table_S, Table_T)
names(to_loop) <- c("Group_R","Group_S","Group_T")

common_genera<-NULL
for( i in 1:length(to_loop)){
  # i = 1
  x <- to_loop[[ i ]]
  increasing_rows_logical <- x[,1] < x[,2]  &  x[,2] < x[,3]  &  x[,3] <= x[,4]   # NB: only in the last comparison being equal is allowed ...
  x <- x[ increasing_rows_logical , ]
  x$Genera <- Taxa.genus.update[ row.names(x) , "Genus" ]
  x <- x[ !grepl("NA_ o NA", x$Genera) , ]
  x$Genera <- gsub("aceae$","ac.", x$Genera)
  x$Genera <- gsub("ales$","al.", x$Genera)
  x$Genera <- gsub("Candidatus","Ca.", x$Genera)
  common_genera <- c(common_genera, x$Genera)
  
  to_plot <- melt(x, id.vars = "Genera" )
  ggplot( data=to_plot , aes(x =variable, y= value , color= Genera ) ) +
    theme_bw() +
    scale_color_manual(values = fill_color_thrends ) +
    theme(   legend.text = element_text ( size = 8 , face="italic"),
             legend.position="bottom",
             legend.margin = margin(-15 ,35, 1 ,5),
             legend.spacing.y = unit(0.1,"cm"),
             legend.key.height = unit(0.28, "cm"),
             legend.key.width = unit(0.25, "cm"),
             panel.grid.major.y = element_line(colour="darkgray", linewidth = 0.075),
             panel.grid.minor.y = element_line(colour="gray", linewidth = 0.03),
    ) +
    scale_x_discrete( expand = c(0,0.05) ) +
    scale_y_continuous( breaks = seq(0,max(to_plot$value)+0.5,0.5) ) +
    geom_point(size=2, alpha=0.95) +
    # geom_point(size=0.5) +
    geom_line( aes(group=Genera, color=Genera) , linewidth=0.035 , show.legend = F ) +
    guides(color=guide_legend(nrow=6)) +
    labs( x= "",y="Percent abundance" , fill="" , color="" )
  
  ggsave(filename = paste("Results/Gradual_increases_along_samples/",names(to_loop)[i],".png"),
         width = 5.5, height = 4, dpi= 300 )
}



### Common ones only ...
common_genera <- gsub("ac.","aceae", common_genera, fixed = T )
common_genera <- gsub("al.","ales", common_genera, fixed = T )
common_in_two <- common_genera[duplicated(common_genera)]
common_in_three <- common_in_two[duplicated(common_in_two)]

con <- file("Results/Gradual_increases_along_samples/Spiegazione_E_dettagli.txt")
sink(con, append = TRUE)
# NB: THE FOLLOWING LINES ARE WRITTEN IN ITALIAN TO DIRECTLY TALK WITH MY COLLABORATOURS
cat("Sono stati cercati tutti i batteri con abbondanza percentuale crescente lungo i tre campioni nei tre gruppi")
cat("\n\nNel dettaglio, questi batteri qui sono risultati in almeno 2 gruppi ...", fill = T)
cat(common_in_two, fill = T)
cat("\n\nQuesti invece sono risultati in comune tra tutti e tre i gruppi ...", fill = T)
cat(common_in_three)
sink()
close(con)
suppressWarnings(rm(con))




##################### R and PACKAGES VERSION #########################

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
package$otherPkgs$mixOmics[c(1,4)]
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
