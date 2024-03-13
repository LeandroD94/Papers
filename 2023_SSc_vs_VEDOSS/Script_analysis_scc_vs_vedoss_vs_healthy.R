##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
library("ggplot2")
library("vegan")
library("ggpubr")
library("ggh4x") #
library("ggvenn") #
# library("egg")
library("dendextend")
library("DESeq2")
library("mixOmics") #
library("xlsx")  
library("dplyr")
library("stringr")
library("qiime2R")
}

{dir.create("Data_check")
dir.create("Data_check/PCoA_test")
dir.create("Results")
dir.create("Results/Abundances")
dir.create("Results/Abundances/Ratio_Firmi_Bacteroi")
dir.create("Results/Hierarchical_clustering")
dir.create("Results/Beta_div")
dir.create("Results/DA_DESeq2/")
dir.create("Results/PICRUST2_LEFSE_results")
dir.create("Results/PLS_DA")
}

options(scipen = 100) # disable scientific annotation



####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample

Metadata <- as.data.frame(read.csv("metadata.R.csv"))
row.names(Metadata)<-str_split(Metadata$FASTQ_code, pattern = "_", simplify = T)[,1] # to match phyloseq codes
head(Metadata)
original_length<-length(Metadata$FASTQ_code[!is.na(Metadata$FASTQ_code)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ_code[!is.na(Metadata$FASTQ_code)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,original_names,sample)

sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels=c("SSc","VEDOSS","Healthy")) # decides the order in plots
sample_names(data)<-sample_data(data)$Sample_name
head(sample_data(data))

# save.image("data_scleroderm_2022_pre_filters.RData")


#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.01% to remove noises/contaminants, too conservative but also safe cutoff, see  PMID: 31164452,  PMID:23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_001_cutoff.csv")

filtered<-as.data.frame(tax_table(filter_taxa(data.genus.temp, function(x) mean(x) <= 0.005, TRUE)))[["Genus"]]
filtered
filtered<-filtered[filtered!="uncultured" & !is.na(filtered)] # to avoid the removal of other uncultured genera
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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Unassigned_domain_checking_pre_filtering.csv", row.names = T, quote = F)
### now filtering out every Unassigned ASV
data<- subset_taxa(data, Kingdom != "d__Eukaryota")
head(tax_table(data))
# write.csv2(tax_table(data), file="Data_check/Every_filtered_ASV_and_taxonomic_assignment.csv", row.names = T)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


### removing chloroplast and mitochondria
if( "Chloroplast" %in% as.data.frame(tax_table(data))[["Genus"]] | "Mitochondria" %in% as.data.frame(tax_table(data))[["Genus"]] ){
  cat("\nRemoving Chloroplast and/or Mitochondria\n\n")
  data<-subset_taxa(data, ! Genus %in% c("Chloroplast","Mitochondria"))
} # NO Chloroplast or Mitochondria!


proof1<- "Marker of the filtering, it is required for the script"

# save.image("data_phyloseq_after_filters.RData")


############################ RAREFACTION ANALYSIS ################################

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
r<-rarecurve(as(t(otu_table(data)),"matrix"), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)


####################### PREPARATION OF THE DATA #######################

if(! "proof1" %in% ls()){
  stop("\n Wait! Did you perform the filtering steps??? \n\n")
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
assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
}
assigned
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)
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

write.csv2(cbind(as(otu_table(unfiltered_data),"matrix"),as(tax_table(unfiltered_data),"matrix")),file="Data_check/OTU_Unfiltered_abundances_and_tax.csv",quote=F)

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

### choosing colors  (see grDevices::colors() )
fill_color_5<-c("coral","springgreen3","gold3","firebrick3","deepskyblue3","darkslategray3") # "others" will be setted as the last one
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","violet","deepskyblue2", "darkslategray3") # "others" will be setted as the last one

### TOP 5 Phyla
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
  facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 12.2 ),
        legend.position="bottom",
        #legend.margin = margin(0,0,0,0),
  ) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Percentual abundance", 
       title = "Five most abundant phyla", 
       caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Abundances/TOP_5_phyla.png",width=8,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top))[["Phylum"]]))

rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others)


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
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) +
  facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") +
  theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_10) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10), 
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text ( size = 11.8 ),
        legend.position="bottom", 
        legend.margin = margin(3,0,0,0)
  ) + 
  guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Percentual abundance",
       title = "Ten most abundant genera",
       caption = " 'Others' includes every genus below rank 10 ")
ggsave(file="Results/Abundances/TOP_10_genera.png",width=8,height=5,dpi=300)
dev.off()

# means of TOP10 genera
write.xlsx(file = "Results/Abundances/TOP_10_genera_Average_abundances.xlsx", row.names = F,
           cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top),1,mean)), "genus"= as.data.frame(tax_table(prune.dat_top))[["Genus"]]))

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))



#################### RATIO FIRMICUTES/BACTEROIDES ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-cbind.data.frame(otu_table(data_fb),tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) )
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

ratio_fb<-ratio_fb[Metadata$Sample_name, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Metadata[,c("Sample_name","Condition")])

# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Condition, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Condition=="VEDOSS","Ratios_Mean"]<-Ratios_Mean["VEDOSS"]
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_Mean"]<-Ratios_Mean["Healthy"]
ratio_fb[ratio_fb$Condition=="SSc","Ratios_Mean"]<-Ratios_Mean["SSc"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Condition, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Condition=="VEDOSS","Ratios_st_err"]<-Ratios_st["VEDOSS"]/sqrt(length(which(ratio_fb$Condition=="VEDOSS")))
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_st_err"]<-Ratios_st["Healthy"]/sqrt(length(which(ratio_fb$Condition=="Healthy")))
ratio_fb[ratio_fb$Condition=="SSc","Ratios_st_err"]<-Ratios_st["SSc"]/sqrt(length(which(ratio_fb$Condition=="SSc")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/Abundances/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio.csv", ratio_fb)

ratio_fb$Condition<-factor(ratio_fb$Condition, levels = c("Healthy","VEDOSS","SSc"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-kruskal.test(Ratio ~ Condition, data=ratio_fb)
p_val<-round(res$p.value, digits = 2)
p_val

DT<-FSA::dunnTest(Ratio ~ Condition, data=ratio_fb, method="bh")
DT

# exporting the results
con <- file("Results/Abundances/Ratio_Firmi_Bacteroi/Statistical_tests_FB_ratios.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Kruskal Wallis test \n", fill=T)
cat("Ratio~Condition :   V=",res$statistic, "p-value=", p_val, "\n",fill=T) # da adattare al kruskal
cat("\n\n Dunnet test (p-value adjusted through BH)\n", fill=T)
print(DT$res)
cat("\n\n\n Shapiro Wilk test p-value: ",distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=Sample_name, fill=Condition, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Healthy"="chartreuse",
                               "VEDOSS"="coral",
                               "SSc"="red3")) +
  facet_grid2(.~Condition, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Condition"))+
  scale_y_sqrt(breaks=c(1,5,seq(10,80,10)) ) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Kruskal Wallis p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi/Ratio_barplot.png",width=8,height=4, dpi=300) 
dev.off()



##### BOXPLOT
ggplot(data=ratio_fb, aes(x=Condition, fill=Condition,y=Ratio)) +
  theme_classic2(base_size =9) +
  scale_fill_manual(values = c("VEDOSS"="coral","Healthy"="chartreuse","SSc"="red3")) +
  scale_color_manual(values = c("VEDOSS"="coral","Healthy"="chartreuse2","SSc"="red3")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(0.5,1,1.5,3,5,7.5,seq(10,80,10)) ) +
  geom_boxplot(width=0.9, size= 0.4, alpha= 0.05, outlier.alpha = 0, show.legend = F) +
  geom_point(position = position_jitterdodge(seed = 1994, jitter.width = 0.8),
             aes(color=Condition), size= 0.7, alpha= 0.65, show.legend = F) + 
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6),
  panel.grid.major.y = element_line(size=0.3)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", caption = paste("Kruskal Wallis p-value:",p_val) )
ggsave(file="Results/Abundances/Ratio_Firmi_Bacteroi/Ratios_FB_boxplot_median_IQR.png",width=4,height=3.5, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))



######################## HIERARCHICAL CLUSTERING ###################

suppressWarnings(rm(c))

#euclidean
{c<-hclust(dist(t(sqrt(otu_table(data.prop)))))
  c<-as.dendrogram(c)
  color_table<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
  colnames(color_table)<-c("Gruppo","Colore")
  colors <- gsub("Healthy","chartreuse",color_table$Colore) #green
  colors <- gsub("VEDOSS","coral",colors) #orange
  colors <- gsub("SSc","red3",colors) # red
  labels_colors(c) <- colors[order.dendrogram(c)] # but sort them based on their order in dendrogram
}
png(file="Results/Hierarchical_clustering/Hierarchical_cluster_Hellinger_ASVs.png",width=3000,height=1800, res=300)
par(mar=c(8,4,4,3), cex.lab=0.4, cex.main=1.4, cex.sub=1.3)
plot(c, main="Community structure computed with Hellinger distance \n (euclidean on sqrt proportional ASVs)", 
     sub="Healthy = Green     VEDOSS = Orange     SSc = Red")
dev.off()

# Bray Curtis
{c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method="bray"))
  c<-as.dendrogram(c)
  color_table<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
  colnames(color_table)<-c("Gruppo","Colore")
  colors <- gsub("Healthy","chartreuse",color_table$Colore) #green
  colors <- gsub("VEDOSS","coral",colors) #orange
  colors <- gsub("SSc","red3",colors) #blue
  labels_colors(c) <- colors[order.dendrogram(c)] # but sort them based on their order in dendrogram
}
png(file="Results/Hierarchical_clustering/Hierarchical_cluster_Bray_Curtis_sqrt_prop_ASVs.png",width=3000,height=1800, res=300)
par(mar=c(8,4,3,4), cex.lab=0.4, cex.main=1.4, cex.sub=1.3)
plot(c, cex.lab= 0.2, main="Community structure computed with Bray-Curtis distance \n on sqrt proportional ASVs", 
     sub="Healthy = Green     VEDOSS = Orange     SSc = Red")
dev.off()

suppressWarnings(rm(c,color_table,labels_colors))


########################## ALFA DIVERSITY ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
# moreover, the alpha diversity is calculated genera (DADA2 overestimates the ASV richness, see doi.org/10.1186/s40168-020-00856-3 )

{pAlpha<-plot_richness(data.genus, measures=c("Shannon", "Observed"), x="Condition", color="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
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
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) +
  theme_classic2(base_size = 11) + 
  scale_color_manual(values = c("VEDOSS"="coral","Healthy"="chartreuse","SSc"="red3")) +
  labs(x="Condition",
       title="Alpha diversity in Healthy, VEDOSS and SSc") +
  guides(fill="none", color="none") +
  stat_compare_means(aes(group = Condition, 
                         label = sprintf("p = %.1g", as.numeric(..p.format..)) # g allows scientific annotation
                         ),
                     method = "kruskal.test", label.x= 0.8, size=3.2, label.y.npc = "top", vjust=-0.35, hjust=-0.2) +
  theme(panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.2) , 
        axis.text.x= element_text(angle=28, vjust=1, hjust=1, size=10),
        axis.title.x = element_text(vjust = -1)
        )
ggsave(file="Results/Alpha_diversity/Alfa_diversity_with_Kruskal.png", width = 5,height =4.2, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data)
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
kruskal.test(Obser_value$value~factor)


# pairwise comparisons
con <- file("Results/Alpha_diversity/Alpha_diversity_groups_pairwise_comparisons.txt")
sink(con, append = TRUE)
cat("Dunn test (Kruskal-Wallis post Hoc) on Alpha diversity Condition sub groups \n p-value corrections done with Benjamini-Hochberg method \n \n", fill=TRUE)

Filter_value<-filter(alphadt, variable=="Observed richness")
a<-FSA::dunnTest(value ~ Condition, data=Filter_value, method="bh")
cat("pairwise observed", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Shannon")
a<-FSA::dunnTest(value ~ Condition, data=Filter_value, method="bh")
cat("\n \n pairwise Shannon", "\n", fill=TRUE)
a
rm(a)

Filter_value<-filter(alphadt, variable=="Evenness")
a<-FSA::dunnTest(value ~ Condition, data=Filter_value, method="bh")
cat("\n \n pairwise Evenness", "\n", fill=TRUE)
a
rm(a)

sink()
close(con)
rm(con, pAlpha, alphadt,H, ev, obs, Obser_value, New_data, factor)


######################## BETA DIVERSITY #######################

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
  perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
  perm_ASV$aov.tab$`Pr(>F)`[1]
  perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
  perm_g<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
  perm_g$aov.tab$`Pr(>F)`[1]
  perm_g_H<-perm_g # needed later for plot
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
  perm_f<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
  perm_o<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
  perm_c<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
}

{sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
  perm_p<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
}

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta


# PLUS: checking it with Bray too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="bray")
perm_ASV$aov.tab
# write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Beta_div/Beta_divers_permanova_BRAY_on_ASV.csv",quote=F,row.names = T)


### pairwise beta diversity
# devtools install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sample_ASV<-as.data.frame(t(sqrt(ASV.genus.prop))) # samples has to be on rows --> t
pair_ASV<- pairwise.adonis(sample_ASV, factors=metadata$Condition, p.adjust.m = "BH", sim.method="euclidean", perm = 9999)
pair_ASV

# exporting beta diversity
suppressWarnings(rm(con))
con<-file("Results/Beta_div/Beta_divers_general_and_pairwise_between_groups.txt")
sink(con, append = TRUE)
cat("General beta diversity on Hellinger \n")
beta
cat("\n \n", fill=TRUE)
cat("Pairwise beta diversity on Hellinger (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
pair_ASV
sink()
close(con)

rm(beta, pair_ASV, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)



# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on GENERA
{BC.dist<-vegan::vegdist(t(sqrt(ASV.genus.prop)), distance="euclidean")
  disper<-vegan::betadisper(BC.dist,metadata$Condition)
  disp_ASV<-vegan::permutest(disper, pairwise=TRUE, permutations=9999)
  disp_ASV$tab
  a<-as.data.frame(disp_ASV$pairwise$permuted)
  colnames(a)<-c("permuted_p_value")
  a$padj_BH<-p.adjust(a$permuted_p_value, method = "BH")
}
a

#export dispersion
suppressWarnings(rm(con))
con<-file("Results/Beta_div/Beta_dispersion_General_and_Pairwise_between_groups.txt")
sink(con, append=TRUE)
cat("General beta dispersion on Hellinger \n")
disp_ASV$tab
cat("\n \n", fill=TRUE)
cat("Pairwise beta dispersion (correction Benjamini-Hochberg on p-value) \n \n", fill=TRUE)
a
sink()
close(con)

rm(disp_ASV,a, con)


####################### PCoA BETA DIV #########################

# on genera
data.prop.labels<-data.genus.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2)  +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)", color="Condition",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  geom_point(size=3, alpha= 0.5) +
  theme_classic(base_size = 10.5) +
  theme(title=element_text(size=10),
        legend.title =element_text(size=10.5)) +
  stat_ellipse(size=0.15)  +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera_points.png", width = 5.2, height = 4.2, dpi=300)

# Genera in 3D --> it is just a bunch of points, it is not readable nor useful

### THE FOLLOWING IS A CHECK: VEDOSS and SSc are actually different?
data.prop.labels<-subset_samples(data.genus.prop, Condition!="Healthy")
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("VEDOSS"="coral","SSc"="red3")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2)  +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance\n (euclidean on Hellinger transformed genera)",
       color="Condition",
       x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption = "This plot is to further check the difference between VEDOSS and SSc")
ggsave(file="Results/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 6, height = 4.5, dpi=300)


##################### VEEN DIAGRAM ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.1% abundance at least in two sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.1, 1, 0)) # if more than 0.1 (then it is not a noise!) --> "a point"
who<-who[!rowSums(who)>6,] # more than 6 "points" --> at least in 7 samples (10% of whole dataset)
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
who<-who[!is.na(who)] # this is to avoid the NA, because there are also "good" NA
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Healthy<-subset_samples(data.venn, Condition=="Healthy")
Healthy<-as.character(tax_table(prune_taxa(taxa_sums(Healthy)>0, Healthy))[,"Genus"])

VEDOSS<-subset_samples(data.venn, Condition=="VEDOSS")
VEDOSS<-as.character(tax_table(prune_taxa(taxa_sums(VEDOSS)>0, VEDOSS))[,"Genus"])

SSc<-subset_samples(data.venn, Condition=="SSc")
SSc<-as.character(tax_table(prune_taxa(taxa_sums(SSc)>0, SSc))[,"Genus"])


ONLY_IN_VEDOSS<- VEDOSS[! VEDOSS %in% SSc & ! VEDOSS %in% Healthy]
ONLY_IN_VEDOSS<- paste(ONLY_IN_VEDOSS, collapse = ", ")
head(ONLY_IN_VEDOSS)

ONLY_IN_Healthy<- Healthy[! Healthy %in% SSc & ! Healthy %in% VEDOSS]
ONLY_IN_Healthy<- paste(ONLY_IN_Healthy, collapse = ", ")
head(ONLY_IN_Healthy)

ONLY_IN_SSc<- SSc[! SSc %in% Healthy & ! SSc %in% VEDOSS]
ONLY_IN_SSc<- paste(ONLY_IN_SSc, collapse = ", ")
head(ONLY_IN_SSc)

IN_HEALTHY_AND_SSc_BUT_NOT_VEDOSS <- Healthy[Healthy %in% SSc & ! Healthy %in% VEDOSS ] 
head(IN_HEALTHY_AND_SSc_BUT_NOT_VEDOSS)

IN_SSc_AND_VEDOSS_BUT_NOT_HEALTHY <- SSc[SSc %in% VEDOSS & ! SSc %in% Healthy ] 
head(IN_SSc_AND_VEDOSS_BUT_NOT_HEALTHY)


con<-file("Results/Beta_div/Venn_Diagr_exclusive_bacteria.txt")
{sink(con, append=TRUE)
cat("ONLY IN SSc", fill=TRUE)
cat(ONLY_IN_SSc)
cat("\n\nONLY IN Healthy", fill=TRUE)
cat(ONLY_IN_Healthy)
cat("\n\nONLY IN VEDOSS", fill=TRUE)
cat(ONLY_IN_VEDOSS)
cat("\n\n\nIn SSc and VEDOSS but not in Healthy ", fill=TRUE)
cat(IN_SSc_AND_VEDOSS_BUT_NOT_HEALTHY)
cat("\n\n\nIn Healthy and SSc but not in VEDOSS ", fill=TRUE)
cat(IN_HEALTHY_AND_SSc_BUT_NOT_VEDOSS)
cat("\n\n\n NB: the bacteria considered for those intersections are only those with a relative abundance over 0.1% at least in seven sample (10% of whole dataset) to avoid the 'noisier' ones")
sink()
}
close(con)


x<-list(Healthy=Healthy,VEDOSS=VEDOSS,SSc=SSc)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("chartreuse","coral","red3")) +
  theme(plot.title = element_text(size=9), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram \n including only genera with minimal abundance > 0.1% \n at least in seven samples (10% of whole dataset)")
ggsave(filename = "Results/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()

rm(ONLY_IN_SSc, ONLY_IN_Healthy, ONLY_IN_VEDOSS, IN_HEALTHY_AND_SSc_BUT_NOT_VEDOSS,
   IN_SSc_AND_VEDOSS_BUT_NOT_HEALTHY, x, con, VEDOSS, Healthy, SSc, data.venn, who)



##################### DA WITH DESEQ2 #######################

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
auto_DESeq_groups<-function(data_pruned,ranks,nume,deno,outfile,fct=1,pt=0.05) {
  DE_results<-NULL
  for (rank in ranks) {
    cat(" WORKING ON",rank,"\n")
    if(rank == "OTU") {
      ori<-data_pruned
      d<-phyloseq_to_deseq2(data_pruned, ~Condition)
    } else {
      ori<-tax_glom(data_pruned, taxrank = rank, NArm = F)
      d<-phyloseq_to_deseq2(ori, ~Condition)
    }
    # DE<-DESeq(d, test="LRT",reduced= ~ 1)
    DE<-DESeq(d)
    for (i in 1:length(deno)) {
      cat("Test of",nume[i],"/",deno[i],":")
      res<-results(DE, contrast=c("Condition", nume[i], deno[i]), test="Wald")
      #res<-results(DE, contrast=c("Condition", nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      #res<-lfcShrink(DE, type="normal", contrast=c("Condition",nume[i], deno[i]), test="Wald", lfcThreshold = fct, alpha= pt)
      sel<-!is.na(res$padj) & res$padj<=pt & abs(res$log2FoldChange)>=fct
      rnum<-sum(ifelse(sel,1,0))
      cat(rnum,"\n")
      if (rnum>0) {
        res_s<-res[sel,]
        ann<-tax_table(ori)[sel,]
        res_s_ann<-cbind( "num"=nume[i],"den"=deno[i],"Rank"=rank,res[sel,],data.frame(as(tax_table(ori)[sel,],"matrix")),row.names(tax_table(ann)) )
        #			cat(paste0(rank,":",nume[i],"/",dena[i],"\n"),file=outfile,append=T)
        # write.table(as.data.frame(res_s_ann),sep="\t",file=outfile,append=T,col.names=F)
        DE_results<-rbind(DE_results, res_s_ann)
      }
    }
  }
  DE_results<<-as.data.frame(DE_results) # superassign into global environment
}

nume<-c("VEDOSS","Healthy","VEDOSS")
deno<-c("Healthy","SSc","SSc")                              
ranks<-c("Phylum","Class","Order","Family","Genus")
options(scipen = 100)
auto_DESeq_groups(data_pruned,ranks,nume,deno,"DE_results.txt",1.5,0.01)   #(automatically computes all DESeq2 analyses)
# NB: fold change threshold is at least 1.5!
DE_results<-DE_results[,colnames(DE_results)!="Species"]
colnames(DE_results)<- c("num","denom","Rank","BaseMean","log2FoldChange","lfcSE","stat","pvalue","p-adj","Domain","Phylum","Class","Order","Family","Genus","ASV")
DE_results<-DE_results[DE_results$BaseMean>100, ]
system(" echo 'The threshold of log2(FC) of the algorithm DESeq2 has been set to be more than 1.5 , \n\n Moreover, every result under the arbitrary threshold Basemean=100 has been removed (AFTER the p-value adjustment) in order to avoid the noisiest results' > Results/DA_DESeq2/NB_results_are_filtered.txt")
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

tabella_tot<-rbind.data.frame(tabella_g,tabella_f,tabella_c,tabella_o)

ggplot(tabella_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  geom_boxplot(width=0.8) + theme_bw( ) + 
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  theme(axis.text.x = element_text(angle = 35, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=10),
        legend.margin=margin(-15, 0, 0, 0), legend.position="bottom") +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa", y="Proportional Abundance %", fill="Condition", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,2),seq(22,max(tabella_tot$Abundance),3)))
ggsave(filename = "Results/DA_DESeq2/Results_DeSeq2_Condition.png", width = 13, height = 7, dpi=300)
dev.off()

# removing redundant results (NB: there are even more names below because initially has been tested the "permissive" logFC>1)
Redund <- c("Bacteroidia","Acidaminococcales","Bacteroidaceae","Rhodospirillales","Tannerellaceae",
            "Streptococcaceae","Peptostreptococcales-Tissierellales",
            "Erysipelatoclostridiaceae","Desulfovibrionales","Burkholderiales")
tabella_tot2<-subset(tabella_tot, ! Bacteria %in% Redund)
tabella_tot2<-tabella_tot2[ ! (tabella_tot2$Bacteria=="[Eubacterium]_coprostanoligenes" & tabella_tot2$Taxa=="Families"), ]
tabella_tot2<-tabella_tot2[ ! (tabella_tot2$Bacteria=="Clostridia_UCG-014" & tabella_tot2$Taxa=="Orders"), ]
tabella_tot2$Bacteria<-gsub("_group","",tabella_tot2$Bacteria)

ggplot(tabella_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), 
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  scale_color_manual(values=c("Healthy"="chartreuse2","VEDOSS"="coral","SSc"="red3")) +
  geom_boxplot(width=0.7, size= 0.4, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.7,jitter.width = 0.32),
             aes(color=Condition), size= 0.62, alpha= 0.6) +  
  theme_classic2(base_size = 11) +
  theme(strip.text.x=element_text(size=12, colour="black"),
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=9.5),
        axis.text.y = element_text(size=8.3),
        legend.margin=margin(-26, 0, 0, 0), 
        legend.position="bottom",
        plot.title= element_text(size=12.5) ,
        plot.margin = margin(5,1,5,28),
        legend.box.margin = margin(5,0,-3,0),
        legend.key.size=unit(0.65,"cm"), 
        legend.text=element_text(size=12),
        panel.grid.minor.y= element_blank(),
        panel.grid.major.y = element_line(size=0.10, color="gray")
  ) + 
  guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa (every result)",
       y="Percentual Abundance",
       fill="Condition", x="") +
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 1, 2.5, 5, 8, seq(12,22,4),seq(25,max(tabella_tot$Abundance)+5,3)))
ggsave(filename = "Results/DA_DESeq2/Every_result_no_redundants2.png", width = 10, height = 6, dpi=300)
dev.off()


################# now plotting result in pairwise (group vs group)

##### VEDOSS vs Healthy
tabella_tot3<-subset(tabella_tot2, Condition != "SSc")
Healthy_vs_VEDOSS<-DE_results[DE_results$num=="VEDOSS" & DE_results$denom=="Healthy", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% Healthy_vs_VEDOSS [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, Healthy_vs_VEDOSS[Healthy_vs_VEDOSS$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral")) +
  scale_color_manual(values=c("Healthy"="chartreuse2","VEDOSS"="coral")) +
  geom_boxplot(width=0.9, size= 0.4, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Condition), size= 0.8, alpha= 0.65) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=11.7,colour="black"),
        legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=8.8),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between Healthy and VEDOSS", y="Percentual Abundance", fill="Condition", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,2),seq(22,max(tabella_tot$Abundance)+4,3)))
ggsave(filename = "Results/DA_DESeq2/Results_only_Healthy_vs_VEDOSS_no_redundants.png", width = 7.5, height = 6, dpi=300)



#####  Healthy vs SSc
tabella_tot3<-subset(tabella_tot2, Condition != "VEDOSS")
SSc_vs_Healthy<-DE_results[DE_results$num=="Healthy" & DE_results$denom=="SSc", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% SSc_vs_Healthy [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, SSc_vs_Healthy[SSc_vs_Healthy$Rank==rank, rank])
}
tot_Bacteria<-gsub("_group","",tot_Bacteria) # if "group" then it must be removed also from "tabella_tot" then also in tot_bacteria
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

tabella_tot3$Bacteria<-gsub("[","",tabella_tot3$Bacteria, fixed=T)
tabella_tot3$Bacteria<-gsub("]","",tabella_tot3$Bacteria, fixed=T)
ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("Healthy"="chartreuse","SSc"="red3")) +
  scale_color_manual(values=c("Healthy"="chartreuse2","SSc"="red3")) +
  geom_boxplot(width=0.9, size= 0.4, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Condition), size= 0.8, alpha= 0.65) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=12.5,colour="black")) + 
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=9.3),
        plot.margin = margin(5,1,5,50),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa between Healthy and SSc", y="Percentual Abundance", fill="Condition", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,2)))
ggsave(filename = "Results/DA_DESeq2/Results_only_Healthy_vs_SSc_DeSeq2_no redundants.png", width = 8, height = 6.5, dpi=300)



##### VEDOSS vs SSc
tabella_tot3<-subset(tabella_tot2, Condition != "Healthy")
SSc_vs_VEDOSS<-DE_results[DE_results$num=="VEDOSS" & DE_results$denom=="SSc", ]

tabella_tot3<-subset(tabella_tot3, OTU %in% SSc_vs_VEDOSS [,"ASV"])
tot_Bacteria<-NULL # just a further checking, in case of same ASV in higher levels
for(rank in c("Phylum","Class","Order","Family","Genus")){
  tot_Bacteria<-c(tot_Bacteria, SSc_vs_VEDOSS[SSc_vs_VEDOSS$Rank==rank, rank])
}
tabella_tot3<-subset(tabella_tot3, Bacteria %in% tot_Bacteria)

ggplot(tabella_tot3, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")),
             scales = "free_x", space="free") +
  scale_fill_manual(values=c("VEDOSS"="coral","SSc"="red3")) +
  scale_color_manual(values=c("VEDOSS"="coral","SSc"="red3")) +
  geom_boxplot(width=0.9, size= 0.4, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Condition), size= 0.8, alpha= 0.65) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  theme(legend.margin=margin(-10, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=9),
        plot.margin = margin(5,2,5,10),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray")
        ) +
  guides( fill=guide_legend(nrow=1), color="none" ) +
  labs(title= "Differently abundant Taxa between VEDOSS and SSc", y="Percentual Abundance", fill="Condition", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,2),seq(22,max(tabella_tot$Abundance)+3,3)))
ggsave(filename = "Results/DA_DESeq2/Results_only_VEDOSS_vs_SSc_DeSeq2_no redundants.png", width = 8, height = 6.5, dpi=300)



####################### CIRCOPLOT FOR DESeq2 ######################

DE_for_Circo<-DE_results

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

# removing redundants
for( t in c("Phylum","Class","Order","Family")){
  DE_for_Circo<-DE_for_Circo[ ! (DE_for_Circo[[t]] %in% Redund  & DE_for_Circo$Rank==t) , ]
}

### plotting VEDOSS vs Healthy
# re-formatting
temp<- DE_for_Circo[DE_for_Circo$num=="VEDOSS" & DE_for_Circo$denom=="Healthy", ]
temp<- temp [ ,! colnames(temp) %in% c("num","denom","ASV","Rank")]
temp<- cbind.data.frame(temp, Taxon=DE_for_Circo[DE_for_Circo$num=="VEDOSS" & DE_for_Circo$denom=="Healthy", "Rank"] ) # Rank (Taxon) as last column
# adjusting labels
temp$Genus <- gsub("-","_",temp$Genus, fixed=T)
temp$Genus <- gsub("[","",temp$Genus, fixed=T)
temp$Genus <- gsub("]","",temp$Genus, fixed=T)
temp$Genus <- gsub("_group","",temp$Genus, fixed=T)
circo<-CircoTax2(temp,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/DA_DESeq2/VEDOSS_vs_Healthy_CIRCOPLOT.png",width=1800,height=1800, res = 300)
circo + labs(title="VEDOSS vs Healthy \n DA taxa log2foldchange\n")
dev.off()


### plotting SSc vs Healthy
# re-formatting for Circoplot
temp<- DE_for_Circo[DE_for_Circo$num=="Healthy" & DE_for_Circo$denom=="SSc", ]
temp<- temp [ ,! colnames(temp) %in% c("num","denom","ASV","Rank")]
temp<- cbind.data.frame(temp, Taxon=DE_for_Circo[DE_for_Circo$num=="Healthy" & DE_for_Circo$denom=="SSc", "Rank"] ) # Rank (Taxon) as last column
# another redundant to remove
temp<- temp[! (temp$Family=="[Eubacterium]_coprostanoligenes_group" & temp$Taxon=="Family") , ]
# adjusting labels
temp$Genus <- gsub("-","_",temp$Genus, fixed=T)
temp$Genus <- gsub("[","",temp$Genus, fixed=T)
temp$Genus <- gsub("]","",temp$Genus, fixed=T)
temp$Genus <- gsub("_group","",temp$Genus, fixed=T)
circo<-CircoTax2(temp,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/DA_DESeq2/SSc_vs_Healthy_CIRCOPLOT.png",width=1800,height=1800, res = 300)
circo + labs(title="SSc vs Healthy \n DA taxa log2foldchange\n")
dev.off()


### plotting VEDOSS vs SSc
# re-formatting for Circoplot
temp<- DE_for_Circo[DE_for_Circo$num=="VEDOSS" & DE_for_Circo$denom=="SSc", ]
temp<- temp [ ,! colnames(temp) %in% c("num","denom","ASV","Rank")]
temp<- cbind.data.frame(temp, Taxon=DE_for_Circo[DE_for_Circo$num=="VEDOSS" & DE_for_Circo$denom=="SSc", "Rank"] ) # Rank (Taxon) as last column
# adjusting labels
temp$Genus <- gsub("-","_",temp$Genus, fixed=T)
temp$Genus <- gsub("[","",temp$Genus, fixed=T)
temp$Genus <- gsub("]","",temp$Genus, fixed=T)
temp$Genus <- gsub("_group","",temp$Genus, fixed=T)
circo<-CircoTax2(temp,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")
png(file="Results/DA_DESeq2/VEDOSS_vs_SSc_CIRCOPLOT.png",width=1800,height=1800, res = 300)
circo + labs(title="VEDOSS vs SSc \n DA taxa log2foldchange\n")
dev.off()


################# PLS-DA and sPLS-DA ###########################

# selecting three sample for each group (as test dataset)
{set.seed(1994)
  Healthy_train<-sample(sample_data(data.genus.prop)[sample_data(data.genus.prop)[["Condition"]]=="Healthy", ] [["Sample_name"]] )[1:3]
}
{set.seed(1994)
  VEDOSS_train<-sample(sample_data(data.genus.prop)[sample_data(data.genus.prop)[["Condition"]]=="VEDOSS", ] [["Sample_name"]] )[1:3]
}
{set.seed(1994)
  SSc_train<-sample(sample_data(data.genus.prop)[sample_data(data.genus.prop)[["Condition"]]=="SSc", ] [["Sample_name"]] )[1:3]
}

# subsetting test and train datasets
data.selected<-subset_samples(data.genus.prop, ! Sample_name %in% c(Healthy_train, VEDOSS_train, SSc_train) )
metadata_PLSDA<-as(sample_data(data.selected),"data.frame")
data.train<-subset_samples(data.genus.prop, Sample_name %in% c(Healthy_train, VEDOSS_train, SSc_train) )
metadata_train<-as(sample_data(data.train),"data.frame")

# starting PLSDA
X <- as.matrix(t(otu_table(data.selected)))
Y <- as.factor(metadata_PLSDA$Condition)

my.splsda <- splsda(X, Y, ncomp = 5) # test from 1 to 5 components

gc()
set.seed(1)
perf.splsda <- perf(my.splsda, validation = "Mfold", # to assest the optimal number of comp
                    folds = 10, nrepeat = 200, cpus = 8,
                    progressBar = TRUE, auc = TRUE)

perf.splsda$choice.ncomp
comp<-perf.splsda$choice.ncomp[4] # centroid distance (2)

png(filename = "Results/PLS_DA/Error for each component PLS-DA.png", width = 2500, height = 1800, res=300)
plot(perf.splsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
title(main="Classification rate error \n for each component used to compute PLS-DA")
dev.off()

my.splsda <- splsda(X, Y, ncomp = comp)

n= 2 # choose the component to use depending on clustering efficiency in plot
png(filename = paste("Results/PLS_DA/normal_PLS_DA_comp_1_and_", n,".png", sep=""), width = 2400, height = 1600, res=300)
plotIndiv(my.splsda, comp = c(1,n), col = c("red","coral","green"), cex = 1.8,
          group = metadata_PLSDA$Condition, ind.names = TRUE, ellipse = TRUE, legend = TRUE, 
          title = paste0("PLSDA among Healthy, VEDOSS and SSc \n on scaled proportional genus data \n (computed with ", comp, " components, plotted on comp 1 and ", n,")"))
dev.off()


################ sPLS-DA for variable selection

possible.pool<- c(3,5,8,10,12,15,18,20) # for test.keepX
gc()
set.seed(1)
my.tuned.splsda <- tune.splsda(X, Y, ncomp = comp, validation = 'Mfold',
                               folds = 10, nrepeat = 200, cpus = 8, # use repeated cross-validation
                               dist = 'centroids.dist', test.keepX =  possible.pool, 
                               measure = "BER",  # use balanced error rate of dist measure
                               progressBar = T)

my.tuned.splsda$choice.ncomp$ncomp # 2
png(filename = paste("Results/PLS_DA/Error_of_splsda_depending_from_computing",comp,"components.png", sep="_"), width = 1500, height = 2000, res=300)
plot(my.tuned.splsda, col = color.jet(comp))
dev.off()
optimal.comp<-my.tuned.splsda$choice.ncomp$ncomp

optimal.keepX <- my.tuned.splsda$choice.keepX[1:optimal.comp] # 15 and 15
optimal.keepX

final.splsda<-splsda(X,Y, ncomp=optimal.comp, keepX = optimal.keepX)
n=2 # choose which component plot depending on plot
# png(filename = paste("Results/PLS_DA/sparse_PLS_DA_comp_1_and",n, sep="_"), width = 2500, height = 1500, res=300)
# plotIndiv(final.splsda, comp = c(1,n), group = metadata_PLSDA2$Condition, ind.names = TRUE, legend = TRUE, ellipse = TRUE,
#           title = paste0("sparse PLSDA between Healthy and VEDOSS \n on scaled proportional genus data \n (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"))
# dev.off()
final.splsda$prop_expl_var
ggplot(mapping = aes(x=final.splsda$variates$X[,"comp1"], y=final.splsda$variates$X[,"comp2"], color=final.splsda$Y)) +
  scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  geom_point(size=2, alpha= 0.9) +
  geom_point(size=3.5, alpha= 0.5) +
  theme_classic(base_size = 10.5) +
  theme(title=element_text(size=10),
        legend.title =element_text(size=10.5)) +
  stat_ellipse(size=0.15)  +  # geom_text(mapping = aes(label=row.names(final.splsda$X)), color="black", size=2.3) +
  labs(title="sparse PLS-DA among Healthy, VEDOSS and SSC", 
       subtitle=paste0("on scaled proportional genus data \n (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"),
       caption = paste0("selected ",optimal.keepX[1]," genera on LC1 and ",optimal.keepX[2]," genera on LC2"),
       color="Condition",
       x=paste("LC1: ",round(final.splsda$prop_expl_var$X[1]*100,digits = 2),"% variation"),
       y=paste("LC2: ",round(final.splsda$prop_expl_var$X[2]*100,digits = 2),"% variation"))
ggsave(filename = paste("Results/PLS_DA/sparsePLS-DA on comp 1 and",n,".png"), width = 5.2, height = 4.3, dpi=300) 


################ plotting loadings

loadings<-as.data.frame(final.splsda$loadings$X)
identical(row.names(Taxa.genus.update[row.names(loadings),]),row.names(loadings))
loadings$Genus<-Taxa.genus.update[row.names(loadings),"Genus"]

a<-loadings[loadings$comp1!=0,]
{tax_a<-Taxa.genus.update
  tax_a$ASV<- row.names(tax_a)
  tax_a<-subset(tax_a, ASV %in% row.names(a))
  tax_a<-tax_a[row.names(a),]
  identical(tax_a$ASV,row.names(a))
  a$Genus<-tax_a$Genus
}

a$comp1<-as.numeric(a$comp1)
a<-a[order(abs(a$comp1)),]
a$Genus <- gsub("uncultured_ f R","uncultured_ f\nR",a$Genus) # for the margin
a$Genus <- factor(a$Genus, levels = a$Genus) # otherwise it would be re-ordered in plot
ggplot(mapping=aes(x=a$Genus,y=a$comp1)) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, size = 11, hjust = 1, vjust = 1)) +
  labs(x="Selected genera on component 1", 
       y="loadings", title = "Loadings of selected genera for component 1 in sPLSDA", 
       caption="\n a dotted line is plot at 50% of max loadings value") +
  geom_hline(yintercept = max(a$comp1)/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  #geom_hline(yintercept = min(a$comp1)/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = "Results/PLS_DA/Loadings of choosen genera for sPLSDA comp 1.png", width = 14, height = 8, dpi=300)

b<-loadings[loadings[[n]]!=0,]
{tax_b<-Taxa.genus.update
  tax_b$ASV<- row.names(tax_b)
  tax_b<-subset(tax_b, ASV %in% row.names(b))
  tax_b<-tax_b[row.names(b),]
  identical(tax_b$ASV,row.names(b))
  b$Genus<-tax_b$Genus
}

b[,n]<-as.numeric(b[,n])
b<-b[order( abs(b[,n])) ,]
b$Genus <- gsub("_group","",b$Genus)
b$Genus <- factor(b$Genus, levels = b$Genus)
ggplot(mapping=aes(x=b$Genus,y=b[,n])) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, size = 11, hjust = 1, vjust = 1)) +
  labs(x=paste("Selected genera on component",n), 
       y="loadings", title = paste("Loadings of selected genera for component",n,"in sPLSDA"), 
       caption="\n dotted lines are plot at 50% of both min and max loadings value") +
  geom_hline(yintercept = max(b[,n])/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  geom_hline(yintercept = min(b[,n])/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = paste("Results/PLS_DA/Loadings_of_choosen_genera_for_sPLSDA_comp",n,".png"), width = 14, height = 8, dpi=300)


############### testing the sPLS-DA model

X.test <- as.matrix(t(otu_table(data.train)))
Y.test <- as.factor(metadata_train$Condition)
predict.splsda <- predict(final.splsda, newdata =as.matrix(X.test), dist = "all")

# evaluating the prediction accuracy
predict.comp <- predict.splsda$class$centroids.dist[,n]
predict.comp
table(factor(predict.comp, levels = levels(Y)), as.factor(Y.test))
# now evaluating the prediction accuracy using ONLY the first component
predict.comp1 <- predict.splsda$class$centroids.dist[,1]
predict.comp1
table(factor(predict.comp1, levels = levels(Y)), as.factor(Y.test))


##### exporting settings and values of sPLSDA

suppressWarnings(rm(con))
con<-file("Results/PLS_DA/Settings_and_results_sPLSDA.txt")
sink(con)
cat("Number of training samples and their original number of genera", fill = TRUE)
cat(dim(X), fill = TRUE)
cat("\n Samples casually selected as test (and then discarded from training data set)", fill=TRUE)
cat(row.names(X.test), fill = TRUE)
cat("\n number of components suggested by perf function for normal PLSDA (centroid dist)", fill=TRUE)
cat(perf.splsda$choice.ncomp[3],fill=TRUE)
cat("\n number of components chosen for normal PLSDA (centroid dist)", fill=TRUE)
cat(comp, fill=TRUE)
cat("\n number of components suggested by tune.splsda function for sPLSDA", fill=TRUE)
cat(my.tuned.splsda$choice.ncomp$ncomp, fill=TRUE)
cat("\n number of components chosen for sPLSDA", fill=TRUE)
cat(optimal.comp, fill = TRUE)
cat("\n number of genera selected by tune.splsda function for each chosen component", fill=TRUE)
cat(optimal.keepX, fill=TRUE)
cat("\n Possible number of genera that could be selected by function tune.splsda", fill=TRUE)
cat(possible.pool, fill = TRUE)
cat("\n \n \n ### testing the prediction efficiency of sPLSDA through confusion matrix using ONLY the first component \n", fill=TRUE)
table(factor(predict.comp1, levels = levels(Y)), Y.test)
cat("\n \n ### testing the prediction efficiency of sPLSDA through confusion matrix component 1 and",n,"\n", fill=TRUE)
table(factor(predict.comp, levels = levels(Y)), Y.test)
cat("\n \n type of distance used: centroid distance")
sink()
close(con)


suppressWarnings( rm(a,b,con,tax_a,tax_b, colors_train, loadings, metadata_train, predict.comp,predict.comp1,train,X,Y,X.test,Y.test,predict.splsda,comp,n,test,optimal.comp,optimal.keepX) )
gc()



############# PLOTTING THE LEFSE RESULT ##################

a <- read.delim("QIIME/PICRUST2_LEFSE/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("QIIME/PICRUST2_LEFSE/Result_LEFSE.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.csv2(Significative_functions_LEFSE, file="Results/PICRUST2_LEFSE_results/Significantly_different_pathways_METACYC.csv", quote = F, row.names = F)
write.xlsx(Significative_functions_LEFSE, file="Results/PICRUST2_LEFSE_results/Significantly_different_pathways_METACYC.xlsx", col.names = T, showNA = F, row.names = F)

# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}

# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" , alpha= 0.5) +
  geom_bar(stat = "identity", width = 0.75, alpha= 0.9) +
  labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = 0)+
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  theme_classic(base_size = 15) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 18)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_every_result.png", width = 8, height = 5, dpi = 300)

# plotting only results over 3
Significative_functions_LEFSE_3<-Significative_functions_LEFSE[abs(Significative_functions_LEFSE$logLDA_score)>3,]
ggplot(data=Significative_functions_LEFSE_3, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" , alpha= 0.5) +
  geom_bar(stat = "identity", width = 0.95, alpha= 0.9) +
  labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = 0)+
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  theme_classic(base_size = 15) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 18)) +
  theme(legend.position = "bottom")
ggsave(filename = "Results/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_METACYC_only_over_3.png", width = 10, height = 4, dpi = 300)

system("echo 'Those are general results of a Kruskal Wallis filtered through a subsequent Wilcoxon test: each pathway (of a group) is differently abundant compared to BOTH of the other groups' > Results/PICRUST2_LEFSE_results/NB.txt")

rm(Significative_functions_LEFSE, a)



###################### SERUM FFAs QUANTITIES PLOTS #########################

### Serum
dir.create("Results/Serum_FFA")

serum<-read.xlsx(file="Serum_FFAs.xlsx", sheetIndex = 1 )
serum<-reshape::melt(serum)
serum$variable<-gsub("X2.","2-",serum$variable)
serum$variable<-gsub("Iso","iso",serum$variable)
serum$variable<-factor(serum$variable, levels=unique(serum$variable))
serum$Group<-factor(serum$Group, levels = c("Healthy","VEDOSS","SSc"))

ggplot(data=serum, aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  scale_color_manual(values=c("Healthy"="chartreuse3","VEDOSS"="coral","SSc"="red3")) +
  geom_boxplot(width=0.5, size= 0.25, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.5, jitter.width = 0.21),
             aes(color=Group), size= 0.08, alpha= 0.5) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=9.5), 
        axis.text.y = element_text(size=7.5),
        plot.margin = margin(5,1,3,1),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray"),
        panel.grid.major.x = element_line(size=0.05, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1), color="none" ) +
  labs(y="Raw quantity", fill="Condition", x="") + 
  scale_x_discrete(expand=c(-1, 1)) + 
  scale_y_sqrt(breaks=c(5,10,25,seq(0,600, 50)))
ggsave(file="Results/Serum_FFA/Serum_EVERY_FFA_boxplot.png", width=7.2, height = 5.8, dpi=300)


###### loop for each FFA

for(i in unique(serum$variable)){
  target<-serum[serum$variable==i , ]
  
  # statistics
  test <- dunn.test::dunn.test(target$value, target$Group, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Group)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Group", y="value", fill="Group", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Group), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Raw quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Results/Serum_FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }

  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Group)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Group", y="value", fill="Group", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Group), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Raw quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Results/Serum_FFA/Serum_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
    
}



##################### STOOL FFAs QUANTITIES PLOTS #########################

### Stool
dir.create("Results/Stool_FFA")

Stool<-as.data.frame(read.csv(file="Fecal_FFAs.csv"))
Stool<-Stool[ , !is.na(Stool[1, ]) ] # some value is available only for VEDOSS group, then unusefull here! (NB: the first row is a SSc sample --> NA for those columns!)
Stool<-reshape::melt(Stool)
Stool$variable<-gsub("X2.","2-",Stool$variable)
Stool$variable<-gsub("Iso","iso",Stool$variable)
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Group<-factor(Stool$Group, levels = c("Healthy","VEDOSS","SSc"))
Stool$value<-Stool$value*100 # percentual

ggplot(data=Stool, aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
  scale_color_manual(values=c("Healthy"="chartreuse3","VEDOSS"="coral","SSc"="red3")) +
  geom_boxplot(width=0.5, size= 0.25, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.5, jitter.width = 0.21),
             aes(color=Group), size= 0.08, alpha= 0.5) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=12,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=9.5), 
        axis.text.y = element_text(size=7.5),
        plot.margin = margin(5,1,3,1),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray"),
        panel.grid.major.x = element_line(size=0.05, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1), color="none" ) +
  labs(y="Percentual quantity", fill="Condition", x="") + 
  scale_x_discrete(expand=c(-1, 1)) + 
  scale_y_sqrt(breaks=c( 1, seq(0,100,5)) )
ggsave(file="Results/Stool_FFA/Stool_EVERY_FFA_boxplot.png", width=7.2, height = 5.8, dpi=300)


###### loop for each FFA

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  # statistics
  test <- dunn.test::dunn.test(target$value, target$Group, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Group)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Group", y="value", fill="Group", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Group), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentual quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Results/Stool_FFA/Stool_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Group)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Group", y="value", fill="Group", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      scale_color_manual(values=c("Healthy"="chartreuse","VEDOSS"="coral","SSc"="red3")) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Group), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentual quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("Results/Stool_FFA/Stool_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}




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
package$otherPkgs$mixOmics[c(1,4)]
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
rm(con)
