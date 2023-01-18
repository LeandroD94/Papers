##################### PREPARING THE ENVIRONMENT ##############

{
library("phyloseq")
library("ggplot2")
library("ggh4x")
library("ggpubr")
library("vegan")
library("dendextend")
library("DESeq2")
library("ggh4x") 
library("egg")
library("qiime2R")
library("dplyr")
library("xlsx")
}

dir.create("Results")

####################### IMPORTING DATA #####################

data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample
sample<-substring(sample, first = 10, last= 14) # strtrim to cut until X
sample<-gsub("-A","",sample)
sample<-gsub("-","",sample)
sample_names(data)<-sample # update

Metadata <- as.data.frame(readxl::read_excel("Metadata_mio.xlsx"))
row.names(Metadata)<-Metadata$IGA_ID
head(Metadata)
original_length<-length(Metadata$IGA_ID[!is.na(Metadata$IGA_ID)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$IGA_ID)),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))
write.csv2(cbind(original_names,Metadata), file="Metadata_with_FASTQ.csv",quote=F,row.names = F)

sample_data(data)$Time<-factor(sample_data(data)$Time, levels = c("PRE","POST"))
rm(original_length,original_names,sample)


##### adding informations about BMI and obesity

BMI<-read.csv(file = "Metadata_obesity.csv", row.names = 1)
head(BMI) # if BMI >30 then they are obese, otherwise they are overweight

length(row.names(BMI)) # 40
sample_length<-length(unique(sample_data(data)$Sample_number)) #36
BMI<-BMI[unique(sample_data(data)$Sample_number), ]
new_length<-length(BMI$Obese_T0[!is.na(BMI$Obese_T0)])
identical(new_length, length(unique(sample_data(data)$Sample_number))) # TRUE

BMI$Sample_number<-row.names(BMI)
temp<-dplyr::left_join(as(sample_data(data),"data.frame"),BMI)
row.names(temp)<-temp$IGA_ID
sample_data(data)<-temp
rm(temp)

head(sample_data(data))

rm(new_length,original_length, sample_length)

# <-- save here!
# save.image("data.RData")

############ CHECKING UNASSIGNED IN PROKARYOTE KINGDOM ######################

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
options(scipen=100)
write.csv2(e[,colnames(e)!="Kingdom"], file="QIIME/Unassigned_domain_checking.csv", row.names = T, quote = F)

rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)

######################### PREPARATION OF DATA #######################

setwd("Results")

{data.phy = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.fam = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
}

{data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
data.phy.prop <- transform_sample_counts(data.phy, function(ASV) ASV/sum(ASV)*100)
data.class.prop <- transform_sample_counts(data.class, function(ASV) ASV/sum(ASV)*100)
data.order.prop <- transform_sample_counts(data.order, function(ASV) ASV/sum(ASV)*100)
data.fam.prop <- transform_sample_counts(data.fam, function(ASV) ASV/sum(ASV)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
}

{Taxa.genus<-as.data.frame(tax_table(data.genus))
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

{taxa_temp<-Taxa.fam
for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
for( x in 1: length(which(is.na(taxa_temp$Family))) ) {
  taxa_temp$Family[ which(is.na(taxa_temp$Family))[1] ]<-paste("NA_ o",taxa_temp[which(is.na(taxa_temp$Family))[1],"Order"])}
for( x in 1: length(which(taxa_temp=="NA_ o NA")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="NA_ o NA")[1] ]<-paste("NA_ c",taxa_temp[which(taxa_temp$Family=="NA_ o NA")[1],"Class"])}
for( x in 1: length(which(duplicated(taxa_temp$Family[taxa_temp$Family=="NA_ c NA"]))) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="NA_ c NA")[1] ]<-paste("NA_ c NA",x+1) }
Taxa.fam.update<-taxa_temp
}

rm(taxa_temp)

####################### % ASSIGNED IN SILVA #########################

{a<-cbind(length(Taxa.genus$Genus),length(which(!is.na(Taxa.genus$Genus))),length(which(!is.na(Taxa.genus$Genus)))/length(Taxa.genus$Genus),"Genus")
b<-cbind(length(Taxa.fam$Family),length(which(!is.na(Taxa.fam$Family))),length(which(!is.na(Taxa.fam$Family)))/length(Taxa.fam$Family),"Family")
c<-cbind(length(Taxa.order$Order),length(which(!is.na(Taxa.order$Order))),length(which(!is.na(Taxa.order$Order)))/length(Taxa.order$Order),"Order")
d<-cbind(length(Taxa.class$Class),length(which(!is.na(Taxa.class$Class))),length(which(!is.na(Taxa.class$Class)))/length(Taxa.class$Class),"Class")
e<-cbind(length(Taxa.phy$Phylum),length(which(!is.na(Taxa.phy$Phylum))),length(which(!is.na(Taxa.phy$Phylum)))/length(Taxa.phy$Phylum),"Phylum")
assigned<-rbind.data.frame(a,b,c,d,e)
colnames(assigned)<-c("Total","Assigned","%","Taxa")
}
assigned
write.csv2(assigned,file="../QIIME/Percentuali_taxa_assegn_in_database_tass.csv",row.names = F, quote = F)
rm(a,b,c,d,e,assigned)

######################## GOOD'S COVERAGE ESTIMATOR #########################

filter<-prune_taxa(taxa_sums(data)==1, data)
length(which(taxa_sums(data)==1)) # if zero there are no singletons
{n<-as.data.frame(otu_table(filter))
  N<-as.data.frame(otu_table(data))
  G<-1-(colSums(n)/colSums(N))
}
con<-file("../QIIME/Percentuale_di_singletons_Good's_coverage.txt")
sink(con)
cat("GOOD'S COVERAGE ESTIMATOR \n", fill=TRUE)
cat("1-(n/N) for each sample, where n is number of singletons and N is Total ASV \n \n", fill=TRUE)
G
sink()
close(con)
rm(con, filter, n, N)

########################### COUNTS EXPORT ##########################################

dir.create("Abundances")

dir.create("Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100) # disable scientific annotation

dir.create("Abundances/TSS_counts")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/TSS_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Abundances/TSS_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Abundances/TSS_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Abundances/TSS_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/TSS_counts/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Abundances/TSS_counts/counts_genus.xlsx",row.names = T, showNA = F)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Abundances/TSS_counts/counts_phylum.xlsx",row.names = T, showNA = F)
}


#################### RATIO FIRMICUTES/BACTEROIDES ###################

suppressWarnings(dir.create("Abundances/Ratio_Firmi_Bacteroi"))

suppressWarnings(rm(data_fb, ratio_fb))
data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
head(tax_table(data_fb)) # The Firmicutes are the first row
ratio_fb<-cbind.data.frame(otu_table(data_fb),tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[1,])/as.vector(ratio_fb[2,]) )
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

ratio_fb<-ratio_fb[Metadata$IGA_ID, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb,Metadata[,c("Sample_number","Condition","Time")])

### placebo
ratio_fb_placebo<-ratio_fb[ratio_fb$Condition=="Placebo", ]
# adding group (PRE and POST) means
Ratios_Mean<-tapply(ratio_fb_placebo$Ratio, ratio_fb_placebo$Time, mean)
ratio_fb_placebo$Ratios_Mean_Condition_Time<-rep("placeholder")
ratio_fb_placebo[ratio_fb_placebo$Time=="PRE","Ratios_Mean_Condition_Time"]<-Ratios_Mean["PRE"]
ratio_fb_placebo[ratio_fb_placebo$Time=="POST","Ratios_Mean_Condition_Time"]<-Ratios_Mean["POST"]
# adding group (PRE and POST) standard error
Ratios_sd<-tapply(ratio_fb_placebo$Ratio, ratio_fb_placebo$Time, sd)
ratio_fb_placebo$Ratios_sterr_Condition_Time<-rep("placeholder")
ratio_fb_placebo[ratio_fb_placebo$Time=="PRE","Ratios_sterr_Condition_Time"]<-Ratios_sd["PRE"]/sqrt(length(which(ratio_fb_placebo$Time=="PRE")))
ratio_fb_placebo[ratio_fb_placebo$Time=="POST","Ratios_sterr_Condition_Time"]<-Ratios_sd["POST"]/sqrt(length(which(ratio_fb_placebo$Time=="POST")))

### probiotic
ratio_fb_probiotic<-ratio_fb[ratio_fb$Condition=="Probiotic", ]
# adding group (PRE and POST) means
Ratios_Mean<-tapply(ratio_fb_probiotic$Ratio, ratio_fb_probiotic$Time, mean)
ratio_fb_probiotic$Ratios_Mean_Condition_Time<-rep("temp")
ratio_fb_probiotic[ratio_fb_probiotic$Time=="PRE","Ratios_Mean_Condition_Time"]<-Ratios_Mean["PRE"]
ratio_fb_probiotic[ratio_fb_probiotic$Time=="POST","Ratios_Mean_Condition_Time"]<-Ratios_Mean["POST"]
# adding group (PRE and POST) sd error
Ratios_sd<-tapply(ratio_fb_probiotic$Ratio, ratio_fb_probiotic$Time, sd)
ratio_fb_probiotic$Ratios_sterr_Condition_Time<-rep("placeholder")
ratio_fb_probiotic[ratio_fb_probiotic$Time=="PRE","Ratios_sterr_Condition_Time"]<-Ratios_sd["PRE"]/sqrt(length(which(ratio_fb_probiotic$Time=="PRE")))
ratio_fb_probiotic[ratio_fb_probiotic$Time=="POST","Ratios_sterr_Condition_Time"]<-Ratios_sd["POST"]/sqrt(length(which(ratio_fb_probiotic$Time=="POST")))

### joining
Firmi_Bacter_Ratio<-rbind.data.frame(ratio_fb_placebo,ratio_fb_probiotic)
head(Firmi_Bacter_Ratio, n=2)
write.csv2(file="Abundances/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio.csv", Firmi_Bacter_Ratio)
write.xlsx(Firmi_Bacter_Ratio, file="Abundances/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio.xlsx", showNA = F)

Firmi_Bacter_Ratio$Sample_number_modified<-gsub("SYN_","",Firmi_Bacter_Ratio$Sample_number)
Firmi_Bacter_Ratio$Time<-factor(Firmi_Bacter_Ratio$Time, levels = c("PRE","POST"))
Firmi_Bacter_Ratio$Ratios_Mean_Condition_Time<-round(as.numeric(Firmi_Bacter_Ratio$Ratios_Mean_Condition_Time),digits = 2)
Firmi_Bacter_Ratio$Ratios_sterr_Condition_Time<-round(as.numeric(Firmi_Bacter_Ratio$Ratios_sterr_Condition_Time),digits = 2)


### BAR PLOT
ggplot(data=Firmi_Bacter_Ratio, aes(x=Time, fill=Time,y=Ratio)) +
  theme_classic(base_size =12) + 
  scale_fill_manual(values = c("PRE"="coral2","POST"="coral4")) +
  facet_grid2(.~Condition+Sample_number_modified, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_sqrt(breaks=c(1,5,seq(10,80,10)) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio")
ggsave(file="Abundances/Ratio_Firmi_Bacteroi/Ratio_barplot.png",width=8,height=4, dpi=300) 
dev.off()

### LINE PLOT (paired sample)
ggplot(data=Firmi_Bacter_Ratio, aes(x=Time, color=Time,y=Ratio)) +
  theme_bw(base_size =13) + 
  scale_color_manual(values = c("PRE"="coral2","POST"="coral4")) +
  facet_grid2(.~Condition+Sample_number_modified, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_sqrt(breaks=c(1,5,seq(10,80,10)) ) +
  geom_point(size=1.3, position=position_dodge(1)) +
  geom_line(aes(group="Sample_number"), color="gray",size=0.5) +
  guides(color="none") +
  theme(panel.grid.major.x = element_blank() )+
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio")
ggsave(file="Abundances/Ratio_Firmi_Bacteroi/Ratio_line_plot.png",width=10,height=3.5, dpi=300) 
dev.off()


### JITTER PLOT (with means too)
set.seed(2)
ggplot(data=Firmi_Bacter_Ratio, aes(x=Time, color=Time,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("PRE"="coral2","POST"="coral4")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_sqrt(breaks=c(1,5,seq(10,80,10)) ) +
  geom_errorbar(aes(ymin = Ratios_Mean_Condition_Time, # to add the mean line
                    ymax = Ratios_Mean_Condition_Time),
                size=0.35,width = 0.35, color= "black") +
  geom_errorbar(aes(ymax= Ratios_Mean_Condition_Time + Ratios_sterr_Condition_Time, # to add the standard deviation
                    ymin= ifelse(Ratios_Mean_Condition_Time - Ratios_sterr_Condition_Time < 0, # the if else is needed to avoid geom bar below zero
                             0, Ratios_Mean_Condition_Time - Ratios_sterr_Condition_Time)),
                width=.1, color="black", size= 0.15) +
  geom_jitter(width = 0.25, size=0.5, show.legend = F) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio")
ggsave(file="Abundances/Ratio_Firmi_Bacteroi/Ratios_jitter_mean.png",width=4,height=3, dpi=300) 
dev.off()


#### STATISTICAL TEST

shapiro.test(Firmi_Bacter_Ratio[Firmi_Bacter_Ratio$Condition=="Placebo","Ratio"])
hist(Firmi_Bacter_Ratio[Firmi_Bacter_Ratio$Condition=="Placebo","Ratio"])
# This is clearly not a Gaussian distribution 

Firmi_Bacter_Ratio<-Firmi_Bacter_Ratio[order(Firmi_Bacter_Ratio$Sample_number,Firmi_Bacter_Ratio$Time),] # reordering the samples for the pair test
# in placebo
res_pla<-wilcox.test(Ratio ~ Time, paired=T, data=Firmi_Bacter_Ratio[Firmi_Bacter_Ratio$Condition=="Placebo",])
round(res_pla$p.value, digits = 2)
# in treated
res_pro<-wilcox.test(Ratio ~ Time, paired=T, data=Firmi_Bacter_Ratio[Firmi_Bacter_Ratio$Condition=="Probiotic",])
round(res_pro$p.value, digits = 2)
# exporting the results
con <- file("Abundances/Ratio_Firmi_Bacteroi/Wilcoxon_paired.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Wilcoxon test (paired)\n", fill=T)
cat("Ratio~Time in placebo:   V=",res_pla$statistic, "p-value=", round(res_pla$p.value, digits = 2), "\n",fill=T)
cat("Ratio~Time in probiotic treatment:   V=",res_pro$statistic, "p-value=", round(res_pro$p.value, digits = 2))
sink()
close(con)
rm(con, res_pla,res_pro,ratio_fb, Ratios_Mean,ratio_fb_placebo, ratio_fb_probiotic)

###################### ABUNDANCES BAR PLOT ##########################

dir.create("Abundances")
setwd("Abundances")

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
}
tabella$Sample_number<-gsub("SYN_","",tabella$Sample_number)
tabella$Time<-factor(tabella$Time, levels = c("PRE","POST"))

ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Phylum)) + theme_classic(base_size =14) + 
  facet_grid2(Condition+Time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP5_phyla_abundances.png",width=6,height=12, dpi=300) 
dev.off()

ggplot(data=tabella, aes(x=Abundance, y=factor(Time, levels = c("POST","PRE")), fill=Phylum)) + 
  theme_classic(base_size =14) + 
  facet_grid2(Condition+Sample_number~., scales = "free", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.5,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant phyla", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP5_phyla_abundances_paired.png",width=6,height=12, dpi=300) 
dev.off()

# TOP 5 Generi
suppressWarnings(rm(top, others, tabella))
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
}
ggplot(data=tabella, aes(x=Abundance, y=Sample, fill=Genus)) + theme_classic(base_size =14) + 
  facet_grid2(Condition+Time~., scales = "free", space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.01,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=11), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="TOP5_genera_abundances.png",width=6,height=12, dpi=300) 
dev.off()

tabella$Sample_number<-gsub("SYN_","",tabella$Sample_number)
ggplot(data=tabella, aes(x=Abundance, 
                         y=factor(Time, levels = c("POST","PRE")), 
                         fill=Genus)) + 
  theme_classic(base_size =14) + 
  facet_grid2(Condition+Sample_number~., scales = "free", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.y = unit(2,"pt"))+
  scale_y_discrete (expand = c(0.5,0) ) +
  geom_bar(stat="identity", position="stack", width = 1) +
  theme(axis.text.y=element_text(angle=0, vjust=0.5, size=10), legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 13 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(y="Patients", x="Relative abundance", title = "Five most abundant genera", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="TOP5_genera_abundances_paired.png",width=6,height=12, dpi=300) 
dev.off()

setwd("..")

###################### HIERARCHICAL CLUSTERING #######################

dir.create("Hierarchical_clustering")
setwd("Hierarchical_clustering")

#################### PLACEBO 

rm(data.temp)
data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))
data.temp <- subset_samples(data.prop, Condition=="Placebo")

tabella_colore<-as.data.frame(cbind(as.character(sample_data(data.temp)$Time),as.character(sample_data(data.temp)$Time)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("PRE","coral",tabella_colore$Colore) #green
colors <- gsub("POST","chartreuse",colors) #orange

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.temp)))))
if(identical(c$labels,sample_names(data.temp))){
  c$labels<-sample_data(data.temp)$Sample_number
}
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Hierarchical_cluster_PLACEBO_Euclidean_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure of placebo group \nusing Euclidean distance on sqrt proportional ASVs",
     sub="PRE = orange     POST = green")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.temp))),method = "bray"))
if(identical(c$labels,sample_names(data.temp))){
  c$labels<-sample_data(data.temp)$Sample_number
}
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Hierarchical_cluster_PLACEBO_Bray_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure of placebo group \nusing Bray-Curtis distance on sqrt proportional ASVs",
     sub="PRE = orange     POST = green")
dev.off()

#################### PROBIOTIC 

rm(data.temp)
data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV))
data.temp <- subset_samples(data.prop, Condition=="Probiotic")

tabella_colore<-as.data.frame(cbind(as.character(sample_data(data.temp)$Time),as.character(sample_data(data.temp)$Time)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("PRE","coral",tabella_colore$Colore) #green
colors <- gsub("POST","chartreuse",colors) #orange

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.temp)))))
if(identical(c$labels,sample_names(data.temp))){
  c$labels<-sample_data(data.temp)$Sample_number
}
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Hierarchical_cluster_PROBIOTIC_Euclidean_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure of probiotic group \nusing Euclidean distance on sqrt proportional ASVs",
     sub="PRE = orange     POST = green")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.temp))),method = "bray"))
if(identical(c$labels,sample_names(data.temp))){
  c$labels<-sample_data(data.temp)$Sample_number
}
c<-as.dendrogram(c)
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Hierarchical_cluster_PROBIOTIC_Bray_sqrt_proportional_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure of probiotic group \nusing Bray-Curtis distance on sqrt proportional ASVs",
     sub="PRE = orange     POST = green")
dev.off()

setwd("..")




################## ALFA DIVERSITY PRE VS POST PLACEBO ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

rm(data.temp)
data.temp<-subset_samples(data, Condition=="Placebo")
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
ev<-H
ev$value<-(H$value)/log((obs$value))
ev$variable<-rep("Evenness")
# updating and ordering samples for pairwise wilcoxon
New_data<-rbind.data.frame(obs,H,ev)
head(New_data)
New_data<-New_data[order(New_data$Time, New_data$Sample),]
pAlpha$data<-New_data
}
pAlpha + theme_bw() + 
  geom_line(aes(group = pAlpha$data$Sample_number), col="grey", size=0.15) +
  labs(x="Time", title="Alpha diversity before and after placebo treatment") +
  guides(fill=FALSE, color=FALSE) +
  theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired = T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4)
ggsave(file="Alfa_diversity_PLACEBO_paired_wilcoxon.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
Obser_value<-Obser_value[order(Obser_value$Time, Obser_value$Sample),]
factor<-Obser_value$Time
factor
Obser_value$Sample # to re-check if they have the same order
wilcox.test(Obser_value$value~factor, paired=T)

################### BRAY CURTIS PLACEBO PRE VS POST #######################

dir.create("Bray_Curtis_Beta_diversity")
setwd("Bray_Curtis_Beta_diversity")

dir.create("Placebo_pre_vs_post")
setwd("Placebo_pre_vs_post")

rm(list=ls(pattern="data.temp"))
data.temp<-subset_samples(data, Condition=="Placebo")
head(sample_data(data.temp))

{data.temp.phy = tax_glom(data.temp, taxrank = "Phylum", NArm = F)
data.temp.class = tax_glom(data.temp, taxrank = "Class", NArm = F)
data.temp.order = tax_glom(data.temp, taxrank = "Order", NArm = F)
data.temp.fam = tax_glom(data.temp, taxrank = "Family", NArm = F)
data.temp.genus = tax_glom(data.temp, taxrank = "Genus", NArm = F)

data.temp.prop <- transform_sample_counts(data.temp, function(ASV) ASV/sum(ASV)*100)
data.temp.phy.prop <- transform_sample_counts(data.temp.phy, function(ASV) ASV/sum(ASV)*100)
data.temp.class.prop <- transform_sample_counts(data.temp.class, function(ASV) ASV/sum(ASV)*100)
data.temp.order.prop <- transform_sample_counts(data.temp.order, function(ASV) ASV/sum(ASV)*100)
data.temp.fam.prop <- transform_sample_counts(data.temp.fam, function(ASV) ASV/sum(ASV)*100)
data.temp.genus.prop <- transform_sample_counts(data.temp.genus, function(ASV) ASV/sum(ASV)*100)
}

rm(ASV.prop)
{ASV.prop<-as.data.frame(otu_table(data.temp.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.temp.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.temp.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.temp.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.temp.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.temp.phy.prop))
}

################## PERMANOVA

metadata<-as(sample_data(data.temp.prop),"data.frame")
library(vegan)

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")
perm_ASV_Bray<-perm_ASV$aov.tab # needed later for plot
perm_ASV_Bray

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")
perm_g_Bray<-perm_g$aov.tab # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Beta_diversity_permanova_PLACEBO_Bray_sqrt_prop_ASV.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
disper<-vegan::betadisper(BC.dist,metadata$Time)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Beta_dispersion_permanova_PLACEBO_Bray_sqrt_prop_ASV.csv",quote=F,row.names = T)


############################### PCOA

# on ASV
data.prop.labels<-data.temp.prop
sample_data(data.prop.labels)$Sample_number<-gsub("SYN_","",sample_data(data.prop.labels)$Sample_number) # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance in placebo group \n computed on sqrt proportional ASV", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_Bray$`Pr(>F)`[1]))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Placebo_on_ASV.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance in placebo group \n computed on sqrt proportional ASV", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Placebo_ASV_no_ellipse.png", width = 8, height = 6, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Bray-Curtis distance in placebo group \n computed on sqrt proportional ASV", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Placebo_on_ASV_points.png", width = 8, height = 6, dpi=300)

# again but on genera
{data.prop.labels<-data.temp.genus.prop
  sample_data(data.prop.labels)$Sample_number<-gsub("SYN_","",sample_data(data.prop.labels)$Sample_number) # if it's needed to change names in plot
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance in placebo group \n computed on sqrt proportional genera", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_Bray$`Pr(>F)`[1]))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Placebo_on_genera.png", width = 8, height = 6, dpi=300)


rm(list=ls(pattern="data.temp"))
setwd("../..")


################## ALFA DIVERSITY PRE VS POST PROBIOTIC ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data, Condition=="Probiotic")
pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Time")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Time, New_data$Sample),]
  pAlpha$data<-New_data
}
pAlpha + theme_bw() + 
  geom_line(aes(group = pAlpha$data$Sample_number),col="gray", size=0.15) +
  labs(x="Time", title="Alpha diversity before and after probiotic treatment") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=90, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = Time), label="p.format", method = "wilcox.test", paired = T, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4)
ggsave(file="Alfa_diversity_PROBIOTIC_paired_wilcoxon.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
library(dplyr)
Obser_value<-filter(alphadt, variable=="Observed")
Obser_value<-Obser_value[order(Obser_value$Time, Obser_value$Sample),]
factor<-Obser_value$Time
factor
Obser_value$Sample # to re-check if they have the same order
wilcox.test(Obser_value$value~factor, paired=T)


################### BRAY CURTIS PROBIOTIC PRE VS POST #######################

dir.create("Bray_Curtis_Beta_diversity")
setwd("Bray_Curtis_Beta_diversity")

system("echo 'Permanova does not consider the samples as pairs, there is not a corrispective analysis in pair' > Nota_riguardo_il_permanova.txt")

dir.create("Probiotic_pre_vs_post")
setwd("Probiotic_pre_vs_post")

rm(list=ls(pattern="data.temp"))
data.temp<-subset_samples(data, Condition=="Probiotic")
head(sample_data(data.temp))

{data.temp.phy = tax_glom(data.temp, taxrank = "Phylum", NArm = F)
  data.temp.class = tax_glom(data.temp, taxrank = "Class", NArm = F)
  data.temp.order = tax_glom(data.temp, taxrank = "Order", NArm = F)
  data.temp.fam = tax_glom(data.temp, taxrank = "Family", NArm = F)
  data.temp.genus = tax_glom(data.temp, taxrank = "Genus", NArm = F)
  
  data.temp.prop <- transform_sample_counts(data.temp, function(ASV) ASV/sum(ASV)*100)
  data.temp.phy.prop <- transform_sample_counts(data.temp.phy, function(ASV) ASV/sum(ASV)*100)
  data.temp.class.prop <- transform_sample_counts(data.temp.class, function(ASV) ASV/sum(ASV)*100)
  data.temp.order.prop <- transform_sample_counts(data.temp.order, function(ASV) ASV/sum(ASV)*100)
  data.temp.fam.prop <- transform_sample_counts(data.temp.fam, function(ASV) ASV/sum(ASV)*100)
  data.temp.genus.prop <- transform_sample_counts(data.temp.genus, function(ASV) ASV/sum(ASV)*100)
}

rm(ASV.prop)
{ASV.prop<-as.data.frame(otu_table(data.temp.prop))
  ASV.genus.prop<-as.data.frame(otu_table(data.temp.genus.prop))
  ASV.fam.prop<-as.data.frame(otu_table(data.temp.fam.prop))
  ASV.class.prop<-as.data.frame(otu_table(data.temp.class.prop))
  ASV.order.prop<-as.data.frame(otu_table(data.temp.order.prop))
  ASV.phy.prop<-as.data.frame(otu_table(data.temp.phy.prop))
}

################## PERMANOVA

metadata<-as(sample_data(data.temp.prop),"data.frame")

sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")
perm_ASV_Bray<-perm_ASV$aov.tab # needed later for plot
perm_ASV_Bray

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")
perm_g_Bray<-perm_g$aov.tab # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Time, data=metadata, permutations = 9999, method="bray")


beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Beta_diversity_permanova_PROBIOTIC_Bray_sqrt_prop_ASV.csv",quote=F,row.names = T)

# Perform an ANOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="bray")
disper<-vegan::betadisper(BC.dist,metadata$Time)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Beta_dispersion_permanova_PROBIOTIC_Bray_sqrt_prop_ASV.csv",quote=F,row.names = T)

############################### PCOA

# on ASV
data.prop.labels<-data.temp.prop
sample_data(data.prop.labels)$Sample_number<-gsub("SYN_","",sample_data(data.prop.labels)$Sample_number) # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance in probiotic group \n computed on sqrt proportional ASV", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_Bray$`Pr(>F)`[1]))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Probiotic_on_ASV.png", width = 8, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance in probiotic group \n computed on sqrt proportional ASV", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Probiotic_ASV_no_ellipse.png", width = 8, height = 6, dpi=300)
# without names neither ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) +
  labs(title="PCoA with Bray-Curtis distance in probiotic group \n computed on sqrt proportional ASV", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Probiotic_on_ASV_points.png", width = 8, height = 6, dpi=300)

# again but on genera
{data.prop.labels<-data.temp.genus.prop
  sample_data(data.prop.labels)$Sample_number<-gsub("SYN_","",sample_data(data.prop.labels)$Sample_number) # if it's needed to change names in plot
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Time") +
  scale_color_manual(values=c("PRE"="coral","POST"="chartreuse")) +
  geom_point(size=3) + theme_classic(base_size = 14) +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +  geom_text(aes(label=sample_data(data.sqrt_prop)$Sample_number), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Bray-Curtis distance in probiotic group \n computed on sqrt proportional genera", color="Time", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_Bray$`Pr(>F)`[1]))
ggsave(file="PCoA_Beta_div_Bray_Curtis_Probiotic_on_genera.png", width = 8, height = 6, dpi=300)

rm(list=ls(pattern="data.temp"))
setwd("../..")

######################### PCoA ON OTHER FACTORS #############################

dir.create("Bray_Curtis_Beta_diversity")
setwd("Bray_Curtis_Beta_diversity")

dir.create("Other_factors")
setwd("Other_factors")

######################## PRE TREATMENT

rm(list=ls(pattern="data.temp"))
data.temp<-subset_samples(data, Time=="PRE")
data.temp.prop <- transform_sample_counts(data.temp, function(ASV) ASV/sum(ASV)*100)
head(sample_data(data.temp))

# on ASV
data.prop.labels<-data.temp.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Smoke", shape = "Condition") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("No"="blue","Ex"="grey","Yes"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA with Bray-Curtis distance before the treatment \n computed on sqrt proportional ASV", color="Smoke", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Smoke_PRE_TREATMENT.png", width = 8, height = 6, dpi=300)

plot_ordination(data.sqrt_prop, ordBC, color = "Gender", shape = "Condition") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("M"="blue","F"="pink")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA with Bray-Curtis distance before the treatment \n computed on sqrt proportional ASV", color="Gender", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Gender_PRE_TREATMENT.png", width = 8, height = 6, dpi=300)

######################## POST TREATMENT

rm(list=ls(pattern="data.temp"))
data.temp<-subset_samples(data, Time=="POST")
data.temp.prop <- transform_sample_counts(data.temp, function(ASV) ASV/sum(ASV)*100)
head(sample_data(data.temp))

# on ASV
data.prop.labels<-data.temp.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Smoke", shape = "Condition") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("No"="blue","Ex"="grey","Yes"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA with Bray-Curtis distance after the treatment \n computed on sqrt proportional ASV", color="Smoke", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Smoke_POST_TREATMENT.png", width = 8, height = 6, dpi=300)

plot_ordination(data.sqrt_prop, ordBC, color = "Gender", shape = "Condition") +
  geom_line(aes(group=Sample_number), col="gray", size= 0.15)+
  scale_color_manual(values=c("M"="blue","F"="pink")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA with Bray-Curtis distance after the treatment \n computed on sqrt proportional ASV", color="Gender", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Gender_POST_TREATMENT.png", width = 8, height = 6, dpi=300)

system("echo 'Fattori distribuiti troppo poco omogeneamente... ma aggiunte queste PCoA giusto come controllo\n\nInoltre, se manca l elisse in qualche gruppo è dovuta a mancanza di sufficienti campioni per il calcolo dell intervallo di confidenza '> Nota_bene.txt")

setwd("../..")

################# DA WITH DESEQ2 _ IN PLACEBO (PAIR) #############################

# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

suppressWarnings(rm(data.temp, data_pruned, data.genus_pruned, data_pruned.genus.prop))

data.temp<-subset_samples(data, Condition=="Placebo")
data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
{ data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
Taxa.genus_pruned<-taxa_temp
rm(taxa_temp)
}
# adding informations to uncultured families
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
{for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
  taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
Taxa.fam_pruned<-taxa_temp
rm(taxa_temp)
}
# other taxonomic vocabularies
{ Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

dir.create("Diff_abundance_DESeq2_PLACEBO")

# genus ##############################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  # View(res_genus)
  rm(res)
  write.csv2(res_genus, file="Diff_abundance_DESeq2_PLACEBO/Analisi diff dei generi_ratio Placebo_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)
  tax_table(data_pruned.genus.prop)<-as.matrix(Taxa.genus_pruned)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  # tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family ##############################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_fam<-as.data.frame(res)
  res_fam$ASV<-row.names(res_fam)
  Taxa.fam_pruned$ASV<-row.names(Taxa.fam_pruned)
  res_fam<-dplyr::left_join(res_fam, Taxa.fam_pruned, by="ASV")
  res_fam$Kingdom<-NULL
  res_fam$Genus<-NULL
  res_fam$Species<-NULL
  # View(res_fam)
  rm(res)
  write.csv2(res_fam, file="Diff_abundance_DESeq2_PLACEBO/Analisi diff delle famiglie con DeSeq2 _ ratio Placebo_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_fam$Family
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.fam.prop, Family %in% target)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
}

# class ###################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_C<-as.data.frame(res)
  res_C$ASV<-row.names(res_C)
  Taxa.class_pruned$ASV<-row.names(Taxa.class_pruned)
  res_C<-dplyr::left_join(res_C, Taxa.class_pruned, by="ASV")
  res_C$Kingdom<-NULL
  res_C$Genus<-NULL
  res_C$Species<-NULL
  res_C$Order<-NULL
  res_C$Family<-NULL
  View(res_C)
  rm(res)
  write.csv2(res_C, file="Diff_abundance_DESeq2_PLACEBO/Analisi diff delle classi_ ratio Placebo_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_C$Class
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.class.prop, Class %in% target)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
}

# order ############################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_O<-as.data.frame(res)
  res_O$ASV<-row.names(res_O)
  Taxa.order_pruned$ASV<-row.names(Taxa.order_pruned)
  res_O<-dplyr::left_join(res_O, Taxa.order_pruned, by="ASV")
  res_O$Kingdom<-NULL
  res_O$Genus<-NULL
  res_O$Species<-NULL
  res_O$Family<-NULL
  View(res_O)
  rm(res)
  write.csv2(res_O, file="Diff_abundance_DESeq2_PLACEBO/Analisi diff degli ordini_ ratio Placebo_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_O$Order
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.order.prop, Order %in% target)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
}

# phylum  ########################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_phy<-as.data.frame(res)
  res_phy$ASV<-row.names(res_phy)
  Taxa.phy_pruned$ASV<-row.names(Taxa.phy_pruned)
  res_phy<-dplyr::left_join(res_phy, Taxa.phy_pruned, by="ASV")
  res_phy$Kingdom<-NULL
  res_phy$Genus<-NULL
  res_phy$Species<-NULL
  res_phy$Class<-NULL
  res_phy$Order<-NULL
  res_phy$Family<-NULL
  View(res_phy)
  rm(res)
  write.csv2(res_phy, file="Diff_abundance_DESeq2_PLACEBO/Analisi diff dei phyla_ratio Placebo_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_phy$Phylum
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.phy.prop, Phylum %in% target)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
}


tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_f$Taxa<-"Families"
tabella_f[,c("Phylum","Order","Class")]<-NULL
colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"
# 
# tabella_o$Taxa<-"Orders"
# tabella_o[,c("Phylum","Class")]<-NULL
# colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"

# tabella_c$Taxa<-"Classes"
# tabella_c[,"Phylum"]<-NULL
# colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"
# 
# tabella_p$Taxa<-"Phyla"
# colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"

# building segment plot basics
tabella_g$ASV<-NULL # no more needed, now the colnames have to match
tabella_tot<-rbind.data.frame(tabella_g,tabella_f)
# to prevent the ggplot reorder (1)
tabella_tot<-tabella_tot[order(match(tabella_tot$Time,levels(tabella_tot$Time))) , ] # reordering based on the level specified order, important for the order in X axis

tabella_tot$Xaxis<-paste0(tabella_tot$Bacteria, tabella_tot$Time)
# to prevent the ggplot reorder (2)
tabella_tot$Xaxis<-factor(tabella_tot$Xaxis, levels = unique(tabella_tot$Xaxis) )
levels(tabella_tot$Xaxis)

plot_1<-ggplot(tabella_tot, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 12) +
  scale_fill_manual(values = c(values=c("PRE"="coral","POST"="chartreuse"))) +
  facet_wrap2(nrow=1,factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, labeller = labeller(group = label_wrap_gen(width = 34)),scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time), size=3) +
  geom_line(aes(group=Sample_number), size=0.5) +
  theme(strip.text.x=element_text(size=13,colour="black"), strip.switch.pad.wrap = unit(10,"line")  ) + 
  theme(axis.text.y = element_text(size=11))+
  scale_x_discrete(labels=rep(unique(levels(tabella_tot$Time))),expand=c(0,0.5)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(1,max(tabella_tot$Abundance)+1,1))) +
  theme(plot.title= element_text(size=14) ,legend.key.size=unit(0.7,"cm"), legend.text=element_text(size=13)) +
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
# scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))
plot_1 +
  labs(title= "Differently abundant taxa in placebo group", y="Proportional Abundance", fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_PLACEBO/Plot_genera_DeSeq2_Time.png", width = 10, height = 6, dpi=300)
dev.off()

suppressWarnings(rm(plot_1, tabella_g, tabella_f, res_genus, res_fam, res))


################# DA WITH DESEQ2 _ IN PROBIOTIC (PAIR) #############################

# Trimming under sum of 10 (see DeSeq2 tutorial) and preparing new data (other ASV may be selected after glomming)

rm(data.temp, data_pruned, data.genus_pruned, data_pruned.genus.prop)

data.temp<-subset_samples(data, Condition=="Probiotic")
data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
unique(sample_data(data_pruned)$Cond_time)
{ data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
  Taxa.genus_pruned<-taxa_temp
  rm(taxa_temp)
}
# adding informations to uncultured families
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
{for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  Taxa.fam_pruned<-taxa_temp
  rm(taxa_temp)
}
# other taxonomic vocabularies
{ Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
  Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
  Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

dir.create("Diff_abundance_DESeq2_PROBIOTIC")


# genus ##############################

suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  # View(res_genus)
  rm(res)
  write.csv2(res_genus, file="Diff_abundance_DESeq2_PROBIOTIC/Analisi diff dei generi_ratio Probiotic_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family ##############################
suppressWarnings(suppressWarnings(rm(res, DE, target)))
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_fam<-as.data.frame(res)
  res_fam$ASV<-row.names(res_fam)
  Taxa.fam_pruned$ASV<-row.names(Taxa.fam_pruned)
  res_fam<-dplyr::left_join(res_fam, Taxa.fam_pruned, by="ASV")
  res_fam$Kingdom<-NULL
  res_fam$Genus<-NULL
  res_fam$Species<-NULL
  # View(res_fam)
  # rm(res)
  # write.csv2(res_fam, file="Diff_abundance_DESeq2_PROBIOTIC/Analisi diff delle famiglie con DeSeq2 _ ratio Probiotic_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_fam$Family
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.fam.prop, Family %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_f<-tabella
}

# class ###################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_C<-as.data.frame(res)
  res_C$ASV<-row.names(res_C)
  Taxa.class_pruned$ASV<-row.names(Taxa.class_pruned)
  res_C<-dplyr::left_join(res_C, Taxa.class_pruned, by="ASV")
  res_C$Kingdom<-NULL
  res_C$Genus<-NULL
  res_C$Species<-NULL
  res_C$Order<-NULL
  res_C$Family<-NULL
  View(res_C)
  rm(res)
  write.csv2(res_C, file="Diff_abundance_DESeq2_PROBIOTIC/Analisi diff delle classi_ ratio Probiotic_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_C$Class
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.class.prop, Class %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_c<-tabella
}

# order ############################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_O<-as.data.frame(res)
  res_O$ASV<-row.names(res_O)
  Taxa.order_pruned$ASV<-row.names(Taxa.order_pruned)
  res_O<-dplyr::left_join(res_O, Taxa.order_pruned, by="ASV")
  res_O$Kingdom<-NULL
  res_O$Genus<-NULL
  res_O$Species<-NULL
  res_O$Family<-NULL
  # View(res_O)
  rm(res)
  write.csv2(res_O, file="Diff_abundance_DESeq2_PROBIOTIC/Analisi diff degli ordini_ ratio Probiotic_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_O$Order
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.order.prop, Order %in% target)
  tabella<-psmelt(target)
  #tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_o<-tabella
}

# phylum  ########################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~Sample_number + Time)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Time", "PRE", "POST"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

if(length(res$log2FoldChange)>0){
  res_phy<-as.data.frame(res)
  res_phy$ASV<-row.names(res_phy)
  Taxa.phy_pruned$ASV<-row.names(Taxa.phy_pruned)
  res_phy<-dplyr::left_join(res_phy, Taxa.phy_pruned, by="ASV")
  res_phy$Kingdom<-NULL
  res_phy$Genus<-NULL
  res_phy$Species<-NULL
  res_phy$Class<-NULL
  res_phy$Order<-NULL
  res_phy$Family<-NULL
  View(res_phy)
  rm(res)
  write.csv2(res_phy, file="Diff_abundance_DESeq2_PROBIOTIC/Analisi diff dei phyla_ratio Probiotic_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_phy$Phylum
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.phy.prop, Phylum %in% target)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_p<-tabella
}

tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_f$Taxa<-"Families"
tabella_f[,c("Phylum","Order","Class")]<-NULL
colnames(tabella_f)[colnames(tabella_f)=="Family"]<-"Bacteria"

#
### The resulter order is clearly a redundant (same result as family, same ASV but collapsed at the order levele)
#
# tabella_o$Taxa<-"Orders"
# tabella_o[,c("Phylum","Class")]<-NULL
# colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"

# tabella_c$Taxa<-"Classes"
# tabella_c[,"Phylum"]<-NULL
# colnames(tabella_c)[colnames(tabella_c)=="Class"]<-"Bacteria"
# 
# tabella_p$Taxa<-"Phyla"
# colnames(tabella_p)[colnames(tabella_p)=="Phylum"]<-"Bacteria"

# THE FAMILY IS A REDUNDANT OF THE ROMBUTSIE GENERA
# tabella_2<-rbind.data.frame(tabella_f)#,tabella_o)

# building segment plot basics
tabella_g<-tabella_g[order(match(tabella_g$Time,levels(tabella_g$Time))),]
tabella_g$Xaxis<-paste0(tabella_g$Bacteria, tabella_g$Time)
# tabella_2$Xaxis<-paste0(tabella_2$Bacteria, tabella_2$Time)
# to prevent the ggplot reorder
tabella_g$Xaxis<-factor(tabella_g$Xaxis, levels = unique(tabella_g$Xaxis) )
# tabella_2$Xaxis<-factor(tabella_2$Xaxis, levels = unique(tabella_2$Xaxis) )

tabella_g$Bacteria<-gsub("Hafnia-Obesumbacterium", "Hafnia\nObesumbacterium", tabella_g$Bacteria)
plot_1<-ggplot(tabella_g, aes(x= Xaxis, y=Abundance, fill=Time)) + 
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c(values=c("PRE"="coral","POST"="chartreuse"))) +
  facet_wrap2(nrow=1,
              factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, 
              labeller = labeller(group = label_wrap_gen(width = 34)),scales = "free", 
              strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Time), size=3) +
  geom_line(aes(group=Sample_number), size=0.5) +
  theme(strip.text.x=element_text(size=13,colour="black"),
        strip.switch.pad.wrap = unit(10,"line")  ) + 
  theme(axis.text.y = element_text(size=12))+
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,2,4,
                        seq(6,max(tabella_g$Abundance)+1.5,2))) +
  scale_x_discrete(labels=unique(levels(tabella_g$Time)),expand=c(0,0.5)) +
  theme(plot.title= element_text(size=14) ,legend.key.size=unit(0.7,"cm"), 
        legend.text=element_text(size=13)) +
  theme(panel.grid.minor.y= element_blank()) +
  theme(legend.margin=margin(0, 0, 0, 0), legend.position = "none")
# scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))
plot_1 +
  labs(title= "Differently abundant genera in probiotic group", y="Proportional Abundance", 
       fill="Time", x="")
ggsave(filename = "Diff_abundance_DESeq2_PROBIOTIC/Plot_genera_DeSeq2_Time.png", width = 10, height = 6, dpi=300)
dev.off()


### THE FAMILIES ARE ONLY REDUNDANT RESULTS
# 
# Redund<-c("Hafniaceae") # to select redundants (same ASV, same results at different levels)
# tabella_2<-subset(tabella_2, ! Bacteria %in% Redund)
# plot_2<-ggplot(tabella_2, aes(x= Xaxis, y=Abundance, fill=Time)) + theme_bw(base_size = 12) +
#   scale_fill_manual(values = c(values=c("PRE"="coral","POST"="chartreuse"))) +
#   facet_wrap2(nrow=1,factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera"))~Bacteria, scales = "free", strip=strip_nested(size = "variable", bleed = T), drop = TRUE) +
#   geom_point(aes(color=Time), size=2.5) +
#   geom_line(aes(group=Sample_number)) +
#   theme(strip.text.x=element_text(size=13,colour="black"), axis.text.y = element_text(size=12)) +
#   scale_x_discrete(labels=rep(unique(levels(tabella_2$Time))), expand=c(0,0.5)) +
#   scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(tabella_2$Abundance),2))) +
#   theme(plot.title= element_text(size=14) ,legend.key.size=unit(0.7,"cm"), legend.text=element_text(size=13)) +
#   theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
#   theme(legend.margin=margin(0, 0, 0, 0), legend.position="none") 
# #scale_y_sqrt(breaks=c(0.1, 0.25,0.75,seq(0,8,0.5)))
# plot_2 +
#   labs(title= "Differently abundant families in probiotic group", y="Proportional Abundance", fill="Time", x="")
# 
# ggsave(filename = "Diff_abundance_DESeq2_PROBIOTIC/Plot_2_DESeq2_Time_no redundants.png", width = 9, height = 6, dpi=300)
# dev.off()
# 
# # now a unique plot
# head_plot<-plot_1 + # refining first plot with the title
#   labs(title= "Differently abundant taxa in probiotic group", y="Proportional abundance", x="") + 
#   theme(plot.title = element_text(size=18))
# png(filename = "Diff_abundance_DESeq2_PROBIOTIC/Plot_TOTAL_DESeq2_Time_no redundants.png", width = 2900, height = 2600, res=300)
# grid.arrange(head_plot,
#              plot_2 + labs(title= "", y="Proportional abundance", x=""))
# dev.off()

system('touch Diff_abundance_DESeq2_PROBIOTIC/The_other_results_are_only_redundants')

rm(tabella_2, plot_2, plot_1, head_plot, tabella_g)

############## PLUS: CONFRONTING PLACEBO AND PROBIOTIC ##################

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data, Time=="POST")
dir.create("PLUS_placebo_vs_Treated")
setwd("PLUS_placebo_vs_Treated")

######################## ALPHA DIVERSITY

pAlpha<-plot_richness(data.temp, measures=c("Shannon", "Observed"), x="Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H<-dplyr::filter(pAlpha$data, variable=="Shannon")
obs<-dplyr::filter(pAlpha$data, variable=="Observed")
# adding evanness
{identical(H$Sample_ID, obs$Sample_ID) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  # updating and ordering samples for pairwise wilcoxon
  New_data<-rbind.data.frame(obs,H,ev)
  head(New_data)
  New_data<-New_data[order(New_data$Time, New_data$Sample),]
  pAlpha$data<-New_data
}
pAlpha + theme_bw() + geom_boxplot(width=0.5) +
  labs(x="Treatment", title="Alpha diversity between not treated and treated") +
  guides(fill=FALSE, color=FALSE) + theme(axis.text.x= element_text(angle=30, vjust=1, hjust=1, size=11)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", paired = F, label.x= 1.5, size=3.5, label.y.npc = "top", vjust=-0.5, hjust=0.4)
ggsave(file="Alfa_diversity_Mann_Withney.png", width = 6,height =6, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed")
Obser_value<-Obser_value[order(Obser_value$Condition, Obser_value$Sample),]
factor<-Obser_value$Condition
factor
wilcox.test(Obser_value$value~factor, paired=F)

##################### BETA DIVERSITY

data.temp.prop <- transform_sample_counts(data.temp, function(ASV) ASV/sum(ASV)*100)

data.prop.labels<-data.temp.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "bray")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variance explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Placebo"="coral","Probiotic"="deepskyblue")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() + 
  labs(title="PCoA with Bray-Curtis distance \n computed on sqrt proportional ASV", color="", x=paste("PC1: ",eigval[1],"% variance"), y=paste("PC2: ",eigval[2],"% variance"))
ggsave(file="PCoA_Beta_div_Bray_Placebo_vs_Treated_on_ASV.png", width = 8, height = 6, dpi=300)

##################### DA ANALYSIS

suppressWarnings(rm(data_pruned, data.genus_pruned, data_pruned_genus.prop))

data_pruned<- prune_taxa(taxa_sums(data.temp) > 10, data.temp)
{ data.genus_pruned<- tax_glom(data_pruned, taxrank = "Genus", NArm = F)
  data.fam_pruned<- tax_glom(data_pruned, taxrank = "Family", NArm = F)
  data.class_pruned<- tax_glom(data_pruned, taxrank = "Class", NArm = F)
  data.order_pruned<- tax_glom(data_pruned, taxrank = "Order", NArm = F)
  data.phy_pruned<- tax_glom(data_pruned, taxrank = "Phylum", NArm = F)
  data_pruned.phy.prop <- transform_sample_counts(data.phy_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.class.prop <- transform_sample_counts(data.class_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.order.prop <- transform_sample_counts(data.order_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.fam.prop <- transform_sample_counts(data.fam_pruned, function(ASV) ASV/sum(ASV)*100)
  data_pruned.genus.prop <- transform_sample_counts(data.genus_pruned, function(ASV) ASV/sum(ASV)*100)
}
# adding informations to uncultured genera
taxa_temp<-as.data.frame(tax_table(data_pruned.genus.prop))
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
  Taxa.genus_pruned<-taxa_temp
  rm(taxa_temp)
}
# adding informations to uncultured families
taxa_temp<-as.data.frame(tax_table(data_pruned.fam.prop))
{for( x in 1: length(which(taxa_temp$Family=="uncultured")) ) {
  taxa_temp$Family[which(taxa_temp$Family=="uncultured")[1]]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Family=="uncultured")[1],"Order"])}
  for( x in 1: length(which(taxa_temp=="uncultured_ o uncultured")) ) {
    taxa_temp$Family[ which(taxa_temp$Family=="uncultured_ o uncultured")[1] ]<-paste("uncultured_ c",taxa_temp[which(taxa_temp$Family=="uncultured_ o uncultured")[1],"Class"])}
  Taxa.fam_pruned<-taxa_temp
  rm(taxa_temp)
}
# other taxonomic vocabularies
{ Taxa.class_pruned<-as.data.frame(tax_table(data.class_pruned))
  Taxa.order_pruned<-as.data.frame(tax_table(data.order_pruned))
  Taxa.phy_pruned<-as.data.frame(tax_table(data.phy_pruned))
}

# genus ##############################
suppressWarnings(rm(res, DE, target, res_genus))
{DEseq_data<-phyloseq_to_deseq2(data.genus_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Placebo", "Probiotic"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05) & abs(res$log2FoldChange)>1,]
  res
}

if(length(res$log2FoldChange)>0){
  res_genus<-as.data.frame(res)
  res_genus$ASV<-row.names(res_genus)
  Taxa.genus_pruned$ASV<-row.names(Taxa.genus_pruned)
  res_genus<-left_join(res_genus, Taxa.genus_pruned, by="ASV")
  res_genus$Kingdom<-NULL
  res_genus$Species<-NULL
  # View(res_genus)
  rm(res)
  write.csv2(res_genus, file="Analisi diff dei generi_ratio Probiotic_PRE vs POST.csv", row.names = F, quote=F, na = "")
  # box plot
  target<-res_genus$Genus
  target<-na.omit(target)
  target<-subset_taxa(data_pruned.genus.prop, Genus %in% target)
  tabella<-psmelt(target)
  tabella$Abundance<-sqrt(tabella$Abundance)
  tabella_g<-tabella
}

# family ##############################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.fam_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Placebo", "Probiotic"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

# class ###################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.class_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Placebo", "Probiotic"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

# order ############################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.order_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Placebo", "Probiotic"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

# phylum  ########################################
suppressWarnings(rm(res, DE, target))
{DEseq_data<-phyloseq_to_deseq2(data.phy_pruned, ~ Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Placebo", "Probiotic"))
  res = res[order(res$padj, na.last=NA), ]
  res<-res[(res$padj < 0.05), ]
  res<-res[abs(res$log2FoldChange)>1, ]
  res
}

suppressWarnings(rm(res, DE, target))

tabella_g$Taxa<-"Genera"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

ggplot(tabella_g, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phyla","Classes","Orders","Families","Genera")), scales = "free_x", space="free") +
  scale_fill_manual(values=c("Probiotic"="deepskyblue","Placebo"="coral")) +
  geom_boxplot(width=0.8) + theme_bw( ) + 
  theme(strip.text.x=element_text(size=10,colour="black")) + 
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,axis.text.x = element_text(angle = 30, vjust=1, hjust=1, size=10), axis.text.y = element_text(size=10)) + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  #scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(tabella_g$Abundance),2))) +
  theme(plot.title= element_text(size=12) ,legend.key.size=unit(0.7,"cm"), legend.text=element_text(size=8)) +
  theme(panel.grid.minor.y= element_blank()) + guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant genera between post placebo and post probiotic (treated) group", y="Proportional Abundance", fill="Condition", x="")
ggsave(filename = "Total_abund_differences_DeSeq2_Condition.png", width = 8, height = 6, dpi=300)
dev.off()

setwd("..")

##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("R_version_and_packages.txt") # create this file and connect it to R through "con" object 
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

sink() # restore STR OUTPUT to R console
close(con)
rm(con)
