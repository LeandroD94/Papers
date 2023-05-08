##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
  library("ggplot2")
  library("vegan")
  library("ggpubr")
  library("ggh4x")
  library("egg")
  library("reshape")
  library("dendextend")
  library("DESeq2")
  library("ggvenn")
  library("Hmisc")
  library("stringr")
  library("xlsx")  
  library("dplyr")
  library("qiime2R")
}

{dir.create("Data_check")
  dir.create("Data_check/PCoA_test")
  dir.create("Results")
}

Subsets<-c("Saliva","Stool") # to create each folder
for(x in Subsets){
  dir.create(paste0('Results/',x))
  dir.create(paste0('Results/',x,'/Abundances'))
  dir.create(paste0('Results/',x,'/Abundances/Ratio_Firmi_Bacteroi'))
  dir.create(paste0('Results/',x,'/Hierarchical_clustering'))
  dir.create(paste0('Results/',x,'/Beta_div'))
  dir.create(paste0('Results/',x,'/DA_DESeq2/'))
  dir.create(paste0('Results/',x,'/PICRUST2_LEFSE_results'))
}

suppressWarnings(rm(Subsets,x))


options(scipen = 100) # disable scientific annotation



######### CIRCOPLOT function
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




####################### IMPORTING DATA #####################

# devtools::install_github("jbisanz/qiime2R")
data<-qza_to_phyloseq(features="QIIME/table.qza", taxonomy="QIIME/taxonomy.qza", tree = "QIIME/rooted-tree.qza")
# changing names
sample<-sample_names(data)
original_names<-sample
sample

if (!require('stringr')) {install.packages('stringr')}
{ sample[stringr::str_starts(original_names,pattern = "ID1548")] <- substring(original_names[stringr::str_starts(original_names,pattern = "ID1548")],first = 19)
sample[stringr::str_starts(original_names,pattern = "ID1642")] <- substring(original_names[stringr::str_starts(original_names,pattern = "ID1642")],first = 19)
sample[stringr::str_starts(original_names,pattern = "ID2117")] <- substring(original_names[stringr::str_starts(original_names,pattern = "ID2117")],first = 11, last=14)
sample
}

sample_names(data)<-sample # update

Metadata <- as.data.frame(read.table("Metadata.csv",sep=",",header = T))
row.names(Metadata)<-Metadata$FASTQ # column with FASTQ/SAMPLE name
head(Metadata)
Metadata<- Metadata[Metadata$FASTQ!="NO_FASTQ", ] # lost sample
original_length<-length(Metadata$FASTQ[!is.na(Metadata$FASTQ)])
Metadata<-Metadata[sample, ]
identical(as.numeric(length(Metadata$FASTQ[!is.na(Metadata$FASTQ)])),as.numeric(original_length))

sample_data(data)<-Metadata
identical(as.numeric(length(sample_names(data))),as.numeric(length(original_names)))

rm(original_length,sample)

sample_data(data)$Condition<-factor(sample_data(data)$Condition, levels=c("Healthy","ICC")) # decides the order in plots
head(sample_data(data))


# save.image("data_pre_noise_filtering.RData")




#################### FILTERING NOISES FROM DATA SET ####################

if(! "proof1" %in% ls()){
  unfiltered_data<-data
}

suppressWarnings(rm(data.genus.temp))
data.genus.temp<-tax_glom(unfiltered_data, taxrank = "Genus", NArm = F)
write.csv2(cbind(otu_table(data.genus.temp),tax_table(data.genus.temp)), file="Data_check/Raw_ASV_Table_pre_filtering.csv", row.names = T)


###### cutting under 0.005% to remove noises/contaminants, too conservative but also safe cutoff, see   23202435   and   DOI: 10.1128/mSystems.00290-19
data.genus.temp<-transform_sample_counts(data.genus.temp, function(x) (x/sum(x))*100)
filtered<-taxa_names(filter_taxa(data.genus.temp, function(x) max(x) <= 0.005, TRUE))
write.csv( cbind(as.data.frame(tax_table(data.genus.temp))[filtered, c("Phylum","Family","Genus")], as.data.frame(otu_table(data.genus.temp))[filtered, ] ), 
           file="Data_check/Filtered_genus_under_MAX_0005_cutoff.csv")

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
write.csv2(e[,colnames(e)!="Kingdom"], file="Data_check/Inferred_domains_proportions.csv", row.names = T, quote = F)
rm(a,b,c,c_a,c_b,d,e,total,Unass,Unass.prop,x)


proof1<- "Marker of the filtering, it is required for the script"


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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Sampling_type") +
  scale_color_manual(values=c("ICC"="coral","Healthy"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= Subject), 
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
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Sampling_type") +
  scale_color_manual(values=c("ICC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= Subject), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

png(filename = "Data_check/PCoA_test/PCoA_BRAY_test.png", width = 3600, height = 1800, res=300)
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
p1<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Sampling_type") +
  scale_color_manual(values=c("ICC"="coral","Healthy"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= Subject), 
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
p2<-plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape="Sampling_type") +
  scale_color_manual(values=c("ICC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= Subject), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

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
p1<-plot_ordination(data.prop.labels, ordBC, color = "Condition", shape="Sampling_type") +
  scale_color_manual(values=c("ICC"="coral","Healthy"="chartreuse")) +
  guides(color="none", shape="none") +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= Subject), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="PCoA wUnifrac (on proportional ASV)\n\n UNfiltered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

##### filtered wUnifrac
suppressWarnings(rm(data.prop.labels, data.sqrt_prop))
data.prop.labels<-data.prop
#{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
{DistBC = phyloseq::distance(data.prop.labels, method = "wunifrac")
  ordBC = ordinate(data.prop.labels, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
p2<-plot_ordination(data.prop.labels, ordBC, color = "Condition", shape="Sampling_type") +
  scale_color_manual(values=c("ICC"="coral","Healthy"="chartreuse")) +
  geom_point(size=3) + theme_bw(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label= Subject), 
            color="black", size=2.5, show.legend = FALSE) +
  labs(title="\n \n filtered data", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))

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
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off()
rm(r)


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
assigned$`%`<-round(as.numeric(assigned$`%`)*100, 2)
write.csv2(assigned,file="Data_check/Percentual_of_taxa_assigned_in_database_after_filters.csv",row.names = F, quote = F)

suppressWarnings(rm(a,b,c,d,e,assigned,data.phy,data.class,data.order,data.fam.data.genus))


########################## GENERAL PCoA ############################

if(! "proof1" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

if(length(unique(sample_data(data)[["Sampling_type"]])) == 1 ){ # It is a control, in case the script has been restarted
  data<-all_data # this object is created during the subsettings
}


{ data.genus = tax_glom(data, taxrank = "Genus", NArm = F)
  data.genus.prop <- transform_sample_counts(data.genus, function(ASV) ASV/sum(ASV)*100)
  data.prop <- transform_sample_counts(data, function(ASV) ASV/sum(ASV)*100)
}


####### on ASV
data.prop.labels<-data.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape = "Sampling_type") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_line(aes(group=Subject), color= "grey", size= 0.15) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same sample")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape = "Sampling_type") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3, alpha=0.68) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  geom_line(aes(group=Subject), color= "grey", size= 0.15) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same sample")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 9, height = 6, dpi=300)

rm(data.sqrt_prop, eigval, DistBC, ordBC,data.prop.labels)


####### again but on Genera
data.prop.labels<-data.genus.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape = "Sampling_type") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_line(aes(group=Subject), color= "grey", size= 0.15) +
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed Genera)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same sample")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera.png", width = 9.5, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition", shape = "Sampling_type") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3, alpha= 0.68) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  geom_line(aes(group=Subject), color= "grey", size= 0.15) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed Genera)", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"),
       caption="the lines are connecting the same sample")
ggsave(file="Results/General_Beta_div/PCoA_Beta_diversity_Hellinger_on_Genera_points.png", width = 9.5, height = 6, dpi=300)


#################### STARTING THE ANALYSIS OF THE SALIVA SUBSET #######################

to_reset<-ls() [! ls() %in% c("data","Metadata","unfiltered_data","contam_data","proof1","all_data", "CircoTax2")]
to_reset<-to_reset [! to_reset %in% c( ls(pattern = "Genera_Saliva.splsda") , ls(pattern="Genera_Saliva.DESEQ2")) ] # to mantain DA and splsda results
rm(list = to_reset)


if(! "proof1" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

if(length(unique(sample_data(data)[["Sampling_type"]])) > 1 ){   # it works only with the whole dataset
  all_data<-data
  rm(data)
}

data<-subset_samples(all_data, Sampling_type=="Saliva")
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

Meta_Saliva<-Metadata[Metadata$Sampling_type=="Saliva",]
  

########################### COUNTS EXPORT Saliva ##########################################


dir.create("Results/Saliva/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Saliva/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Saliva/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Saliva/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Saliva/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Saliva/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Saliva/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/Saliva/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Saliva/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT Saliva ##########################

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
ggplot(data=tabella, aes(x=Subject, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (Saliva)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Saliva/Abundances/TOP_5_phyla_Saliva.png",width=9,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Saliva/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
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
ggplot(data=tabella, aes(x=Subject, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera (Saliva)", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Saliva/Abundances/TOP_5_genera_Saliva.png",width=9,height=5,dpi=300)
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
ggplot(data=tabella, aes(x=Subject, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Saliva/Abundances/TOP_8_genera_Saliva.png",width=9,height=5,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


#################### RATIO FIRMICUTES/BACTEROIDES Saliva ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_Saliva)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_Saliva[,c("FASTQ","Condition")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Condition, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Condition=="ICC","Ratios_Mean"]<-Ratios_Mean["ICC"]
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_Mean"]<-Ratios_Mean["Healthy"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Condition, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Condition=="ICC","Ratios_st_err"]<-Ratios_st["ICC"]/sqrt(length(which(ratio_fb$Condition=="ICC")))
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_st_err"]<-Ratios_st["Healthy"]/sqrt(length(which(ratio_fb$Condition=="Healthy")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/Saliva/Abundances/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Saliva.csv", ratio_fb)

ratio_fb$Condition<-factor(ratio_fb$Condition, levels = c("ICC","Healthy"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-wilcox.test(Ratio ~ Condition, paired=F, data=ratio_fb)
p_val<-round(res$p.value, digits = 3)
p_val
# exporting the results
con <- file("Results/Saliva/Abundances/Ratio_Firmi_Bacteroi/Mann_Whi_Wilcoxon_Saliva.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Condition (Saliva) :   V=",res$statistic, "p-value=", p_val, "\n",fill=T)
cat("\n Shapiro Wilk test p-value: ",distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ, fill=Condition, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Healthy"="chartreuse",
                               "ICC"="coral")) +
  facet_grid2(.~Condition, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Condition"))+
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, seq(2,10,0.5))) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Mann Whitney p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/Saliva/Abundances/Ratio_Firmi_Bacteroi/Ratio_barplot_Saliva.png",width=6,height=4.5, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Condition, color=Condition,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("ICC"="coral","Healthy"="chartreuse")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, seq(2,10,0.5))) +
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
        axis.text.y=element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio in Saliva",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/Saliva/Abundances/Ratio_Firmi_Bacteroi/Ratios_jitter_with_mean_and_STerr_Saliva.png",width=4,height=3.5, dpi=300) 
dev.off()


##### BOXPLOT
ggplot(data=ratio_fb, aes(x=Condition, fill=Condition,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_fill_manual(values = c("ICC"="coral","Healthy"="chartreuse")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, seq(2,10,0.5))) +
  geom_boxplot(show.legend = F, color="black", size=0.2, width=0.5) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio in Saliva", caption = paste("Mann Withney p-value:",p_val) )
ggsave(file="Results/Saliva/Abundances/Ratio_Firmi_Bacteroi/Ratios_FB_boxplot_median_IQR.png",width=4,height=3.5, dpi=300) 
dev.off()


suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))



##################### SETTING THE GROUP COLORS ######################

# same function down there will search the colors from here
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("Healthy","chartreuse",tabella_colore$Colore) #green
colors <- gsub("ICC","coral",colors) #orange


#################### HIERARCHICAL CLUSTERING Saliva ###################

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.prop))))) # then Hellinger
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/Saliva/Hierarchical_clustering/Hierarchical_cluster_Hellinger_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure of Saliva subset using Hellinger distance",
     sub="Healthy = green     ICC = orange")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/Saliva/Hierarchical_clustering/Hierarchical_cluster_Bray_sqrt_prop_normalized_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure of Saliva subset using\n Bray-Curtis distance on sqrt proportional ASVs",
     sub="Healthy = green     ICC = orange")
dev.off()

suppressWarnings(rm(c,color_table,labels_colors))


########################## ALFA DIVERSITY Saliva ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition", color = "Condition")
pAlpha
# plot_richness( ) compute diversity like estimate_diversity( )
H <- pAlpha$data[pAlpha$data$variable=="Shannon", ]
obs <- pAlpha$data[pAlpha$data$variable=="Observed", ]
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

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha= 0) + theme_classic2() +  # alpha 0 overwrites the original black color
  scale_color_manual(values = c("Healthy"="chartreuse", "ICC"="coral")) +
  labs(x="Condition", title="Alpha diversity between ICC and Healthy (Saliva)") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 1.4, size=3.5,
                     label.y.npc = "top", vjust=-0.5, hjust=0.3) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Saliva/Alfa_diversity_with_Mann_Withn_Wilcox_Saliva.png", width = 6,height =5.3, dpi=300)


# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY  Saliva #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[1]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[1]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Saliva/Beta_div/Beta_diversity_permanova_Hellinger_Saliva.csv",quote=F,row.names = T)


### PLUS: checking it with Jaccard too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="jaccard")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Saliva/Beta_div/Beta_divers_permanova_JACCARD_on_ASV.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Condition)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Saliva/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS Saliva #####################

# on ASV
data.prop.labels<-data.prop
# sample_names(data.prop.labels)<-gsub(" ","",sample_names(data.prop.labels))       # if it's needed to change names in plot
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on Saliva subset", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Saliva/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_Saliva.png", width = 9, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on Saliva subset",
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
ggsave(file="Results/Saliva/Beta_div/PCoA_Beta_diversity_Hellinger_ASV_no_ellipse.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3, alpha=0.5) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on Saliva subset", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
ggsave(file="Results/Saliva/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 9, height = 6, dpi=300)

rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels)


### again but on genera
# Hellinger
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3, alpha=0.5) + theme_classic(base_size = 13.4) + stat_ellipse(size=0.2) + 
  #geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=1.5, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on Saliva subset", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Saliva/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 9, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))



#################### VEEN DIAGRAM Saliva ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.5% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.5, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Healthy<-subset_samples(data.venn, Condition=="Healthy")
Healthy<-as.character(tax_table(prune_taxa(taxa_sums(Healthy)>0, Healthy))[,"Genus"])

ICC<-subset_samples(data.venn, Condition=="ICC")
ICC<-as.character(tax_table(prune_taxa(taxa_sums(ICC)>0, ICC))[,"Genus"])


ONLY_IN_ICC<- ICC[! ICC %in% Healthy]
ONLY_IN_ICC<- paste(ONLY_IN_ICC, collapse = ", ")
head(ONLY_IN_ICC)

ONLY_IN_Healthy<- Healthy[! Healthy %in% ICC]
ONLY_IN_Healthy<- paste(ONLY_IN_Healthy, collapse = ", ")
head(ONLY_IN_Healthy)


x<-list(Healthy=Healthy,ICC=ICC)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("chartreuse","coral")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram of Saliva subset \n(only genera with minimal abundance > 0.5% \n at least in three samples)")
ggsave(filename = "Results/Saliva/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()


suppressWarnings(rm(ONLY_IN_ICC, ONLY_IN_Healthy, x, con, ICC, Healthy, data.venn, who))



################### DA WITH DESEQ2 Saliva #####################

if(! "proof1" %in% ls() ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Healthy", "ICC"))
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
    # write.csv2(r, file=paste0("Results/Saliva/DA_DESeq2/DA_",t,"_ratio_Healthy_vs_ICC_Saliva.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Saliva/DA_DESeq2/Every_result_DESeq2_Saliva.csv", row.names = F)
write.xlsx(Res_tot, file="Results/Saliva/DA_DESeq2/Every_result_DESeq2_Saliva.xlsx", showNA = F, col.names = T)

ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 50, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance),2))) +
  labs(title= "Differently abundant Taxa (Saliva)", y="Proportional Abundance", 
       fill="Condition", x="")
ggsave(filename = "Results/Saliva/DA_DESeq2/DA_Condition_every_result_Saliva.png", width = 12, height = 8.5, dpi=300)
dev.off()

Redund<- c("Flavobacteriaceae","Saccharimonadales","Saccharimonadia","Patescibacteria")
Table_tot2<-subset(Table_tot, ! Bacteria %in% Redund) # to remove redundant results
DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  scale_fill_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 30, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot$Abundance+2),2))) +
  labs(title= "Differently abundant Taxa (Saliva)", y="Proportional Abundance", 
       fill="Condition", x="")

DE_plot
ggsave(filename = "Results/Saliva/DA_DESeq2/DA_Condition_WITHOUT_redundants_Saliva.png", width = 12, height = 8, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Saliva/DA_DESeq2/NB.txt ")

# for further comparisons ...
Genera.DESEQ2_Saliva<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])



##### CIRCOPLOT TO PLOT THE LOG2FC

### preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
### removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
Res_tot2$Genus <- gsub("Saccharimonadaceae", "       Saccharimonadaceae", Res_tot2$Genus) # more space for the plot

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12, fc_col=2,sort="no")

png(file="Results/Saliva/DA_DESeq2/CircoTax_plot.png",width=1800,height=1800, res = 300)
circo + labs(title="ICC vs Healthy \n Saliva DA taxa log2foldchange\n")
dev.off()

circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))

# diff abundances + logfoldchange
png(file="Results/Saliva/DA_DESeq2/Final_plot.png",width=4560, height=2150, res = 300)
ggarrange(DE_plot, circo_no_m, nrow = 1, widths= c(1,0.58), labels = c("A)", "B)"))
dev.off()

rm(circo, circo_no_m)


################# CORRELATIONS WITH FFA VALUES (Saliva) ######################

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data.genus.prop, Condition != "Healthy") # many variables have the same value for each Healthy subjects --> not usable 
prune.dat_target <- subset_taxa(data.temp, Genus %in% Genera.DESEQ2_Saliva)
sample_names(prune.dat_target) <- sample_data(data.temp)$Subject
target <-as.data.frame(otu_table(prune.dat_target))
target_tax<-as.data.frame(tax_table(prune.dat_target))
identical(row.names(target_tax),row.names(target))
row.names(target)<-target_tax$Genus
colnames(target)
rm(target_tax)
target<-as.data.frame(t(target))

SCFA <- read.table("FFA_values.tsv", header=T, sep="\t", dec=",")
SCFA[ , !colnames(SCFA)%in%c("Sample_Code","Condition")] <- apply(SCFA[, !colnames(SCFA)%in%c("Sample_Code","Condition")], 2, as.numeric)

# same order with the phyloseq data
row.names(SCFA)<-SCFA$Sample_Code
SCFA<-SCFA[row.names(target), ]
length(SCFA$Sample_Code[! is.na(SCFA$Sample_Code)]) == 8  # TRUE --> any sample cutted avay (there are FFA values regarding 8 ICC)
SCFA <- SCFA[! is.na(SCFA$Sample_Code), ]

target<-target[row.names(SCFA), ] # to re-subset according to missing subjects among SCFA values


# Computing the correlations...
x<-cbind.data.frame(target,SCFA)
{x<-as.matrix(x[ ! colnames(x) %in% c("Sample_Code","Condition")])
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("DA_Genera","FFA","Corr","pvalue")
}
corr<-subset(correlation, DA_Genera %in% colnames(target))
corr<-subset(corr, FFA %in% colnames(SCFA))
colnames(corr)<-c("DA_Genera","SCFA","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$DA_Genera, y = corr$SCFA, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "Serum DA SCFA", x= "Saliva DA Genera", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Saliva/DA_genera_vs_SCFA_using_only_ICC.csv.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="Results/Saliva/DA_genera_vs_SCFA_using_only_ICC.csv", quote = F, na="", row.names = F)

corr2<-corr[corr$SCFA!="X2.MethylButyrric", ]
ggplot(corr2, aes(x = DA_Genera, y = SCFA, fill = Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "Serum DA SCFA", x= "Saliva DA Genera", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Saliva/DA_genera_vs_SCFA_using_only_ICC_WITHOUT_MethButyrric.csv.png", dpi=300, width = 7, height = 7)


################### ANALYSING LEFSE RESULTS SALIVA ###########################

suppressWarnings(rm(a))
a <- read.delim("QIIME/PICRUST2_LEFSE/picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("QIIME/PICRUST2_LEFSE/Result_LEFSE_SALIVA.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="Results/Saliva/PICRUST2_LEFSE_results/Significantly_different_pathways_KEGG_Stool.xlsx", col.names = T, showNA = F, row.names = F)


# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="ICC"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="ICC"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), 
                                               fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), 
            hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Healthy",1,0), size=3.2) +
  scale_fill_manual(values=c("ICC"="coral", "Healthy"="chartreuse")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15))+
  theme(legend.position = "bottom")
ggsave(filename = "Results/Saliva/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_KEGG_every_result.png", width = 8, height = 7.2, dpi = 300)

#plotting only results over 3
Significative_functions_LEFSE_3<-Significative_functions_LEFSE[abs(Significative_functions_LEFSE$logLDA_score)>3,]
ggplot(data=Significative_functions_LEFSE_3, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black") + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= -0.7, label=Pathway)) +
  scale_fill_manual(values=c("ICC"="coral", "Healthy"="Chartreuse")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15))+
  theme(legend.position = "bottom")
ggsave(filename = "Results/Saliva/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_KEGG_over_3.png", width = 7, height = 2.5, dpi = 300)

rm(Significative_functions_LEFSE, a)


#################### STARTING THE ANALYSIS OF THE STOOL SUBSET #######################

to_reset<-ls() [! ls() %in% c("data","Metadata","unfiltered_data","contam_data","proof1","all_data","CircoTax2")]
to_reset<-to_reset [! to_reset %in% c( ls(pattern = "Genera.splsda") , ls(pattern="Genera.DESEQ2")) ] # to mantain DA and splsda results
rm(list = to_reset)


if(! "proof1" %in% ls()){
  cat("\n Wait! Did you perform the filtering step??? \n\n")
  Sys.sleep(2)
}

if(length(unique(sample_data(data)[["Sampling_type"]])) > 1 ){
  all_data<-data
  rm(data)
}

data<-subset_samples(all_data, Sampling_type=="Stool")
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

Meta_Stool<-Metadata[Metadata$Sampling_type=="Stool",]


########################### COUNTS EXPORT Stool ##########################################

dir.create("Results/Stool/Abundances/Raw_counts")
{write.csv2(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="Results/Stool/Abundances/Raw_counts/counts_otu.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.phy),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Stool/Abundances/Raw_counts/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Stool/Abundances/Raw_counts/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Stool/Abundances/Raw_counts/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Stool/Abundances/Raw_counts/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Stool/Abundances/Raw_counts/counts_genus.csv",quote=F)
}

options(scipen = 100)
dir.create("Results/Stool/Abundances/Relative_abundances")
{write.csv2(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_phylum.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.class.prop),"matrix"),as(tax_table(data.class),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_class.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.order.prop),"matrix"),as(tax_table(data.order),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_order.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.fam.prop),"matrix"),as(tax_table(data.fam),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_family.csv",quote=F)
  write.csv2(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_genus.csv",quote=F)
  write.xlsx(cbind(as(otu_table(data.genus.prop),"matrix"),as(tax_table(data.genus),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_genus.xlsx",showNA = F, col.names = T)
  write.xlsx(cbind(as(otu_table(data.phy.prop),"matrix"),as(tax_table(data.phy),"matrix")),file="Results/Stool/Abundances/Relative_abundances/counts_phylum.xlsx",showNA = F, col.names = T)
}


###################### ABUNDANCES BAR PLOT Stool ##########################

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
ggplot(data=tabella, aes(x=Subject, y=Abundance, fill=Phylum)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Patients", y="Relative abundance", title = "Five most abundant phyla (Stool)", caption = " 'Others' includes every phylum below rank 5 ")
ggsave(file="Results/Stool/Abundances/TOP_5_phyla_Stool.png",width=9,height=5, dpi=300)
dev.off()

# means of TOP5 phyla
write.xlsx(file = "Results/Stool/Abundances/TOP5_phyla_average_abundance.xlsx", row.names = F,
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
# tabella$Condition<-gsub("ICC"," ",tabella$Condition)    # if it is needed to rename
ggplot(data=tabella, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_5) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + 
  labs(x="Patients", y="Relative abundance", title = "Five most abundant genera (Stool)", caption = " 'Others' includes every genus below rank 5 ")
ggsave(file="Results/Stool/Abundances/TOP_5_genera_Stool.png",width=9,height=5,dpi=300)
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
ggplot(data=tabella, aes(x=Subject, y=Abundance, fill=Genus)) + facet_grid(cols= vars(Condition),scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size =14) +
  scale_fill_manual(values=fill_color_8) +
  theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) +
  theme(legend.position="bottom", legend.margin = margin(3,0,0,0)) + guides(fill=guide_legend(nrow=3)) + 
  labs(x="Patients", y="Relative abundance", title = "Eigth most abundant genera", caption = " 'Others' includes every genus below rank 8 ")
ggsave(file="Results/Stool/Abundances/TOP_8_genera_Stool.png",width=9,height=5,dpi=300)
dev.off()

suppressWarnings(rm(tabella,prune.data.others, prune.dat_top, tabella_top, tabella_others, top, others))


#################### RATIO FIRMICUTES/BACTEROIDES Stool ###################

suppressWarnings(rm(data_fb, ratio_fb))

data_fb<-subset_taxa(data.phy.prop, Phylum %in% c("Bacteroidota","Firmicutes"))
F_index<-grep("Firmicutes",tax_table(data_fb)[,"Phylum"])  # to annotate which one is the first row
B_index<-grep("Bacteroidota",tax_table(data_fb)[,"Phylum"])
ratio_fb<-otu_table(data_fb)
ratio_fb <- rbind.data.frame(ratio_fb, as.vector(ratio_fb[F_index,])/as.vector(ratio_fb[B_index,]) ) # to compute the ratio
row.names(ratio_fb)<-c(as.vector(tax_table(data_fb)[,"Phylum"]) , "Ratio")
ratio_fb<-t(ratio_fb)

name_match<- row.names(Meta_Stool)
ratio_fb<-ratio_fb[name_match, ] # same order
identical(length(row.names(ratio_fb)),length(sample_names(data_fb))) # TRUE
ratio_fb<-cbind.data.frame(ratio_fb, Meta_Stool[,c("FASTQ","Condition")])


# mean
Ratios_Mean<-tapply(ratio_fb$Ratio, ratio_fb$Condition, mean)
ratio_fb$Ratios_Mean<-rep("temp")
ratio_fb[ratio_fb$Condition=="ICC","Ratios_Mean"]<-Ratios_Mean["ICC"]
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_Mean"]<-Ratios_Mean["Healthy"]
ratio_fb$Ratios_Mean<-as.numeric(ratio_fb$Ratios_Mean)
# st err
Ratios_st<-tapply(ratio_fb$Ratio, ratio_fb$Condition, sd)
ratio_fb$Ratios_st_err<-rep("temp")
ratio_fb[ratio_fb$Condition=="ICC","Ratios_st_err"]<-Ratios_st["ICC"]/sqrt(length(which(ratio_fb$Condition=="ICC")))
ratio_fb[ratio_fb$Condition=="Healthy","Ratios_st_err"]<-Ratios_st["Healthy"]/sqrt(length(which(ratio_fb$Condition=="Healthy")))
ratio_fb$Ratios_st_err<-as.numeric(ratio_fb$Ratios_st_err)

head(ratio_fb, n=2)
write.csv2(file="Results/Stool/Abundances/Ratio_Firmi_Bacteroi/Firmi_Bacter_Ratio_Stool.csv", ratio_fb)

ratio_fb$Condition<-factor(ratio_fb$Condition, levels = c("ICC","Healthy"))


####### STATISTICAL TEST

distr<-shapiro.test(ratio_fb[,"Ratio"])
distr
hist(ratio_fb[,"Ratio"])

res<-wilcox.test(Ratio ~ Condition, paired=F, data=ratio_fb)
p_val<-round(res$p.value, digits = 3)
p_val

# exporting the results
con <- file("Results/Stool/Abundances/Ratio_Firmi_Bacteroi/Mann_Whi_Wilcoxon_Stool.txt") # create this file and connect it to R through "con" object 
sink(con, append = TRUE) # redirect STR ERROR to "con"
cat("Mann_Whitney_Wilcoxon test \n", fill=T)
cat("Ratio~Condition (Stool) :   V=",res$statistic, "p-value=", p_val, "\n",fill=T)
cat("\n Shapiro Wilk test p-value: ",distr$p.value)
sink()
close(con)


######## BAR PLOT
ggplot(data=ratio_fb, aes(x=FASTQ, fill=Condition, y=Ratio)) +
  theme_bw(base_size =12) + 
  scale_fill_manual(values = c("Healthy"="chartreuse",
                               "ICC"="coral")) +
  facet_grid2(.~Condition, scales = "free_x", 
              space="free", strip = strip_nested(size="constant"))+
  theme(panel.spacing.x = unit(2,"pt"))+
  scale_x_discrete (expand = c(0.01,0) ) +
  geom_line(aes(y= ratio_fb$Ratios_Mean, group="Condition"))+
  scale_y_continuous(breaks=c(0, 1, 3, 5, seq(2,20,2))) +
  geom_bar(stat="identity", position="stack", width = 0.8) +
  guides(fill="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=8)) +
  labs(x="Patients", y="Firmicutes/Bacteroidetes Ratio",
       subtitle = paste("Mann Whitney p-value:",p_val),
       caption = "the line is the average ratio of the group")
ggsave(file="Results/Stool/Abundances/Ratio_Firmi_Bacteroi/Ratio_barplot_Stool.png",width=6,height=4.5, dpi=300) 
dev.off()


####### JITTER PLOT (with means and str err)
set.seed(2)
ggplot(data=ratio_fb, aes(x=Condition, color=Condition,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_color_manual(values = c("ICC"="coral","Healthy"="chartreuse")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_continuous(breaks=c(0, 1, 3, 5, seq(2,20,2))) +
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
        axis.text.y=element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio", title = "Firmicutes/Bacteroidetes Ratio in Stool",
       caption = paste("Mann Whitney p-value:",p_val) )
ggsave(file="Results/Stool/Abundances/Ratio_Firmi_Bacteroi/Ratios_jitter_with_mean_and_STerr_Stool.png",width=4,height=3.5, dpi=300) 
dev.off()


##### BOXPLOT
ggplot(data=ratio_fb, aes(x=Condition, fill=Condition,y=Ratio)) +
  theme_bw(base_size =9) +
  scale_fill_manual(values = c("ICC"="coral","Healthy"="chartreuse")) +
  facet_grid(.~Condition, scales = "free_x", space="free")+
  theme(panel.spacing.x = unit(20,"pt"))+
  scale_x_discrete (expand = c(0.2,0.3) ) +
  scale_y_continuous(breaks=c(0, 1, 3, 5, seq(2,20,2))) +
  geom_boxplot(show.legend = F, color="black", size=0.2, width=0.5) +
  guides(color="none") +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, size=8),
        axis.text.y=element_text(size=6)) +
  labs(x="", y="Firmicutes/Bacteroidetes Ratio in Stool", caption = paste("Mann Withney p-value:",p_val) )
ggsave(file="Results/Stool/Abundances/Ratio_Firmi_Bacteroi/Ratios_FB_boxplot_median_IQR.png",width=4,height=3.5, dpi=300) 
dev.off()

suppressWarnings(rm(con, res, ratio_fb, Ratios_Mean, Ratios_st_err,data_fb, p_val))


##################### SETTING THE GROUP COLORS ######################

# same function down there will search the colors from here
tabella_colore<-as.data.frame(cbind(as.character(sample_data(data)$Condition),as.character(sample_data(data)$Condition)))
colnames(tabella_colore)<-c("Gruppo","Colore")
colors <- gsub("Healthy","chartreuse",tabella_colore$Colore) #green
colors <- gsub("ICC","coral",colors) #orange


#################### HIERARCHICAL CLUSTERING Stool ###################

# euclidean
c<-hclust(dist(t(sqrt(otu_table(data.prop))))) # then Hellinger
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/Stool/Hierarchical_clustering/Hierarchical_cluster_Hellinger_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.4, cex.sub=1.3)
plot(c,main="Community structure of Stool subset using Hellinger distance",
     sub="Healthy = green     ICC = orange")
dev.off()

# Bray Curtis
c<-hclust(vegan::vegdist(t(sqrt(otu_table(data.prop))),method = "bray"))
c<-as.dendrogram(c)
labels_colors(c) <- colors[order.dendrogram(c)]
png(file="Results/Stool/Hierarchical_clustering/Hierarchical_cluster_Bray_sqrt_prop_normalized_ASV.png",width=2500,height=1800, res=300)
par(mar=c(6,4,5,4), cex.lab=1, cex.main=1.35, cex.sub=1.3)
plot(c,main="Community structure of Stool subset using\n Bray-Curtis distance on sqrt proportional ASVs",
     sub="Healthy = green     ICC = orange")
dev.off()

suppressWarnings(rm(c,color_table,labels_colors))


########################## ALFA DIVERSITY Stool ############################

# no normalization according to phyloseq https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

pAlpha<-plot_richness(data, measures=c("Shannon", "Observed"), x="Condition", color="Condition")
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

pAlpha + geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0) + theme_classic2() + 
  scale_color_manual(values = c("Healthy"="chartreuse", "ICC"="coral")) +
  labs(x="Condition", title="Alpha diversity between ICC and Healthy (Stool)") +
  guides(fill="none", color="none") + theme(axis.text.x= element_text(angle=0, vjust=0, hjust=0.5, size=10)) +
  stat_compare_means(aes(group = Condition), label="p.format", method = "wilcox.test", label.x= 1.4, size=3.5, 
                     label.y.npc = "top", vjust=-0.5, hjust=0.3) +
  theme(panel.grid.major.y = element_line(size=0.4), panel.grid.minor.y = element_line(size=0.25) , 
        axis.title.x = element_text(vjust = -1))
ggsave(file="Results/Stool/Alfa_diversity_with_Mann_Withn_Wilcox_Stool.png", width = 6,height =5.3, dpi=300)

# just to test the plotted p value
alphadt<- as.data.frame(pAlpha$data) # the column "value" contains the alpha value
Obser_value<-filter(alphadt, variable=="Observed richness")
factor<-Obser_value$Condition
wilcox.test(Obser_value$value~factor)

rm(pAlpha, alphadt, H, ev, obs, Obser_value, New_data, factor)


##################### BETA DIVERSITY  Stool #######################

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
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
perm_ASV$aov.tab$`Pr(>F)`[1]
perm_ASV_H<-perm_ASV$aov.tab$`Pr(>F)`[1] # needed later for plot

sample_OTU<-as.data.frame(t(sqrt(ASV.genus.prop)))
perm_g<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")
perm_g$aov.tab$`Pr(>F)`[1]
perm_g_H<-perm_g$aov.tab$`Pr(>F)`[1] 

sample_OTU<-as.data.frame(t(sqrt(ASV.fam.prop)))
perm_f<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.class.prop)))
perm_o<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.order.prop)))
perm_c<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

sample_OTU<-as.data.frame(t(sqrt(ASV.phy.prop)))
perm_p<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="euclidean")

beta<-rbind(perm_ASV$aov.tab[1,],perm_g$aov.tab[1,],perm_f$aov.tab[1,],perm_o$aov.tab[1,],perm_c$aov.tab[1,],perm_p$aov.tab[1,])
row.names(beta)<-c("Raw_ASV","Genera","Families","Orders","Classes","Phyla")
beta
write.csv2(beta, file="Results/Stool/Beta_div/Beta_diversity_permanova_Hellinger_Stool.csv",quote=F,row.names = T)


### PLUS: checking it with Jaccard too
suppressWarnings(rm(sample_OTU,perm_ASV))
sample_OTU<-as.data.frame(t(sqrt(ASV.prop))) # samples has to be on rows --> t
perm_ASV<- vegan::adonis(sample_OTU ~Condition, data=metadata, permutations = 9999, method="jaccard")
write.csv2(as.data.frame(perm_ASV$aov.tab), file="Results/Stool/Beta_div/Beta_divers_permanova_JACCARD_on_ASV.csv",quote=F,row.names = T)


# Perform an ANCOVA-like test to determine if the variances differ by groups --> betadispersion test on distance
# on ASV
BC.dist<-vegan::vegdist(t(sqrt(ASV.prop)), distance="euclidean")
disper<-vegan::betadisper(BC.dist,metadata$Condition)
disp_ASV<-vegan::permutest(disper, permutations=9999)
disp_ASV
write.csv2(disp_ASV$tab, file="Results/Stool/Beta_div/Beta_dispersion_permanova_Helling.csv",quote=F,row.names = T)

rm(beta, perm_g, perm_f, perm_o, perm_c, perm_p, perm_ASV)


########################### PCoA BRAY CURTIS Stool #####################

# on ASV
data.prop.labels<-data.prop
{data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues # to get eigen values of every PC
  eigval<- round((eigval/sum(eigval))*100, 1) # to get variation explained by every PC
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on Stool subset", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"))
ggsave(file="Results/Stool/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_Stool.png", width = 9, height = 6, dpi=300)
# without ellipses
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3) + theme_classic(base_size = 14) + 
  geom_text(aes(label=sample_names(data.sqrt_prop)), color="black", size=3, show.legend = FALSE) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on Stool subset",
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
ggsave(file="Results/Stool/Beta_div/PCoA_Beta_diversity_Hellinger_ASV_no_ellipse.png", width = 9, height = 6, dpi=300)
# without names
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3, alpha=0.5) + theme_classic(base_size = 14) + stat_ellipse(size=0.2) +
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed ASV)\n computed on Stool subset", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_ASV_H))
ggsave(file="Results/Stool/Beta_div/PCoA_Beta_diversity_Hellinger_on_ASV_points.png", width = 9, height = 6, dpi=300)

rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels)


### again but on genera
# Hellinger
{data.prop.labels<-data.genus.prop
  data.sqrt_prop<-transform_sample_counts(data.prop.labels, sqrt) # square root of proportion
  DistBC = phyloseq::distance(data.sqrt_prop, method = "euclidean")
  ordBC = ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)
  eigval<-ordBC$values$Eigenvalues
  eigval<- round((eigval/sum(eigval))*100, 1)
}
plot_ordination(data.sqrt_prop, ordBC, color = "Condition") +
  scale_color_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  geom_point(size=3, alpha=0.5) + theme_classic(base_size = 13.4) + stat_ellipse(size=0.2) + 
  labs(title="PCoA with Hellinger distance (euclidean on Hellinger transformed genera)\n computed on Stool subset", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation"), subtitle = paste("PERMANOVA Pr(>F) =",perm_g_H))
ggsave(file="Results/Stool/Beta_div/PCoA_Beta_diversity_Hellinger_on_genera.png", width = 9, height = 6, dpi=300)

suppressWarnings(rm(data.sqrt_prop_perm, eigval, DistBC, ordBC,data.prop.labels))


#################### VEEN DIAGRAM STOOL ##########################

data.genus.temp<-data.genus.prop
tax_table(data.genus.temp)<-as.matrix(Taxa.genus.update)
data.venn<-data.genus.temp

### abundance AND prevalence filter (0.5% abundance at least in 3 sample)
who<-as.data.frame(otu_table(data.venn))
who<-apply(who, MARGIN=2, function(x) ifelse(x>0.5, 1, 0)) # if more than 0.5 --> "a point"
who<-who[!rowSums(who)>2,] # more than 2 "points" --> at least in 3 samples
who<-as.vector(tax_table(data.venn)[row.names(who),"Genus"]) # to translate in ASV code (of genus collapsed object)
data.venn<-subset_taxa(data.venn, ! Genus %in% who)


Healthy<-subset_samples(data.venn, Condition=="Healthy")
Healthy<-as.character(tax_table(prune_taxa(taxa_sums(Healthy)>0, Healthy))[,"Genus"])

ICC<-subset_samples(data.venn, Condition=="ICC")
ICC<-as.character(tax_table(prune_taxa(taxa_sums(ICC)>0, ICC))[,"Genus"])


ONLY_IN_ICC<- ICC[! ICC %in% Healthy]
ONLY_IN_ICC<- paste(ONLY_IN_ICC, collapse = ", ")
head(ONLY_IN_ICC)

ONLY_IN_Healthy<- Healthy[! Healthy %in% ICC]
ONLY_IN_Healthy<- paste(ONLY_IN_Healthy, collapse = ", ")
head(ONLY_IN_Healthy)


x<-list(Healthy=Healthy,ICC=ICC)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, show_percentage = F,
       fill_color = c("chartreuse","coral")) +
  theme(plot.title = element_text(size=10), plot.caption = element_text(size=7) ) +
  labs(title = "Venn Diagram of Stool subset \n(only genera with minimal abundance > 0.5% \n at least in three samples)",
       caption=paste("Taxa only in Healthy subjects:",ONLY_IN_Healthy) )
ggsave(filename = "Results/Stool/Beta_div/Venn_Diagramm.png", width = 4, height = 4, dpi=300, bg = "white")
dev.off()


suppressWarnings(rm(ONLY_IN_ICC, ONLY_IN_Healthy, x, con, ICC, Healthy, data.venn, who))



################### DA WITH DESEQ2 Stool #####################

if(! "proof1" %in% ls() ){
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
  DEseq_data<-phyloseq_to_deseq2(d, ~Condition)
  DE<-DESeq(DEseq_data)
  res<-results(DE, contrast= c("Condition", "Healthy", "ICC"))
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
    # write.csv2(r, file=paste0("Results/Stool/DA_DESeq2/DA_",t,"_ratio_Healthy_vs_ICC_Stool.csv"), row.names = F, quote=F, na = "")
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
write.csv(Res_tot, file="Results/Stool/DA_DESeq2/Every_result_DESeq2_Stool.csv", row.names = F)
write.xlsx(Res_tot, file="Results/Stool/DA_DESeq2/Every_result_DESeq2_Stool.xlsx", showNA = F, col.names = T)


ggplot(Table_tot, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  scale_fill_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=12), plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank() ) +   
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,36,2), seq(38,max(Table_tot$Abundance+3),3))) +
  labs(title= "Differently abundant Taxa (Stool)", y="Proportional Abundance", 
       fill="Condition", x="")
ggsave(filename = "Results/Stool/DA_DESeq2/DA_Condition_every_result_Stool.png", width = 15, height = 8, dpi=300)
dev.off()

# removing redundants
Table_tot2<-Table_tot[!(Table_tot$Bacteria=="Clostridia_vadinBB60_group" & Table_tot$Taxa %in% c("Family","Order","Class") ), ]
Redund<- c("NA_ f Erysipelotrichaceae","Selenomonadaceae", "Tannerellaceae","Streptococcaceae","Bacteroidia")
Table_tot2<-subset(Table_tot2, ! Bacteria %in% Redund) # to remove redundant results
DE_plot<-ggplot(Table_tot2, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")), scales = "free_x", space="free") +
  geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
  scale_fill_manual(values=c("Healthy"="chartreuse","ICC"="coral")) +
  theme(strip.text.x=element_text(size=14,colour="black")) + 
  guides( fill=guide_legend(nrow=1) ) +
  theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
        legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 30, vjust=1, hjust=1, size=12), 
        axis.text.y = element_text(size=10), 
        plot.title= element_text(size=18),
        panel.grid.minor.y= element_blank(),
        plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
  scale_x_discrete(expand=c(-0.2, 1)) +
  scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,36,2), seq(38,max(Table_tot$Abundance+3),3))) +
  labs(title= "Differently abundant Taxa (Stool)", y="Proportional Abundance", 
       fill="Condition", x="")

DE_plot
ggsave(filename = "Results/Stool/DA_DESeq2/DA_Condition_WITHOUT_redundants_Stool.png", width = 12, height = 8.5, dpi=300)
dev.off()

system(" echo 'Every result under the arbitrary threshold of basemean=50 has been removed in order to avoid the most noisy results' > Results/Stool/DA_DESeq2/NB.txt ")

# for further comparisons ...
Genera.DESEQ2_Stool<-unique(Table_tot[Table_tot$Taxa=="Genus","Bacteria"])



##### CIRCOPLOT TO PLOT THE LOG2FC

### preparing the object
Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
### removing the redundant results
Res_tot2<-Res_tot2[!(Res_tot2$Genus %in% "NA_ f Erysipelotrichaceae" & Res_tot2$Taxon %in% "Genus"), ] # almost absent
Res_tot2<-Res_tot2[!(Res_tot2$Family %in% Redund & Res_tot2$Taxon %in% "Family"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Order %in% Redund & Res_tot2$Taxon %in% "Order"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Class %in% Redund & Res_tot2$Taxon %in% "Class"), ]
Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% Redund & Res_tot2$Taxon %in% "Phylum"), ]
Res_tot2$Genus<-gsub("CAG-352"," CAG_352",Res_tot2$Genus)

# plotting
circo<-CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),tax_col=8:12,fc_col=2,sort="no")

png(file="Results/Stool/DA_DESeq2/CircoTax_plot.png",width=1800,height=1800, res = 300)
circo + labs(title="ICC vs Healthy \n Stool DA taxa log2foldchange\n")
dev.off()

circo_no_m <- circo + theme(plot.margin= unit(c(0,0,0,0.2), unit="cm"))

# diff abundances + logfoldchange
png(file="Results/Stool/DA_DESeq2/Final_plot.png",width=4560, height=2400, res = 300)
ggarrange(DE_plot, circo_no_m, nrow = 1, widths= c(1,0.58), labels = c("A)", "B)"))
dev.off()

rm(circo, circo_no_m)


################# CORRELATIONS WITH FFA VALUES (Stool) ######################

suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data.genus.prop, Condition != "Healthy") # many variables have the same value for each Healthy subjects --> not usable 
prune.dat_target <- subset_taxa(data.temp, Genus %in% Genera.DESEQ2_Stool)
sample_names(prune.dat_target) <- sample_data(data.temp)$Subject
target <-as.data.frame(otu_table(prune.dat_target))
target_tax<-as.data.frame(tax_table(prune.dat_target))
identical(row.names(target_tax),row.names(target))
row.names(target)<-target_tax$Genus
colnames(target)
rm(target_tax)
target<-as.data.frame(t(target))

SCFA <- read.table("FFA_values.tsv", header=T, sep="\t", dec=",")
SCFA[ , !colnames(SCFA)%in%c("Sample_Code","Condition")] <- apply(SCFA[, !colnames(SCFA)%in%c("Sample_Code","Condition")], 2, as.numeric)

# same order with the phyloseq data
row.names(SCFA)<-SCFA$Sample_Code
SCFA<-SCFA[row.names(target), ]
length(SCFA$Sample_Code[! is.na(SCFA$Sample_Code)]) == 7  # TRUE --> any sample cutted avay (there are FFA values regarding 8 ICC but no stool for one of them)
SCFA <- SCFA[! is.na(SCFA$Sample_Code), ]

target<-target[row.names(SCFA), ] # to re-subset according to missing subjects among SCFA values


# Computing the correlations...
x<-cbind.data.frame(target,SCFA)
{x<-as.matrix(x[ ! colnames(x) %in% c("Sample_Code","Condition")])
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  identical(correlation_corr[,1:2],correlation_pvalue[,1:2])
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("DA_Genera","FFA","Corr","pvalue")
}
corr<-subset(correlation, DA_Genera %in% colnames(target))
corr<-subset(corr, FFA %in% colnames(SCFA))
colnames(corr)<-c("DA_Genera","SCFA","Corr","pvalue")
corr$padj<-p.adjust(corr$pvalue, method = "holm")
corr$Sign<-corr$padj
corr$Sign[corr$Sign < 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""
ggplot(corr, aes(x = corr$DA_Genera, y = corr$SCFA, fill = corr$Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "Serum DA SCFA", x= "Stool DA Genera", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Stool/DA_genera_vs_SCFA_using_only_ICC.csv.png", dpi=300, width = 7, height = 7)
write.csv2(corr, file="Results/Stool/DA_genera_vs_SCFA_using_only_ICC.csv", quote = F, na="", row.names = F)

corr2<-corr[corr$SCFA!="X2.MethylButyrric", ]
corr2<-corr2[! corr2$DA_Genera %in% c("Alloprevotella","Mitsuokella"), ]
ggplot(corr2, aes(x = DA_Genera, y = SCFA, fill = Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle = -25, hjust = 0, size= 11)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =6) +
  labs(title = "", y= "Serum DA SCFA", x= "Stool DA Genera", 
       caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(1, 'cm'), legend.title = element_text(size=12))
ggsave(file="Results/Stool/DA_genera_vs_SCFA_using_only_ICC_WITHOUT_NAs.csv.png", dpi=300, width = 7, height = 7)


################### ANALYSING LEFSE RESULTS STOOL ###########################

suppressWarnings(rm(a))
a <- read.delim("QIIME/PICRUST2_LEFSE/picrust2/KEGG_pathways/path_abun_unstrat_descriptions.tsv.gz") 
Descriptions<-a[,c("pathway","description")]

Significative_functions_LEFSE<- read.delim("QIIME/PICRUST2_LEFSE/Result_LEFSE_STOOL.res", header=FALSE)
head(Significative_functions_LEFSE, n=4)
head(Descriptions$pathway, n=4) # just a further checking of the row matching (different text format but same information)
colnames(Significative_functions_LEFSE)<-c("Pathway","Log_highest_mean","Class_with_highest_mean","logLDA_score","p-value")
Significative_functions_LEFSE$Pathway<-Descriptions$description
Significative_functions_LEFSE$Pathway_ID<-Descriptions$pathway
head(Significative_functions_LEFSE, n=4)

# if sensed as significative then the logLDA is displayed
Significative_functions_LEFSE<-Significative_functions_LEFSE[!is.na(Significative_functions_LEFSE$logLDA_score),]
write.xlsx(Significative_functions_LEFSE, file="Results/Stool/PICRUST2_LEFSE_results/Significantly_different_pathways_KEGG_Stool.xlsx", col.names = T, showNA = F, row.names = F)


# modifing names for the plot)
Significative_functions_LEFSE$Pathway<-gsub("&beta;-","", Significative_functions_LEFSE$Pathway, fixed = T)
{Significative_functions_LEFSE$Pathway<-paste("",Significative_functions_LEFSE$Pathway,"") # needed to distance text from lines
  Significative_functions_LEFSE<-Significative_functions_LEFSE[order(abs(as.numeric(Significative_functions_LEFSE$logLDA_score))), ] # order based on the effect size
  Significative_functions_LEFSE$Pathway<-factor(Significative_functions_LEFSE$Pathway, levels = Significative_functions_LEFSE$Pathway) # to prevent alphabetical re-sorting
}
# inverting the values of a group to make a simmetric plot
Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="ICC"]<- -Significative_functions_LEFSE$logLDA_score[Significative_functions_LEFSE$Class_with_highest_mean=="ICC"]


# plotting every result
ggplot(data=Significative_functions_LEFSE, aes(y=Pathway, x=as.numeric(logLDA_score), 
                                               fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), 
            hjust = ifelse(Significative_functions_LEFSE$Class_with_highest_mean=="Healthy",1,0), size=3.2) +
  scale_fill_manual(values=c("ICC"="coral", "Healthy"="chartreuse")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15))+
  theme(legend.position = "bottom")
ggsave(filename = "Results/Stool/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_KEGG_every_result.png", width = 8, height = 6, dpi = 300)

#plotting only results over 3
Significative_functions_LEFSE_3<-Significative_functions_LEFSE[abs(Significative_functions_LEFSE$logLDA_score)>3,]
ggplot(data=Significative_functions_LEFSE_3, aes(y=Pathway, x=as.numeric(logLDA_score), fill=Class_with_highest_mean)) +
  geom_bar(stat = "identity", width = 1, color="black" ) + labs(x="log LDA score", y="", fill="") +
  geom_text(aes(x= 0, label=Pathway), hjust = ifelse(Significative_functions_LEFSE_3$Class_with_highest_mean=="Healthy",1,0), size=3.2) +
  scale_fill_manual(values=c("ICC"="coral", "Healthy"="Chartreuse")) +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(-3,-2,-1, 1, 2, 3))+
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.3)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 15))+
  theme(legend.position = "bottom")
ggsave(filename = "Results/Stool/PICRUST2_LEFSE_results/PICRUST2_LEFSE_plot_diff_KEGG_over_3.png", width = 7, height = 2.5, dpi = 300)

rm(Significative_functions_LEFSE, a)


#################### BOX PLOT FOR EACH FFA QUANTITIES ####################

SCFA <- read.table("FFA_values.tsv", header=T, sep="\t", dec=",")
colnames(SCFA) <- gsub("X2","2-",colnames(SCFA))
colnames(SCFA) <- gsub(".","",colnames(SCFA), fixed=T)
head(SCFA, n=2)
table_SCFA<-reshape::melt(SCFA)

ggplot(data=table_SCFA, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("Healthy"="chartreuse","ICC"="coral")) +
  theme_bw(base_size = 15) +
  theme( axis.text.x = element_text(angle=-20, size = 9, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7.5),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.1, 0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          100, 250, 500)) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.55, size=0.3, outlier.size = 1) + 
  labs(title="FFA concentrations", x="", fill="", y="Concentration (µM/L)")
ggsave(filename = "Results/FFA_values_con_ZOOM_asse_y_ZOOM.png", height=4.5, width = 8.5, dpi=300)

# without 2-MethBut (same value for every subject in each group)
table_SCFA2<-table_SCFA[table_SCFA$variable!="2-MethylButyrric", ]
ggplot(data=table_SCFA2, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("Healthy"="chartreuse","ICC"="coral")) +
  theme_bw(base_size = 15) +
  theme( axis.text.x = element_text(angle=-20, size = 9, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7.5),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.1, 0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          100, 250, 500)) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.55, size=0.3, outlier.size = 1) + 
  labs(title="FFA concentrations", x="", fill="", y="Concentration (µM/L)")
ggsave(filename = "Results/FFA_values_con_ZOOM_asse_y_ZOOM_without_2MethButirr.png", height=4.5, width = 8.5, dpi=300)


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
package$otherPkgs$pca3d[c(1,4)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$mixOmics[c(1,4)]
cat("\n \n \nEvery package: \n", fill=TRUE)
print(package$otherPkgs)

sink()
close(con)
suppressWarnings(rm(con))
