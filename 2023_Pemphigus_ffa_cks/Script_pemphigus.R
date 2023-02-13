###################### PREPARING DATA #########################

{library(Hmisc)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(mixOmics)
library(ecodist)
library(reshape2)
library(readxl)
}

options(scipen=100)

CKS <- read_excel("mio_CKS.xlsx")
FFA <- read_excel("mio_FFAs.xlsx") # pre treatment
FFA_post <- read_excel("mio_FFAs_post_treatment.xlsx") # post treatment


for(x in c("IFN-gamma","IL17A","IL-10","IL21","IL-5")) {
  CKS[ , x]<-as.numeric(CKS[[x]]) }
rm(x)

CKS$Condition<-factor(CKS$Condition, levels = c("Healthy","PEMFIGO"))

FFA[! colnames(FFA) %in% c("Subject","Condition")]<-apply(FFA[!colnames(FFA)%in% c("Subject","Condition")], MARGIN = 2, as.numeric)
FFA$Condition<-factor(FFA$Condition, levels = c("Healthy","PEMFIGO"))

SCFA<-FFA[,1:9] # from acetic to valeric 
MCFA<-FFA[,c(1,2,10:15)] # to dodecanoic
LCFA<-FFA[,c(1,2,16:18)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

identical(colnames(FFA_post[3:18]), colnames(FFA[3:18])) # TRUE
SCFA_post<-FFA_post[,1:9]
MCFA_post<-FFA_post[,c(1,2,10:15)]
LCFA_post<-FFA_post[,c(1,2,16:18)]


PDAI <- read.csv("mio_PEMFIGO_PDAI.csv")
# loading info regarding PDAI
row.names(PDAI)<- PDAI$Subject
# subsetting and reordering the samples
PDAI_CKS<-PDAI[PDAI$Subject %in% CKS$Subject, ] 
PDAI_CKS<-PDAI_CKS[CKS$Subject %in% PDAI$Subject, ] 
identical(PDAI_CKS$Subject, CKS$Subject[CKS$Subject %in% PDAI_CKS$Subject])
PDAI_CKS<-PDAI_CKS[! PDAI_CKS$PDAI_gravity=="NOT_KNOWN",]
PDAI_FFA<-PDAI[FFA$Subject %in% PDAI$Subject,]
identical(PDAI_FFA$Subject, FFA$Subject[FFA$Subject %in% PDAI_FFA$Subject])
PDAI_FFA<-PDAI_FFA[! PDAI_FFA$PDAI_gravity=="NOT_KNOWN",]


dir.create("Results")
setwd("Results")


######################### NORMALIZATIONS ##############################

CKS_raw<-CKS
SCFA_raw<-SCFA
MCFA_raw<-MCFA
LCFA_raw<- LCFA

# #### normalization through relative abundance in each sample, e.g. https://www.nature.com/articles/s41598-021-81420-3
# CKS[ , c("IFN-gamma","IL17A","IL-10","IL21","IL-5")]<-CKS[ , c("IFN-gamma","IL17A","IL-10","IL21","IL-5")]/rowSums(CKS[ , c("IFN-gamma","IL17A","IL-10","IL21","IL-5")])
# rowSums(CKS[ , c("IFN-gamma","IL17A","IL-10","IL21","IL-5")])
# 
# #### normalization of FFA respect total SCFA, MCFA or LCFA
# SCFA[ !colnames(SCFA)%in% c("Subject","Condition")]<-SCFA[!colnames(SCFA)%in% c("Subject","Condition")]/rowSums(SCFA[!colnames(SCFA)%in% c("Subject","Condition")])
# rowSums(SCFA[!colnames(SCFA)%in% c("Subject","Condition")])
# 
# MCFA[ !colnames(MCFA)%in% c("Subject","Condition")]<-MCFA[!colnames(MCFA)%in% c("Subject","Condition")]/rowSums(MCFA[!colnames(MCFA)%in% c("Subject","Condition")])
# rowSums(MCFA[!colnames(MCFA)%in% c("Subject","Condition")])
# 
# LCFA<- LCFA[ , !colnames(LCFA)=="Octadecanoic"] # it has a NA in sample HS58 
# LCFA[ !colnames(LCFA)%in% c("Subject","Condition")]<- LCFA[!colnames(LCFA)%in% c("Subject","Condition")]/rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])
# rowSums(LCFA[!colnames(LCFA)%in% c("Subject","Condition")])


####################### MANN WITHNEY TESTs ####################

 
dir.create("CKS")

# CKS raw
results<-NULL
for(x in c("IFN-gamma","IL17A","IL-10","IL21","IL-5")) {
  test<-wilcox.test(CKS_raw[[x]]~CKS$Condition)
  wilc<-cbind(paste(x),test$statistic,test$p.value)
  results<-rbind.data.frame(results,wilc)
}
colnames(results)<-c("Cytokine","W-statistic","p.value")
results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
row.names(results)<-NULL
results
write.csv2(results, file="CKS/Mann_Whitney_raw_CKS.csv", row.names = F)

table <- melt(data = CKS_raw)
ggplot(data=table, mapping=aes(x=variable, y=value, fill=Condition)) +
  geom_boxplot() + labs(title="Raw cytokines quantities", x="", y="Quantity") +
  theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05)) +
  scale_fill_manual(values = c("Healthy"="deepskyblue","PEMFIGO"="coral"))
ggsave(filename = "CKS/raw_CKS_quantity.png", height=3, width = 6, dpi=300)

table2<- table[! c(table$Subject=="19_267" & table$variable=="IL21"), ]
ggplot(data=table2, mapping=aes(x=variable, y=value, fill=Condition)) +
  geom_boxplot() + scale_y_continuous(breaks = seq(0, 40, 5)) +
  theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05)) +
  scale_fill_manual(values = c("Healthy"="deepskyblue","PEMFIGO"="coral")) +
  labs(title="Raw cytokines quantities", x="", y="Quantity", subtitle = "without the IL-21 measurement of subject 19_267")
ggsave(filename = "CKS/raw_CKS_quantity_without_THE_outlier.png", height=3, width = 6, dpi=300)


# # CKS relative abundance
# results<-NULL
# for(x in c("IFN-gamma","IL17A","IL-10","IL21","IL-5")) {
#   test<-wilcox.test(CKS[[x]]~CKS$Condition)
#   wilc<-cbind(paste(x),test$statistic,test$p.value)
#   results<-rbind.data.frame(results,wilc)
# }
# colnames(results)<-c("Cytokine","W-statistic","p.value")
# results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
# row.names(results)<-NULL
# results
# write.csv2(results, file="Mann_Whitney_relative_abund_CKS.csv", row.names = F)
# 
# table <- melt(data = CKS)
# ggplot(data=table, mapping=aes(x=variable, y=value, fill=Condition)) +
#   geom_boxplot() + labs(title="Raw cytokines quantities", x="", y="Relative quantity") +
#   scale_y_continuous(breaks = seq(0, 1, 0.20))
# ggsave(filename = "relative_CKS_quantity.png", height=3, width = 6, dpi=300)




dir.create("FFA")

##### FFA raw abundances

list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Subject","Condition")]) {
    test<-wilcox.test(temp[[x]]~temp$Condition)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit.csv"), row.names = F)
  
  table <- melt(data = temp)
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Condition)) +
    scale_fill_manual(values = c("Healthy"="deepskyblue","PEMFIGO"="coral")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05)) +
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance.png"), height=3, width = 6, dpi=300)
}


########################### CORRELATIONS ##########################

dir.create("Correlations")


######### BETWEEN FFA AND IL

# cutting same subjects from data sets

common_FFA <- as.data.frame(FFA[FFA$Subject %in% CKS$Subject, ])
row.names(common_FFA)<- common_FFA$Subject
common_FFA <- common_FFA[CKS$Subject,]
identical(common_FFA$Subject,CKS$Subject) # TRUE

head(common_FFA, n=3)

SCFA_common<-common_FFA[,3:9] # from acetic to valeric 
MCFA_common<-common_FFA[,c(10:15)]
LCFA_common<-common_FFA[,c(16:18)]  

### starting correlation

rm(list)
list<-list("SCFA"=SCFA_common, "MCFA"=MCFA_common, "LCFA"=LCFA_common)
             #1                 #2                  #3

for(y in 1:length(list)){
  rm(x, correlation, data_corr)
  x<-cbind.data.frame(list[[y]], CKS_raw)
  
  PEM<-x[x$Condition=="PEMFIGO", ! colnames(x) %in% c("Subject","Condition") ] # cutting only values of PEMFIGO subjects
  Healthy<-x[x$Condition=="Healthy", ! colnames(x) %in% c("Subject","Condition") ] # cutting only values of Healthy subjects
  
  # sub loop to correlate values of every subject, PEM subjects or Healthy subjects
  sub_list<- list("every"=x[, ! colnames(x) %in% c("Subject","Condition")], "PEMFIGO"=PEM, "Healthy"=Healthy)
                  #1                                                          #2            #3           
  for(z in 1:length(sub_list)) {
  r<-rcorr(as.matrix(sub_list[[z]]), type = "spearman")
  data_corr<-as.data.frame(as.table(r$r))
  data_pvalue<-as.data.frame(as.table(r$P))
  correlation<-cbind(data_corr,data_pvalue[,3])
  colnames(correlation)<-c(paste(names(list)[y]),"CKS","Corr","pvalue")
  correlation<-correlation[ correlation[[1]] %in% colnames(list[[y]]), ]
  correlation<-correlation[ correlation$CKS %in% colnames(CKS_raw), ]
  correlation$'padj(BH)'<-p.adjust(correlation$pvalue, method = "BH")
  correlation$Sign<-correlation$`padj(BH)`
  correlation$Sign[correlation$Sign<0.05]<-"*"
  correlation$Sign[correlation$Sign>0.05]<-" "
  write.csv2(correlation, file = paste0("Correlations/Spearm_",names(list)[y]," vs CKS_in ",names(sub_list)[z]," subject.csv"), row.names=F, quote=F)
  
  # heatmap
  ggplot(correlation, aes(x = correlation[[1]], y = correlation[[2]], fill = Corr)) +
    geom_tile(color = "white", lwd = 0.5, linetype = 1) +
    scale_fill_gradient2(low="blue", mid = "white", high="red") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 9)) +
    guides(fill= guide_colourbar(title ="rho")) +
    geom_text(aes(label= Sign), color= "white", size =4.5) +
    labs(title = "Spearman correlation between",
         subtitle = paste0(names(list)[y]," and cytokines raw abundances\n in ",names(sub_list)[z]," subject"),
         y= "Cytokines", x= paste(names(list)[y]), caption = "p-adj lower than 0.05 are displayed through * symbol") +
    theme(legend.text = element_text(size=10), legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(0.6, 'cm'), legend.title = element_text(size=10))
  ggsave(filename = paste0("Correlations/Heatmap_",paste(names(list)[y]," vs CKS_in ",names(sub_list)[z],"subject.png")), width = 5, height = 4, dpi= 300 )
  }
}

#### following lines are made just for the manual test of the significative result 
# y= 2 # MCFA
# z= 1 # every
# # --> return to x creation and then...
# cor.test(x$IsoHexanoic,CKS_raw$IL21, method = "spearman") # OK!

rm(correlation,x,z,y)



########### BETWEEN PDAI AND IL

CKS_PDAI<-as.data.frame(CKS)
rownames(CKS_PDAI)<-CKS_PDAI$Subject
CKS_PDAI<-CKS_PDAI[rownames(PDAI_CKS), ]
identical(CKS_PDAI$Subject,PDAI_CKS$Subject)
x<-cbind(PDAI_CKS, CKS_PDAI)   # the order of the samples has been already checked in the first paragraph
x<-x[, ! colnames(x) %in% c("Subject","PDAI_gravity","Condition")]

r<-rcorr(as.matrix(x), type = "spearman")
data_corr<-as.data.frame(as.table(r$r))
data_pvalue<-as.data.frame(as.table(r$P))
correlation<-cbind(data_corr,data_pvalue[,3])
colnames(correlation)<- c("PDAI","CKS","corr","pvalue")
correlation<-correlation[ correlation[[1]] =="PDAI_pre", ]
correlation<-correlation[ correlation[[2]] %in% colnames(CKS_PDAI), ]
correlation$'padj(BH)'<-p.adjust(correlation$pvalue, method = "BH")
correlation$Sign<-correlation$`padj(BH)`
correlation$Sign[correlation$Sign<0.05]<-"*"
correlation$Sign[correlation$Sign>0.05]<-" "

write.csv2(correlation, file = "Correlations/Spearm_PDAI_vs_CKS_.csv", row.names=F, quote=F)

# heatmap
ggplot(correlation, aes(x = correlation[[1]], y = correlation[[2]], fill = corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle = -20, hjust = 0, size= 9)) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "white", size =4.5) +
  labs(title = "Spearman correlation between",
       subtitle = "PDAI and cytokines raw abundances",
       y= "Cytokines", x= "PDAI", caption = "p-adj lower than 0.05 are displayed through * symbol") +
  theme(legend.text = element_text(size=10), legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.6, 'cm'), legend.title = element_text(size=10))
ggsave(filename = "Correlations/Heatmap_PDAI_vs_IL_subject.png", width = 5, height = 4, dpi= 300 )



########################## PCoA ##############################

dir.create("PCoA")


############### CKS

dir.create("PCoA/based_on_CKS")

CKS_PCoA<-as.data.frame(CKS_raw)
row.names(CKS_PCoA)<-CKS_PCoA$Subject
dist<-vegdist(CKS_PCoA[, !colnames(CKS_PCoA) %in% c("Subject","Condition")], method = "bray") # NB: sample on rows
coord<- pco(dist)

ggplot(data=as.data.frame(coord$vectors), aes(x=V1, y=V2, color=CKS_PCoA$Condition) ) + 
  geom_point(size=3) + theme_classic() +  stat_ellipse() +
  geom_text(aes(label=row.names(CKS_PCoA)),size=2, color="black")+
  labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
       y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
       color="Condition", title = "PCoA computed on CKS abundances\n of healthy and PEMFIGO subjects\n with Bray Curtis index")
ggsave(file="PCoA/based_on_CKS/PCoA_BRAY_CKS_Healthy_PEMFIGO.png", height = 5, width = 6, dpi=300)              

# checking for significant dispersion
permutest(betadisper(dist, CKS_PCoA$Condition), permutations = 9999)
system("echo 'Despite appearences there is not a significant difference \nin dispersion between Healthy and PEMFIGO' > PCoA/based_on_CKS/Dispersion.txt")

# only PERMFIGO (for PDAI)
rm(dist,coord)
CKS_PCoA<-CKS_PCoA[CKS_PCoA$Condition =="PEMFIGO",]
dist<-vegdist(CKS_PCoA[, !colnames(CKS_PCoA) %in% c("Subject","Condition")], method = "bray")
coord<- pco(dist)

ggplot(data=as.data.frame(coord$vectors), aes(x=V1, y=V2, color=PDAI_CKS$PDAI_gravity) ) + 
  geom_point(size=3) + theme_classic() +  stat_ellipse() +
  scale_color_manual(values=c("Low_PDAI"="coral", "High_PDAI"="red4")) +
  geom_text(aes(label=row.names(PDAI_CKS)),size=2, color="black")+
  labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
       y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
       color="Condition", title = "PCoA computed on CKS abundances\n of PEMFIGO subjects (pre treatment)\n with Bray Curtis index")
ggsave(file="PCoA/based_on_CKS/PCoA_BRAY_CKS_PDAI.png", height = 5, width = 6, dpi=300)              


########### on FFA subsets

dir.create("PCoA/based_on_FFA")

list_PCoA<-list("SCFA"=SCFA_raw,"MCFA"=MCFA_raw,"LCFA"=LCFA_raw)
                  #1               #2            #3

p_values<-NULL
for(x in 1:length(list_PCoA)){
  rm(data, coord)
  data<-as.data.frame(list_PCoA[[x]])
  row.names(data)<-data$Subject
  if(x==3) {data<-data[!is.na(data$Octadecanoic),] } # there is a NA in LCFA
  dist<-vegdist(data[, !colnames(data) %in% c("Subject","Condition")], method = "bray") # NB: sample on rows
  coord<- pco(dist)
  
  ggplot(data=as.data.frame(coord$vectors), aes(x=V1, y=V2, color=data$Condition) ) + 
    geom_point(size=3) + theme_classic() +  stat_ellipse() +
    geom_text(aes(label=row.names(data)),size=2, color="black")+
    labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
         y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
         color="Condition", title = paste("PCoA computed on",names(list_PCoA)[x],"abundances\n of healthy and PEMFIGO subjects\n with Bray Curtis index"))
  ggsave(file=paste("PCoA/based_on_FFA/PCoA_BRAY_",names(list_PCoA)[x],"Healthy_PEMFIGO.png", sep="_"), height = 5, width = 6, dpi=300)              
  
  # checking for significant dispersion
  rm(p_value_diver, p_value_disp)
  p_value_diver<-adonis(dist ~ data$Condition, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Condition), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
 
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
  
  # only PERMFIGO (for PDAI)
  rm(dist,coord)
  data<-data[data$Condition =="PEMFIGO",]
  data<-data[data$Subject %in% PDAI_FFA$Subject, ] # there are subjects with unknown PDAI
  dist<-vegdist(data[, !colnames(data) %in% c("Subject","Condition")], method = "bray")
  coord<- pco(dist)
  
  ggplot(data=as.data.frame(coord$vectors), aes(x=V1, y=V2, color=PDAI_FFA$PDAI_gravity) ) + 
    geom_point(size=3) + theme_classic() +  stat_ellipse() +
    scale_color_manual(values=c("Low_PDAI"="coral", "High_PDAI"="red4")) +
    geom_text(aes(label=row.names(PDAI_FFA)),size=2, color="black")+
    labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
         y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
         color="Condition", title = paste("PCoA computed on",names(list_PCoA)[x],"abundances\n of PEMFIGO subjects\n with Bray Curtis index"))
  ggsave(file=paste("PCoA/based_on_FFA/PCoA_BRAY_",names(list_PCoA)[x],"PEMFIGO_PDAI.png", sep="_"), height = 5, width = 6, dpi=300)              

}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="PCoA/based_on_FFA/p_values_BRAY_Healthy_vs_Pemfigo.csv", quote = F,row.names = T)



############# PCoA on FFAs pre and post Treatment ############

dir.create("PCoA/FFAs_pre_post_Treatment")

colnames(FFA)[1:4]
# "Subject"   "Condition" "Acetic"    "Propionic"

colnames(FFA_post)[1:4]
# "Sample_pair_pre_treat" "Subject"   "Acetic"  "Propionic" 

list_PCoA<-list("SCFA"=SCFA_post,"MCFA"=MCFA_post,"LCFA"=LCFA_post)
list_PCoA_pre<-list("SCFA"=SCFA_raw,"MCFA"=MCFA_raw,"LCFA"=LCFA_raw)
#1               #2            #3

for(x in 1:length(list_PCoA)){
  rm(data, temp, temp_pre, dist, coord)
  data<-as.data.frame(list_PCoA[[x]])
  # creting a unique matrix
  temp_pre<-as.data.frame(list_PCoA_pre[[x]])
  temp_pre$Subject_ID<-temp_pre[["Subject"]]
  temp<-cbind.data.frame(data[,"Subject"],rep("Treated"),data[, ! colnames(data) %in% c("Sample_pair_pre_treat", "Subject")])
  colnames(temp)[1:2]<-c("Subject","Condition")
  temp$Subject_ID<-data[["Sample_pair_pre_treat"]]
  temp<-as.data.frame(rbind(temp_pre, temp))
  # identical(length(which(duplicated(temp$Subject_ID))),length(data$Subject)) # just a check
  
  row.names(temp)<-temp$Subject
  if(x==3) {temp<-temp[!is.na(temp$Octadecanoic),] } # there is a NA in LCFA
  dist<-vegdist(temp[, !colnames(temp) %in% c("Subject","Condition","Subject_ID")], method = "bray") # NB: sample on rows
  coord<- pco(dist)
  
  ggplot(data=as.data.frame(coord$vectors), aes(x=V1, y=V2, color=temp$Condition) ) + 
    geom_point(size=3) + theme_classic() +  stat_ellipse() +
    geom_text(aes(label=row.names(temp)),size=1.8, color="black") +
    geom_line(aes(group=temp$Subject_ID), size=0.2, color="grey") +
    labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
         y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
         color="Condition", caption = "lines connect the same subject pre and post treatment",
         title = paste("PCoA computed on",names(list_PCoA)[x],"abundances\n with Bray Curtis index"))
  
  ggsave(file=paste("PCoA/FFAs_pre_post_Treatment/PCoA_BRAY",names(list_PCoA)[x],"treatment.png", sep="_"), height = 6, width = 7, dpi=300)              
  
  # checking for significant dispersion
  rm(p_value_diver)
  p_value_diver<- pairwise.adonis(dist, factors=temp$Condition, p.adjust.m = "BH", sim.method="bray", perm = 9999)
  write.csv2(p_value_diver, file=paste("PCoA/FFAs_pre_post_Treatment/Beta_diversity_",names(list_PCoA)[x],"healthy_PEM_treat.csv"), quote=F, row.names=F)
}


################ PLS-DA on FFA and Condition ##############

dir.create("PLS_DA")

FFA_PLS<-as.data.frame(FFA[, ! colnames(FFA) %in% c("SCFA_raw_total","MCFA_raw_total","LCFA_raw_total") ])
row.names(FFA_PLS)<-FFA_PLS$Subject
head(FFA_PLS,n=2)
data.selected<-FFA_PLS[ , ! colnames(FFA_PLS) %in% c("Subject","Condition")]
metadata<-FFA_PLS[ , colnames(FFA_PLS) %in% c("Subject","Condition")]

X <- as.matrix(data.selected)
Y <- as.factor(metadata$Condition)

my.splsda <- splsda(X, Y, ncomp = 10)

gc()
perf.splsda <- perf(my.splsda, validation = "Mfold", # to assest the optimal number of comp
                    folds = 10, nrepeat = 5*(dim(X)[1]), cpus = 7,
                    progressBar = TRUE, auc = TRUE)

perf.splsda$choice.ncomp
# comp<-perf.splsda$choice.ncomp[3] # centroid distance
comp<-3 # it should be 1 but choosing 3 (less error than 2) in order to plot in 2D

png(filename = "PLS_DA/Error for each component PLS-DA", width = 2500, height = 1800, res=300)
plot(perf.splsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
title(main="Classification rate error \n for each component used to compute PLS-DA")
dev.off()

my.splsda <- splsda(X, Y, ncomp = comp)

n= 2 # choose the component to use depending on clustering efficiency in plot
png(filename = paste("PLS_DA/normal_PLS_DA_comp_1_and", n, sep="_"), width = 2500, height = 1500, res=300)
plotIndiv(my.splsda, comp = c(1,n),
          group = metadata$Condition, ind.names = TRUE, ellipse = TRUE, legend = TRUE, 
          title = paste0("PLSDA between Healthy and PEMFIGO \n on FFAs raw abundances \n computed with ", comp, " components, plotted on comp 1 and ", n,")"))
dev.off()

################ sPLS-DA for variable selection

# cutting out same sample for test (otherwise clear overfitting during the test)
set.seed(1)
train <- sample(1:nrow(X), 30) # select a huge number of sample
train
test <- setdiff(1:nrow(X), train)

X.test <- X[test,]
X <- X[train, ]
Y.test <- Y[test]
Y<- Y[train]
metadata2<-metadata[row.names(X),]

# 30 sample --> 10 folds divide samples equaly among k groups

possible.pool<- seq(3,round(dim(X)[2]),1) # for test.keepX
gc()
set.seed(1)
my.tuned.splsda <- tune.splsda(X, Y, ncomp = comp, validation = 'Mfold',
                               folds = 10, nrepeat = 4*(dim(X)[1]), cpus = 7, # use repeated cross-validation
                               dist = 'centroids.dist', test.keepX =  possible.pool, # testing from 10 to 1/5 MAX number of variables
                               measure = "BER", progressBar = T) # use balanced error rate of dist measure

my.tuned.splsda$choice.ncomp$ncomp
png(filename = paste("PLS_DA/Error_of_splsda_depending_of_computing_components",n, sep="_"), width = 1500, height = 2000, res=300)
plot(my.tuned.splsda, col = color.jet(comp))
dev.off()
# optimal.comp <- my.tuned.splsda$choice.ncomp$ncomp
optimal.comp<-2

optimal.keepX <- my.tuned.splsda$choice.keepX[1:optimal.comp]
optimal.keepX

final.splsda<-splsda(X,Y, ncomp=optimal.comp, keepX = optimal.keepX)
n=2 # choose which component plot depending on plot
ggplot(mapping = aes(x=final.splsda$variates$X[,"comp1"], y=final.splsda$variates$X[,"comp2"], color=final.splsda$Y)) +
  geom_point(size=3) + theme_classic(base_size = 14) + stat_ellipse() +
  geom_text(mapping = aes(label=row.names(final.splsda$X)), color="black", size=2) +
  scale_colour_manual(values = c("Healthy" = "deepskyblue", "PEMFIGO" = "coral"))+
  labs(title="sparse PLS-DA between Healthy and PEMFIGO subjects", 
       subtitle=paste0(" on raw FFAs abundances \n  (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"),
       caption = paste0("selected ",optimal.keepX[1]," FFAs on LC1 and ",optimal.keepX[2]," FFAs on LC2"),
       color="Condition", x=paste("X-Var1: ",round(final.splsda$prop_expl_var$X[1]*100,digits = 2),"% variance"),
       y=paste("LC2: ",round(final.splsda$prop_expl_var$X[2]*100,digits = 2),"% variance"))
ggsave(filename = paste("PLS_DA/sPLS-DA on comp 1 and",n,".png"), width = 6.5, height = 5, dpi=300) 
dev.off()

################ plotting loadings

loadings<-as.data.frame(final.splsda$loadings$X)
loadings$FFA<-row.names(loadings)

a<-loadings[loadings$comp1!=0,]
a$comp1<-as.numeric(a$comp1)
a<-a[order(abs(a$comp1)),]
a$FFA <- factor(a$FFA, levels = a$FFA) # otherwise it would be re-ordered in plot
ggplot(mapping=aes(x=a$FFA,y=a$comp1)) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, size = 11, hjust = 1, vjust = 1)) +
  labs(x="Selected FFAs on component 1", 
       y="loadings", title = "Loadings of selected FFAs for component 1 in sPLSDA", 
       caption="\n a dotted line is plot at 50% of the lower loading value") +
  # geom_hline(yintercept = max(a$comp1)/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  geom_hline(yintercept = min(a$comp1)/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = "PLS_DA/Loadings of choosen FFAs for sPLSDA comp 1.png", width = 14, height = 8, dpi=300)

b<-loadings[loadings[[n]]!=0,]
b[,n]<-as.numeric(b[,n])
b<-b[order( abs(b[,n])) ,]
b$FFA <- factor(b$FFA, levels = b$FFA)
ggplot(mapping=aes(x=b$FFA,y=b[,n])) + geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, size = 11, hjust = 1, vjust = 1)) +
  labs(x=paste("Selected FFAs on component",n), 
       y="loadings", title = paste("Loadings of selected FFAs for component",n,"in sPLSDA"), 
       caption="\n a dotted line is plot at 50% of the lower loading value") +
  #geom_hline(yintercept = max(b[,n])/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, size= 1) +
  geom_hline(yintercept = min(b[,n])/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = paste("PLS_DA/Loadings_of_choosen_FFAs_for_sPLSDA_comp",n,".png"), width = 14, height = 8, dpi=300)

############### testing the sPLS-DA model

predict.splsda <- predict(final.splsda, newdata =as.matrix(X.test), dist = "all")

# evaluating the prediction accuracy
predict.comp <- predict.splsda$class$centroids.dist[,n]
predict.comp
table(factor(predict.comp, levels = levels(Y)), as.factor(Y.test))
# now evaluating the prediction accuracy using ONLY the first component
predict.comp1 <- predict.splsda$class$centroids.dist[,1]
predict.comp1
table(factor(predict.comp1, levels = levels(Y)), as.factor(Y.test))

##### exporting impostations and values of sPLSDA

rm(con)
con<-file("Valori_parametri_risultati_prediz_di_sPLSDA.txt")
sink(con)
cat("Number of training samples and original number of FFAs", fill = TRUE)
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
cat("\n number of FFAs selected by tune.splsda function for each chosen component", fill=TRUE)
cat(optimal.keepX, fill=TRUE)
cat("\n Possible number of FFAs that could be selected by function tune.splsda", fill=TRUE)
cat(possible.pool, fill = TRUE)
cat("\n \n \n ### testing the prediction efficiency of sPLSDA through confusion matrix using ONLY the first component \n", fill=TRUE)
table(factor(predict.comp1, levels = levels(Y)), Y.test)
cat("\n", length(which(predict.comp1!="")),"on", length(predict.comp1), "samples predicted" )
cat("\n \n ### testing the prediction efficiency of sPLSDA through confusion matrix component 1 and",n,"\n", fill=TRUE)
table(factor(predict.comp, levels = levels(Y)), Y.test)
cat("\n", length(which(predict.comp!="")),"on", length(predict.comp), "samples predicted" )
cat("\n \n type of distance used: centroid distance")
sink()

rm(a,b,con, loadings, metadata, metadata2, predict.comp,predict.comp1,train,X,Y,X.test,Y.test,my.tuned.splsda,predict.splsda,comp,n,test,optimal.comp,optimal.keepX)
gc()

setwd("..")