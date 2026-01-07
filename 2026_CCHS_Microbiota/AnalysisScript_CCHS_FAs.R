##################### PREPARING THE ENVIRONMENT #################

{ library("phyloseq")
  library("ggplot2")
  library("ggpubr")
  library("ggh4x") #
  library("vegan")
  library("ecodist")
  library("mixOmics")
  library("Hmisc")
}

{
  dir.create("Results")
  dir.create("Results/FFAs")

}

options(scipen = 100) # disable scientific annotation




##################### IMPORTING THE METADATA #########################

Metadata <- as.data.frame(read.delim(file="Metadata_CCHS.txt", sep="\t", header = T))
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
# head(Metadata)

Metadata <- Metadata[ ! Metadata$Sample_name %in% c("C1_bis","C17G","H2_bis","H5","HP4_bis","HP9G" , "C5","C10","C11" ) , ]




##################### IMPORTING THE FFA TABLE #########################

Original_Stool<-as.data.frame(read.csv2(file="Fecal_FFAs.csv"))  # first column (row names) is samples
samples_FFA<-Original_Stool[,1]  # it will be lose with "apply"
Original_Stool<-Original_Stool[,-1]
# transforming in percentages within samples
Original_Stool<-apply(Original_Stool, MARGIN=2, FUN= as.numeric )
totals<-rowSums(Original_Stool)
Original_Stool<-as.data.frame( (Original_Stool/totals)*100 )
Original_Stool$Samples<-samples_FFA

# Metadata used specifically for the FFAs analysis
new_meta <- Metadata[ !Metadata$Sample_name %in% c("C1_bis","C17G","H2_bis","H5","HP4_bis","HP9G" , "C5","C10","C11" ) , ]
row.names(new_meta) <- gsub("_bis","", new_meta$Sample_name)

Original_Stool <- Original_Stool[ Original_Stool$Samples%in%row.names(new_meta), ]
if(length(Original_Stool$Samples)!=32){  # NB: two samples were not enough to analyse also their FFAs --> 34 -2 = 32 FFA expected
  stop("Something went wrong during the import of the FFAs table...")
} else { print("OK ;-) ") }




########################  BOXPLOTS AND MANN WHITNEY FFA and ACAs #########################

Stool <- Original_Stool

# reshaping for ggplot2
Stool<-reshape::melt(Stool, id="Samples")
Stool$variable<-gsub("X2.","2-",Stool$variable)
Stool$variable<-gsub("Iso","iso",Stool$variable)
# adding short medium long infos
Stool$SML<-rep("temp")
Stool$SML[Stool$variable %in% c("Acetic","Propionic","isoButyric","Butyric","isoValeric","Valeric","2-ethylButyric")] <- "SCFA"
Stool$SML[Stool$variable %in% c("Hexanoic","isoHexanoic","Heptanoic","Octanoic","Nonanoic","Decanoic","Dodecanoic")] <- "MCFA"
Stool$SML[Stool$variable %in% c("Tetradecanoic","Hexadecanoic","Octadecanoic")] <- "LCFA"
Stool$SML[Stool$variable %in% c("Benzoic","Phenylacetyc","Phenylpropionic")] <- "ACA"
# final adjustments
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Group <- ifelse( grepl("^C",Stool$Samples) , "CCHS", "Healthy" )
Stool$Group<-factor(Stool$Group, levels = c("Healthy","CCHS"))
Stool$SML<-factor(Stool$SML, levels = c("SCFA","MCFA","LCFA","ACA"))
Stool$Age_groups<- new_meta[Stool$Samples,"Age_3groups"]
Stool$Age<-new_meta[Stool$Samples,"Age"]
# plotting
ggplot(data=Stool, aes(x=variable, y=value, fill=Group)) +
  scale_fill_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2"  )) +
  scale_color_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2" )) +
  facet_grid(~ SML, scales = "free_x", space = "free_x") +
  geom_boxplot(width=0.65, linewidth= 0.2, alpha= 0.1, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.65, jitter.width = 0.22),
             aes(color=Group), size= 0.2, alpha= 0.85) +  
  theme_classic(base_size = 11) + 
  theme(strip.text.x=element_text(size=11.5,colour="black")) + 
  theme(legend.margin=margin(-20, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1, size=9), 
        axis.text.y = element_text(size=7.5),
        plot.margin = margin(3,1,2,1),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=13),
        panel.grid.major.y = element_line(linewidth=0.10, color="darkgray"),
        panel.grid.major.x = element_line(linewidth=0.04, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1), color="none" ) +
  labs(y="Percentual quantity", fill="", x="") + 
  scale_x_discrete(expand=c(-1, 1)) + 
  scale_y_sqrt(breaks=c( 1, seq(5,100,5)) , expand=c(0,0.01) )
ggsave(file="Results/FFAs/Stool_EVERY_FFA_boxplot.png", width=5.5, height = 4.5, dpi=300)
ggsave(file="Results/FFAs/Stool_EVERY_FFA_boxplot_LARGE.png", width=7, height = 4, dpi=450)




#################### METABOLITES STATISTICS (Wilcoxon) ###############################

dir.create("Results/FFAs/Classic_Wilcoxon")

table_pvalue<-NULL

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  test <- wilcox.test( target$value ~ target$Group, paired=F  )
  p_value_value <- test$p.value
  p_value_value <- round(p_value_value, 3)
  p_value_value[p_value_value==0]<-"0.001"
  
  table_pvalue<- rbind.data.frame(table_pvalue, cbind(i,p_value_value))
}
table_pvalue$p_adj<- round( p.adjust(table_pvalue$p_value_value, method = "BH") , 2 )
table_pvalue$signif<- ifelse(table_pvalue$p_adj<0.05,"*","")
colnames(table_pvalue)<-c("variable","MannW_p_value","p_adj","signif")

i<-"Benzoic"
# trying after removing the healthy outlier
target<-Stool[Stool$variable==i , ]
target<-target[target$Samples!="HP10", ]

test <- wilcox.test( target$value ~ target$Group, paired=F  )
p_value_value <- test$p.value
p_value_value <- round(p_value_value, 3)
p_value_value[p_value_value==0]<-"0.001"

new_row<- cbind( variable = paste(i,"without outliers"), 
                 MannW_p_value = p_value_value ,
                 p_adj= table_pvalue[table_pvalue$variable=="Benzoic", "p_adj" ],   # using the same p adj as the other benzoic row, to no inflate the adjustment penalities
                 signif = table_pvalue[table_pvalue$variable=="Benzoic", "signif" ]
)
table_pvalue<- rbind.data.frame(table_pvalue, new_row )

write.csv2(table_pvalue , file="Results/FFAs/Classic_Wilcoxon/MannWithney_FFAs.csv" , quote = F, row.names = F)




##### PLOTTING SIGNIFICATIVES

for(i in unique(table_pvalue$variable)){
  
  target<-Stool[Stool$variable==i , ]
  target_ages <- target[, "Age"]
  
  if( i == "Benzoic without outliers" ) {
    target <- Stool[Stool$variable=="Benzoic" & Stool$Samples!="HP10", ]
    target_ages <-  target[, "Age"] # a "NULL" is here overwritten 
  }
  
  target_pvalue<- table_pvalue[table_pvalue$variable==i, ]
  # plot
  ggboxplot(data=target, x="Group", y="value", fill="Group", 
            width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
    scale_fill_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2"  )) +
    scale_color_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2" )) +
    geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
               aes(color=Group), size= 1, alpha= 0.5) +  
    theme_classic2(base_size = 7) + 
    #geom_text(aes(label=target_ages), color="black", size=2.5, show.legend = FALSE) +
    theme(strip.text.x=element_text(size=12,colour="black"),
          axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
          axis.text.y = element_text(size=5),
          plot.margin = margin(5,2,0,2),
          plot.title= element_text(size=8, hjust = 0.5, vjust=1.8) ,
          plot.subtitle = element_text(size=6, vjust=1) ,
          legend.key.size=unit(0.8,"cm"), 
          legend.text=element_text(size=7),
          panel.grid.major.y = element_line(linewidth =0.12, color="gray"),
          panel.grid.minor.y = element_line(linewidth =0.03, color="gray"),
          panel.grid.major.x = element_blank()
    ) +
    guides( color="none", fill="none" ) +
    labs(y="Percentual quantity", x="",
         title = i,
         subtitle = paste0("p-value: ", target_pvalue$MannW_p_value, " (adjusted :",target_pvalue$p_adj , ")" )
    ) + 
    scale_x_discrete(expand=c(0.2, 0))
  
  ggsave(file=paste0("Results/FFAs/Classic_Wilcoxon/",i,"_with_Mann_Whitney.png"), width = 2, height = 2.3, dpi=300)
}




######################### METABOLITES STATISTICS (LM WITH AGE CONTINUOUS) ###############################

dir.create("Results/FFAs/LM_with_age_CONTINUOUS")

table_pvalue<-NULL

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  model_tested <- lm( target$value ~ target$Group + target$Age)  # for the plot  (DO NOT USE RANKS HERE)
  ranked_model_tested <- lm( rank(target$value) ~ target$Group + target$Age)  # for the statistics
  coef_model <- coef( model_tested )  # base intercept, slope and covariate-intercept-diff
  test <- summary( ranked_model_tested )
  p_value_value <- test$coefficients[2,4]  # Condition CCHS' p-value , after age correction
  p_value_value <- round(p_value_value, 3)
  p_value_value[p_value_value==0]<-"0.001"
  age_pvalue <- round(test$coefficients[3,4] , 3)
  
  table_pvalue<- rbind.data.frame(table_pvalue, 
                                  cbind(i,coef_model[1],coef_model[2],coef_model[3],p_value_value,age_pvalue)
                                  )
}
table_pvalue$p_adj<- round( p.adjust(table_pvalue$p_value_value, method = "BH") , 2 )
table_pvalue$signif<- ifelse(table_pvalue$p_value_value<=0.05,"*","")
table_pvalue$signif_adh<- ifelse(table_pvalue$p_adj<=0.05,"*","")
colnames(table_pvalue)<-c("variable",
                          "intercept","Age_intercept_diff","Cond_slope",
                          "p_value_COND_after_AGE","p_value_AGE","p_adj_COND","signif","signif_adj"
                          )

write.csv2(table_pvalue , file="Results/FFAs/LM_with_age_CONTINUOUS/LINEAR_MODEL_INCLUDING_age_as_covar_CONTINUOUS.csv" , quote = F, row.names = F)



##### PLOTTING SIGNIFICATIVES

plot_these_regardless<- c("Butyric","Valeric","Phenylpropionic","Octanoic","Nonanoic")
for(i in c( unique(table_pvalue$variable[table_pvalue$signif=="*"] ) , plot_these_regardless ) ){   # yes, I want to plot "Butyric" acid regarless of the p-value

  target<-Stool[Stool$variable==i , ]
  target_ages <- target[, "Age"]

  if( i == "Benzoic without outliers" ) {
    target <- Stool[Stool$variable=="Benzoic" & Stool$Samples!="HP10", ]
    target_ages <-  target[, "Age"] # a "NULL" is here overwritten
  }

  target_pvalue<- table_pvalue[table_pvalue$variable==i, ]
  # plot
  ggplot(data=target, aes(y=value, x=Age, fill=Group) ) +
    geom_abline(intercept= as.numeric(target_pvalue[["intercept"]]) , 
                slope= as.numeric(target_pvalue[["Cond_slope"]]) , 
                colour="deepskyblue2" , linetype= 2 
                ) +
    geom_abline(intercept= as.numeric(target_pvalue[["intercept"]])+ as.numeric(target_pvalue[["Age_intercept_diff"]]) ,
                slope= as.numeric(target_pvalue[["Cond_slope"]]),
                colour="coral", linetype= 2
                ) +
    scale_fill_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2"  )) +
    scale_color_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2" )) +
    geom_point(aes(color=Group), size= 1.3, alpha= 0.5) +
    geom_point(aes(color=Group), size= 0.45, alpha= 0.9) +
    theme_classic2(base_size = 7) +
    theme(strip.text.x=element_text(size=12,colour="black"),
          axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=5),
          axis.text.y = element_text(size=5),
          plot.margin = margin(5,2,0,2),
          plot.title= element_text(size=8, hjust = 0.5, vjust=1.8) ,
          plot.subtitle = element_text(size=6, vjust=1) ,
          legend.key.size=unit(0.8,"cm"),
          legend.text=element_text(size=7),
          panel.grid.major.y = element_line(linewidth =0.12, color="gray"),
          panel.grid.minor.y = element_line(linewidth =0.03, color="gray"),
          panel.grid.major.x = element_blank()
    ) +
    guides( color="none", fill="none" ) +
    labs(y="Percentual quantity", x="Age",
         title = i,
         subtitle = paste0("p-value in corrected model: ", target_pvalue$p_value_COND_after_AGE, " (adjusted :",target_pvalue$p_adj_COND , ")" )
    ) +
    scale_x_continuous(breaks = seq(10, max(target$Age)+1, 5))

  ggsave(file=paste0("Results/FFAs/LM_with_age_CONTINUOUS/",i,"_with_AgeContinuous.png"), width = 2.3, height = 2, dpi=300)
}




######################### METABOLITES STATISTICS (LM WITH AGE GROUPS) ###############################

dir.create("Results/FFAs/LM_with_age_Groups")

table_pvalue<-NULL

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  model_tested <- lm( target$value ~ target$Group + target$Age)  # for the plot  (DO NOT USE RANKS HERE, MOREOVER AGE HAS TO BE CONTINUOUS)
  ranked_model_tested <- lm( rank(target$value) ~ target$Group + target$Age_groups)  # for the statistics
  coef_model <- coef( model_tested )  # base intercept, slope and covariate-intercept-diff
  test <- summary( ranked_model_tested )
  p_value_value <- test$coefficients[2,4]  # Condition CCHS' p-value , after age correction
  p_value_value <- round(p_value_value, 3)
  p_value_value[p_value_value==0]<-"0.001"
  age_pvalue <- round(test$coefficients[3,4] , 3)
  
  table_pvalue<- rbind.data.frame(table_pvalue, 
                                  cbind(i,coef_model[1],coef_model[2],coef_model[3],p_value_value,age_pvalue)
  )
}
table_pvalue$p_adj<- round( p.adjust(table_pvalue$p_value_value, method = "BH") , 2 )
table_pvalue$signif<- ifelse(table_pvalue$p_value_value<=0.05,"*","")
table_pvalue$signif_adh<- ifelse(table_pvalue$p_adj<=0.05,"*","")
colnames(table_pvalue)<-c("variable",
                          "intercept","Age_intercept_diff","Cond_slope",
                          "p_value_COND_after_AGE","p_value_AGE","p_adj_COND","signif","signif_adj"
)

write.csv2(table_pvalue , file="Results/FFAs/LM_with_age_Groups/LINEAR_MODEL_INCLUDING_age_as_cofact_GROUPS.csv" , quote = F, row.names = F)



##### PLOTTING SIGNIFICATIVES

plot_these_regardless<- c("Butyric","Valeric","Phenylpropionic","Octanoic","Nonanoic")
for(i in c( unique(table_pvalue$variable[table_pvalue$signif=="*"] ) , plot_these_regardless ) ){   # yes, I want to plot "Butyric" acid regarless of the p-value
  
  target<-Stool[Stool$variable==i , ]
  target_ages <- target[, "Age"]
  
  if( i == "Benzoic without outliers" ) {
    target <- Stool[Stool$variable=="Benzoic" & Stool$Samples!="HP10", ]
    target_ages <-  target[, "Age"] # a "NULL" is here overwritten
  }
  
  target_pvalue<- table_pvalue[table_pvalue$variable==i, ]
  # plot
  ggplot(data=target, aes(y=value, x=Age, fill=Group) ) +
    geom_abline(intercept= as.numeric(target_pvalue[["intercept"]]) , 
                slope= as.numeric(target_pvalue[["Cond_slope"]]) , 
                colour="deepskyblue2" , linetype= 2 
    ) +
    geom_abline(intercept= as.numeric(target_pvalue[["intercept"]])+ as.numeric(target_pvalue[["Age_intercept_diff"]]) ,
                slope= as.numeric(target_pvalue[["Cond_slope"]]),
                colour="coral", linetype= 2
    ) +
    scale_fill_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2"  )) +
    scale_color_manual(values=c( "CCHS"="coral","Healthy"="deepskyblue2" )) +
    geom_point(aes(color=Group), size= 1.3, alpha= 0.5) +
    geom_point(aes(color=Group), size= 0.45, alpha= 0.9) +
    theme_classic2(base_size = 7) +
    theme(strip.text.x=element_text(size=12,colour="black"),
          axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=5),
          axis.text.y = element_text(size=5),
          plot.margin = margin(5,2,0,2),
          plot.title= element_text(size=8, hjust = 0.5, vjust=1.8) ,
          plot.subtitle = element_text(size=6, vjust=1) ,
          legend.key.size=unit(0.8,"cm"),
          legend.text=element_text(size=7),
          panel.grid.major.y = element_line(linewidth =0.12, color="gray"),
          panel.grid.minor.y = element_line(linewidth =0.03, color="gray"),
          panel.grid.major.x = element_blank()
    ) +
    guides( color="none", fill="none" ) +
    labs(y="Percentual quantity", x="Age",
         title = i,
         subtitle = paste0("p-value in corrected model: ", target_pvalue$p_value_COND_after_AGE, " (adjusted :",target_pvalue$p_adj_COND , ")" )
    ) +
    scale_x_continuous(breaks = seq(10, max(target$Age)+1, 5))
  
  ggsave(file=paste0("Results/FFAs/LM_with_age_Groups/",i,"_with_AgeGROUPS.png"), width = 2.3, height = 2, dpi=300)
}




################## TEST PCoA METABOL ###################

dir.create("Results/FFAs/PLS_DA_on_FFAs")

Metab_PCoA <-as.data.frame(Original_Stool)   # NB: already in proportions, see the transformation in import sections
row.names(Metab_PCoA)<-Original_Stool$Samples



### Bray Curtis
dist<-vegdist(Metab_PCoA[, !colnames(Metab_PCoA) %in% c("Samples")], method = "bray") # NB: sample on rows
coord<- pco(dist)
which_color <- ifelse ( grepl("C",row.names(Metab_PCoA)) , "coral", "deepskyblue2" )
which_shape <- ifelse ( grepl("C",row.names(Metab_PCoA)) , 1 , 3 )

new_names <- new_meta[ row.names(Metab_PCoA) , "Sample_Family" ]
age_here <- new_meta[ row.names(Metab_PCoA) , "Age_groups" ]

ggplot(data=as.data.frame(coord$vectors), aes(x=X1, y=X2,
                                              color= which_color,
                                              shape= age_here) ) + 
  geom_point(size=5 , alpha= 0.5) + 
  geom_point(size=3) + 
  geom_point(size=3) + 
  theme_classic() +
  stat_ellipse() +
  geom_text(aes(label= new_names ),size=2, color="black") +
  labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
     y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
     color="Condition", shape="Age",
     title = "PCoA computed on metab prop abundances \nwith Bray Curtis index") +
  guides(color="none")
ggsave(file="Results/FFAs/PCoA_metabolites.png", height = 4, width = 5.2, dpi=300)       



### Again, on Euclidean distance
dist<-vegdist(Metab_PCoA[, !colnames(Metab_PCoA) %in% c("Samples")], method = "euclidean") # NB: sample on rows
coord<- pco(dist)
which_color <- ifelse ( grepl("C",row.names(Metab_PCoA)) , "coral", "deepskyblue2" )
which_shape <- ifelse ( grepl("C",row.names(Metab_PCoA)) , 1 , 3 )

new_names <- new_meta[ row.names(Metab_PCoA) , "Sample_Family" ]
new_names <- gsub("_bis", "", new_names)
age_here <- new_meta[ row.names(Metab_PCoA) , "Age_groups" ]

ggplot(data=as.data.frame(coord$vectors), aes(x=X1, y=X2,
                                              color= which_color,
                                              shape= age_here) ) + 
  geom_point(size=5 , alpha= 0.5) + 
  geom_point(size=3) + 
  geom_point(size=3) + 
  theme_classic(base_size = 8) +
  stat_ellipse( linewidth= 0.11 ) +
  geom_text(aes(label= new_names ),size=2, color="black") +
  labs(x=paste("PC1:",round(coord[["values"]][1]/sum(coord[["values"]])*100,digits=2),"%"),
       y=paste("PC2:",round(coord[["values"]][2]/sum(coord[["values"]])*100,digits=2),"%"),
       color="Condition", shape="Age",
       title = "PCoA computed on metab prop abundances \nwith Euclidean distance index") +
  guides(color="none")
ggsave(file="Results/FFAs/PCoA_metabolites_Euclidean.png", height = 3.2, width = 3.5, dpi=300)       




####################### sPLS-DA on FFAs ###########################

FFA_PLS<-as.data.frame(Original_Stool)
row.names(FFA_PLS)<-FFA_PLS$Samples
head(FFA_PLS,n=2)
data.selected<-FFA_PLS[ , ! colnames(FFA_PLS) %in% "Samples"]
metadata_FFA<-FFA_PLS[ , colnames(FFA_PLS) %in% "Samples" , drop = F]
metadata_FFA$Condition <- ifelse( grepl("^H",metadata_FFA$Samples), "Healthy", "CCHS" )

X <- as.matrix(data.selected)
Y <- as.factor(metadata_FFA$Condition)


### classic PLS DA
my.splsda <- splsda(X, Y, ncomp = 5)

gc()
set.seed(1994)

perf.splsda <- perf(my.splsda, validation = "Mfold", # to assest the optimal number of comp
                    folds = 10, nrepeat = 100 , cpus = 7,
                    progressBar = TRUE, auc = TRUE)

perf.splsda$choice.ncomp
# comp<-perf.splsda$choice.ncomp[3] # centroid distance   (4, but the difference is very small, then let's set the value to 2)
comp<-2 

png(filename = "Results/FFAs/PLS_DA_on_FFAs/Error for each component PLS-DA.png", width = 2500, height = 1800, res=300)
plot(perf.splsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
title(main="Classification rate error \n for each component used to compute PLS-DA")
dev.off()

my.splsda <- splsda(X, Y, ncomp = comp)

n= 2 # choose the component to use depending on clustering efficiency in plot
png(filename = paste("Results/FFAs/PLS_DA_on_FFAs/normal_PLS_DA_comp_1_and", n, ".png", sep="_"), width = 2500, height = 1500, res=300)
plotIndiv(my.splsda, comp = c(1,n),
          group = metadata_FFA$Condition, ind.names = TRUE, ellipse = TRUE, legend = TRUE, 
          title = paste0("PLSDA between HC and PV \n on FFAs raw abundances \n computed with ", comp, " components, plotted on comp 1 and ", n,")"))

dev.off()



################## sPLS-DA


# cutting out same sample for testing the model
set.seed(1994)
train <- sample(1:nrow(X), 27) # selecting a huge number of  training samples
train
test <- setdiff(1:nrow(X), train)

X.test <- X[test,]
X <- X[train, ]
Y.test <- Y[test]
Y<- Y[train]
metadata2<-metadata_FFA[row.names(X),]

# about 30 sample (27) --> 10 folds divide samples equally among k groups... yet they would be only 3 sample in each subgroups to test the data on!
# if k-fold is set to 5 (which is still a conventionally often used "default" value), there would be about 5 subgroups with 5 samples each --> This should be better...

#possible.pool<- seq(1,round(dim(X)[2]),1) # for test.keepX
possible.pool<- c(1,seq(2,8,2)) # up to max 8
gc()
set.seed(1994)
my.tuned.splsda <- tune.splsda(X, Y, ncomp = comp, validation = 'Mfold',
                               # folds = 10,
                               folds = 5, # see above
                               nrepeat = 100, cpus = 8, # use repeated cross-validation
                               dist = 'centroids.dist', test.keepX =  possible.pool, # testing from 10 to 1/5 MAX number of variables
                               measure = "BER",
                               progressBar = T) # use balanced error rate of dist measure

my.tuned.splsda$choice.ncomp$ncomp  # 1 is enough ... but I need the second component to plot a PCoA 
png(filename = paste("Results/FFAs/PLS_DA_on_FFAs/Error_of_splsda_depending_of_computing_components",n,".png" , sep="_"), width = 1500, height = 2000, res=300)
plot(my.tuned.splsda, col = color.jet(comp))
dev.off()
# optimal.comp <- my.tuned.splsda$choice.ncomp$ncomp
optimal.comp<-2   # read above

optimal.keepX <- my.tuned.splsda$choice.keepX[1:optimal.comp]
optimal.keepX  # 8 and 8

final.splsda<-splsda(X,Y, ncomp=optimal.comp, keepX = optimal.keepX)
n=2 # choose which component plot depending on plot
ggplot(mapping = aes(x=final.splsda$variates$X[,"comp1"], y=final.splsda$variates$X[,"comp2"], color=final.splsda$Y)) +
  geom_point(size=4, alpha = 0.6) +
  geom_point(size=3) + theme_classic(base_size = 13) + stat_ellipse() +
  geom_text(mapping = aes(label=row.names(final.splsda$X)), color="black", size=2) +
  scale_colour_manual(values = c("Healthy" = "deepskyblue", "CCHS" = "coral"))+
  labs(title="sparse PLS-DA on 28 training samples", 
       subtitle=paste0(" on proportional scaled FFAs values \n  (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"),
       caption = paste0("selected ",optimal.keepX[1]," FFAs on LC1 and ",optimal.keepX[2]," FFAs on LC2"),
       color="Condition", x=paste("LC1: ",round(final.splsda$prop_expl_var$X[1]*100,digits = 2),"% variance"),
       y=paste("LC2: ",round(final.splsda$prop_expl_var$X[2]*100,digits = 2),"% variance"))
ggsave(filename = paste("Results/FFAs/PLS_DA_on_FFAs/sPLS-DA on comp 1 and",n,".png"), width = 6.5, height = 5, dpi=300) 
dev.off()
# without names
ggplot(mapping = aes(x=final.splsda$variates$X[,"comp1"], y=final.splsda$variates$X[,"comp2"], color=final.splsda$Y)) +
  geom_point(size=2.38, alpha= 1) +
  geom_point(size=3.8, alpha= 0.65) +
  theme_classic(base_size = 11) + stat_ellipse(linewidth=0.18) +
  #theme( legend.text = element_text(size=7.5) ) +
  scale_colour_manual(values = c("Healthy" = "deepskyblue", "CCHS" = "coral"))+
  labs(title="sparse PLS-DA", 
       subtitle=paste0(" on proportional scaled FFAs values \n  (computed with ", optimal.comp, " components, plotted on comp 1 and ", n,")"),
       caption = paste0("selected ",optimal.keepX[1]," FFAs on LC1 and ",optimal.keepX[2]," FFAs on LC2"),
       color="Condition",
       x=paste("LC1: ",round(final.splsda$prop_expl_var$X[1]*100,digits = 2),"% variance"),
       y=paste("LC2: ",round(final.splsda$prop_expl_var$X[2]*100,digits = 2),"% variance"))
ggsave(filename = paste("Results/FFAs/PLS_DA_on_FFAs/sPLS-DA on comp 1 and",n,"_points_only.png"), width = 4.8, height = 4, dpi=300) 
dev.off()



################ plotting loadings

loadings<-as.data.frame(final.splsda$loadings$X)
loadings$FFA<-row.names(loadings)

a<-loadings[loadings$comp1!=0,]
a$comp1<-as.numeric(a$comp1)
a<-a[order(abs(a$comp1)),]
a$FFA <- gsub("X2","2-", a$FFA)
a$FFA <- factor(a$FFA, levels = a$FFA) # otherwise it would be re-ordered in plot
ggplot(mapping=aes(x=a$FFA,y=a$comp1)) +
  geom_bar(stat = "identity", width = 0.6) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, size = 11, hjust = 1, vjust = 1),
        plot.title = element_text(size= 8 ),
        plot.caption = element_text(size= 5.5 )
        ) +
  labs(x="Selected FFAs on component 1", 
       y="loadings", title = "Loadings of selected FFAs for component 1 in sPLSDA", 
       caption="\n a dotted line is plot at 50% of the highest loading value") +
  geom_hline(yintercept = max(a$comp1)/2, colour="red", linetype="longdash") +
  geom_hline(yintercept = 0, linewidth= 1) +
  # geom_hline(yintercept = min(a$comp1)/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,0.5,1), "cm"))
ggsave(filename = "Results/FFAs/PLS_DA_on_FFAs/Loadings of choosen FFAs for sPLSDA comp 1.png", width = 6.25, height = 4, dpi=300)

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
  geom_hline(yintercept = 0, linewidth= 1) +
  #geom_hline(yintercept = min(b[,n])/2, colour="red", linetype="longdash") +
  theme(plot.margin = unit(c(1,0.5,1,1.8), "cm"))
ggsave(filename = paste("Results/FFAs/PLS_DA_on_FFAs/Loadings_of_choosen_FFAs_for_sPLSDA_comp",n,".png"), width = 8, height = 5, dpi=300)




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


##### exporting settings and values of sPLSDA

suppressWarnings(rm(con))
con<-file("Results/FFAs/PLS_DA_on_FFAs/Valori_parametri_risultati_prediz_di_sPLSDA.txt")
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

rm(a,b,con, loadings, metadata_FFA, metadata2, predict.comp,predict.comp1,train,X,Y,X.test,Y.test,my.tuned.splsda,predict.splsda,comp,n,test,optimal.comp,optimal.keepX)
gc()




################# CORRELATIONS BETWEEN BACTERIA AND FFA  (EVERY SAMPLE) ######################

dir.create("Results/FFAs/FFAs_correlations_EVERYsample")

load("data_prepared_after_filters_CCHS.RData")  # this is saved from the 16S scripts


suppressWarnings(rm(data.temp))
data.temp<-subset_samples(data.genus.prop, Sample_name %in% new_meta$Sample_name )
tax_table(data.temp) <- as.matrix( Taxa.genus.update )
#prune.dat_target <- subset_taxa(data.temp, Genus %in% c("Lactobacillus","Prevotella","Streptococcus","Eubacterium ventriosum group"))
tax_table(data.temp)[,"Genus"]<- gsub("[", "", tax_table(data.temp)[,"Genus"], fixed=T)
tax_table(data.temp)[,"Genus"]<- gsub("]", "", tax_table(data.temp)[,"Genus"], fixed=T)
tax_table(data.temp)[,"Genus"]<- gsub("_group", "", tax_table(data.temp)[,"Genus"], fixed=T)
tax_table(data.temp)[,"Genus"]<- gsub("_Incertae_Sedis", "", tax_table(data.temp)[,"Genus"])
prune.dat_target <- filter_taxa( data.temp ,function(x) mean(x)> 0.1, prune = T)
sample_names(prune.dat_target) <- sample_data(data.temp)$Sample_name
target <-as.data.frame(otu_table(prune.dat_target))
target_tax<-as.data.frame(tax_table(prune.dat_target))
identical(row.names(target_tax),row.names(target))
row.names(target)<-target_tax$Genus
colnames(target) <- gsub("_bis","", colnames(target) )
# Removing the sample without FFA (not enough biomass to analyse both)
target <- target[ , colnames(target)%in%Original_Stool$Samples  ]
target<-as.data.frame(t(target))
rm(target_tax)

# which_ones <- c("Acetic","Propionic","isoButyric","Butyric","isoValeric","Valeric","2-ethylButyric")   # SCFA
# which_ones <- c(which_ones, "Benzoic", "Phenylacetyc" ,"Phenylpropionic")  # adding also the aromatic carboxylic acids
# which_ones <- c(which_ones, c("Hexanoic","isoHexanoic","Heptanoic","Octanoic","Nonanoic","Decanoic","Dodecanoic") )   # few MCFA are significant, then they are included ...
# FFA <- Original_Stool[ Original_Stool$Samples %in% row.names(new_meta) , colnames(Original_Stool) %in% which_ones ]
# row.names(FFA)<-Original_Stool$Samples

FFA <- Original_Stool[ Original_Stool$Samples %in% row.names(new_meta) , ]
row.names(FFA)<-FFA$Samples
FFA$Samples<-NULL
  
# same order with the phyloseq data
FFA<-FFA[row.names(target), ]
length(FFA[[1]]) # 32 --> ok

target<-target[row.names(FFA), ] # to re-subset according to missing subjects among FFA values


# Computing the correlations...
x<-cbind.data.frame(target,FFA)
{x<-as.matrix(x[ ! colnames(x) %in% "Condition" ])
  r<-rcorr(x, type = "spearman")
  correlation_corr<-as.data.frame(as.table(r$r))
  correlation_pvalue<-as.data.frame(as.table(r$P))
  if(!  identical(correlation_corr[,1:2],correlation_pvalue[,1:2]) ){
    stop("Check here ...")
  }
  correlation<-cbind(correlation_corr,correlation_pvalue[,3])
  colnames(correlation)<-c("Genera","FFA","Corr","pvalue")
}
corr<-correlation[ correlation$Genera %in% colnames(target), ]
corr<-corr[ corr$FFA %in% colnames(FFA), ]
# colnames(corr)<-c("Genera","FFA","Corr","pvalue")
corr$padj<-round( p.adjust(corr$pvalue, method = "BH") ,2 )
corr$pvalue<-round( corr$pvalue ,2 )
corr$Sign<-corr$padj
corr$Sign[corr$Sign <= 0.05]<-"*"
corr$Sign[corr$Sign > 0.05]<-""

#corr2 <- corr[corr$Genera %in% corr$Genera[corr$Sign=="*"], ]
# corr2$Genera <- gsub("uncultured_ ","uncult_", corr2$Genera)
# corr2$Genera <- gsub("coprostanoligenes","coprostan.", corr2$Genera)
# ggplot(corr2, aes(y = Genera, x = FFA, fill = Corr)) +
#   geom_tile(color = "white", lwd = 0.5, linetype = 1) +
#   scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
#   theme_bw(base_size=12) +
#   theme(axis.text.x=element_text(angle = -40, hjust = 0, size= 10.5)) +
#   guides(fill= guide_colourbar(title ="rho")) +
#   geom_text(aes(label= Sign), color= "white", size =8) +
#   labs( y= "Genera signif. correlated", x= "FFA") +
#   theme(legend.text = element_text(size=9), legend.key.height= unit(1.8, 'cm'),
#         legend.key.width= unit(0.9, 'cm'), legend.title = element_text(size=10))
# ggsave(file="Results/FFAs/FFAs_correlations_EVERYsample/OnlySignificantlyCorrelatedGenera_Heatmap.png",
#        dpi=300, width = 5.8, height = 5)

write.csv2(corr[order(corr$padj, decreasing = F), ] ,
           file="Results/FFAs/FFAs_correlations_EVERYsample/SpearmanCorr_Genera_vs_FFA.csv", quote = F, na="", row.names = F)


# better version
corr2<-corr
corr2$Genera<-factor(corr2$Genera, levels = sort(unique(as.character(corr2$Genera)))) # following alphabetical order
original_levels<-gsub(".","-",levels(corr2$FFA), fixed=T)
corr2$FFA<-gsub(".","-",corr2$FFA, fixed=T)
corr2$FFA<-factor(corr2$FFA, levels = original_levels)
levels(corr2$FFA) <- gsub("^X","", levels(corr2$FFA))
# corr2$SignUNadj <- ifelse(corr2$pvalue<=0.05, "*", "")
ggplot(corr2, aes(x = Genera, y = FFA, fill = Corr)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1)) +
  theme_bw(base_size=13) +
  guides(fill= guide_colourbar(title ="rho")) +
  geom_text(aes(label= Sign), color= "gray35", size =7.25, vjust=0.81) +
  geom_text(aes(label= Sign), color= "gray90", size =4, vjust=0.81) +
  # geom_text(aes(label= SignUNadj), color= "gray35", size =7.25, vjust=0.81) +
  # geom_text(aes(label= SignUNadj), color= "gray90", size =4, vjust=0.81) +
  theme(axis.text.x=element_text(angle = 28, hjust = 1, vjust=1.08, size= 7),
        axis.text.y=element_text(size= 6.3, hjust=1.03),
        #axis.text.y=element_text(size= 8.5),
        legend.text = element_text(size=9.5, hjust=0.05),
        legend.key.height= unit(2.64, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.title = element_text(size=14.2, vjust=0.25),
        legend.margin = margin(0,0.5,0,-5),
        plot.margin = margin(1,0.5,-16.8,0.5)
        ) +
  scale_x_discrete(expand = c(0,0.01)) +
  labs(y= "", x= "Genera with average % abund > 0.1" #,
       #caption= "\n adjusted p-value (BH) lower than 0.05 are displayed through * sign"
       )  +
coord_flip()   # NB: flips everything!
#ggsave(filename = "Results/FFAs/FFAs_correlations_EVERYsample/Heatmap_corr.png", width = 8.5, height = 4.5, dpi=300 )
ggsave(filename = "Results/FFAs/FFAs_correlations_EVERYsample/Heatmap_corr.png", height = 5.95, width = 6, dpi=300 )



### Single plots
to_plot_list<- corr[corr$Sign=="*", ]

pdf(file="Results/FFAs/FFAs_correlations_EVERYsample/Check_on_signif_linear_corr.pdf", width = 6.2, height = 6)
par( mfrow= c(3,3) )

for (i in 1:length(to_plot_list$Sign)){
  Signif_G <- as.character(to_plot_list[i,"Genera"])
  Signif_S <- as.character(to_plot_list[i,"FFA"])
  FFA_values <- FFA[ , Signif_S , drop=F]
  Gen_abund <- otu_table(subset_taxa ( prune.dat_target , Genus == Signif_G ) )
  colnames(Gen_abund) <- gsub("_bis","", colnames(Gen_abund))
  Gen_abund <- Gen_abund[ ,row.names(FFA_values), drop=F]  # same order
  FFA_values <- round( as.numeric(FFA_values[,1]) , 1)
  Gen_abund <- round( as.numeric(Gen_abund) , 2)
  
  plot(Gen_abund, FFA_values, 
       xlab= paste0 (Signif_G , "%") ,
       ylab= paste0 (Signif_S ,"%") 
  )
  abline( lm(FFA_values~Gen_abund) , col=2, lwd= 2)
}

dev.off()


