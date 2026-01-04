################ PREPARING THE ENVIRONMENT ###############

{
  library("ggplot2")
  library("reshape2")
  library("ggh4x")
  library("DESeq2")
}

options(scipen=100)

suppressWarnings( dir.create("RNA_PHA_Results") )
suppressWarnings( dir.create("RNA_PHA_Results/Differentially_expressed") )

fill_color_30_BIS<-c("chartreuse","deeppink", "firebrick3",
                     "orange","darkmagenta", 
                     "darkblue","grey60", "grey95",
                     "bisque","cyan","yellow3","brown4",
                     "springgreen2", "grey20", 
                     "deepskyblue2","darkgreen","orange3",
                     "darkorchid", "lightblue1" , "violet", "blue",
                     "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                     "lightgrey", "red","darkslategray3")


load(file="Tabelle_RNA_PHA_after_filters.RData")




############## *** DEFINING THE PHA RELATED GENES *** #################

# PHA related genes (see paper introduction for references)
only_PHAgenes_IDs<-dictionary [ grepl("^pha", dictionary$Gene ) |
                                  grepl("^phb", dictionary$Gene ) |
                                  grepl("^acs", dictionary$Gene ) |   # Acetyl synthetase, see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^alkK", dictionary$Gene ) |   # see DOI: 10.1080/21655979.2025.2458363
                                  grepl("^croR", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^ccr", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^ecm", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^fad", dictionary$Gene ) |   # DOI 10.1080/21655979.2025.2458363 (fatty acid degradation)
                                  grepl("^mcd", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^mch", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^mcm", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^mce", dictionary$Gene ) |   # methylmalonyl-CoA epimerase, see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^epi", dictionary$Gene ) |   # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("[P-p]hasin", dictionary$Description )  |
                                  grepl("3-hydroxycarboxylic", dictionary$Description , fixed=T )  |
                                  grepl("3-hydroxyvalerate", dictionary$Description , fixed=T )  |
                                  grepl("3-enoyl-CoA hydratase", dictionary$Description , fixed=T )  |
                                  grepl("3-ketoacyl-CoA reductase", dictionary$Description , fixed=T )  |
                                  grepl("3-hydroxyacyl-CoA", dictionary$Description , fixed=T )  |
                                  grepl("3-hydroxyacyl-ACP thioester", dictionary$Description , fixed=T )  |  # 10.1080/21655979.2025.2458363 
                                  grepl("cetoacetyl-CoA thiolase", dictionary$Description , fixed=T )  |    # 10.1080/21655979.2025.2458363 
                                  grepl("Ketoacyl-CoA thiolase", dictionary$Description , fixed=T )  |
                                  grepl("ketoacyl-CoA thiolase", dictionary$Description , fixed=T )  |
                                  grepl("2,3-trans-enoyl-CoA", dictionary$Description , fixed=T )  |      # 10.1080/21655979.2025.2458363 
                                  grepl("[P-p]olyhydroxyalkanoate", dictionary$Description )  |
                                  grepl("[P-p]olyhydroxybutyr", dictionary$Description )  |
                                  grepl("-hydroxyalkanoate", dictionary$Description , fixed=T)  |
                                  grepl("hydroxyalkanoic acid", dictionary$Description , fixed=T)  |
                                  grepl("Poly(3-hydroxyalkanoate)", dictionary$Description , fixed=T)  |
                                  grepl("poly(3-hydroxyalkanoate)", dictionary$Description, fixed=T )  |
                                  grepl("Poly(3-hydroxybutyrate)", dictionary$Description )  |
                                  grepl("Polyhydroxyalkanoic", dictionary$Description )  |
                                  grepl("[P-p]oly(R)-hydroxyalkanoic", dictionary$Description )  |
                                  grepl("-hydroxybutyr", dictionary$Description , fixed=T )  |
                                  grepl("hydroxybutanoic", dictionary$Description , fixed=T )  |   # product of degradation, see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("Acetyl-CoA acetyltransf", dictionary$Description , fixed=T )  | # simil phbA
                                  grepl("Acetoacetyl-CoA reduct", dictionary$Description , fixed=T )  |  # simil phbA
                                  grepl("acetoacetyl-CoA reduct", dictionary$Description , fixed=T )  |  # simil phbA
                                  grepl("^Crotonoyl-CoA", dictionary$Description , fixed=T )  | #  # NOT METHYL-! see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl(" Crotonoyl-CoA", dictionary$Description , fixed=T )  | #  # NOT METHYL-! see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl(" crotonoyl-CoA", dictionary$Description , fixed=T )  | # NOT METHYL-! see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^[E-e]thylmalonyl", dictionary$Description )  | # see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("Erythro-3-methylmalyl-CoA", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("erythro-3-methylmalyl-CoA", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("ketothiolase", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("Mesaconyl-CoA", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("mesaconyl-CoA", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("Methylsuccinyl-CoA", dictionary$Description , fixed=T )  | #  with no "-" after methyl, see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("-methylsuccinyl-CoA", dictionary$Description , fixed=T )  | #  with no "-" after methyl, see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("^Mce related", dictionary$Description , fixed=T )  | #  methylmalonyl-CoA epimerase, see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("Propionyl-CoA synth", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("propionyl-CoA synth", dictionary$Description , fixed=T )  | #  see fig 7 of paper DOI 10.1016/j.jece.2025.116522
                                  grepl("PHB depol", dictionary$Description , fixed=T )  |  #
                                  grepl("PHA/PHB", dictionary$Description , fixed=T )  |  #
                                  grepl("PHB/PHA", dictionary$Description , fixed=T )  |  #
                                  grepl("PHV", dictionary$Description , fixed=T )  ,
]


### Searching for other pha genes using UniprotKB ...
## downloading results as table from https://www.uniprot.org/uniprotkb?query=Polyhydroxyalkanoate 
# system("wget -r https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28Polyhydroxyalkanoate%29 -O pha_results_UniProtKB.gz ")
# system("gzip -d pha_results_UniProtKB.gz ")
# file.rename("pha_results_UniProtKB","pha_results_UniProtKB.tsv")
UniProt_KB_IDs<-read.delim("pha_results_UniProtKB.tsv", header = T, sep = "\t")
UniProt_KB_IDs$Gene.Names <- UniProt_KB_IDs$Gene.Names
UniProt_KB_IDs<- UniProt_KB_IDs[ !duplicated(UniProt_KB_IDs$Gene.Names) , ]

### removing genes which are too *in*directly involved...
UniProt_KB_IDs <- UniProt_KB_IDs[!UniProt_KB_IDs$Gene.Names%in%c("thiG","ychF","cat","ilvD","ccrB","dna") , ] 
UniProt_KB_IDs <- UniProt_KB_IDs[!grepl("DNA|Dna ", UniProt_KB_IDs$Protein.names) , ]
only_PHAgenes_IDs <- only_PHAgenes_IDs[!only_PHAgenes_IDs$Gene%in%c("thiG","ychF","cat","ilvD","ccrB","dna") , ] 
only_PHAgenes_IDs <- only_PHAgenes_IDs[!grepl("4-hydroxyb", only_PHAgenes_IDs$Description) , ]
only_PHAgenes_IDs <- only_PHAgenes_IDs[!grepl("^phaAB", only_PHAgenes_IDs$Gene) , ]  # NOT involved with bioplastics

only_PHAgenes_IDs <- only_PHAgenes_IDs[!grepl("Adineta|Rotaria|Diploscapter|Carchesium|Stentor|Nippostrongylus",only_PHAgenes_IDs$Taxon) , ]

### adding the "others" pha related genes to the list
UniProt_Plus <- UniProt_KB_IDs [ UniProt_KB_IDs$Entry %in% dictionary$Entry , "Entry" ]
UniProt_Plus <- UniProt_Plus [! UniProt_Plus %in% only_PHAgenes_IDs$Entry ]
length(UniProt_Plus)   # other 26 entries retrieved using UniProt
rm(UniProt_KB_IDs)


PHA_IDs <- c( only_PHAgenes_IDs$Entry , UniProt_Plus )




######################## COMPUTING THE PROPORTIONS ##############################

if(! "proof_filter" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
}

### Proportions computed after the filters and with the unmapped
prop_with_unmap <- apply(table_abs, MARGIN = 2, function(x) x/sum(x))
prop_with_unmap <- prop_with_unmap * 100

### Proportions computed after the filters BUT without the unmapped (for statistical tests or alike)
prop_with_NO_unmap <- apply(table_abs[! row.names(table_abs) %in% "UNMAPPED", ], MARGIN = 2, function(x) x/sum(x))
prop_with_NO_unmap <- prop_with_NO_unmap * 100

### Gene proportions
prop_GENE_with_unmap <- apply(table_gene, MARGIN = 2, function(x) x/sum(x))
prop_GENE_with_unmap <- as.data.frame(prop_GENE_with_unmap) * 100




############ SIGNIF DIFFERENCES S vs F in PHA REACTORS  (TRANSCRIPTS) ################

#        CPM         %          
# 0.5:1.000.000 = x:100 --> x= 0.00005 
keep  <- rowSums(prop_with_unmap> 0.00005) >= 2   
# reference about the threshold: PMID 27508061 (in CPM) ... the ">=2" is because the samples for groups are 3 --> at least in 2!
# cts   <- table_abs[rowSums(table_abs)>10, ] 
# cts   <- cts[row.names(cts) %in% names(which(keep)), ]
cts   <- table_abs[names(which(keep)),  ]       # reference: DESeq2 base tutorial
cts   <- as.matrix( cts[ ! row.names(cts) %in% "UNMAPPED", ] )
cat("testing", length(cts[,1]) , "on" , length(table_abs[,1]) )
coldata<-metadata
row.names(coldata)<-coldata$Sample_name
coldata$Time_point<-paste0("character",coldata$Time_point)  # to avoid being used as numeric (use "Experiment day" for a more precise flow of time indicator)
coldata$Time_point<- as.factor(coldata$Time_point)
coldata$Reactor<- as.factor(coldata$Reactor)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Time_point+Reactor)

DE<-DESeq(dds)  # computing (this will take at least a minute... have a coffie!)

all_res<-results(DE, contrast= c("Reactor", "R1", "R2"))   # R1 is the base level --> compared to R2 (+ foldchange is higher in R1)
all_res <- all_res[order(all_res$log2FoldChange, na.last=NA, decreasing = T), ]
# hist(all_res$baseMean, breaks = 500, xlim = c(0,250))    # zoom
all_res <-all_res[all_res$padj < 0.05 & !is.na(all_res$padj) & abs(all_res$log2FoldChange)>0.5 & all_res$baseMean> 50 , ]




################### FOCUSING ON THE MOST DIFFERENT TRANSCRIPTS ###################

res_pos <-all_res[1:25, ]   # already ordered ---> positives ---> increased in R1
res_neg<-all_res[(length(all_res$log2FoldChange)-24):(length(all_res$log2FoldChange)), ] #  increased in R2


##### plotting increased functions in R1 ...
prop<-prop_with_unmap[row.names(res_pos), ]
colnames(prop)<- paste0(metadata$Reactor, "_", metadata$Experiment_day) # , "th")
to_plot <- melt(prop)
to_plot$Reactor<-as.character(to_plot$Var2)
to_plot$Reactor<-gsub("_.*","", to_plot$Reactor)
to_plot$Reactor<-gsub("R1","1/1", to_plot$Reactor)
to_plot$Reactor<-gsub("R2","1/2", to_plot$Reactor)
to_plot$Reactor<-factor(to_plot$Reactor, levels = c("1/1","1/2"))
to_plot$Time<-as.character(to_plot$Var2)
to_plot$Time<-gsub(".*_","", to_plot$Time)
to_plot$Xaxis<-paste0(to_plot$Var1,to_plot$Reactor)  # Var1 is the Protein
ggplot(to_plot, aes(x= Reactor, y=value, fill=Reactor)) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("1/1"="coral","1/2"="deepskyblue")) +
  facet_wrap2( .~Var1,
               labeller = labeller(group = label_wrap_gen(width = 34)),
               scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Reactor), size=1.25, alpha=1) +
  geom_point(aes(color=Reactor), size=2.5, alpha=0.4) +
  geom_text(aes(label=Time), size=2.1, color="gray20") +
  geom_line(aes(group=Time), linewidth=0.1, alpha=0.5) +
  scale_x_discrete(labels=unique(levels(to_plot$Reactor)), expand=c(0,0.5)) +
  theme(strip.text.x=element_text(size=7.5,colour="black"),
        strip.switch.pad.wrap = unit(6,"line") ,
        axis.text.y = element_text(size=5.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title= element_text(size=8),
        panel.grid.minor.y= element_blank(),
        plot.margin=margin(1, 1, 1, 1, unit="pt"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.margin=margin(-6, 0, -2, 0),
        legend.position = "bottom") +
  labs(title="Twentyfive most significant increased transcripts in reactor R1",
       y="Percent counts",
       color="F/F ratio", fill="F/F ratio"
  )
ggsave(filename = "RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_RNA_in_R1.png", width = 6.2, height = 5.2, dpi=300)
dev.off()

# saving the info about the results too
dic_res<-dictionary
row.names(dic_res)<-dic_res$Entry
dic_res$Entry<-NULL
table_res<-cbind.data.frame(res_pos[ , c(2,4,6)], dic_res[row.names(res_pos), ])
write.csv2(file="RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_RNA_in_R1_infos.csv", table_res, row.names=T, quote=F)



###### plotting increased functions in R2 ...
prop<-prop_with_unmap[row.names(res_neg), ]
colnames(prop)<- paste0(metadata$Reactor, "_", metadata$Experiment_day) #, "th")
to_plot <- melt(prop)
to_plot$Reactor<-as.character(to_plot$Var2)
to_plot$Reactor<-gsub("_.*","", to_plot$Reactor)
to_plot$Reactor<-gsub("R1","1/1", to_plot$Reactor)
to_plot$Reactor<-gsub("R2","1/2", to_plot$Reactor)
to_plot$Reactor<-factor(to_plot$Reactor, levels = c("1/1","1/2"))
to_plot$Time<-as.character(to_plot$Var2)
to_plot$Time<-gsub(".*_","", to_plot$Time)
to_plot$Xaxis<-paste0(to_plot$Var1,to_plot$Reactor)  # Var1 is the Protein
ggplot(to_plot, aes(x= Reactor, y=value, fill=Reactor)) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("1/1"="coral","1/2"="deepskyblue")) +
  facet_wrap2( .~Var1,
               labeller = labeller(group = label_wrap_gen(width = 34)),
               scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Reactor), size=1.25, alpha=1) +
  geom_point(aes(color=Reactor), size=2.5, alpha=0.4) +
  geom_text(aes(label=Time), size=2.1, color="gray20" ) +
  geom_line(aes(group=Time), linewidth=0.1, alpha=0.5) +
  scale_x_discrete(labels=unique(levels(to_plot$Reactor)), expand=c(0,0.5)) +
  theme(strip.text.x=element_text(size=7.5,colour="black"),
        strip.switch.pad.wrap = unit(6,"line") ,
        axis.text.y = element_text(size=5.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title= element_text(size=8),
        panel.grid.minor.y= element_blank(),
        plot.margin=margin(1, 1, 1, 1, unit="pt"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.margin=margin(-6, 0, -2, 0),
        legend.position = "bottom") +
  labs(title="Twentyfive most significant increased trascripts in reactor R2",
       y="Percent counts",
       color="F/F ratio", fill="F/F ratio"
  )
ggsave(filename = "RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_RNA_in_R2.png", width = 6.2, height = 5.2, dpi=300)
dev.off()

# saving the info about the results too
dic_res<-dictionary
row.names(dic_res)<-dic_res$Entry
dic_res$Entry<-NULL
table_res<-cbind.data.frame(res_neg[ , c(2,4,6)], dic_res[row.names(res_neg), ])
write.csv2(file="RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_RNA_in_R2_infos.csv", table_res, row.names=T, quote=F)


suppressWarnings( rm(res, res_neg, res_pos, dic_res, table_res, to_plot, coldata, cts) )




############# FOCUSING ON PHA RELATED DIFF mRNAs ##############

res_PHA <- as.data.frame(all_res)
res_PHA <- res_PHA[ row.names(res_PHA)%in%PHA_IDs , ]

prop<-prop_with_unmap[row.names(res_PHA), ]
colnames(prop)<- paste0(metadata$Reactor, "_", metadata$Experiment_day) #, "th")
to_plot <- melt(prop)
to_plot$Reactor<-as.character(to_plot$Var2)
to_plot$Reactor<-gsub("_.*","", to_plot$Reactor)
to_plot$Reactor<-gsub("R1","1/1", to_plot$Reactor)
to_plot$Reactor<-gsub("R2","1/2", to_plot$Reactor)
to_plot$Reactor<-factor(to_plot$Reactor, levels = c("1/1","1/2"))
to_plot$Time<-as.character(to_plot$Var2)
to_plot$Time<-gsub(".*_","", to_plot$Time)
to_plot$Xaxis<-paste0(to_plot$Var1,"_",to_plot$Reactor)  # Var1 is the mRNA
# to_plot$value<-round( to_plot$value, 5 )   # better use the package scale in scaleYcontinuous

# using descriptions in the plot ...
selected_infos <- dictionary[ dictionary$Entry%in% row.names(res_PHA) , ]
selected_infos <- selected_infos[! duplicated(selected_infos$Entry) , ]
row.names(selected_infos) <- selected_infos$Entry
duplications <- selected_infos$Description[duplicated(selected_infos$Description) & duplicated(selected_infos$Taxon)] # NB: only a duplication will be True
duplications <- selected_infos$Description %in% duplications # now every duplication is True
selected_infos[  duplications , "Description"] <- paste0( selected_infos[duplications,"Description"] , "\n", selected_infos[duplications,"Entry"] )

selected_infos <- selected_infos[to_plot$Var1 , ]  # ordered repetition according to the table

to_plot$infos<-paste0( selected_infos$Description,"\n(", selected_infos$Taxon,")" )
# unique(to_plot$infos)
to_plot$infos <- gsub("Acetoacetyl-CoA reduct PhbB2","Acetoacetyl-CoA reduct\nPhbB2", to_plot$infos , fixed = T)
to_plot$infos <- gsub("Acyl-coenzyme A dehydrog","Acyl-CoA dehydrogenase", to_plot$infos , fixed = T)
to_plot$infos <- gsub("class I poly(R)-hydroxyalkanoic acid synth","PHA synthetase class I", to_plot$infos , fixed = T )
to_plot$infos <- gsub("class I poly(R)-hydroxyal. acid synth","PHA synthetase class I", to_plot$infos , fixed = T )
to_plot$infos <- gsub("Class I poly(R)-hydroxyalkanoic acid synth","PHA synthetase class I", to_plot$infos , fixed = T )
to_plot$infos <- gsub("3-hydroxybutyryl-CoA","3HB-CoA", to_plot$infos , fixed=T)
# to_plot$infos <- gsub("dehydrog","dehydrogen", to_plot$infos , fixed=T)
to_plot$infos <- gsub("Poly(3-hydroxybutyrate)","P3HB", to_plot$infos , fixed=T)
to_plot$infos <- gsub("Poly(3-hydroxybutyrate)","P3HB", to_plot$infos, fixed = T )
to_plot$infos <- gsub("Poly (3-hydroxybutyrate)","P3HB", to_plot$infos , fixed = T )
to_plot$infos <- gsub("3-hydroxybutyrate","3HB", to_plot$infos , fixed = T )
to_plot$infos <- gsub("[P-p]olyhydroxyalkanoate","PHA", to_plot$infos )
# to_plot$infos <- gsub("poly(R)-hydroxyalkanoic", "poly(R)-hydroxyal.", to_plot$infos, fixed=T)
to_plot$infos <- gsub("depol\n(","depolimerase\n(", to_plot$infos , fixed=T )
to_plot$infos <- gsub("fam prot","family", to_plot$infos )
to_plot$infos <- gsub("^phasin","Phasin", to_plot$infos )
to_plot$infos <- gsub("Phasin prot","Phasin", to_plot$infos )
to_plot$infos <- gsub("PHA depol, intracellular","Intracellular PHA depolim", to_plot$infos )
to_plot$infos <- gsub("3HB dehydrog","3HB dehydrogenase", to_plot$infos )
to_plot$infos <- gsub("3HB-CoA dehydrog","3HB-CoA dehydrogenase", to_plot$infos )
to_plot$infos <- gsub("Granule-associated prot (PHASIN)","Phasin family", to_plot$infos , fixed = T )
to_plot$infos <- gsub("Oceanicaulis sp. HLUCCA04","Oceanicaulis HLUCCA04", to_plot$infos )
to_plot$infos <- gsub("Neomegalonema perideroedes","Neomegalonema perider.", to_plot$infos )
to_plot$infos <- gsub("uncult bact","uncultured bacterium", to_plot$infos )
to_plot$infos <- gsub("bact)","bacterium)", to_plot$infos )
to_plot$infos <- gsub("Parazoarcus communis SWub3 = DSM 12120", "Parazoarcus communis", to_plot$infos, fixed=T)
to_plot$infos <- gsub("Paracoccus aestuariivivens", "Paracoccus aestuariiviv.", to_plot$infos, fixed=T)
to_plot$infos <- gsub("phenylacetica B4P", "phenylacetica", to_plot$infos, fixed=T)
# unique(to_plot$infos)

ggplot(to_plot, aes(x= Reactor, y=value, fill=Reactor)) +
  theme_classic(base_size = 11) +
  scale_color_manual(values = c("1/1"="coral","1/2"="deepskyblue2")) +
  facet_wrap2( .~infos, nrow = 8,
               labeller = labeller(group = label_wrap_gen(width = 34)),
               scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Reactor), size=0.85, alpha=1) +
  geom_point(aes(color=Reactor), size=1.8, alpha=0.4) +
  geom_text(aes(label=Time), size=1.25, color="gray15",
            hjust=ifelse(to_plot$Reactor=="1/1",1.35,-0.35)
  ) +
  geom_line(aes(group=Time), linewidth=0.1, alpha=0.85, color="gray45") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001) ) +
  scale_x_discrete(labels=unique(levels(to_plot$Reactor)), expand=c(0,0.5)) +
  theme(strip.text.x=element_text(size=5.2,colour="black", 
                                  lineheight = 0.8,
                                  margin = margin(2,0,2,0, "pt")
                                  ),
        strip.switch.pad.wrap = unit(6,"line") ,
        axis.text.y = element_text(size=3.85, hjust = 0.6),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth =0.25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title= element_text(size=5),
        panel.grid.minor.y= element_blank(),
        # panel.border = element_rect(linewidth =0.1, fill = NA ),
        strip.background = element_rect(color = "black", size = 0.65),        
        panel.spacing.x = unit(0.15, "line"),
        plot.margin=margin(1, 1, 1, 1, unit="pt"),
        legend.title = element_text(size=12.5),
        legend.text = element_text(size=12),
        legend.margin=margin(-8, 0, -2, 0),
        legend.position = "bottom") +
  labs(
    # title="PHA related signif diff transcripts",
    y="Percent counts",
    color="F/F ratio", fill="F/F ratio"
  ) +
  guides ( color="none", fill="none")
ggsave(filename = "RNA_PHA_Results/Differentially_expressed/DESeq2_transcripts_PHArelated.png", width = 4.9, height = 5.85, dpi=300)
dev.off()

# saving the info about the results too
dic_res<-dictionary
row.names(dic_res)<-dic_res$Entry
dic_res$Entry<-NULL
table_res<-cbind.data.frame(   dic_res[row.names(res_PHA), ]   ,  round(res_PHA[ , c(2,4,6)],3)   )
table_res$padj[table_res$padj==0.000] <- "<0.001"
table_res$Description <- gsub(" prot","", table_res$Description,fixed=T)
table_res$Description <- gsub(" fam"," family protein", table_res$Description,fixed=T)
table_res$Description <- gsub("synth","synthetase", table_res$Description,fixed=T)
table_res$Description <- gsub("Granule-associated (PHASIN)","Phasin", table_res$Description,fixed=T)
table_res$Description <- gsub("^phasin","Phasin", table_res$Description)
table_res$Description <- gsub("[P-p]olyhydroxyalkanoate","PHA", table_res$Description)
table_res$Description <- gsub("PHA depol, intracellular","PHA depolimerase (intracellular)", table_res$Description,fixed=T)
table_res$Description <- gsub(" reduct"," reductase", table_res$Description,fixed=T)
table_res$Taxon <- gsub("uncult ","uncultured ", table_res$Taxon)
table_res$Taxon <- gsub("= ","", table_res$Taxon, fixed=T)
write.csv2(file="RNA_PHA_Results/Differentially_expressed/DESeq2_transcripts_PHArelated_infos.csv", table_res, row.names=T, quote=F)

suppressWarnings( rm(res, res_PHA, res_pos, dic_res, table_res, to_plot, coldata, cts) )




############ SIGNIF DIFFERENCES S vs F in PHA REACTORS  (GENES) ################

#        CPM         %          
# 0.5:1.000.000 = x:100 --> x= 0.00005 
keep  <- rowSums(prop_GENE_with_unmap> 0.00005) >= 2   
# reference about the threshold: PMID 27508061 (in CPM) ... the ">=2" is because the samples for groups are 3 --> at least in 2!
# cts   <- table_abs[rowSums(table_abs)>10, ] 
# cts   <- cts[row.names(cts) %in% names(which(keep)), ]
cts   <- table_gene[names(which(keep)),  ]       # reference: DESeq2 base tutorial
cts   <- as.matrix( cts[ ! row.names(cts) %in% "UNMAPPED", ] )
cat("testing", length(cts[,1]) , "on" , length(table_abs[,1]) )
coldata<-metadata
row.names(coldata)<-coldata$Sample_name
coldata$Time_point<-paste0("character",coldata$Time_point)  # to avoid being used as numeric (use "Experiment day" for a more precise flow of time indicator)
coldata$Time_point<- as.factor(coldata$Time_point)
coldata$Reactor<- as.factor(coldata$Reactor)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Time_point+Reactor)

DE<-DESeq(dds)  # computing (this will take at least a minute... have a coffie!)

all_res<-results(DE, contrast= c("Reactor", "R1", "R2"))   # R1 is the base level --> compared to R2 (+ foldchange is higher in R1)
all_res <- all_res[order(all_res$log2FoldChange, na.last=NA, decreasing = T), ]
# hist(all_res$baseMean, breaks = 500, xlim = c(0,250))    # zoom
all_res <-all_res[all_res$padj < 0.05 & !is.na(all_res$padj) & abs(all_res$log2FoldChange)>0.5 & all_res$baseMean> 50 , ]




################### FOCUSING ON THE MOST DIFFERENT GENES ###################

res_pos <-all_res[1:25, ]   # already ordered ---> positives ---> increased in R1
res_neg<-all_res[(length(all_res$log2FoldChange)-24):(length(all_res$log2FoldChange)), ] #  increased in R2


##### plotting increased functions in R1 ...
prop<-prop_GENE_with_unmap[row.names(res_pos), ]
colnames(prop)<- paste0(metadata$Reactor, "_", metadata$Experiment_day) # , "th")
prop$gene <- row.names(prop)
to_plot <- melt(prop)
to_plot$Reactor<-as.character(to_plot$variable)
to_plot$Reactor<-gsub("_.*","", to_plot$Reactor)
to_plot$Reactor<-gsub("R1","1/1", to_plot$Reactor)
to_plot$Reactor<-gsub("R2","1/2", to_plot$Reactor)
to_plot$Reactor<-factor(to_plot$Reactor, levels = c("1/1","1/2"))
to_plot$Time<-as.character(to_plot$variable)
to_plot$Time<-gsub(".*_","", to_plot$Time)
to_plot$Xaxis<-paste0(to_plot$Var1,to_plot$Reactor)  # Var1 is the Protein
ggplot(to_plot, aes(x= Reactor, y=value, fill=Reactor)) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("1/1"="coral","1/2"="deepskyblue")) +
  facet_wrap2( .~gene,
               labeller = labeller(group = label_wrap_gen(width = 34)),
               scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Reactor), size=1.25, alpha=1) +
  geom_point(aes(color=Reactor), size=2.5, alpha=0.4) +
  geom_text(aes(label=Time), size=2.1, color="gray20") +
  geom_line(aes(group=Time), linewidth=0.1, alpha=0.5) +
  scale_x_discrete(labels=unique(levels(to_plot$Reactor)), expand=c(0,0.5)) +
  theme(strip.text.x=element_text(size=7.5,colour="black"),
        strip.switch.pad.wrap = unit(6,"line") ,
        axis.text.y = element_text(size=5.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title= element_text(size=8),
        panel.grid.minor.y= element_blank(),
        plot.margin=margin(1, 1, 1, 1, unit="pt"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.margin=margin(-6, 0, -2, 0),
        legend.position = "bottom") +
  labs(title="Twentyfive most significant increased genes in reactor R1",
       y="Percent counts",
       color="F/F ratio", fill="F/F ratio"
  )
ggsave(filename = "RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_GENES_in_R1.png", width = 6.2, height = 5.2, dpi=300)
dev.off()

# saving the info about the results too
dic_res<-dictionary[ dictionary$Gene %in% prop$gene , ]
row.names(dic_res)<-dic_res$Gene
dic_res$Gene<-NULL
table_res<-cbind.data.frame(res_pos[ , c(2,4,6)], dic_res[row.names(res_pos), ])
write.csv2(file="RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_GENES_in_R1_infos.csv", table_res, row.names=T, quote=F)



###### plotting increased functions in R2 ...
prop<-prop_GENE_with_unmap[row.names(res_neg), ]
colnames(prop)<- paste0(metadata$Reactor, "_", metadata$Experiment_day) # , "th")
prop$gene <- row.names(prop)
to_plot <- melt(prop)
to_plot$Reactor<-as.character(to_plot$variable)
to_plot$Reactor<-gsub("_.*","", to_plot$Reactor)
to_plot$Reactor<-gsub("R1","1/1", to_plot$Reactor)
to_plot$Reactor<-gsub("R2","1/2", to_plot$Reactor)
to_plot$Reactor<-factor(to_plot$Reactor, levels = c("1/1","1/2"))
to_plot$Time<-as.character(to_plot$variable)
to_plot$Time<-gsub(".*_","", to_plot$Time)
to_plot$Xaxis<-paste0(to_plot$Var1,to_plot$Reactor)  # Var1 is the Protein
ggplot(to_plot, aes(x= Reactor, y=value, fill=Reactor)) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("1/1"="coral","1/2"="deepskyblue")) +
  facet_wrap2( .~gene,
               labeller = labeller(group = label_wrap_gen(width = 34)),
               scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Reactor), size=1.25, alpha=1) +
  geom_point(aes(color=Reactor), size=2.5, alpha=0.4) +
  geom_text(aes(label=Time), size=2.1, color="gray20") +
  geom_line(aes(group=Time), linewidth=0.1, alpha=0.5) +
  scale_x_discrete(labels=unique(levels(to_plot$Reactor)), expand=c(0,0.5)) +
  theme(strip.text.x=element_text(size=7.5,colour="black"),
        strip.switch.pad.wrap = unit(6,"line") ,
        axis.text.y = element_text(size=5.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title= element_text(size=8),
        panel.grid.minor.y= element_blank(),
        plot.margin=margin(1, 1, 1, 1, unit="pt"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.margin=margin(-6, 0, -2, 0),
        legend.position = "bottom") +
  labs(title="Twentyfive most significant increased genes in reactor R2",
       y="Percent counts",
       color="F/F ratio", fill="F/F ratio"
  )
ggsave(filename = "RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_GENES_in_R2.png", width = 6.2, height = 5.2, dpi=300)
dev.off()

# saving the info about the results too
dic_res<-dictionary[ dictionary$Gene %in% prop$gene , ]
row.names(dic_res)<-dic_res$Gene
dic_res$Gene<-NULL
table_res<-cbind.data.frame(res_neg[ , c(2,4,6)], dic_res[row.names(res_neg), ])
write.csv2(file="RNA_PHA_Results/Differentially_expressed/DESeq2_Twentyfive_most_increased_RNA_in_R2_infos.csv", table_res, row.names=T, quote=F)

suppressWarnings( rm(res, res_neg, res_pos, dic_res, table_res, to_plot, coldata, cts) )




############# FOCUSING ON PHA RELATED DIFF GENES ##############

these_genes <- dictionary[ dictionary$Entry %in% PHA_IDs , ]
these_genes <- these_genes[these_genes$Gene!="" , ]  # PHA related *BUT* not linked to any gene --> Impossible to merge to other transcripts
these_genes <- these_genes[!these_genes$Gene%in%c("thiG","ychF","cat","ilvD","ccrB") , ] # not related

res_PHA <- as.data.frame(all_res)
res_PHA <- res_PHA[ row.names(res_PHA)%in%these_genes$Gene , ]

prop<-prop_GENE_with_unmap[row.names(res_PHA), ]
colnames(prop)<- paste0(metadata$Reactor, "_", metadata$Experiment_day) #, "th")
prop$gene <- row.names(prop)
to_plot <- melt(prop)
to_plot$Reactor<-as.character(to_plot$variable)
to_plot$Reactor<-gsub("_.*","", to_plot$Reactor)
to_plot$Reactor<-gsub("R1","1/1", to_plot$Reactor)
to_plot$Reactor<-gsub("R2","1/2", to_plot$Reactor)
to_plot$Reactor<-factor(to_plot$Reactor, levels = c("1/1","1/2"))
to_plot$Time<-as.character(to_plot$variable)
to_plot$Time<-gsub(".*_","", to_plot$Time)
to_plot$Xaxis<-paste0(to_plot$Var1,to_plot$Reactor)  # Var1 is the Protein
to_plot$value<-round( to_plot$value, 4 )

# using descriptions in the plot ...
selected_infos <- these_genes[ these_genes$Gene%in% row.names(res_PHA) , ]
selected_infos <- selected_infos[!duplicated(selected_infos$Gene) , ]
row.names(selected_infos) <- selected_infos$Gene
selected_infos <- selected_infos[to_plot$gene , ]  # ordered repetition according to the table
to_plot$infos<-paste0( selected_infos$Description,"\n(", selected_infos$Gene,")" )
# unique(to_plot$infos)
to_plot$infos <- gsub("3-hydroxybutyryl-CoA","3HB-CoA", to_plot$infos , fixed=T)
to_plot$infos <- gsub("dehydrog","dehydrogen", to_plot$infos , fixed=T)
to_plot$infos <- gsub("Poly(3-hydroxybutyrate)","P3HB", to_plot$infos , fixed=T)
to_plot$infos <- gsub("Poly(3-hydroxybutyrate)","P3HB", to_plot$infos, fixed = T )
to_plot$infos <- gsub("Poly (3-hydroxybutyrate)","P3HB", to_plot$infos , fixed = T )
to_plot$infos <- gsub("Polyhydroxyalkanoate","PHA", to_plot$infos )
to_plot$infos <- gsub("depol\n(","depolimerase\n(", to_plot$infos , fixed=T )
to_plot$infos <- gsub("depol, intracellular","intracellular depol", to_plot$infos )
to_plot$infos <- gsub("fam prot","family", to_plot$infos )
to_plot$infos <- gsub("Phasin prot","Phasin", to_plot$infos )

ggplot(to_plot, aes(x= Reactor, y=value, fill=Reactor)) +
  theme_classic(base_size = 11) +
  scale_color_manual(values = c("1/1"="coral","1/2"="deepskyblue")) +
  facet_wrap2( .~infos, nrow = 6,
               labeller = labeller(group = label_wrap_gen(width = 34)),
               scales = "free", strip=strip_nested(size = "variable", bleed = T),drop = TRUE) +
  geom_point(aes(color=Reactor), size=0.45, alpha=1) +
  geom_point(aes(color=Reactor), size=2, alpha=0.4) +
  geom_text(aes(label=Time), size=1.4, color="gray20",
            hjust=ifelse(to_plot$Reactor=="1/1",1.25,-0.25)
  ) +
  geom_line(aes(group=Time), linewidth=0.11, alpha=0.85, color="gray45") +
  scale_x_discrete(labels=unique(levels(to_plot$Reactor)), expand=c(0,0.5)) +
  theme(strip.text.x=element_text(size=6.7,colour="black", lineheight =0.76 ),
        strip.switch.pad.wrap = unit(6,"line") ,
        axis.text.y = element_text(size=4.3),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title= element_text(size=5),
        panel.grid.minor.y= element_blank(),
        panel.grid.major.y= element_line(color="gray", linewidth = 0.075),
        plot.margin=margin(1, 1, 1, 1, unit="pt"),
        legend.title = element_text(size=12.5),
        legend.text = element_text(size=12),
        legend.margin=margin(-10, 0, -2, 0),
        legend.position = "bottom") +
  labs(
    # title="PHA related signif diff transcripts",
    y="Percent counts (whole dataset)",
    color="F/F ratio", fill="F/F ratio"
  )
ggsave(filename = "RNA_PHA_Results/Differentially_expressed/DESeq2_GENES_PHA_related.png", width = 4.2, height = 5, dpi=300)
dev.off()

suppressWarnings( rm(res, res_PHA, res_pos, dic_res, table_res, to_plot, coldata, cts) )

