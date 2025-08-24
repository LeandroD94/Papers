# 26/07/2025

### Script part 1: preparing the required objs and checking the rarefactions



##################### PREPARING THE ENVIRONMENT #################

library("phyloseq")
library("vegan")
options(scipen = 100)


{ 
  dir.create("Data_check")
  dir.create("Results")
  dir.create("Results/Abundances")
  dir.create("Results/Comparison_with_particular_AGS")
  dir.create("Results/Alpha_diversity")
  dir.create("Results/PCoA_with_Every_microbes")
  dir.create("Results/PCoA_with_Prokaryotes_only")
  dir.create("Results/PCoA_with_Eukar_only")
  dir.create("Results/PCoA_with_Viruses_only")
  #dir.create("Results/Full_vs_Pilot_vs_Lab")
}

options(scipen = 100) # disable scientific annotation



### Colors for the stacked bar plots 
# choosing colors --> (see grDevices::colors() )
fill_color_5<-c("magenta3", "gold2", "firebrick3","springgreen2","deepskyblue2", "darkslategray3") 
fill_color_10<-c("springgreen2","darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet",  "darkslategray3")
fill_color_15<-c("brown3","springgreen2","wheat","darkmagenta","coral","yellow3","magenta","pink3", "blue2","firebrick2","gold","gray","chartreuse3","violet", "deepskyblue2","darkslategray3")
fill_color_19<-c("darkblue","brown4","springgreen2","wheat","lightcoral","coral","yellow3","darkmagenta","pink3", "blue","firebrick3","gray","gold","darkgreen","violet", "deepskyblue2","wheat3","red","chartreuse3","darkslategray3")
fill_color_25<-c("springgreen2","deeppink","darkmagenta",
                 "darkblue","grey50",
                 "bisque2","cyan","yellow3","brown",
                 "green", "grey20",
                 "deepskyblue2","lightgreen","orange2",
                 "violet", "blue",
                 "yellow", "pink3","yellow4", "chocolate4",
                 "coral2","darkgreen",
                 "lightgrey", "red","darkslategray3")
fill_color_30<-c("wheat3","deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "aquamarine2",
                 "bisque2","cyan","yellow3","brown",
                 "springgreen4", "grey20", "firebrick3",
                 "deepskyblue2","lightgreen","orange2",
                 "darkorchid", "lightblue1" , "violet", "blue1",
                 "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                 "lightgrey", "red","darkslategray3")
fill_color_30_BIS<-c("chartreuse","deeppink", "firebrick3",
                     "darkcyan","darkmagenta", 
                     "darkblue","grey55", "aquamarine3",
                     "bisque","cyan","yellow3","brown4",
                     "springgreen4", "grey20", 
                     "deepskyblue4","lightgreen","orange2",
                     "darkorchid", "lightblue1" , "violet", "blue",
                     "yellow", "pink3","yellow4", "chocolate4","coral2","darkgreen",
                     "lightgrey", "red","darkslategray3")




##### IMPORTING THE CLASSIFICATION DATA

load(file="Phyloseq_objs_from_Kaiju_classif.RData")

# removing AS samples
data <- prune_samples(!sample_names(data) %in% c("1623906", "SRR22378311","SRR24975224") , data)
data.genus <- prune_samples(!sample_names(data.genus) %in% c("1623906", "SRR22378311","SRR24975224") , data.genus)



###### BUILDING THE SAMPLE DATA
Metadata <- as.data.frame(read.table(file = "metadata.txt", header = T, sep="\t"))
# Three samples from project PRJNA783874 (8) had anomalous issues during their download ... Removing them
Metadata <- Metadata[! Metadata$FASTQ_ID %in% c("SRR17089439", "SRR17089440", "SRR17089441" ) , ]
row.names(Metadata)<-Metadata$FASTQ_ID # column with FASTQ/SAMPLE name
#head(Metadata)
originals<- Metadata$FASTQ_ID
original_length<-length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])
#Metadata$FASTQ_ID [! Metadata$FASTQ_ID %in% sample_names(data) ]
#sample_names(data) [ ! sample_names(data) %in% Metadata$FASTQ_ID ]
Metadata<-Metadata[sample_names(data), ]

identical(as.numeric(length(Metadata$FASTQ_ID[!is.na(Metadata$FASTQ_ID)])),as.numeric(original_length))

sample_data(data) <- Metadata
sample_names(data)<-sample_data(data)$Sample_name
sample_names(data)<-factor(sample_data(data)$Sample_name, levels=unique(sample_data(data)$Sample_name))
sample_data(data)$Reactor_number<-factor(sample_data(data)$Reactor_number,
                                         levels = seq(min(sample_data(data)$Reactor_number),
                                                      max(sample_data(data)$Reactor_number),1 )
)
sample_data(data)$Project_grid<-paste0("Proj_",sample_data(data)$Project_number)
sample_data(data)$Project_grid<-factor( sample_data(data)$Project_grid,
                                        levels = paste0( "Proj_", seq(seq(min(sample_data(data)$Project_number),max(sample_data(data)$Project_number),1) ) )
)
sample_data(data)$Project_number<-factor(sample_data(data)$Project_number,
                                         levels = seq(min(sample_data(data)$Project_number),max(sample_data(data)$Project_number),1)
)
sample_data(data)$Reactor_scale<-factor(sample_data(data)$Reactor_scale,
                                        levels = c("Full","Pilot","Lab")
)
sample_data(data)$Influent<-factor(sample_data(data)$Influent,
                                   levels = c("Real","Synthetic")
)

sample_names(data.genus)<-sample_data(data)$Sample_name
sample_data(data.genus)<-sample_data(data)



### solving two imperfections of the genera taxonomy
temp <- tax_table(data.genus)[ , "Genus"]
temp[temp=="Candidatus_Propionivibrio"] <- "Propionivibrio"  # it is not Candidatus at genus level, only few species...
tax_table(data.genus)[ , "Genus"] <- temp
data.genus <- tax_glom(data.genus, taxrank = "Genus", NArm=F )
rm(temp)




################## DEFINING THE KEY ORGANISMS GROUPS #################

### FUNCTIONAL GROUPS

# common PAO GAO list (PMID:38524765 ,  PMID: 22827168,  PMID: 15774634,  PMID: 33404187, MIDAS website)
PAO_list<-c("Candidatus Phosphoribacter","Ca Phosphoribacter","Phosphoribacter",
            "Tetrasphaera", "Nostocoides", # Synonyms
            "Candidatus Accumulimonas","Accumulimonas","Ca Accumulimonas",
            "Candidatus Accumulibacter","Accumulibacter","Ca Accumulibacter",
            "Candidatus Microthrix","Microthrix","Ca Microthrix", "Candidatus Neomicrothrix",
            "Candidatus Dechloromonas","Dechloromonas",
            "Malikia",
            "Quatrionicoccus",
            "Beggiatoa",
            "Gemmatimonas",
            "Friedmaniella",
            "Tessaracoccus",
            "Azonexus",
            "Knoellia", # https://doi.org/10.1016/j.ese.2024.100387
            "Phycicoccus", # https://doi.org/10.1016/j.ese.2024.100387
            "Lampropedia", # https://www.sciencedirect.com/science/article/pii/S004313549600351X
            # "Pseudomonas",   # too aspecific!
            "Candidatus Accumulimonas","Accumulimonas","Ca Accumulimonas",
            "Microlunatis","Microlunatus",
            "Candidatus Lutibacillus","Lutibacillus","Ca Lutibacillus")
GAO_list<-c("Candidatus Competibacter","Candidatus Contendobacter","Candidatus Proximibacter",
            "Ca Competibacter","Ca Contendobacter","Ca Proximibacter",
            "Competibacter","Contendobacter",
            # "Proximibacter", # may also be a PAO!
            "Defluviicoccus","Micropruina",
            "Candidatus Propionivibrio", "Propionivibrio")

# AOBs derive from from MIDAS, also in other articles I can't find other names...
AOB_list <- c("Nitrosomonas", "Nitrosococcus",   # synonymous, see MIDAS
              "Nitrosospira", "Nitrosolobus", "Nitrosovibrio",   # synonymous, see MIDAS
              "966-1",     # https://www.sciencedirect.com/science/article/pii/S0045653523014649
              "Ellin6067" )   # https://www.sciencedirect.com/science/article/pii/S0960852423008970
AOB_list <- c(AOB_list, "Candidatus_Nitrososphaera" , "Candidatus Nitrososphaera", "Nitrososphaera", "Nitrosocaldus", "Nitrosotalea")  # adding also AOA, from PMID: 24559743
NOB_list <- c("Nitrospira","Candidatus_Nitrospira", "Candidatus Nitrospira", "Nitrospirae","Nitrospirales","Nitrospiraceae",  # alcuni da MIDAS, ma tutti elencati in DOI: 10.1016/S0076-6879(11)86005-2
              # few Nitrospira names are quite similar, these are synonims or, at least, they still are of Nitrospirales order which is still known for AO activity, see https://www.nature.com/articles/s41396-023-01518-6
              "Nitrospina","Nitrospinae","Nitrospinaceae","Nitrospinota",
              "Nitrobacter",
              "Nitrococcus",
              "Nitrotoga", "Candidatus Nitrotoga", "Candidatus_Nitrotoga")

# PHA accumulators
table_PHA<-read.table(file="PHA_Accumulators_maggio2025.tsv", fill = T, header = T, sep="\t")




### SESSILE CILIATES 

sessilida_list <- c("Astylozoon",   # from NCBI taxonomy
                    "Ambiphrya",
                    "Apocarchesium",
                    "Campanella",
                    "Carchesium",
                    "Epicarchesium",
                    "Cothurnia",
                    "Epistylis",
                    "Hastatella",
                    "Lagenophrys",
                    "Myoschiston",
                    "Opisthostyla",
                    "Orborhabdostyla",
                    "Opercularia",
                    "Ophrydium",
                    "Parapiosoma",
                    "Platycola",
                    "Pseudepistylis",
                    "Pseudovorticella",
                    "Planeticovorticella",
                    "Platycola",
                    "Pyxidium",
                    "Rhabdostyla",
                    "Sessilida sp.",
                    "Telotrochidium",
                    "Thuricola",
                    "Usconophrys",
                    "Vaginicola",
                    "Vorticella",
                    "Vorticellides",
                    "Zoothamnium",
                    "Zoothamnopsis"
)
# not from sessilida order but still sessile ciliates (e.g. see http://web.tiscali.it/ifts/protozoi_ciliati_sessili_i.htm )
sessilida_list <- c(sessilida_list, 
                    "Stentor",   # Ciliato dell'ordine Heterotrichida, appartenente alla Classe Polyhymenophora. A differenzia dei Peritrichi, non presenta lo stadio di telotroco; esso infatti. pur essendo forma sessile, è in grado di spostarsi all'occorrenza per colonizzare nuovi substrati. La forma è quella di una tromba ma. se disturbato, può contrarsi assumendo la forma di una palla
                    "Folliculina",   # Heterotrichida, vive in una teca chiamata lorica, stesso ordine di Stentor
                    "Prostentor"     # Heterotrichida, simile a stentor ma meno comune
)




################## SAVING THE OBJECT READY FOR THE ANALYSIS ####################

save.image("Data_prepared_for_analysis.RData")




############################ RAREFACTION ANALYSES ################################

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
    if(length(x) < 2) {
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
  # legend("bottomright",paste(sat,"saturated samples"),bty="n")
}



steps <- 1000   # This requires a lot of time...

png(file="Data_check/Rarefaction_curve_PROKARYOTES_GENERA.png",width=2700,height=2100, res=300)
r<-rarecurve(t(as(otu_table( subset_taxa(data.genus,Domain%in%c("Bacteria","Archaea")) ) ,"matrix") ) ,
             step=steps,label=F, col="gray50",
             ylab = "Prokaryota Genera", xlab= "Related reads amount")
evalslopes(r,sample_data(data.genus)$Project_number,lim=0.001,cex=1.4)
dev.off()


png(file="Data_check/Rarefaction_curve_PROKARYOTES_SPECIES.png",width=2700,height=2100, res=300)
r<-rarecurve(t(as(otu_table( subset_taxa(data,Domain%in%c("Bacteria","Archaea")) ) ,"matrix") ) ,
             step=steps,label=F, col="gray50",
             ylab = "Prokaryota Species", xlab= "Related reads amount")
evalslopes(r,sample_data(data)$Project_number,lim=0.001,cex=1.4)
dev.off()



# WHAT IF AN ABUNDANCE FILTER IS APPLIED ...
png(file="Data_check/Rarefaction_curve_PROKARYOTES_SPECIES_AFTER_FILTERS.png",width=2700,height=2100, res=300)
data.prop<-transform_sample_counts(data, function(x) (x/sum(x))*100 )
who_over <- tax_table(data.prop)[ rowMeans( otu_table(data.prop) )>0.0001 , "Species"]
r<-rarecurve(t(as(otu_table( subset_taxa(data, Domain%in%c("Bacteria","Archaea") & Species%in%who_over ) ) ,"matrix") ) ,
             step=steps,label=F, ylab = "Prokaryota Species (filtered 0.0001 mean)", xlab= "Related reads amount")
evalslopes(r,sample_data(data)$Project_number,lim=0.001,cex=1.4)
dev.off()
write.csv2( file="Data_check/Prokaryote_species_with_mean_under_00001.csv",
            tax_table( subset_taxa(data, Domain%in%c("Bacteria","Archaea") & !Species%in%who_over )) [,"Domain_species"] ,
            quote=F , row.names = F )
rm(who_over)



steps <- 100   # these clades are much less populated ...

png(file="Data_check/Rarefaction_curve_EUKARYOTA_GENERA.png",width=2700,height=2100, res=300)
r<-rarecurve(t(as(otu_table( subset_taxa(data.genus,Domain=="Eukaryota") ) ,"matrix") ) ,
             step=steps,label=F, ylab = "Eukaryota Genera", xlab= "Related reads amount")
evalslopes(r,sample_data(data.genus)$Project_number,lim=0.001,cex=1.4)
dev.off()


# WHAT IF AN ABUNDANCE FILTER IS APPLIED ...
png(file="Data_check/Rarefaction_curve_EUKARYOTA_GENERA_AFTER_FILTERS.png",width=2700,height=2100, res=300)
data.genus.prop<-transform_sample_counts(data.genus, function(x) (x/sum(x))*100 )
who_over <- as.character( tax_table(data.genus)[ rowMeans( otu_table(data.genus.prop) )>0.0001 , "Genus"] )
r<-rarecurve(t(as(otu_table( subset_taxa(data.genus, Domain=="Eukaryota" & Genus%in%who_over ) ) ,"matrix") ) ,
             step=steps,label=F, ylab = "Eukaryota Genera\n(filtered 0.0001 mean)", xlab= "Related reads amount")
evalslopes(r,sample_data(data.genus)$Project_number,lim=0.001,cex=1.4)
dev.off()
write.csv2( file="Data_check/Eukaryota_genera_with_mean_under_00001.csv",
            cbind( otu_table( subset_taxa(data, Domain=="Eukaryota" & !Genus%in%who_over )),
                   tax_table( subset_taxa(data, Domain=="Eukaryota" & !Genus%in%who_over )) [,"Domain_species"] ),
            quote=F , row.names = F )
rm(who_over)

png(file="Data_check/Rarefaction_curve_VIRUSES_SPECIES.png",width=2700,height=2100, res=300)
r<-rarecurve(t( as(otu_table( subset_taxa(data,Domain=="Viruses") ) ,"matrix") ) ,
             step=steps,label=F, ylab = "Virus Species", xlab= "Related reads amount")
evalslopes(r,sample_data(data)$Project_number,lim=0.001,cex=1.4)
dev.off()




###################### AVERAGE READS NUMBER #########################

options(scipen=100)
table <- read.table(file="FASTQ_check/original_number_of_seqs.txt", sep=":")
table <- table[ grepl("_1",table$V1 )|grepl("_R1_",table$V1) , ]   # only forwards, then doubling the result afterwards
# table$V1 <- gsub("_1.fastq ","",table$V1, fixed = T)
table$V1 <- gsub("_.*","",table$V1)
table$Project <- Metadata[table$V1, "Project_number"]
table <- table[!is.na(table$Project), ]
table <- aggregate( V2~Project , table, mean )
colnames(table)<- c("Project","Average_read_number")
table$Average_read_number <- round(table$Average_read_number, digits = -7)
table$Average_in_millions<-gsub("000000","", table$Average_read_number)
write.csv(table, file="Data_check/Average_Read_for_project.csv", row.names = F, quote = F)
rm(table)




##################### R AND PACKAGES VERSION #########################

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
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))

