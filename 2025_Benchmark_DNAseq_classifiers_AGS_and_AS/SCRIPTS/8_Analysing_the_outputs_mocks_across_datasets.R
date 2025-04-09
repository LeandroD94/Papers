# R 4.3.2


########################### PREPARING THE ENVIRONMENT ########################### 

dir.create("Results_Mock")
dir.create("Results_Mock/Venn_diagrams")
options(scipen=100)

library("vegan")
library("ecodist")
library("ggplot2")
library("reshape")
library("ggvenn")

fill_color_20<-c("blue", "deeppink", "chocolate3", "darkmagenta","bisque2","cyan","brown","yellow2","firebrick3","springgreen4","violet","darkslategray3","wheat3", "gray", "pink3","yellow4","red","darkgreen","lightgrey","coral","grey28") 
fill_color_30<-c("wheat3","deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "aquamarine3",
                 "bisque2","cyan","yellow3","lightcoral",
                 "springgreen4", "firebrick3",
                 "grey15","deepskyblue2","lightgreen","orange2",
                 "darkorchid", "lightblue1" ,"blue1", "violet", 
                 "yellow", "pink3","yellow4", "chocolate4","coral","darkgreen",
                 "lightgrey", "red","darkslategray3")
fill_color_30_species<-c("wheat3","deeppink", "darkcyan","darkmagenta",
                 "darkblue","grey50", "aquamarine3",
                 "bisque2","cyan","yellow3","lightcoral",
                 "springgreen4", "firebrick3",
                 "deepskyblue2","lightgreen","orange2",
                 "darkorchid", "lightblue1" ,"blue1", "violet", 
                 "yellow", "pink3","yellow4", "chocolate4","coral","darkgreen",
                 "lightgrey", "red","darkslategray3", "grey28")




#### custom function to import the tables from the classifiers in tables with the same format (rows=genera, columns=samples)

Import_from_Tot_Seq<- function( input_folder , table_suffix , column_tax , column_abund ) {
  
  # input_folder= folder in which search the objects (include the / in the folder path)
  # table_suffix= suffix which the files (tsv) to import have in common
  # column_tax= the name of the taxonomy column in such files
  # column_abund= abundance column
  
  files_to_import<-dir(input_folder)
  files_to_import<- files_to_import[grep(table_suffix, files_to_import, fixed=T)]
  files_list<-list()
  for( i in 1:length(files_to_import) ) {
    if( grepl("kMetaShot", files_to_import[i]) ){ # kMetaShot output requires a re-formatting to be imported by this output
      raw_table<- as.data.frame(read.table( paste0(input_folder , files_to_import[i]) , sep=",", header = T, dec = ".", quote="") )
      # NB: the a2r table is already at genus level ... but some genus names are cut if they feature a space (e.g. "Candidatus Accumulibacter")
      raw_table <- raw_table[ raw_table$organism_name!="-" , ]  # removing the unclassified
      if( grepl("MAGs", files_to_import[i]) ){  # almost every contigs have confidences too low to be subsetted --> only MAGs
        raw_table <- raw_table[ raw_table$ass2ref!=0 , ]  # if zero then they are the obs below the set conf score
      }
      raw_table$genus <- gsub( "Candidatus ","Candidatus_",raw_table$genus )
      #raw_table[ grepl("Candidatus", raw_table$organism_name) , "organism_name"] <- paste( "Candidatus", raw_table[grepl("Candidatus", raw_table$organism_name), "genus"]  , sep="_")  # 
      raw_table[ , "organism_name"] <- gsub(" .*","", raw_table[ , "organism_name"]) # solving eventual spaces (species level name --> genus level)
      raw_table[ raw_table$organism_name=="Candidatus" , "organism_name"] <- paste( raw_table[raw_table$organism_name=="Candidatus","organism_name"] , raw_table[raw_table$organism_name=="Candidatus","genus"] , sep="_")  # there are few entries which are just called "Candidatus" but have a known tax ID
      raw_table<-raw_table[ , "organism_name", drop=F]
      files_list[[i]] <- cbind.data.frame( "organism_name"= names(table(raw_table)) , "abund"=as.numeric(table(raw_table)) ) # counting each occurrence
    } else {
      # the output of the other programs are already formatted as required by this part of the function
      files_list[[i]]<- as.data.frame(read.table( paste0(input_folder , files_to_import[i]) , sep="\t", header = T, dec = ".", quote="") )
    }
  }
  names(files_list)<-gsub(table_suffix, "", files_to_import)
  # NB: uncleaned file names can cause problem with the unclassified count (due to the row name matching used afterward)
  # joining every report to a temporary unique table for the script (to search for unique bacteria in every assignment and then creating the related row)
  complete_table<-NULL
  for(i in 1:length(files_list)){
    new_sample<-cbind.data.frame(sample=names(files_list)[i], files_list[[i]]) # sample name (repeated along the column) with the related results
    # kraken2 core_nt gives few duplicated entries sometime! Solving it by adding the tax id...
    if( grepl("bracken",table_suffix) ) {
      duplicated_entries <- duplicated(new_sample$name)
      new_sample$name[ duplicated_entries ] <- paste0(new_sample$name[ duplicated_entries ] , "_taxID", new_sample$taxonomy_id[ duplicated_entries ] )
    }
    complete_table<-rbind.data.frame(complete_table, new_sample)
  }
  
  colnames(complete_table)[colnames(complete_table)== column_tax ]<-"Bacterium"
  colnames(complete_table)[colnames(complete_table)==column_abund ]<-"Abundance"  # selecting the estimation by bracken as abundance
  
  
  # from list (each sample is repeated along a column) to feature table (each sample is a column)
  feature_table<-cbind("Temp_column"=unique(complete_table$Bacterium))
  # each bacterium in a row --> every bacterium in this table --> searching across samples --> if absent then 0
  for( s in 1:length(unique(complete_table$sample) )) {  # NB: the numeric index is required to deal with the column name afterward!
    sample<- unique(complete_table$sample)[s] # 1 sample at time
    table_4_that_sample<- complete_table[complete_table$sample==sample , ]
    abundances<-NULL # resets the vector, ready for the next bacterium
    for( b in unique(complete_table$Bacterium)) {
      # scan each bacterium in order: if found in that sample then appends the correspective abundances, if absent paste a 0
      if( b %in% table_4_that_sample$Bacterium ){
        abundances<-c(abundances, table_4_that_sample[ table_4_that_sample$Bacterium==b, "Abundance"])
      } else {
        abundances<-c(abundances, 0 ) # that bacterium is absent --> zero
      }
    }
    feature_table<-cbind.data.frame(feature_table, abundances)
    colnames(feature_table)[s+1] <- sample   # the +1 is for first column, which is temporarily the taxonomy column
  }
  colnames(feature_table)[ colnames(feature_table)=="Temp_column" ] <- "Genus"
  table_obtained <<- feature_table
  print("The imported counts are now stored in the object 'table_obtained' available in your main R scope\n")
}




############## IMPORTING MOCKS ABUNDANCES ACCORDING TO THE VARIOUS PROGRAMS ###################### 


### importing the mock abundances
table_mock <- read.delim(file="Obtained_Counts_mock/true_abudances_of_DNAseq_mock.tsv", sep=" ", col.names= c("Genus", "Known_mock") , header= F )
table_mock$Known_mock <-apply( table_mock [ ,  "Known_mock" , drop=F ] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
table_mock$Known_mock <- round( table_mock [ , "Known_mock" , drop=F ] , 4 )


### import from RiboFrame (full counts)
Import_from_Tot_Seq(input_folder="Obtained_Counts_mock/", table_suffix="_full_count.genus.cnt", column_tax="Name", column_abund="Perc" )
Counts_RiboF_full<-table_obtained
row.names(Counts_RiboF_full)<-Counts_RiboF_full$Genus


### import from RiboFrame (V3-V4)
Import_from_Tot_Seq(input_folder="Obtained_Counts_mock/", table_suffix="_V3_V4.genus.cnt", column_tax="Name", column_abund="Perc" )
Counts_RiboF_V3V4<-table_obtained
row.names(Counts_RiboF_V3V4)<-Counts_RiboF_V3V4$Genus


### import from Kraken2 (BOTH K-SILVA and K-NT_core)
Import_from_Tot_Seq(input_folder="Obtained_Counts_mock/", table_suffix="bracken_abund_Genus", column_tax="name", column_abund="new_est_reads" )
Counts_Kraken<-table_obtained
row.names(Counts_Kraken)<-Counts_Kraken$Genus


### importing MEGAHIT-kMetaShot tables and parsing them manually
Import_from_Tot_Seq(input_folder="Obtained_Counts_mock/", table_suffix="classif_contigs_", column_tax="organism_name", column_abund="abund" )
Counts_Contigs<-table_obtained
# head(Counts_Contigs, n=2)
row.names(Counts_Contigs)<-Counts_Contigs$Genus
Counts_Contigs <- Counts_Contigs[ , !grepl("conf0.4", colnames(Counts_Contigs), fixed=T ) ]  # conf threshold set to 0 (see import function) --> 
colnames(Counts_Contigs)<- gsub("_conf0.2","", colnames(Counts_Contigs))  # not used here


### importing MEGAHIT-MetaBAT2-kMetashot tables and parsing them manually
Import_from_Tot_Seq(input_folder="Obtained_Counts_mock/", table_suffix="classif_MAGs_", column_tax="organism_name", column_abund="abund" )
Counts_MAGs<-table_obtained
# head(Counts_MAGs, n=2)
row.names(Counts_MAGs)<-Counts_MAGs$Genus


### importing Kaiju output and parsing it manually
Kaiju_list <- read.delim(file= "Obtained_Counts_mock/kaiju_classif_summary.tsv", header=T)
Kaiju_list <- Kaiju_list[ , c("file","reads","taxon_name") ]
colnames(Kaiju_list)<- c("sample","Abundance","Bacterium")
# NB: few genera are repeated (those which derive from the species level, see Kaiju script), hence solving them is required for the code below
Kaiju_list<-aggregate(Abundance~sample+Bacterium, Kaiju_list, FUN=sum)
# from list (each sample is repeated along a column) to feature table (each sample is a column)
feature_table<-cbind("Temp_column"=unique(Kaiju_list$Bacterium))
# each bacterium in a row --> every bacterium in this table --> searching across samples --> if absent then 0
for( s in 1:length(unique(Kaiju_list$sample) )) {  # NB: the numeric index is required to deal with the column name afterward!
  sample<- unique(Kaiju_list$sample)[s] # 1 sample at time
  table_4_that_sample<- Kaiju_list[Kaiju_list$sample==sample , ]
  abundances<-NULL # resets the vector, ready for the next bacterium
  for( b in unique(Kaiju_list$Bacterium)) {
    # scan each bacterium in order: if found in that sample then appends the correspective abundances, if absent paste a 0
    if( b %in% table_4_that_sample$Bacterium ){
      abundances<-c(abundances, table_4_that_sample[ table_4_that_sample$Bacterium==b, "Abundance"])
    } else {
      abundances<-c(abundances, 0 ) # that bacterium is absent --> zero
    }
  }
  feature_table<-cbind.data.frame(feature_table, abundances)
  colnames(feature_table)[s+1] <- sample   # the +1 is for first column, which is temporarily the taxonomy column
}
colnames(feature_table)[ colnames(feature_table)=="Temp_column" ] <- "Genus"
Counts_Kaiju <- feature_table
row.names(Counts_Kaiju) <- Counts_Kaiju$Genus
# mainteining only the mock samples (Kaiju outputs are already half parsed in its script)
Counts_Kaiju<- Counts_Kaiju [ , grepl("Genus",colnames(Counts_Kaiju))|grepl("mock",colnames(Counts_Kaiju) )]




################# SEARCHING AND SOLVING TAXONOMIES INCOMPABILITIES ################# 


### translating names between SILVA and NCBI (at least the most abundants)
top_Kraken <- names( sort (rowSums( Counts_Kraken[, colnames(Counts_Kraken)!="Genus" & !grepl("SILVA",colnames(Counts_Kraken))] ), decreasing=TRUE ) ) [1:30]
top_Kraken <-gsub(" ", "_", top_Kraken)
top_Kaiju <- names( sort (rowSums( Counts_Kaiju[, colnames(Counts_Kaiju)!="Genus"] ), decreasing=TRUE ) ) [1:30]
top_Kaiju <-gsub(" ", "_", top_Kaiju)
Kraken_only_SILVA <- Counts_Kraken[, colnames(Counts_Kraken)=="Genus" | grepl("SILVA",colnames(Counts_Kraken))] 
Kraken_only_SILVA <- Kraken_only_SILVA[ rowSums(Kraken_only_SILVA[,colnames(Kraken_only_SILVA)!="Genus"])>0 , ]

unique_names_Kraken <- top_Kraken [!top_Kraken %in% Kraken_only_SILVA$Genus ]
unique_names_Kaiju <- top_Kaiju[!top_Kaiju %in% Kraken_only_SILVA$Genus  & !top_Kaiju%in%unique_names_Kraken ]
# unique_names_Kraken
# all of the Kraken unique names exists in the top Kraken-SILVA, aside from the NCBI genus Circaeaster does not have a genus level in SILVA (it is "incertae sedis" in SILVA)...

# unique_names_Kaiju
# checking these names, I do not see any synonymies to modify


### solving the different names where possible (through gsub), or adding from the species level output in the previous script...
for (x in c( "Counts_Kraken" , "Counts_Kaiju")) {
  that_table <- get(x)   # x is still a character, getting the corresponding object
  ### solving the different names where possible (through gsub), or adding from the species level output in the previous script...
  that_table$Genus <- gsub("Nostocoides","Tetrasphaera", that_table$Genus ) 
  that_table$Genus <- gsub("Candidatus_Moraniibacteriota bacterium","Candidatus_Moranbacteria", that_table$Genus )
  that_table$Genus <- gsub("Candidatus_Moranbacteria.*","Candidatus_Moranbacteria", that_table$Genus )
  that_table$Genus <- gsub("Candidatus_Kaiserbacteria bacterium","Candidatus_Kaiserbacteria", that_table$Genus )
  that_table$Genus <- gsub("Candidatus_Nomurabacteria bacterium","Candidatus_Nomurabacteria", that_table$Genus )
  that_table$Genus <- gsub("Candidatus_Magasanikbacteria bacterium","Candidatus_Magasanikbacteria", that_table$Genus )
  # the SILVA genus SH-PL14 has an NCBI taxID 1632864 according to MIDAS, which corresponds to the "missing" genus of the species "Planctomyces sp. SH-PL14"
  that_table$Genus <- gsub("Planctomyces sp. SH-PL14","SH-PL14", that_table$Genus )
  that_table$Genus <- gsub("Planctomyces_sp._SH-PL14","SH-PL14", that_table$Genus )
  # the SILVA genus IMCC26207 has the NCBI code "1641811" according to MIDAS --> corresponds to the "missing" genus of the species Actinobacteria bacterium IMCC26207
  # the SILVA genus 67-14 corresponds to the "missing" genus of the species Solirubrobacterales bacterium 67-14 on NBCI
  that_table$Genus <- gsub("Solirubrobacterales bacterium","67-14", that_table$Genus )
  that_table$Genus <- gsub("Solirubrobacterales_bacterium","67-14", that_table$Genus )
  # the SILVA genus 37-13 corresponds to the "missing" genus of the species Bacteroidetes bacterium 37-13
  that_table$Genus <- gsub("Bacteroidetes_bacterium_37-13","37-13", that_table$Genus)
  that_table$Genus <- gsub("Bacteroidetes bacterium 37-13","37-13", that_table$Genus)
  # Bdellovibrio, Lentimicrobium and Prosthecobacter are actually also in NCBI, they are just not featured in the Kraken dataset ...
  # but Bdellovibrio and Prosthecobacter are featured in Kaiju dataset !
  # "Pir4 lineage" (of Pirellulaceae) is feautured in both SILVA and MIDAS but it is not in NCBI at all!
  # "SM1A02" does not exist on NCBI either...
  # codes such as "1-20"  are not found in NCBI...
  that_table$Genus <- gsub("Gracilibacteria","Altimarinota", that_table$Genus )
  that_table$Genus <- gsub("Nostocoides","Tetrasphaera", that_table$Genus )
  that_table$Genus <- gsub("Tequatrovirus","phage_T4", that_table$Genus )
  that_table <- aggregate ( .~Genus , that_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Genus level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 
  row.names(that_table)<-that_table$Genus
  assign( x, that_table ) # overwrites the original object
  rm(x,that_table)
}


# RiboFrame does not classify 67-14 as such but as "Incertae_Sedis_f_67_14"
Counts_RiboF_full$Genus <-gsub("Incertae_Sedis_f_67_14","67-14", Counts_RiboF_full$Genus)
Counts_RiboF_V3V4$Genus <-gsub("Incertae_Sedis_f_67_14","67-14", Counts_RiboF_V3V4$Genus)
row.names(Counts_RiboF_full)<-Counts_RiboF_full$Genus
row.names(Counts_RiboF_V3V4)<-Counts_RiboF_V3V4$Genus

# Counts_Kraken$Genus[grepl("67.14",Counts_Kraken$Genus)]

# uncultured observations (as unique genus, even though from different families) only in SILVA ...
# Counts_Kraken[grep("uncultured",Counts_Kraken$Genus) ,]
# Counts_Kaiju$Genus[grep("uncultured",Counts_Kaiju$Genus)]
Counts_Kraken <- Counts_Kraken[!grepl("uncultured",Counts_Kraken$Genus) ,]  # removing these SILVA unique lines
Counts_RiboF_full <- Counts_RiboF_full[ !grepl("uncultured",Counts_RiboF_full$Genus) ,]  # removing these SILVA unique lines
Counts_RiboF_V3V4 <- Counts_RiboF_V3V4[ ! grepl("uncultured",Counts_RiboF_V3V4$Genus) ,]  # removing these SILVA unique lines


### Solving the missing names of the Candidatus genera in kMetaShot contig's a2r output (at least the most abundant) and other small issues ...
for (x in c( "Counts_Contigs" , "Counts_MAGs")) {
  that_table <- get(x)
  that_table<- that_table[ order(rowMeans(that_table[,colnames(that_table)!="Genus"]),decreasing = T) , ]
  # that_table[grepl("Candidatus",that_table$Genus) , ][1 :15 , ]
  that_table$Genus <- gsub("Candidatus_327159","Candidatus_Accumulibacter",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_221279","Candidatus_Competibacter",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_1579476","Candidatus_Chloroploca",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_622681","Candidatus_Planktophila",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_198251","Candidatus_Pelagibacter",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_1400859","Candidatus_Contendobacter",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_84565","Candidatus_Sodalis",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_453161","Candidatus_Nitrotoga",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_33055","Candidatus_Kinetoplastibacterium",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_273135","Candidatus_Cardinium",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_33926","Candidatus_Phytoplasma",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_34019","Candidatus_Liberibacter",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_49082","Candidatus_Arthromitus",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_1679002","Candidatus_Methylopumilus",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_41949","Candidatus_Microthrix",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_1470349","Candidatus_Stoquefichus",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_1826864","Candidatus_Nitrosocosmicus",that_table$Genus)
  that_table$Genus <- gsub("Candidatus_1752737","Candidatus_Moranbacteria",that_table$Genus) # it is not featured ... but just in case !
  that_table$Genus <- gsub("Candidatus_2045217","Candidatus_Moranbacteria",that_table$Genus) # it is not featured ... but just in case !
  that_table$Genus <- gsub("Candidatus_83766","Propionivibrio",that_table$Genus) # NB: misclassification of MAGs approach (NB: no, it is not a Candidatus according to NCBI)
  # now every other crucial taxonomy issue will be solved
  that_table$Genus <- gsub("Solirubrobacter","67-14",that_table$Genus)
  that_table$Genus <- gsub("Gracilibacteria","Altimarinota", that_table$Genus )
  that_table$Genus <- gsub("Nostocoides","Tetrasphaera", that_table$Genus )
  that_table$Genus <- gsub("Tequatrovirus.*","phage_T4", that_table$Genus )
  that_table$Genus <- gsub(" .*","", that_table$Genus) # solving eventual spaces (species level --> genus level)
  that_table$Genus <- gsub("[","", that_table$Genus, fixed = T) # solving eventual spaces (species level --> genus level)
  that_table$Genus <- gsub("]","", that_table$Genus, fixed = T) # solving eventual spaces (species level --> genus level)
  that_table$Genus <- gsub("cf.","Phormidesmis", that_table$Genus, fixed = T) # solving eventual spaces (species level --> genus level)
  # Checked manually: the observation called just "Candidatus" has the same tax ID as Accumulibacter ...
  that_table$Genus <- gsub("Candidatus$","Candidatus_Accumulibacter", that_table$Genus) # solving eventual spaces (species level --> genus level)
  that_table <- aggregate ( .~Genus , that_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Genus level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 
  row.names( that_table )<- that_table$Genus
  assign( x, that_table ) # overwrites the original object
  rm(x,that_table)
}



### Working on the mock names too
table_mock$Genus<-gsub("Solirubrobacterales_bacterium","67-14", table_mock$Genus)
# table_mock$Genus




###################### WORKING ON VALUES: ROUNDING AND FILTERING ###################### 

# if (! "proof1" %in% ls() ){
#   unfilt_Counts_RiboF_V3V4<-Counts_RiboF_V3V4
#   unfilt_Counts_RiboF_full<-Counts_RiboF_full
#   unfilt_Counts_Kraken<-Counts_Kraken
#   unfilt_Counts_Kaiju<-Counts_Kaiju
#   unfilt_Counts_Contigs<-Counts_Contigs
#   unfilt_Counts_MAGs<-Counts_MAGs
# }


# From raw abund to percentages
Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"]<-apply( Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"]<-apply( Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"]<-apply( Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_MAGs[ , !colnames(Counts_MAGs)%in%"Genus"]<-apply( Counts_MAGs[ , !colnames(Counts_MAGs)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
# RiboFrame counts are already in percentages...


# ### FILTERING LOW VALUES (while they are separated datasets) ...
# filter_threshold <- 0.001       # still performed, only for RiboFrame
# # {
# # proof1<-"The samples were already filtered"
# # 
# # Counts_RiboF_V3V4 <- Counts_RiboF_V3V4[ rowMeans( Counts_RiboF_V3V4[ , !colnames(Counts_RiboF_V3V4)%in%"Genus"] ) > filter_threshold , ] 
# # Counts_RiboF_full <- Counts_RiboF_full[ rowMeans( Counts_RiboF_full[ , !colnames(Counts_RiboF_full)%in%"Genus"] ) > filter_threshold , ] 
# # 
# # # in Kraken, the mock was analysed by both SILVA and nt core database --> two dataset in one! The filtering will be performed on two different row means.
# # Counts_Kraken[ rowMeans( Counts_Kraken[, !colnames(Counts_Kraken)=="Genus" & !grepl("SILVA",colnames(Counts_Kraken))] ) < filter_threshold , !colnames(Counts_Kraken)=="Genus" & !grepl("SILVA",colnames(Counts_Kraken)) ] <- 0
# # Counts_Kraken[ rowMeans( Counts_Kraken[, !colnames(Counts_Kraken)=="Genus" & grepl("SILVA",colnames(Counts_Kraken))] ) < filter_threshold , !colnames(Counts_Kraken)=="Genus" & grepl("SILVA",colnames(Counts_Kraken)) ] <- 0
# # Counts_Kraken <- Counts_Kraken[ rowSums(Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"] ) > 0 , ]   # every row with only zeros is discarded
# # 
# # same in RiboF ...
# Counts_RiboF_V3V4[ rowMeans(Counts_RiboF_V3V4[,!colnames(Counts_RiboF_V3V4)%in%c("Genus")]) < filter_threshold , !colnames(Counts_RiboF_V3V4)%in%c("Genus") ] <- 0
# Counts_RiboF_V3V4 <- Counts_RiboF_V3V4[ rowSums(Counts_RiboF_V3V4[ , !colnames(Counts_RiboF_V3V4)%in%"Genus"] ) > 0 , ]   # every row with only zeros is discarded
# Counts_RiboF_full[ rowMeans(Counts_RiboF_full[,!colnames(Counts_RiboF_full)%in%c("Genus")]) < filter_threshold , !colnames(Counts_RiboF_full)%in%c("Genus") ] <- 0
# Counts_RiboF_full <- Counts_RiboF_full[ rowSums(Counts_RiboF_full[ , !colnames(Counts_RiboF_full)%in%"Genus"] ) > 0 , ]   # every row with only zeros is discarded
# # # same also in Kaiju ...
# # Counts_Kaiju[ rowMeans(Counts_Kaiju[,!colnames(Counts_Kaiju)%in%c("Genus")]) < filter_threshold , !colnames(Counts_Kaiju)%in%c("Genus") ] <- 0
# # Counts_Kaiju <- Counts_Kaiju[ rowSums(Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"] ) > 0 , ]   # every row with only zeros is discarded
# # # same also in kMetaShot
# # Counts_Contigs[ rowMeans(Counts_Contigs[,!colnames(Counts_Contigs)%in%c("Genus")]) < filter_threshold , !colnames(Counts_Contigs)%in%c("Genus") ] <- 0
# # Counts_Contigs <- Counts_Contigs[ rowSums(Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"] ) > 0 , ]   # every row with only zeros is discarded
# # 
# # }
# # write.csv2(unfilt_Counts_Kraken[ !unfilt_Counts_Kraken$Genus %in% Counts_Kraken$Genus, ], file = "Results_Mock/Noises_Kraken_filteret_out_abund_threshold.tsv", row.names = F, quote = F)
# # write.csv2(unfilt_Counts_Kaiju[ !unfilt_Counts_Kaiju$Genus %in% Counts_Kaiju$Genus, ], file = "Results_Mock/Noises_Kaiju_filteret_out_abund_threshold.tsv", row.names = F, quote = F)
# 
# 
# # Filtering (while they are separated datasets) ...
# # REMOVING SINGLETONS FROM THE OBJECTS, NB: the unfiltered object is used (--> not transformed and reproducible)
# for (x in c( "Counts_Kraken" , "Counts_Kaiju" , "Counts_Contigs")) {   # NB: Counts_MAGs is not included in the filter due to the excessive low counts !
#   that_table <- get(x)   # x is still a character, getting the corresponding object
#   unfilt_table <- get(paste0("unfilt_",x))
#   sing_doublet <-unfilt_table[ apply(unfilt_table[,colnames(unfilt_table)!="Genus"], MARGIN=1, function(x) max(x)<=1), "Genus"]
#   that_table <-that_table[!that_table$Genus%in% sing_doublet, ]
#   assign( x, that_table ) # overwrites the original object
#   write.csv2(unfilt_table[ unfilt_table$Genus %in% sing_doublet, ], file = paste0("Results_Mock/",x,"filtered_being_singletons.csv"), row.names = F, quote = F)
#   rm(x,that_table)
#   }
# NB: RIBOFRAME IS NOT INCLUDED HERE AS ITS OUTPUT IS ALREADY IN PERCENTAGES -->  USING THE FILTER ABOVE
# NB: Neither Counts_MAGs is included in the filter due to the excessive low counts !



# the UNCLASSIFIED reads are removed from Kaiju (being the only one with this row between its mock counts)
Counts_Kaiju <- Counts_Kaiju[ ! grepl("unclassified", Counts_Kaiju$Genus) , ]
# the "cannot_be_assigned_to_a_(non-viral)_genus" are reads which have been classified, but not at genus level and/or into a non-viral clade... but classified! See https://github.com/bioinformatics-centre/kaiju/issues/138


# Re-computing the percentages after the filters
Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"]<-apply( Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"]<-apply( Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"]<-apply( Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_MAGs[ , !colnames(Counts_MAGs)%in%"Genus"]<-apply( Counts_MAGs[ , !colnames(Counts_MAGs)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_RiboF_V3V4[ , !colnames(Counts_RiboF_V3V4)%in%"Genus"]<-apply( Counts_RiboF_V3V4[ , !colnames(Counts_RiboF_V3V4)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
Counts_RiboF_full[ , !colnames(Counts_RiboF_full)%in%"Genus"]<-apply( Counts_RiboF_full[ , !colnames(Counts_RiboF_full)%in%"Genus"] , MARGIN = 2, FUN= function(x) (x/sum(x))*100 )
# Rounding
Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"]<- round( Counts_Kraken[ , !colnames(Counts_Kraken)%in%"Genus"] , 4 )
Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"]<- round( Counts_Kaiju[ , !colnames(Counts_Kaiju)%in%"Genus"] , 4 )
Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"]<- round( Counts_Contigs[ , !colnames(Counts_Contigs)%in%"Genus"] , 4 )
Counts_MAGs[ , !colnames(Counts_MAGs)%in%"Genus"]<- round( Counts_MAGs[ , !colnames(Counts_MAGs)%in%"Genus"] , 4 )
Counts_RiboF_V3V4[ , !colnames(Counts_RiboF_V3V4)%in%"Genus"]<- round( Counts_RiboF_V3V4[ , !colnames(Counts_RiboF_V3V4)%in%"Genus"] , 4 )
Counts_RiboF_full[ , !colnames(Counts_RiboF_full)%in%"Genus"]<- round( Counts_RiboF_full[ , !colnames(Counts_RiboF_full)%in%"Genus"] , 4 )




##################### PREPARING THE FINAL OBJECTS ##################### 

### adding prefix to each colname
colnames(Counts_Kaiju)[colnames(Counts_Kaiju)!="Genus"] <- gsub("mock_AGS_DNAseq_Evalue","Kaiju_E", colnames(Counts_Kaiju)[colnames(Counts_Kaiju)!="Genus"])
colnames(Counts_Kaiju)[colnames(Counts_Kaiju)!="Genus"] <- gsub("mock_AGS_DNAseq_contigs_Evalue","Kaiju_contigs_E", colnames(Counts_Kaiju)[colnames(Counts_Kaiju)!="Genus"])
colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"] <- gsub("mock_AGS_DNAseq_conf","Kraken_c", colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"])
colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"] <- gsub("_.tsv","", colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"])
colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"] <- gsub(".tsv","", colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"])
colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"] <- gsub("mock_AGS_DNAseq_contigs_conf", "Kraken_contigs_c", colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"])
colnames(Counts_RiboF_full)[colnames(Counts_RiboF_full)!="Genus"] <- gsub("mock_AGS_DNAseq_conf","RiboF_full_c", colnames(Counts_RiboF_full)[colnames(Counts_RiboF_full)!="Genus"])
colnames(Counts_RiboF_V3V4)[colnames(Counts_RiboF_V3V4)!="Genus"] <- gsub("mock_AGS_DNAseq_conf","RiboF_V3V4_c", colnames(Counts_RiboF_V3V4)[colnames(Counts_RiboF_V3V4)!="Genus"])


### parsing between programs (while adding the known "true" abundances of the mock)
files_list<- list( Counts_Kraken, Counts_RiboF_full, Counts_RiboF_V3V4, Counts_Kaiju , Counts_Contigs, Counts_MAGs, table_mock )
names(files_list) <-  c("Counts_Kraken", "Counts_RiboF_full", "Counts_RiboF_V3V4", "Counts_Kaiju" , "Counts_Contigs", "Counts_MAGs", "table_mock" )
every_entry<-NULL
for (x in 1:length(files_list)) {
  every_entry <- c(every_entry, files_list[[x]][["Genus"]])
}
# from list (each set_pipeline is repeated along a column) to feature table (each set_pipeline is a column)
feature_table<-cbind("Temp_column"=unique(every_entry))
# each bacterium in a row --> every bacterium in this table --> searching across datasets
for( s in 1:length(files_list) ) {
  set_pipeline<- names(files_list)[[s]]  # 1 set_pipeline at time
  table_4_that_set_pipeline<- files_list[[s]]
  abundances<-NULL  # resets the object, ready for the next bacterium
  # an empty column to be filled with each sample...
  temp_set_pipeline<- matrix( "temp", length(unique(every_entry)) , length(colnames(table_4_that_set_pipeline)) )
  temp_set_pipeline <- as.data.frame( temp_set_pipeline )
  temp_set_pipeline[ ,1] <- unique(every_entry)
  for( b in unique(every_entry)) {
    # scan each bacterium in order: if found in that set_pipeline then appends the correspective abundances, if absent paste a 0
    if( b %in% table_4_that_set_pipeline$Genus ){
      abundances<- as.numeric(table_4_that_set_pipeline[ table_4_that_set_pipeline$Genus==b, !colnames(table_4_that_set_pipeline)%in%"Genus"  ])   # NB: the entire ROW will be captured this time
    } else {
      abundances<- rep(0, length( colnames(table_4_that_set_pipeline)[!colnames(table_4_that_set_pipeline)%in%"Genus" ]) )  # that bacterium is absent --> a series of zero
    }
    temp_set_pipeline[ temp_set_pipeline$V1 == b , ] <- c(b , abundances)  # overwrite the "temp"s row with actual values
  }
  colnames(temp_set_pipeline)<- colnames(table_4_that_set_pipeline)
  feature_table<-cbind.data.frame(feature_table, temp_set_pipeline)  # this adds the sample to the complete table
}
feature_table <- feature_table[ , ! colnames(feature_table) %in% "Genus" ]
beasts <- feature_table$Temp_column  # this column correspond to another "Genus" column
# NB: "apply" will delete the colnames, then the genera names will be assigned as row.names afterwards
feature_table$Temp_column <- NULL
feature_table<-as.data.frame(apply(feature_table, MARGIN=2, FUN= as.numeric))
row.names(feature_table) <- beasts

# write.table(feature_table, file="Results_Mock/Counts_MOCK_from_every_pipeline.tsv", sep = "\t", quote=F, row.names = T, col.names = T)


### preparing the metadata
metadata <- cbind.data.frame( "Workflow"=colnames(feature_table)[ colnames(feature_table)!="Genus"] )
metadata$Program <- gsub("_.*","", metadata$Workflow )
metadata$Program <- gsub("Known","Mock_true_abund", metadata$Program )
metadata$Program_sub <- metadata$Program
metadata$Program_sub[grepl("_plus",metadata$Workflow)] <- "Kaiju_nr_euk+"
#metadata$Program_sub[grepl("Kaiju_contigs_",metadata$Workflow)] <- "Kaiju_nr_euk+"
metadata$Program_sub[grepl("_SILVA",metadata$Workflow)] <- "Kraken_SILVA"
metadata$Program_sub[grepl("_full_",metadata$Workflow)] <- "RiboF_full"
metadata$Program_sub[grepl("_V3V4_",metadata$Workflow)] <- "RiboF_V3V4"
metadata$Program_sub[grepl("kMetaShot_MAGs",metadata$Workflow)] <- "kMetaShot_MAGs"
metadata$Program_sub[grepl("Kaiju_contigs",metadata$Workflow)] <- "Kaiju_contigs"
metadata$Program_sub[grepl("Kraken_contigs",metadata$Workflow)] <- "Kraken_contigs"
metadata$Program_sub <- gsub( "kMetaShot$","kMetaShot_contigs", metadata$Program_sub, perl = T)
metadata$Settings <- gsub(".*_[c,E]","",metadata$Workflow )
metadata$Settings <- ifelse( grepl("Kaiju",metadata$Program) , paste0("E", metadata$Settings) , paste0("c", metadata$Settings) )
metadata$Settings[grepl("Kaiju_contigs_",metadata$Workflow) ] <- paste0(metadata$Settings[grepl("Kaiju_contigs_",metadata$Workflow) ],"_plus")
metadata$Settings <- gsub("_plus","_nr_euk+",metadata$Settings )
metadata$Settings <- gsub("conf","c",metadata$Settings )
metadata$Settings <- gsub("cKnown_mock","Mock",metadata$Settings )
metadata$Settings <- gsub("EUKA","E0.001_m22_NOeuk",metadata$Settings )
metadata$Settings[grepl("default_",metadata$Workflow)] <- paste0("default_", metadata$Settings[grepl("default_",metadata$Workflow)]) 
metadata$Settings[grepl("metalarge_",metadata$Workflow)] <- paste0("large_", metadata$Settings[grepl("metalarge_",metadata$Workflow)]) 
metadata$Settings[grepl("custom_",metadata$Workflow)] <- paste0("custom_", metadata$Settings[grepl("custom_",metadata$Workflow)]) 
metadata$Settings <- gsub("custom$","custom_no_conf", metadata$Settings)
metadata$Settings <- gsub("ckMetaShot_default","default_no_conf", metadata$Settings)
metadata$Settings <- gsub("ckMetaShot_metalarge","metalarge_no_conf", metadata$Settings)


### cleaning the environment before saving the data
to_maintain <- ls()[ grepl("Counts_", ls() ) | grepl("unfilt_", ls() ) ]
to_maintain <- c(to_maintain, "proof1", "metadata", "feature_table", "table_mock", "fill_color_20", "fill_color_30", "fill_color_30_species")
rm(list= ls()[!ls()%in%to_maintain] )

# save.image( "Results_Mock/data_prepared_for_analysis.RData" )
# write.csv2(metadata, file="Metadata_workflows.csv", quote=F, row.names = F)
write.csv2(feature_table, file="Results_Mock/Abundances_every_classif.csv", quote=F, row.names = T)




##################### ABUNDANCES (TRUE BACTERIA) ########################

true_genera<-table_mock$Genus
judged_table<-feature_table
judged_table$Genus<-row.names(judged_table)
judged_table[! judged_table$Genus %in%table_mock$Genus , "Genus"] <- "misclassified" 
judged_table <- aggregate ( .~Genus , judged_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Genus level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 
judged_table$Genus <- factor(judged_table$Genus , levels = c(table_mock$Genus,"misclassified"))

write.csv2(file = "Results_Mock/Abundances_true_bacteria_no_missclass_table.csv", judged_table[ judged_table$Genus!="misclassified",  ] , row.names = F, quote=F)

judged_table<-melt(judged_table , id.vars = "Genus" )
#colnames(judged_table)

meta_barplot<- metadata
row.names(meta_barplot)<-meta_barplot$Workflow
meta_barplot <- meta_barplot[judged_table$variable , ]
meta_barplot$Program_sub <- gsub("Mock_true_abund","T", meta_barplot$Program_sub)
meta_barplot$Settings[grepl("contigs", meta_barplot$Program_sub) & !grepl("kMetaShot",meta_barplot$Program_sub) ] <- paste0( "contigs_", meta_barplot$Settings[grepl("contigs", meta_barplot$Program_sub) & !grepl("kMetaShot",meta_barplot$Program_sub) ] )
meta_barplot$Settings <- gsub("large_c","metalarge_c",meta_barplot$Settings)
meta_barplot$Program_sub <- gsub("Kaiju_contigs","Kaiju_nr_euk+", meta_barplot$Program_sub)
meta_barplot$Program_sub <- gsub("_contigs","", meta_barplot$Program_sub)   # impossible to fit in the plot !
meta_barplot$Program_sub <- gsub("kMetaShot$","kMetaShot_contigs", meta_barplot$Program_sub)
meta_barplot$Program_sub <- gsub("Kraken","Kraken2", meta_barplot$Program_sub)
meta_barplot$Program_sub <- gsub("Kraken2$","Kraken2_nt core", meta_barplot$Program_sub, perl=T) 
meta_barplot$Program_sub <- gsub("Kaiju$","Kaiju_nr euk", meta_barplot$Program_sub, perl=T) 
meta_barplot$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju_nr euk+", meta_barplot$Program_sub, fixed =T) 
# meta_barplot$Workflow <- gsub("_plus", "", meta_barplot$Workflow )
meta_barplot$Workflow <- gsub("True", "T", meta_barplot$Workflow )
meta_barplot$Program_sub <- factor(meta_barplot$Program_sub, 
                                   levels = c("T", "Kaiju_nr euk", "Kaiju_nr euk+", "Kraken2_nt core", "Kraken2_SILVA",
                                              "RiboF_full", "RiboF_V3V4","kMetaShot_contigs","kMetaShot_MAGs") 
)

table <- cbind.data.frame( judged_table, meta_barplot )
table$Settings <-  gsub("_nr_euk+", "", table$Settings, fixed=T)
levels(table$Program_sub)<-gsub("_","\n", levels(table$Program_sub) )


# plotting
ggplot(data=table, aes(x=Settings, y=value, fill=Genus)) +
  facet_grid( ~ Program_sub, space="free_x", scales = "free_x") +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_20) +
  theme(axis.text.x=element_text(angle=38,
                                 vjust=1,
                                 hjust=1,
                                 size= 5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=5.25),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 7.5 ),
  legend.position="bottom",
  legend.margin = margin(-20 ,5, 2 ,-20),
  plot.margin = margin(0.5,1,2,0.5)) +
  guides(fill=guide_legend(nrow=7)) +
  labs(x="", y="Percentual abundance",
       fill="",
       #title="True reads",
       caption = "'misclassified' includes every misclassified read"
  )
ggsave(file="Results_Mock/Abundances_true_bacteria_with_missclassif.png",width=6,height=4, dpi=300)
dev.off()



# TRUE -NON - BACTERIA
table2<-table[ ! table$Genus %in% "misclassified", ]
table2<-table2[ table2$Genus %in% c("phage_T4", "Homo", "Diploscapter", "Paramecium"), ]
#table2$Genus <- factor(table2$Genus , levels = c("misclassified", table_mock$Genus))

# plotting
ggplot(data=table2, aes(x=Settings, y=value, fill=Genus)) +
  #scale_y_continuous( limits = c(0,2.05)) +
  facet_grid( ~ Program_sub, space="free_x", scales = "free_x") +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_20[17:20]) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 4.8
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=5.4),
  legend.key.height = unit(0.15, "cm"),
  legend.key.width = unit(0.4, "cm"),
  legend.text = element_text ( size = 7.5 ),
  legend.position="bottom",
  legend.margin = margin(-25 ,5, 2 ,-20),
  plot.margin = margin(0.1,1,1,0.1)) +
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Percentual abundance",
       fill="",
       #title="True reads",
       caption = "ZOOM ON TRUE NON-PROKARYOTES"
  )
ggsave(file="Results_Mock/Abundances_true_bacteria_with_missclassif_ZOOM_on_metaz.png",width=6,height=4.3, dpi=300)
dev.off()




##################### TOP ABUNDANCES (misclassification) ########################

true_genera<-table_mock$Genus
judged_table<-feature_table
judged_table$Genus<-row.names(judged_table)
judged_table <- judged_table[ , ! colnames(judged_table) %in% "Known_mock" ]
judged_table <- judged_table[! judged_table$Genus %in%table_mock$Genus , ] 
genus_means <- rowMeans(judged_table[ ,colnames(judged_table)!="Genus" ])
top <- sort( genus_means , decreasing=TRUE) [1:29]
judged_table[! judged_table$Genus %in% names(top) , "Genus" ] <- "Others"
judged_table <- aggregate ( .~Genus , judged_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Genus level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 

judged_table2 <- judged_table
judged_table2 <- judged_table2[ order(rowMeans(judged_table2[ ,colnames(judged_table2)!="Genus" ]), decreasing = T) , ]
write.table( judged_table2, file="Results_Mock/Top_misclassifications_percentages.tsv", sep="\t", row.names = F, quote = F )

judged_table$Genus <- gsub("cannot be assigned to a (non-viral) genus","As_generic_virus",judged_table$Genus, fixed = T)
judged_table$Genus <- gsub("Incertae_Sedis","I.S.",judged_table$Genus)
names(top) <- gsub("cannot be assigned to a (non-viral) genus","As_generic_virus",names(top) , fixed=T )
names(top) <- gsub("Incertae_Sedis","I.S.",names(top) )
judged_table$Genus <- factor(judged_table$Genus , levels = c(names(top),"Others"))
judged_table<-melt(judged_table , id.vars = "Genus" )
#colnames(judged_table)

meta_barplot<- metadata
row.names(meta_barplot)<-meta_barplot$Workflow
meta_barplot <- meta_barplot[judged_table$variable , ]
meta_barplot$Program_sub <- gsub("Mock_true_abund","True", meta_barplot$Program_sub)
meta_barplot$Settings[grepl("contigs", meta_barplot$Program_sub) & !grepl("kMetaShot",meta_barplot$Program_sub) ] <- paste0( "contigs_", meta_barplot$Settings[grepl("contigs", meta_barplot$Program_sub) & !grepl("kMetaShot",meta_barplot$Program_sub) ] )
meta_barplot$Settings <- gsub("large_c","metalarge_c",meta_barplot$Settings)
meta_barplot$Program_sub <- gsub("Kaiju_contigs","Kaiju_nr_euk+", meta_barplot$Program_sub)
meta_barplot$Program_sub <- gsub("_contigs","", meta_barplot$Program_sub)   # impossible to fit in the plot !
meta_barplot$Program_sub <- gsub("kMetaShot$","kMetaShot_contigs", meta_barplot$Program_sub)
meta_barplot$Program_sub <- gsub("Kraken","Kraken2", meta_barplot$Program_sub)
meta_barplot$Program_sub <- gsub("Kraken2$","Kraken2_nt core", meta_barplot$Program_sub, perl=T) 
meta_barplot$Program_sub <- gsub("Kaiju$","Kaiju_nr euk", meta_barplot$Program_sub, perl=T) 
meta_barplot$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju_nr euk+", meta_barplot$Program_sub, fixed =T) 

meta_barplot$Program_sub <- factor(meta_barplot$Program_sub, 
                                   levels = c("Kaiju_nr euk", "Kaiju_nr euk+", "Kraken2_nt core", "Kraken2_SILVA",
                                              "RiboF_full", "RiboF_V3V4","kMetaShot_contigs","kMetaShot_MAGs")
)

table <- cbind.data.frame( judged_table, meta_barplot )
table$Settings <-  gsub("_nr_euk+", "", table$Settings, fixed=T)
levels(table$Program_sub)<-gsub("_","\n", levels(table$Program_sub) )


# plotting
ggplot(data=table, aes(x=Settings, y=value, fill=Genus)) +
  facet_grid( ~ Program_sub, space="free_x", scales = "free_x") +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x=element_text(angle=38,
                                 vjust=1,
                                 hjust=1,
                                 size= 5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=6.5),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.38, "cm"),
  legend.text = element_text ( size = 7.05 ),
  legend.position="bottom",
  legend.margin = margin(-15 ,5, 1 ,-30),
  plot.margin = margin(2,1,2,10)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance",
       fill="",
       # title="Most abundant misclassifications",
       caption="'Others' includes less abundant misclassifications"
  )
ggsave(file="Results_Mock/Top_misclassifications_barplot.png",width=6,height=4, dpi=300)
dev.off()



# 
# ##################### TOP ABUNDANCES (misclassification, ONLY TOP PROGRAMS) ########################
# 
# true_genera<-table_mock$Genus
# judged_table<-feature_table[ , metadata$Workflow[metadata$Program_sub %in% c("Kaiju","Kaiju_nr_euk+","Kraken","kMetaShot","RiboF_full") & metadata$Program_sub != "Kraken_SILVA"]]
# # head( judged_table , n=2)
# judged_table$Genus<-row.names(judged_table)
# judged_table <- judged_table[! judged_table$Genus %in%table_mock$Genus , ] 
# genus_means <- rowMeans(judged_table[ ,colnames(judged_table)!="Genus" ])
# top <- sort( genus_means , decreasing=TRUE) [1:29]
# judged_table[! judged_table$Genus %in% names(top) , "Genus" ] <- "Others"
# judged_table <- aggregate ( .~Genus , judged_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Genus level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 
# 
# judged_table2 <- judged_table
# judged_table2 <- judged_table2[ order(rowMeans(judged_table2[ ,colnames(judged_table2)!="Genus" ]), decreasing = T) , ]
# write.table( judged_table2, file="Results_Mock/Top_misclassifications_ONLY_KAIJU_KRAKEN_kMETA_RIBO_percentages.tsv", sep="\t", row.names = F, quote = F )
# 
# judged_table$Genus <- gsub("cannot be assigned to a (non-viral) genus","As_generic_virus",judged_table$Genus, fixed = T)
# names(top) <- gsub("cannot be assigned to a (non-viral) genus","As_generic_virus",names(top) , fixed=T )
# judged_table$Genus <- factor(judged_table$Genus , levels = c(names(top),"Others"))
# judged_table<-melt(judged_table , id.vars = "Genus" )
# #colnames(judged_table)
# 
# meta_barplot<- metadata[ metadata$Program_sub %in% c("Kaiju","Kaiju_nr_euk+" "Kraken","RiboF_full", "kMetaShot") & metadata$Program_sub != "Kraken_SILVA" , ]
# row.names(meta_barplot)<-meta_barplot$Workflow
# meta_barplot <- meta_barplot[judged_table$variable , ]
# meta_barplot$Program_sub <- gsub("Kaiju_contigs","Kaiju_nr_euk+", meta_barplot$Program_sub)
# meta_barplot$Program_sub <- gsub("kMetaShot$","kMetaShot\ncontigs", meta_barplot$Program_sub)   # impossible to fit in the plot !
# meta_barplot$Program_sub <- gsub("RiboF_full$","RiboFrame\nfull 16S", meta_barplot$Program_sub)   # impossible to fit in the plot !
# meta_barplot$Program_sub <- gsub("Kaiju$","Kaiju_nr euk", meta_barplot$Program_sub, perl=T) 
# meta_barplot$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju_nr euk+", meta_barplot$Program_sub, fixed =T) 
# meta_barplot$Program_sub <- factor(meta_barplot$Program_sub, 
#                                    levels = c( "Kaiju\nnr euk", "Kaiju\nnr euk+", "Kraken","RiboFrame\nfull 16S", "kMetaShot\ncontigs") )
# 
# levels(meta_barplot$Program_sub) <- gsub("Kraken","Kraken2\nnt core", levels(meta_barplot$Program_sub))
# 
# table <- cbind.data.frame( judged_table, meta_barplot )
# table$Settings <-  gsub("_nr_euk+", "", table$Settings, fixed=T)
# 
# # plotting
# ggplot(data=table, aes(x=Settings, y=value, fill=Genus)) +
#   facet_grid( ~ Program_sub, space="free_x", scales = "free_x") +
#   geom_bar( stat="identity", position="stack", na.rm = F) +
#   theme_classic(base_size =8.5) +
#   scale_fill_manual(values=fill_color_30) +
#   theme(axis.text.x=element_text(angle=35,
#                                  vjust=1,
#                                  hjust=1,
#                                  size= 5
#   ),
#   axis.text.y=element_text(size=6.2),
#   axis.title.y = element_text(size=8),
#   axis.title =element_text(size=10),
#   strip.text = element_text(size=6.25),
#   legend.key.height = unit(0.2, "cm"), 
#   legend.key.width = unit(0.38, "cm"),
#   legend.text = element_text ( size = 7.5 ),
#   legend.position="bottom",
#   legend.margin = margin(-8 ,5, 1 ,-30),
#   plot.margin = margin(2,1,2,5)) +
#   guides(fill=guide_legend(nrow=6)) +
#   labs(x="", y="Percentual abundance",
#        fill="",
#        # title="Most abundant misclassifications (top programs only)",
#        caption="'Others' includes less abundant misclassifications"
#   )
# ggsave(file="Results_Mock/Top_misclassifications__ONLY_KAIJU_KRAKEN_kMETA_RiboF_barplot.png",width=6,height=4, dpi=300)
# dev.off()




##################### TOP ABUNDANCES (misclassification, ONLY KAIJU) ########################

true_genera<-table_mock$Genus
judged_table<-feature_table[ , metadata$Workflow[metadata$Program_sub %in% c("Kaiju","Kaiju_nr_euk+") ]]
# head( judged_table , n=2)
judged_table$Genus<-row.names(judged_table)
judged_table <- judged_table[! judged_table$Genus %in%table_mock$Genus , ]
# NB: "the follow obs is removed to allow the visualization of other obs in the barplot ...
judged_table <- judged_table[! judged_table$Genus %in% "cannot be assigned to a (non-viral) genus" , ]

genus_means <- rowMeans(judged_table[ ,colnames(judged_table)!="Genus" ])
top <- sort( genus_means , decreasing=TRUE) [1:29]
judged_table[! judged_table$Genus %in% names(top) , "Genus" ] <- "Others"
judged_table <- aggregate ( .~Genus , judged_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Genus level... because now replicated can be found (Kraken is launched on both SILVA and NCBI)

judged_table2 <- judged_table
judged_table2 <- judged_table2[ order(rowMeans(judged_table2[ ,colnames(judged_table2)!="Genus" ]), decreasing = T) , ]
write.table( judged_table2, file="Results_Mock/Top_misclassifications_ONLY_KAIJU_percentages.tsv", sep="\t", row.names = F, quote = F )

judged_table$Genus <- factor(judged_table$Genus , levels = c(names(top),"Others"))
judged_table<-melt(judged_table , id.vars = "Genus" )
#colnames(judged_table)

meta_barplot<- metadata[ metadata$Program_sub %in% c("Kaiju","Kaiju_nr_euk+") , ]
meta_barplot$Program_sub <- gsub("Kaiju$","Kaiju\nnr euk", meta_barplot$Program_sub , perl=T )
meta_barplot$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju\nnr euk+", meta_barplot$Program_sub, fixed = T)
row.names(meta_barplot)<-meta_barplot$Workflow
meta_barplot <- meta_barplot[judged_table$variable , ]
meta_barplot$Program_sub <- factor(meta_barplot$Program_sub,
                                   levels = c( "Kaiju\nnr euk","Kaiju\nnr euk+" ) )

table <- cbind.data.frame( judged_table, meta_barplot )
table$Settings <-  gsub("_nr_euk+", "", table$Settings, fixed=T)

# plotting
ggplot(data=table, aes(x=Settings, y=value, fill=Genus)) +
  facet_grid( ~ Program_sub, space="free_x", scales = "free_x") +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x=element_text(angle=35,
                                 vjust=1,
                                 hjust=1,
                                 size= 5
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=8),
  legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.35, "cm"),
  legend.text = element_text ( size = 7.5 ),
  legend.position="bottom",
  legend.margin = margin(-8 ,5, 1 ,-30),
  plot.margin = margin(2,1,2,5)) +
  guides(fill=guide_legend(nrow=6)) +
  labs(x="", y="Percentual abundance",
       fill="",
       #title="Most abundant misclassifications (Kaiju only)",
       caption="'Others' includes less abundant misclassifications"
  )
ggsave(file="Results_Mock/Top_misclassifications__ONLY_KAIJU_barplot.png",width=6,height=4, dpi=300)
dev.off()


### checking nr euk + unique misclassific
mis_k_plus<- row.names(feature_table)[feature_table$Kaiju_E0.01_m30_plus>0]
mis_k_plus<- mis_k_plus[!mis_k_plus %in%   row.names(feature_table)[feature_table$Kaiju_E0.01_m30>0] ]
mis_k_plus<- mis_k_plus [ ! mis_k_plus %in% table_mock$Genus ]
write.csv2(feature_table[mis_k_plus, grepl("_plus", colnames(mis_k_plus))], file="Results_Mock/Unique_misclassific_Kaiju_db_plus.csv", row.names = T, quote = F)




##################### PCoA ##################### 

colors_sub_program <- c("Kaiju_nr euk"="darkgreen", "Kaiju_nr euk+"="chartreuse3", "Kaiju_contigs"="chartreuse",  "Kraken2_nt core"="red3", "Kraken2_contigs"="coral3",  "Kraken2_SILVA"="coral", 
                        "kMetaShot_contigs"="violet", "kMetaShot_MAGs"="darkviolet",
                        "RiboFrame_full"="blue2", "RiboFrame_V3V4"="deepskyblue", "Mock_true_abundances"="yellow2")
colors_sub_program2 <- colors_sub_program
names(colors_sub_program2) <- gsub( "_" ," ", names(colors_sub_program) )



### Hellinger (euclidean + sqrt)
this <- sqrt( t(feature_table) )
dist <- vegdist(this, method = "euclidean", na.rm = F)
pcoa <- cmdscale (dist, eig = TRUE)
colnames(pcoa$points) <- c("PC1","PC2")
# ensuring the same order between the PCoA samples in the related obj and the metadata
meta_pcoa<-metadata
row.names(meta_pcoa)<-meta_pcoa$Workflow
meta_pcoa<- meta_pcoa[row.names(pcoa$points) , ] 
meta_pcoa$Program_sub<-gsub("abund","abundances",meta_pcoa$Program_sub)
meta_pcoa$Program_sub<-gsub("RiboF","RiboFrame",meta_pcoa$Program_sub)
meta_pcoa$Program_sub<-gsub("kMetaShot$","kMetaShot_contigs",meta_pcoa$Program_sub)
meta_pcoa$Program_sub <- gsub("Kraken","Kraken2", meta_pcoa$Program_sub)
meta_pcoa$Program_sub <- gsub("Kraken2$","Kraken2_nt core", meta_pcoa$Program_sub, perl=T) 
meta_pcoa$Program_sub <- gsub("Kaiju$","Kaiju_nr euk", meta_pcoa$Program_sub, perl=T) 
meta_pcoa$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju_nr euk+", meta_pcoa$Program_sub, fixed =T) 
meta_pcoa$Program_sub<-gsub("_"," ",meta_pcoa$Program_sub)
meta_pcoa$Settings <- gsub("_nr_euk.*","", meta_pcoa$Settings)
# test plot
# ordiplot (pcoa, display = 'sites', type = 'text',
#           xlab=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
#           ylab=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
# )
# plot
pcoa$points <- as.data.frame(pcoa$points)
ggplot( data= pcoa$points, aes( x= PC1, y=PC2 , color= meta_pcoa$Program_sub )) +
  xlim( c(-5,10) ) + # to include the whole labels
  theme_classic() +
  scale_color_manual( values = colors_sub_program2 ) +
  geom_point(size= 3, alpha= 0.6 ) +
  geom_text( label= gsub("_SILVA","",meta_pcoa$Settings) ,
             show.legend = F , size= 2 , color="black") +
  theme( legend.text = element_text(size=6.5),
         legend.margin = margin( 1,1,1, 1) 
         )+
  labs(color=NULL, 
       #title= "PCoA on processed mock (Hellinger)",
       x=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
       y=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
  )
ggsave(file="Results_Mock/PCoA_Hellinger_mock.png", dpi=300, width = 5, height = 4.5)



### Hellinger again, but magnifing the most effective methods
#subset_table<- feature_table[ , grepl("RiboF_full",colnames(feature_table)) | grepl("Kaiju_",colnames(feature_table)) | grepl("MAGs",colnames(feature_table)) | grepl("mock",colnames(feature_table)) ]
#subset_table<- subset_table[ , ! grepl("SILVA|contigs", colnames(subset_table) ) ]
subset_table<- feature_table[ , ! grepl("SILVA", colnames(feature_table) ) ]
this <- sqrt( t(subset_table ) )
dist <- vegdist(this, method = "euclidean", na.rm = F)
pcoa <- cmdscale (dist, eig = TRUE)
colnames(pcoa$points) <- c("PC1","PC2")
meta_pcoa<-metadata  # ensuring the same order between the PCoA samples in the related obj and the metadata
row.names(meta_pcoa)<-meta_pcoa$Workflow
# this time, the following reordering acts also as sub-setting (only the selected samples are maintained)
meta_pcoa<- meta_pcoa[row.names(pcoa$points) , ] 
meta_pcoa$Program_sub<-gsub("abund","abundances",meta_pcoa$Program_sub)
meta_pcoa$Program_sub<-gsub("kMetaShot$","kMetaShot_contigs",meta_pcoa$Program_sub)
meta_pcoa$Program_sub<-gsub("RiboF_","RiboFrame_",meta_pcoa$Program_sub)
meta_pcoa$Program_sub <- gsub("Kraken","Kraken2", meta_pcoa$Program_sub)
meta_pcoa$Program_sub <- gsub("Kraken2$","Kraken2_nt core", meta_pcoa$Program_sub, perl=T) 
meta_pcoa$Program_sub <- gsub("Kaiju$","Kaiju_nr euk", meta_pcoa$Program_sub, perl=T) 
meta_pcoa$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju_nr euk+", meta_pcoa$Program_sub, fixed =T) 
meta_pcoa$Program_sub<-gsub("_"," ",meta_pcoa$Program_sub)
meta_pcoa$Settings <- gsub("_nr_euk.*","", meta_pcoa$Settings)
pcoa$points <- as.data.frame(pcoa$points)
ggplot( data= pcoa$points, aes( x= PC1, y=PC2 , color= meta_pcoa$Program_sub )) +
  theme_classic() +
  xlim( c(-7.1, 4.2) ) + # to include the whole labels
  scale_color_manual( values = colors_sub_program2[ unique(meta_pcoa$Program_sub ) ] ) +
  geom_point(size= 3.65, alpha= 0.55 ) +
  geom_text( label= gsub("_SILVA","",meta_pcoa$Settings) ,
             show.legend = F , size= 2.2 , color="black") +
  theme( legend.text = element_text(size=6.5 ),
         legend.margin = margin( 1,-2,1, -6),
         legend.spacing.x = unit(0,"pt")
  )+
  labs(color=NULL, 
       # title= "PCoA on processed mock (Hellinger, top pipelines)",
       x=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
       y=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
  )
ggsave(file="Results_Mock/PCoA_Hellinger_mock_BETTER_PIPELINES.png", dpi=300, width = 5, height = 4.6)


### Bray-Curtis
#this <- t(feature_table)
this <- t(subset_table)  # only most effective methods
dist <- vegdist(this, method = "bray", na.rm = F)
pcoa <- cmdscale (dist, eig = TRUE)
colnames(pcoa$points) <- c("PC1","PC2")
# ensuring the same order between the PCoA samples in the related obj and the metadata
meta_pcoa<-metadata
row.names(meta_pcoa)<-meta_pcoa$Workflow
meta_pcoa<- meta_pcoa[row.names(pcoa$points) , ] 
meta_pcoa$Program_sub<-gsub("abund","abundances",meta_pcoa$Program_sub)
meta_pcoa$Program_sub<-gsub("RiboF","RiboFrame",meta_pcoa$Program_sub)
meta_pcoa$Program_sub<-gsub("kMetaShot$","kMetaShot_contigs",meta_pcoa$Program_sub)
meta_pcoa$Program_sub <- gsub("Kraken","Kraken2", meta_pcoa$Program_sub)
meta_pcoa$Program_sub <- gsub("Kraken2$","Kraken2_nt core", meta_pcoa$Program_sub, perl=T) 
meta_pcoa$Program_sub <- gsub("Kaiju$","Kaiju_nr euk", meta_pcoa$Program_sub, perl=T) 
meta_pcoa$Program_sub <- gsub("Kaiju_nr_euk+","Kaiju_nr euk+", meta_pcoa$Program_sub, fixed =T) 
meta_pcoa$Program_sub<-gsub("_"," ",meta_pcoa$Program_sub)
meta_pcoa$Settings <- gsub("_nr_euk.*","", meta_pcoa$Settings)
pcoa$points <- as.data.frame(pcoa$points)
ggplot( data= pcoa$points, aes( x= PC1, y=PC2 , color= meta_pcoa$Program_sub )) +
  xlim( c(-0.55, 0.42) ) + # to include the whole labels
  theme_classic() +
  scale_color_manual( values = colors_sub_program2 ) +
  geom_point(size= 3.65, alpha= 0.55 ) +
  geom_text( label= gsub("_SILVA","",meta_pcoa$Settings) ,
             show.legend = F , size= 2 , color="black") +
  theme( legend.text = element_text(size=6.5 ),
         legend.margin = margin( 1,1,1, 1) 
  )+
  labs(color=NULL, 
       # title= "PCoA on processed mock (Bray-Curtis)",
       x=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
       y=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
  )
ggsave(file="Results_Mock/PCoA_Bray_mock_BETTER_PIPELINES.png", dpi=300, width = 5.2, height = 4.6)


# ### Bray again, but magnifying the most effective methods (Kaiju only)
# subset_table<- feature_table[ , grepl("Kaiju_",colnames(feature_table)) | grepl("mock",colnames(feature_table)) ]
# this <- t( subset_table )
# dist <- vegdist(this, method = "bray", na.rm = F)
# pcoa <- cmdscale (dist, eig = TRUE)
# colnames(pcoa$points) <- c("PC1","PC2")
# meta_pcoa<-metadata  # ensuring the same order between the PCoA samples in the related obj and the metadata
# row.names(meta_pcoa)<-meta_pcoa$Workflow
# # this time, the following reordering acts also as sub-setting (only the selected samples are maintained)
# meta_pcoa<- meta_pcoa[row.names(pcoa$points) , ] 
# meta_pcoa$Program_sub<-gsub("_"," ",meta_pcoa$Program_sub)
# # meta_pcoa$Settings<-gsub("_","\n",meta_pcoa$Settings)
# ggplot( data= pcoa$points, aes( x= PC1, y=PC2 , color= meta_pcoa$Program_sub )) +
#   theme_classic() +
#   scale_color_manual( values = colors_sub_program2 ) +
#   geom_point(size= 3.5, alpha= 0.6 ) +
#   geom_text( label= gsub("_SILVA","",meta_pcoa$Settings) ,
#              #lineheight = .6
#              show.legend = F , size= 2 , 
#              color="black"
#              
#              ) +
#   theme( legend.text = element_text(size=6),
#          legend.margin = margin( 1,1,1, 1) 
#   )+
#   labs(color=NULL, 
#        title= "PCoA on processed mock (Bray, only Kaiju)",
#        x=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
#        y=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
#   )
# ggsave(file="Results_Mock/PCoA_Bray_mock_BETTER_PIPELINES.png", dpi=300, width = 5.5, height = 5)




######################## VEEN DIAGRAMMS ###########################

# venn custom function (exports txt and ggplot computing everything by itself ...)
custom_venn<- function( selection_program, selection_settings1, selection_settings2, colors_venn) {
  GROUP1 <- feature_table[ , metadata$Workflow[metadata$Program_sub %in%selection_program & metadata$Settings==selection_settings1], drop=F  ]
  GROUP2 <- feature_table[ , metadata$Workflow[metadata$Program_sub %in% selection_program & metadata$Settings==selection_settings2], drop=F  ]
  GROUP1 <- row.names( GROUP1[ rowSums(GROUP1)>0 , , drop=F] )
  GROUP2 <- row.names( GROUP2[ rowSums(GROUP2)>0 , , drop=F] )
  
  
  true_genera <- table_mock$Genus
  x<-list( GROUP1, GROUP2, true_genera)
  names(x) <- c( paste0(selection_settings1,"  ") , paste0("  ",selection_settings2)  ,"True")
  
  ggvenn(x, stroke_size = 0.6,
         set_name_size = 4.5,
         show_percentage = F,
         stroke_color = "darkgray",
         text_size = 5.2,
         text_color = "black",
         fill_color = colors_venn ) +
    theme(plot.margin = margin(-5,-15,0,-15)
    )
  ggsave(filename =  paste0("Results_Mock/Venn_diagrams/Exclusive_obs_",selection_program,"___",selection_settings1,"_vs_",selection_settings2,".png"),
         width = 4.5, height = 4.5, dpi=300, bg = "white")
  #dev.off()
  
  
  ### exporting the names ...
  ONLY_IN_GROUP1<- GROUP1[! GROUP1 %in% GROUP2]
  # Now collapsing in one vector for the script  
  ONLY_IN_GROUP1<- paste(ONLY_IN_GROUP1, collapse = ", ")
  #head(ONLY_IN_GROUP1)
  
  ONLY_IN_GROUP2<- GROUP2[! GROUP2 %in% GROUP1]
  ONLY_IN_GROUP2<- paste(ONLY_IN_GROUP2, collapse = ", ")
  #head(ONLY_IN_GROUP2)
  
  ONLY_IN_MOCK<- true_genera[ ! true_genera %in% GROUP1 & ! true_genera %in% GROUP2 ]
  ONLY_IN_MOCK_AND_GROUP1 <- true_genera[ ! true_genera %in% GROUP2 & true_genera %in% GROUP1 ]
  ONLY_IN_MOCK_AND_GROUP2 <- true_genera[ ! true_genera %in% GROUP1 & true_genera %in% GROUP2 ]
  
  con<-file( paste0("Results_Mock/Venn_diagrams/Exclusive_obs_",selection_program,"___",selection_settings1,"_vs_",selection_settings2,".txt") )
  sink(con, append=TRUE)
  cat( paste0("\n\nONLY IN ", selection_settings2, " and not in ", selection_settings1), fill=TRUE)
  cat(ONLY_IN_GROUP2)
  cat( paste0("\n\nONLY IN ", selection_settings1, " and not in ", selection_settings2), fill=TRUE)
  cat(ONLY_IN_GROUP1)
  cat( "\n\nONLY IN THE MOCK (TRUE) and not classified from neither of the settings", fill=TRUE)
  cat(ONLY_IN_MOCK)
  cat( paste0("\n\nONLY IN IN THE MOCK (TRUE) and in ", selection_settings1), fill=TRUE)
  cat(ONLY_IN_MOCK_AND_GROUP1)
  cat( paste0("\n\nONLY IN IN THE MOCK (TRUE) and in ", selection_settings2), fill=TRUE)
  cat(ONLY_IN_MOCK_AND_GROUP2)
  sink()
  close(con)
  
}


# unique(metadata$Program_sub)
# unique(metadata$Settings)
custom_venn(selection_program = "Kraken", selection_settings1 = "c0.45", selection_settings2 = "c0.85",
            colors_venn= c("coral","red4","yellow3") )
custom_venn(selection_program = "Kraken", selection_settings1 = "c0.15", selection_settings2 = "c0.99",
            colors_venn= c("coral","red","yellow3") )
custom_venn(selection_program = "Kaiju", selection_settings1 = "E0.01_m11", selection_settings2 = "E0.00001_m42",
            colors_venn= c("green","darkgreen","yellow3") )
custom_venn(selection_program = "Kaiju", selection_settings1 = "E0.00001_m30", selection_settings2 = "E0.00001_m42",
            colors_venn= c("green","darkgreen","yellow3") )
custom_venn(selection_program = "Kaiju", selection_settings1 = "E0.01_m42", selection_settings2 = "E0.00001_m42",
            colors_venn= c("green","darkgreen","yellow3") )
custom_venn(selection_program = "Kaiju_nr_euk+", selection_settings1 = "E0.00001_m11_nr_euk+", selection_settings2 = "E0.00001_m42_nr_euk+",
            colors_venn= c("green3","darkgreen","yellow3") )
custom_venn(selection_program = "Kaiju_nr_euk+", selection_settings1 = "E0.00001_m30_nr_euk+", selection_settings2 = "E0.00001_m42_nr_euk+",
            colors_venn= c("green3","darkgreen","yellow3") )
custom_venn(selection_program = "RiboF_full", selection_settings1 = "c0.8", selection_settings2 = "c0.99",
            colors_venn= c("deepskyblue","blue","yellow3") )
# custom_venn(selection_program = "kMetaShot_contigs", selection_settings1 = "default_c0.2", selection_settings2 = "large_c0.2",
#             colors_venn= c("pink","violet","yellow3") )
# custom_venn(selection_program = "kMetaShot_contigs", selection_settings1 = "default_c0.4", selection_settings2 = "large_c0.4",
#             colors_venn= c("pink","violet","yellow3") )
# custom_venn(selection_program = "kMetaShot_contigs", selection_settings1 = "default_c0.2", selection_settings2 = "custom_c0.2",
#             colors_venn= c("pink","violet","yellow3") )
custom_venn(selection_program = "kMetaShot_MAGs", selection_settings1 = "default_c0.2", selection_settings2 = "large_c0.2",
            colors_venn= c("magenta","darkviolet","yellow3") )
suppressWarnings( rm(GROUP1, GROUP2, ONLY_IN_GROUP1, ONLY_IN_GROUP2, ONLY_IN_MOCK, ONLY_IN_MOCK_AND_GROUP1, ONLY_IN_GROUP2, con,
                     selection_program, selection_settings1, selection_settings2 ) )




#################### !!! GOING TO SPECIES LEVEL !!! #########################

dir.create("Results_Mock/Species_level")

# Kaiju
Kaiju_s<-read.delim(file="Obtained_Counts_mock/Species_level/kaiju_classif_summary_species_level.tsv", sep="\t")
Kaiju_s<-Kaiju_s[grepl("DNAseq_Evalue0.00001_m42_output.txt", Kaiju_s$file) , c("reads","taxon_name")]
Kaiju_s$taxon_name<-gsub(";","",Kaiju_s$taxon_name)
colnames(Kaiju_s)<-c("Abundance","Bacterium")
Kaiju_s <- Kaiju_s[ order(Kaiju_s$Abundance, decreasing = T) , ]
# Kaiju counts are already aggregated (unique entries)


# Kaiju +
Kaiju_sp<-read.delim(file="Obtained_Counts_mock/Species_level/kaiju_classif_summary_species_level.tsv", sep="\t")
Kaiju_sp<-Kaiju_sp[grepl("DNAseq_Evalue0.00001_m42_output_PLUS.txt", Kaiju_sp$file) , c("reads","taxon_name")]
Kaiju_sp$taxon_name<-gsub(";","",Kaiju_sp$taxon_name)
colnames(Kaiju_sp)<-c("Abundance","Bacterium")
Kaiju_sp <- Kaiju_sp[ order(Kaiju_sp$Abundance, decreasing = T) , ]


# Kraken (0.15)
kraken_s15<-read.delim(file="Obtained_Counts_mock/Species_level/mock_AGS_DNAseq_bracken_abund_Species15.tsv", sep="\t")
kraken_s15<-kraken_s15[ , c("new_est_reads", "name")]
colnames(kraken_s15)<-c("Abundance","Bacterium")
kraken_s15 <- kraken_s15[ order(kraken_s15$Abundance, decreasing = T) , ]

# Kraken (0.99)
kraken_s99<-read.delim(file="Obtained_Counts_mock/Species_level/mock_AGS_DNAseq_bracken_abund_Species99.tsv", sep="\t")
kraken_s99<-kraken_s99[ , c("new_est_reads", "name")]
colnames(kraken_s99)<-c("Abundance","Bacterium")
kraken_s99 <- kraken_s99[ order(kraken_s15$Abundance, decreasing = T) , ]

# kMetaShot
kMeta_s<- read.delim(file="Obtained_Counts_mock/kMetaShot_MAGs_classif_MAGs_metalarge_conf0.2", sep=",")
kMeta_s<-kMeta_s[kMeta_s$species!=0 , c("organism_name"), drop=F]  # NB: if species 0 --> the obs have a confidence threshold below the set value for that level --> the confidence is taken into account
kMeta_s<- cbind.data.frame( "Abundance"=as.numeric(table(kMeta_s)) , "Bacterium"=names(table(kMeta_s)) )
kMeta_s <- kMeta_s[ order(kMeta_s$Abundance, decreasing = T) , ]


# checking the true species abundances in the mock
mock_species<- read.delim(file="Obtained_Counts_mock/Species_level/True_mock_species_level_counts.tsv", sep="\t", header = T)




### Solving taxonomies incompatibilities
Kaiju_s$Bacterium <- gsub(" ","_", Kaiju_s$Bacterium)
Kaiju_sp$Bacterium <- gsub(" ","_", Kaiju_sp$Bacterium)
kMeta_s$Bacterium <- gsub(" ","_", kMeta_s$Bacterium)
kraken_s15$Bacterium <- gsub(" ","_", kraken_s15$Bacterium)
kraken_s99$Bacterium <- gsub(" ","_", kraken_s99$Bacterium)
mock_species$Bacterium <- gsub(" ","_", mock_species$Bacterium)

# Kaiju_s$Bacterium[ ! Kaiju_s$Bacterium %in% mock_species$Bacterium]
# Kaiju_s$Bacterium <- gsub("non ci sono classificaz in T4","Enterobacteria_phage_T4", Kaiju_s$Bacterium )
Kaiju_s$Bacterium <- gsub("Candidatus_Moranbacteria_bacterium.*","Candidatus_Moranbacteria_sp.", Kaiju_s$Bacterium )
Kaiju_s$Bacterium <- gsub("Solirubrobacterales_bacterium_67-14.*","Solirubrobacterales_bacterium_sp.", Kaiju_s$Bacterium )
Kaiju_s$Bacterium <- gsub("_sp..*","_sp.", Kaiju_s$Bacterium )
Kaiju_s$Bacterium <- gsub("[","", Kaiju_s$Bacterium , fixed = T)
Kaiju_s$Bacterium <- gsub("]","", Kaiju_s$Bacterium , fixed = T)
Kaiju_s <- Kaiju_s[ ! Kaiju_s$Bacterium %in% c("unclassified","cannot_be_assigned_to_a_(non-viral)_sp.") , ]
Kaiju_s$Bacterium <-  gsub("Nostocoides","Tetrasphaera", Kaiju_s$Bacterium)

Kaiju_sp$Bacterium <- gsub("Candidatus_Moranbacteria_bacterium.*","Candidatus_Moranbacteria_sp.", Kaiju_sp$Bacterium )
Kaiju_sp$Bacterium <- gsub("Solirubrobacterales_bacterium_67-14.*","Solirubrobacterales_bacterium_sp.", Kaiju_sp$Bacterium )
Kaiju_sp$Bacterium <- gsub("_sp..*","_sp.", Kaiju_sp$Bacterium )
Kaiju_sp$Bacterium <- gsub("[","", Kaiju_sp$Bacterium , fixed = T)
Kaiju_sp$Bacterium <- gsub("]","", Kaiju_sp$Bacterium , fixed = T)
Kaiju_sp <- Kaiju_sp[ ! Kaiju_sp$Bacterium %in% c("unclassified","cannot_be_assigned_to_a_(non-viral)_sp.") , ]
Kaiju_sp$Bacterium <-  gsub("Nostocoides","Tetrasphaera", Kaiju_sp$Bacterium)


# kMeta_s$Bacterium[ ! kMeta_s$Bacterium %in% mock_species$Bacterium]
kMeta_s$Bacterium <- gsub("Candidatus_Competibacter_denitrificans.*","Candidatus_Competibacter_denitrificans", kMeta_s$Bacterium)
kMeta_s$Bacterium <- gsub("Candidatus_Competibacter_phosphatis.*","Candidatus_Competibacter_phosphatis", kMeta_s$Bacterium)
kMeta_s$Bacterium <- gsub("Thauera_phenylacetica_B4P","Thauera_phenylacetica", kMeta_s$Bacterium)
kMeta_s$Bacterium <- gsub("Thauera_butanivorans_NBRC_103042","Thauera_butanivorans", kMeta_s$Bacterium)
kMeta_s$Bacterium <- gsub("_sp..*","_sp.", kMeta_s$Bacterium)

# kraken_s15$Bacterium[ ! kraken_s15$Bacterium %in% mock_species$Bacterium]
kraken_s15$Bacterium <- gsub("Candidatus_Moranbacteria_bacterium.*","Candidatus_Moranbacteria_sp.", kraken_s15$Bacterium )
kraken_s15$Bacterium <- gsub("Solirubrobacterales_bacterium_67-14.*","Solirubrobacterales_bacterium_sp.", kraken_s15$Bacterium )
kraken_s15$Bacterium <- gsub("_sp..*","_sp.", kraken_s15$Bacterium )
kraken_s15$Bacterium <- gsub("[","", kraken_s15$Bacterium , fixed = T)
kraken_s15$Bacterium <- gsub("]","", kraken_s15$Bacterium , fixed = T)
kraken_s15$Bacterium <- gsub("Tequatrovirus_T4","Enterobacteria_phage_T4", kraken_s15$Bacterium )
kraken_s15$Bacterium <- gsub("Nostocoides","Tetrasphaera", kraken_s15$Bacterium )

kraken_s99$Bacterium <- gsub("Candidatus_Moranbacteria_bacterium.*","Candidatus_Moranbacteria_sp.", kraken_s99$Bacterium )
kraken_s99$Bacterium <- gsub("Solirubrobacterales_bacterium_67-14.*","Solirubrobacterales_bacterium_sp.", kraken_s99$Bacterium )
kraken_s99$Bacterium <- gsub("_sp..*","_sp.", kraken_s99$Bacterium )
kraken_s99$Bacterium <- gsub("[","", kraken_s99$Bacterium , fixed = T)
kraken_s99$Bacterium <- gsub("]","", kraken_s99$Bacterium , fixed = T)
kraken_s99$Bacterium <- gsub("Tequatrovirus_T4","Enterobacteria_phage_T4", kraken_s99$Bacterium )
kraken_s99$Bacterium <- gsub("Nostocoides","Tetrasphaera", kraken_s99$Bacterium )

Kaiju_s <- aggregate( . ~ Bacterium , Kaiju_s , FUN=sum )
Kaiju_sp <- aggregate( . ~ Bacterium , Kaiju_sp , FUN=sum )
kMeta_s <- aggregate( . ~ Bacterium , kMeta_s , FUN=sum )
kraken_s15 <- aggregate( . ~ Bacterium , kraken_s15 , FUN=sum )
kraken_s99 <- aggregate( . ~ Bacterium , kraken_s99 , FUN=sum )




### Filtering singletons
# Kaiju_s_discarded<- Kaiju_s[Kaiju_s$Abundance<=1, ]
# Kraken_s_discarded<- kraken_s15[kraken_s15$Abundance<=1, ]
# if( length(which( Kaiju_s_discarded$Bacterium %in% mock_species$Bacterium ) ) >= 1 ) {
#   stop("Wait! There are true mock bacteria among the singletons!!!!")
# }
# if( length(which( kraken_s15$Bacterium %in% mock_species$Bacterium ) ) >= 1 ) {
#   stop("Wait! There are true mock bacteria among the singletons!!!!")
# }
#Kaiju_s_discarded$Bacterium[ Kaiju_s_discarded$Bacterium %in% mock_species$Bacterium ]
#Kraken_s_discarded$Bacterium[ Kraken_s_discarded$Bacterium %in% mock_species$Bacterium ]

# There are true mock bacteria among the Kraken2 singletons (Paramecium sonneborni)... it is not possible to filter them!




### Parsing between programs (while adding the known "true" abundances of the mock)
files_list<- list( kMeta_s, kraken_s15, kraken_s99, Kaiju_s , Kaiju_sp,  mock_species )
names(files_list) <-  c("kMeta_s", "kraken_s15", "kraken_s99", "Kaiju_s" , "Kaiju_sp" , "mock_species" )
every_entry<-NULL
for (x in 1:length(files_list)) {
  every_entry <- c(every_entry, files_list[[x]][["Bacterium"]])
}
# from list (each set_pipeline is repeated along a column) to feature table (each set_pipeline is a column)
feature_table<-cbind("Temp_column"=unique(every_entry))
# each bacterium in a row --> every bacterium in this table --> searching across datasets
for( s in 1:length(files_list) ) {
  set_pipeline<- names(files_list)[[s]]  # 1 set_pipeline at time
  table_4_that_set_pipeline<- files_list[[s]]
  abundances<-NULL  # resets the object, ready for the next bacterium
  # an empty column to be filled with each sample...
  temp_set_pipeline<- matrix( "temp", length(unique(every_entry)) , length(colnames(table_4_that_set_pipeline)) )
  temp_set_pipeline <- as.data.frame( temp_set_pipeline )
  temp_set_pipeline[ ,1] <- unique(every_entry)
  for( b in unique(every_entry) ) {
    # scan each bacterium in order: if found in that set_pipeline then appends the correspective abundances, if absent paste a 0
    if( b %in% table_4_that_set_pipeline$Bacterium ){
      abundances<- as.numeric(table_4_that_set_pipeline[ table_4_that_set_pipeline$Bacterium==b, "Abundance" , drop = T ])   # NB: the entire ROW will be captured this time
    } else {
      abundances<- 0  # that bacterium is absent --> a zero
    }
    temp_set_pipeline[ temp_set_pipeline$V1 == b , ] <- c(b , abundances)  # overwrite the "temp"s row with actual values
  }
  colnames(temp_set_pipeline)<- c("Bacterium", names(files_list)[[s]] )
  feature_table<-cbind.data.frame(feature_table, temp_set_pipeline)  # this adds the sample to the complete table
}
feature_table <- feature_table[ , ! colnames(feature_table) %in% "Bacterium" ]
beasts <- feature_table$Temp_column  # this column correspond to another "Genus" column
# NB: "apply" will delete the colnames, then the genera names will be assigned as row.names afterwards
feature_table$Temp_column <- NULL
feature_table<-as.data.frame(apply(feature_table, MARGIN=2, FUN= as.numeric))
row.names(feature_table) <- beasts



### Transforming in proportions
feature_table <- as.data.frame( apply( feature_table, MARGIN = 2, FUN= function(x) (x/sum(x))*100 ) )
feature_table <- round( feature_table, digits = 4)




##################### ABUNDANCES TOP TRUE SPECIES ########################

# NB: it is not possible to plot the TRUE species because they are too much! I do not have enough colors... 
# true_species<-mock_species$Bacterium
judged_table<-feature_table
judged_table$Bacterium<-row.names(judged_table)
row.names(judged_table)<-NULL
#judged_table[! judged_table$Bacterium %in% mock_species$Bacterium , "Bacterium"] <- "misclassified" 
#TOP<-judged_table[ order(rowMeans(judged_table[,!colnames(judged_table)%in%"Bacterium"]),decreasing=T) , "Bacterium"] [1:20]
TOP<-mock_species$Bacterium[ order(mock_species$Abundance,decreasing=T) ] [1:28]  # most abundant TRUE bacteria
original_vector <- judged_table$Bacterium  # one vector, two subsequents modifications (otherwise the second one would depend on the first one)
judged_table[ ! original_vector %in% TOP , "Bacterium"] <- "Others"
judged_table[ ! original_vector %in% mock_species$Bacterium , "Bacterium"] <- "misclassified"

write.csv2(file = "Results_Mock/Species_level/Abundances_TOP28_TRUE_SPECIES_also_misclassified_in_others__table.csv", judged_table[ !judged_table$Bacterium %in% c("Others","misclassified"),  ] , row.names = F, quote=F)

judged_table <- aggregate ( .~Bacterium , judged_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Bacterium level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 

judged_table$Bacterium <- factor(judged_table$Bacterium , levels = c( judged_table$Bacterium[!judged_table$Bacterium%in%c("Others","misclassified")] ,"Others","misclassified"))
judged_table<-melt(judged_table , id.vars = "Bacterium" )
#colnames(judged_table)

meta_barplot<- judged_table$variable  # column name from the function melt, it corresponds to the original tables names
meta_barplot<- gsub("kMeta_s","kMetaShot\nMAGs", meta_barplot)
#meta_barplot<- gsub("kraken_s15","Kraken2 (conf 0.15)", meta_barplot)
#meta_barplot<- gsub("kraken_s99","Kraken2 (conf 0.99)", meta_barplot)
meta_barplot<- gsub("kraken_s.*","Kraken2\nnt core", meta_barplot)  # unique Kraken2 grid
meta_barplot<- gsub("Kaiju_sp","Kaiju\nnr euk+", meta_barplot)
meta_barplot<- gsub("Kaiju_s","Kaiju\nnr euk", meta_barplot)
meta_barplot<- gsub("mock_species","True", meta_barplot)
#meta_barplot <- factor(meta_barplot , levels = c("True", "Kaiju (E 0.001 m 40)", "Kraken2 (conf 0.15)", "Kraken2 (conf 0.99)", "kMetaShot\nMAGs (metalarge)") )
meta_barplot <- factor(meta_barplot , levels = c("True", "Kaiju\nnr euk", "Kaiju\nnr euk+", "Kraken2\nnt core", "kMetaShot\nMAGs") )

table <- cbind.data.frame( judged_table, meta_barplot )

table$variable <- gsub("_s$","", table$variable)
table$variable <- gsub("kraken_s","c0.", table$variable)
table$variable <- gsub("Kaiju.*","E.00001_m42", table$variable)
table$variable <- gsub("kMeta","metalarge_c0.2", table$variable)

# plotting
ggplot(data=table, aes(x=variable, y=value, fill=Bacterium)) +
  facet_grid( ~ meta_barplot, space="free_x", scales = "free_x") +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30_species) +
  theme(axis.text.x=element_text(angle=20,
                                 vjust=1,
                                 hjust=1,
                                 size= 6
  ),
  axis.text.y=element_text(size=6.5),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=7),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.21, "cm"),
  legend.text = element_text ( size = 6.8 ),
  legend.position="bottom",
  legend.margin = margin(-13 ,5, 1 ,-33),
  plot.margin = margin(1,1,2,10)) +
  guides(fill=guide_legend(nrow=8)) +
  labs(x="", y="Percentual abundance",
       fill="",
       # title="28 Most abundant true species",
       caption = "'Others' includes less abundant true species\n'misclassified' include every misclassified read"
  )
ggsave(file="Results_Mock/Species_level/Abundances_TOP28_TRUE_SPECIES_also_misclassified_in_others.png",width=6,height=4, dpi=300)
dev.off()




################# ABUNDANT TOP misclassified (Species) #################

judged_table<-feature_table
judged_table$Bacterium<-row.names(judged_table)
#row.names(judged_table)<-NULL
judged_table <- judged_table[ , ! colnames(judged_table) %in% "mock_species" ]
judged_table <- judged_table[! judged_table$Bacterium %in%mock_species$Bacterium , ] 
Bacterium_means <- rowMeans(judged_table[ ,colnames(judged_table)!="Bacterium" ])
top <- sort( Bacterium_means , decreasing=TRUE) [1:29]
judged_table[! judged_table$Bacterium %in% names(top) , "Bacterium" ] <- "Others"
judged_table <- aggregate ( .~Bacterium , judged_table , FUN=sum)   # the counts in every column (which are counts of each sample) are summed based on the Bacterium level... because now replicated can be found (Kraken is launched on both SILVA and NCBI) 

judged_table2 <- judged_table
judged_table2 <- judged_table2[ order(rowMeans(judged_table2[ ,colnames(judged_table2)!="Bacterium" ]), decreasing = T) , ]
write.table( judged_table2, file="Results_Mock/Species_level/Top_misclassifications_species_percentages.tsv", sep="\t", row.names = F, quote = F )

# judged_table$Bacterium <- gsub("cannot be assigned to a (non-viral) Bacterium","As_generic_virus",judged_table$Bacterium, fixed = T)
# judged_table$Bacterium <- gsub("Incertae_Sedis","I.S.",judged_table$Bacterium)
# names(top) <- gsub("cannot be assigned to a (non-viral) Bacterium","As_generic_virus",names(top) , fixed=T )
# names(top) <- gsub("Incertae_Sedis","I.S.",names(top) )
judged_table$Bacterium <- factor(judged_table$Bacterium , levels = c(names(top),"Others"))
judged_table<-melt(judged_table , id.vars = "Bacterium" )
#colnames(judged_table)

meta_barplot<- judged_table$variable  # column name from the function melt, it corresponds to the original tables names
meta_barplot<- gsub("kMeta_s","kMetaShot\nMAGs", meta_barplot)
#meta_barplot<- gsub("kraken_s15","Kraken2 (conf 0.15)", meta_barplot)
#meta_barplot<- gsub("kraken_s99","Kraken2 (conf 0.99)", meta_barplot)
meta_barplot<- gsub("kraken_s.*","Kraken2\nnt core", meta_barplot)  # unique Kraken2 grid
meta_barplot<- gsub("Kaiju_sp","Kaiju\nnr euk+", meta_barplot)
meta_barplot<- gsub("Kaiju_s","Kaiju\nnr euk", meta_barplot)
meta_barplot<- gsub("mock_species","True", meta_barplot)
#meta_barplot <- factor(meta_barplot , levels = c("True", "Kaiju (E 0.001 m 40)", "Kraken2 (conf 0.15)", "Kraken2 (conf 0.99)", "kMetaShot\nMAGs (metalarge)") )
meta_barplot <- factor(meta_barplot , levels = c("True", "Kaiju\nnr euk", "Kaiju\nnr euk+", "Kraken2\nnt core", "kMetaShot\nMAGs") )
# meta_barplot <- factor(meta_barplot , levels = c("True", "Kaiju (E 0.001 m 40)", "Kraken (conf 0.15)", "Kraken (conf 0.99)", "kMetaShot\nMAGs (metalarge)") )

table <- cbind.data.frame( judged_table, meta_barplot )


# plotting
ggplot(data=table, aes(x=variable, y=value, fill=Bacterium)) +
  facet_grid( ~ meta_barplot, space="free_x", scales = "free_x") +
  geom_bar( stat="identity", position="stack", na.rm = F) +
  theme_classic(base_size =8.5) +
  scale_fill_manual(values=fill_color_30) +
  theme(axis.text.x=element_text(angle=30,
                                 vjust=1,
                                 hjust=1,
                                 size= 5.4
  ),
  axis.text.y=element_text(size=6.2),
  axis.title.y = element_text(size=8),
  axis.title =element_text(size=10),
  strip.text = element_text(size=5.7),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.text = element_text ( size = 6.9 ),
  legend.position="bottom",
  legend.margin = margin(-16 ,5, 1 ,-32),
  plot.margin = margin(2,1,2,10)) +
  guides(fill=guide_legend(nrow=10)) +
  labs(x="", y="Percentual abundance",
       fill="",
       title="Most abundant misclassifications",
       caption="'Others' includes less abundant misclassifications"
  )
ggsave(file="Results_Mock/Species_level/Top_misclassifications_SPECIES.png",width=6.05,height=4.5, dpi=300)
dev.off()




##################### PCoA (Species) ##################### 

colors_sub_program <- c("Kaiju (E 0.0001 m 42)"="darkgreen", "Kaiju+ (E 0.0001 m 42)"="green", "Kraken (conf 0.15)"="red2", "Kraken (conf 0.99)"="red4", "kMetaShot\nMAGs (metalarge)"="darkviolet",
                        "Mock true abundances"="yellow2")



### Hellinger (euclidean + sqrt)
this <- sqrt( t(feature_table) )
dist <- vegdist(this, method = "euclidean", na.rm = F)
pcoa <- cmdscale (dist, eig = TRUE)
colnames(pcoa$points) <- c("PC1","PC2")
# ensuring the same order between the PCoA samples in the related obj and the metadata
meta_barplot<- row.names(pcoa$points)
meta_barplot<- gsub("kMeta_s","kMetaShot\nMAGs (metalarge)", meta_barplot)
meta_barplot<- gsub("kraken_s15","Kraken (conf 0.15)", meta_barplot)
meta_barplot<- gsub("kraken_s99","Kraken (conf 0.99)", meta_barplot)
meta_barplot<- gsub("Kaiju_sp","Kaiju+ (E 0.0001 m 42)", meta_barplot)
meta_barplot<- gsub("Kaiju_s","Kaiju (E 0.0001 m 42)", meta_barplot)
meta_barplot<- gsub("mock_species","Mock true abundances", meta_barplot)
meta_barplot <- factor(meta_barplot , levels = c("Mock true abundances", "Kaiju (E 0.0001 m 42)", "Kaiju+ (E 0.0001 m 42)", "Kraken (conf 0.15)", "Kraken (conf 0.99)",  "kMetaShot\nMAGs (metalarge)") )

labels_plot <- ifelse( row.names(pcoa$points)=="mock_species", "Mock", "")  # empty if it is not the mock

# test plot
# ordiplot (pcoa, display = 'sites', type = 'text',
#           xlab=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
#           ylab=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
# )

# plot
pcoa$points <- as.data.frame(pcoa$points)
ggplot( data= pcoa$points, aes( x= PC1, y=PC2 , color= meta_barplot )) +
  # xlim( c(-5,10) ) + # to include the whole labels
  theme_classic() +
  scale_color_manual( values = colors_sub_program ) +
  geom_point(size= 4, alpha= 0.65 ) +
  geom_text( label= labels_plot ,
             show.legend = F , size= 2.5 , color="black", alpha= 1) +
  theme( legend.text = element_text(size=6.5),
         legend.margin = margin( 1,1,1, 1) 
  )+
  labs(color=NULL, 
       title= "PCoA on processed mock (Hellinger)",
       x=paste( "PC1 (", round( pcoa$eig[1]/sum(pcoa$eig)*100,2), "%)"),
       y=paste( "PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,2), "%)")
  )
ggsave(file="Results_Mock/Species_level/PCoA_Hellinger_mock_species.png", dpi=300, width = 5, height = 4.5)




##################### VENN DIAGRAMS (Species) #########################


# venn custom function 2 (exports txt and ggplot computing everything by itself ...)
custom_venn_speciesLevel<- function( selection_program1, selection_program2, colors_venn) {
  GROUP1 <- feature_table[ , selection_program1, drop=F  ]
  GROUP2 <- feature_table[ , selection_program2, drop=F  ]
  GROUP1 <- row.names( GROUP1[ rowSums(GROUP1)>0 , , drop=F] )
  GROUP2 <- row.names( GROUP2[ rowSums(GROUP2)>0 , , drop=F] )
  
  
  true_species <- mock_species$Bacterium
  x<-list( GROUP1, GROUP2, true_species)
  names(x) <- c( paste0(selection_program1,"  ") , paste0("  ",selection_program2)  ,"True")
  
  ggvenn(x, stroke_size = 0.6,
         set_name_size = 4.5,
         show_percentage = F,
         stroke_color = "darkgray",
         text_size = 5.2,
         text_color = "black",
         fill_color = colors_venn ) +
    theme(plot.margin = margin(-5,-15,0,-15)
    )
  ggsave(filename =  paste0("Results_Mock/Species_level/Exclusive_obs_",selection_program1,"__vs__",selection_program2,".png"),
         width = 4.5, height = 4.5, dpi=300, bg = "white")
  #dev.off()
  
  
  ### exporting the names ...
  ONLY_IN_GROUP1<- GROUP1[! GROUP1 %in% GROUP2]
  # Now collapsing in one vector for the script  
  ONLY_IN_GROUP1<- paste(ONLY_IN_GROUP1, collapse = ", ")
  #head(ONLY_IN_GROUP1)
  
  ONLY_IN_GROUP2<- GROUP2[! GROUP2 %in% GROUP1]
  ONLY_IN_GROUP2<- paste(ONLY_IN_GROUP2, collapse = ", ")
  #head(ONLY_IN_GROUP2)
  
  ONLY_IN_MOCK<- true_species[ ! true_species %in% GROUP1 & ! true_species %in% GROUP2 ]
  ONLY_IN_MOCK_AND_GROUP1 <- true_species[ ! true_species %in% GROUP2 & true_species %in% GROUP1 ]
  ONLY_IN_MOCK_AND_GROUP2 <- true_species[ ! true_species %in% GROUP1 & true_species %in% GROUP2 ]
  
  con<-file( paste0("Results_Mock/Species_level/Exclusive_obs_",selection_program1,"__vs__",selection_program2,".txt") )
  sink(con, append=TRUE)
  cat( paste0("\n\nONLY IN ", selection_program2, " and not in ", selection_program1), fill=TRUE)
  cat(ONLY_IN_GROUP2)
  cat( paste0("\n\nONLY IN ", selection_program1, " and not in ", selection_program2), fill=TRUE)
  cat(ONLY_IN_GROUP1)
  cat( "\n\nONLY IN THE MOCK (TRUE) and not classified from neither of the settings", fill=TRUE)
  cat(ONLY_IN_MOCK)
  cat( paste0("\n\nONLY IN IN THE MOCK (TRUE) and in ", selection_program1), fill=TRUE)
  cat(ONLY_IN_MOCK_AND_GROUP1)
  cat( paste0("\n\nONLY IN IN THE MOCK (TRUE) and in ", selection_program2), fill=TRUE)
  cat(ONLY_IN_MOCK_AND_GROUP2)
  sink()
  close(con)
  
}


#selection_program1 <- colnames(feature_table)[1]
#selection_program2 <- colnames(feature_table)[2]
custom_venn_speciesLevel(selection_program1 = "Kaiju_s", selection_program2 = "kraken_s15",
            colors_venn= c("darkgreen","red2","yellow3") )
custom_venn_speciesLevel(selection_program1 = "kMeta_s", selection_program2 = "kraken_s99",
                         colors_venn= c("magenta3","red4","yellow3") )
custom_venn_speciesLevel(selection_program1 = "Kaiju_s", selection_program2 = "kMeta_s",
                         colors_venn= c("darkgreen","magenta3","yellow3") )
custom_venn_speciesLevel(selection_program1 = "Kaiju_s", selection_program2 = "Kaiju_sp",
                         colors_venn= c("darkgreen","green2","yellow3") )
custom_venn_speciesLevel(selection_program1 = "Kaiju_sp", selection_program2 = "kMeta_s",
                         colors_venn= c("green2","magenta3","yellow3") )
custom_venn_speciesLevel(selection_program1 = "Kaiju_sp", selection_program2 = "kMeta_s",
                         colors_venn= c("green2","magenta3","yellow3") )
custom_venn_speciesLevel(selection_program1 = "Kaiju_sp", selection_program2 = "kraken_s99",
                         colors_venn= c("green2","red4","yellow3") )


suppressWarnings( rm(GROUP1, GROUP2, ONLY_IN_GROUP1, ONLY_IN_GROUP2, ONLY_IN_MOCK, ONLY_IN_MOCK_AND_GROUP1, ONLY_IN_GROUP2, con,
                     selection_program, selection_settings1, selection_settings2 ) )




##################### R AND PACKAGES VERSION #########################

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("Results_Mock/R_version_and_packages.txt")
sink(con, append = TRUE)

cat(package$R.version$version.string)
cat("   running on", package$running)
cat("\n", "\n", fill=TRUE)
cat("\n", "\n", fill=TRUE)
package$otherPkgs$ggplot2[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$vegan[c(1,3)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$ecodist[1:2]
cat("\n", "\n", fill=TRUE)
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))

