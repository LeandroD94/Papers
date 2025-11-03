################### PREPARING THE ENVIRONMENT ##################

library("reshape2")

options(scipen=100)

table_abs<-read.table("Diamond_output/Abs_counts_from_Diamond.tsv", row.names = 1, header = T)
dictionary<-read.delim("Diamond_output/IDs_infos.tsv", header = T, sep="\t")
metadata<-read.delim("PHA_RNA_Metadata.tsv", header = T, sep="\t")

colnames(table_abs)<-gsub("^X","",colnames(table_abs))



#### Adding unmapped reads
mapping<-read.table("Diamond_output/Diamond_matching_percents.tsv", row.names = 1, header = T)
mapping<-mapping[colnames(table_abs),  ]
unmapped_count <- mapping$INPUT_SEQS_FROM_BOTH - mapping$CONCORD_MATCHES   # these are unmapped but still (theoretically) microbial
# NB: "inputs from both" does not mean summing up R1 and R2, moreover the "Concord_matches" includes also the unique R1 and R2 matches (therefore is the count of every successful not discarded mapping)
table_abs<-rbind(table_abs, UNMAPPED=unmapped_count)

# just a further sanity check ...
# not_micr<-as.data.frame(read.delim("Diamond_output/abs_NOT_microb_counts_from_Diamond.tsv", row.names = 1, header = T, sep="\t"))
# colnames(not_micr)<-gsub("^X","",colnames(not_micr))
# not_micr_maps<-colSums(not_micr[ ,  colnames(table_abs)] ) #  this step both orders and subsets
# mapping [1, "CONCORD_MATCHES"] - not_micr_maps[[1]] == colSums(table_abs[row.names(table_abs)!="UNMAPPED",])[[1]]    # the "not microbial" were removed from the total successful mappings in a last "manual" filter during the alignments parsing



#### Custom names
table_abs <- table_abs[ , metadata$FASTQ_ID] # this re-orders
colnames(table_abs)<-metadata$Sample_name
# table_TPM <- table_TPM[ , metadata$FASTQ_ID]
# colnames(table_TPM)<-metadata$Sample_name




################### FILTERING ABS TABLE ###################

# this is an important step according to PMID 27441086 and PMID 27508061

if ( ! "table_abs_original" %in% ls() ) {
  table_abs_original <- table_abs
}

max_values <- apply( table_abs_original , MARGIN = 1, max )
# discarded <- table_abs_original[ ! row.names(table_abs_original) %in% names(max_values[max_values>2]) , ]
# discarded
table_abs<-table_abs_original[ names(max_values[max_values>2]), ]
remaining_message<- paste( length(table_abs[,1]), "maintained on", length(table_abs_original[,1]) ,"filtered reads -->", round( (length(table_abs[,1]) /length(table_abs_original[,1]))*100 , 3) ,"% remaining")
remaining_message

# table_abs<-table_abs_original # to reset the filter



### also the dictionary has to be re-filtered (to avoid bugs and NA during the re-subsettings)
dictionary <- dictionary[ dictionary$Entry %in% row.names(table_abs) , ]
length(unique(dictionary$Entry))
length(unique(dictionary$Taxon))

proof_filter<-"Singletons and Doubletons have been removed, see the 'remaining message'"




###################### PREPARING THE GENE LEVEL OBJECT ##########################
# this is a slow step, hence it is performed only once, here, and after the filtering step!

# glomming to gene level
table_gene<-cbind.data.frame(table_abs[dictionary$Entry, ], Gene=dictionary$Gene)
table_gene<-rbind.data.frame(table_gene , c(table_abs["UNMAPPED", ], Gene="UNMAPPED" ) )
table_gene <- table_gene[ ! table_gene$Gene=="" , ]  # fragments and proteins without name
table_gene <- table_gene[ ! is.na(table_gene$Gene), ]
table_gene <- table_gene[ grepl("[A-z]",table_gene$Gene), ]  # at least one letter in the gene name, to avoid that rare rows which derives from errors due to gene without name
table_gene <- table_gene[ ! grepl("^[A-Z]$", table_gene$Gene) , ]  # other fragments without name, which had an anomalous gene name (just a letter) during the gene table formation
table_gene<-aggregate(.~Gene, table_gene, FUN=sum)

# head( table_gene[order(table_gene$Gene), ] )

row.names(table_gene)<-table_gene$Gene
table_gene$Gene<-NULL

# dim(table_gene)

dictionary$Pathway <- NULL   # to save space (the file is already too big)




############## DEFINING THE PHA RELATED GENES #################

# PHA related genes (see paper introduction for references)
only_PHAgenes_IDs<-dictionary [ grepl("^pha", dictionary$Gene ) |
                                  grepl("^phb", dictionary$Gene ) |
                                  # grepl("phosphate", dictionary$Description )  |  # such as psi (phosphate stress inducible)
                                  grepl("[P-p]hasin", dictionary$Description )  |
                                  grepl("[P-p]olyhydroxyalkanoate", dictionary$Description )  |
                                  grepl("[P-p]olyhydroxybutyr", dictionary$Description )  |
                                  grepl("-hydroxyalkanoate", dictionary$Description , fixed=T)  |
                                  grepl("Polyhydroxyalkanoic", dictionary$Description )  |
                                  grepl("beta-hydroxybutyrate", dictionary$Description , fixed=T )  |
                                  grepl("Acetyl-CoA acetyltransf", dictionary$Description , fixed=T )  | # simil phbA
                                  grepl("Acetoacetyl-CoA reduct", dictionary$Description , fixed=T )  |  # simil phbA
                                  grepl("PHB depol", dictionary$Description , fixed=T )  |  #
                                  grepl("PHA/PHB", dictionary$Description , fixed=T )  |  #
                                  grepl("PHB/PHA", dictionary$Description , fixed=T )  |  #
                                  grepl("[P-p]oly(R)-hydroxyalkanoic", dictionary$Description )  ,
]


### Searching for other pha genes using UniprotKB 
### downloading results as table from https://www.uniprot.org/uniprotkb?query=Polyhydroxyalkanoate 
# system("wget -r https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28Polyhydroxyalkanoate%29 -O pha_results_UniProtKB.gz ")
# system("gzip -d pha_results_UniProtKB.gz ")
# file.rename("pha_results_UniProtKB","pha_results_UniProtKB.tsv")
UniProt_KB_IDs<-read.delim("pha_results_UniProtKB.tsv", header = T, sep = "\t")
UniProt_KB_IDs$Gene.Names <- UniProt_KB_IDs$Gene.Names
UniProt_KB_IDs<- UniProt_KB_IDs[ !duplicated(UniProt_KB_IDs$Gene.Names) , ]
# UniProt_KB_IDs[grepl("pha", UniProt_KB_IDs)]  # just to check

### adding the "others" pha related genes to the list
UniProt_Plus <- UniProt_KB_IDs [ UniProt_KB_IDs$Entry %in% dictionary$Entry , "Entry" ]
UniProt_Plus <- UniProt_Plus [! UniProt_Plus %in% only_PHAgenes_IDs$Entry ]
length(UniProt_Plus)   # other 27 entries retrieved using UniProt
rm(UniProt_KB_IDs)

PHA_IDs <- c( only_PHAgenes_IDs$Entry , UniProt_Plus )
  



############## SAVING THE FINAL OBJECT TO ANALYZE ... ##############

save( file= "Tabelle_RNA_PHA_after_filters.RData",
      list = c("table_abs", "table_gene", "mapping", "unmapped_count", 
               "dictionary", "metadata", "remaining_message", "proof_filter",
               "PHA_IDs"
               ) 
      )




######## \\\\\\ R AND PACKAGES VERSION , ANALYSIS NOTES \\\\\\ ###############

### if on Windows, change "$otherPkgs" with "$loadedOnly"

package<-sessionInfo()

con <- file("R_version_and_packages.txt")
sink(con, append = TRUE)

cat(package$R.version$version.string)
cat("   running on", package$running)
cat("\n", "\n", fill=TRUE)
package$otherPkgs$ggplot2[1:2]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$vegan[c(1,3)]
cat("\n", "\n", fill=TRUE)
package$otherPkgs$DESeq2[1:2]
cat("\n", "\n", fill=TRUE)
cat("\n", "\n", fill=TRUE)
print("vegan")
packageVersion("vegan")
cat("\n", "\n", fill=TRUE)
print("ecodist")
packageVersion("ecodist")
cat("\n", "\n", fill=TRUE)
cat("\n \n \nEvery package: \n", fill=TRUE)
#print(package$otherPkgs)
print(package$loadedOnly)

sink()
close(con)
suppressWarnings(rm(con))
