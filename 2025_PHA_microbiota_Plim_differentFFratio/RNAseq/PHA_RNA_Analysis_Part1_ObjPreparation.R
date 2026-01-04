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




############# OVERWRITING BIOLOGICALLY IMPOSSIBLE TAXONOMIC ASSEGNATIONS ##########################

# Ca. Accumulibacter and Nitrosomonas cannot grow in this reactors, indeed they are not find in the 16S DNA dataset
# The corresponding abundant transcripts are actually almost identical proteins (SAME FUNCTION) but belonging to other bacteria, usually Thauera or Azorcus --> overwriting

# A0A080MAY5 (supposed Uncharacterized of Accumulibacter)
dictionary[ dictionary$Entry=="A0A080MAY5" , "Entry"] <- "Q5P6Y71"    # 63% identity according to blast (protein db), the best one of those related to bacteria featured in the reactor (see DNA)
dictionary[ dictionary$Entry=="Q5P6Y71" , "Taxon"] <- "Azoarcus sp."  # Should be "Aromatoleum aromaticum", which is actually an old name of Azoarcus
row.names(table_abs)[row.names(table_abs)=="A0A080MAY5"]<-"Q5P6Y71"

# A0A011MM28 (Another supposed uncharacterized of Accumulibacter)
dictionary[ dictionary$Entry=="A0A011MM28" , "Entry"] <- "Q5P9111"    # 63% identity according to blast (protein db), the best one of those related to bacteria featured in the reactor (see DNA)
dictionary[ dictionary$Entry=="Q5P9111" , "Taxon"] <- "Azoarcus sp."  # Should be "Aromatoleum aromaticum", which is actually an old name of Azoarcus
row.names(table_abs)[row.names(table_abs)=="A0A011MM28"]<-"Q5P9111"

# A0A08M108 (supposed ClpA of Accumulibacter)
dictionary[ dictionary$Entry=="A0A963ML30" , "Entry"] <- "A0A974SMY9"    # 92% identity according to blast (protein db)
dictionary[ dictionary$Entry=="A0A974SMY9" , "Taxon"] <- "Azospira sp."  # NB: not Azoarcus (confirmed in LPSN)
row.names(table_abs)[row.names(table_abs)=="A0A963ML30"]<-"A0A974SMY9"

# A0A0F7KJ89 (supposed CspA of Nitrosomonas)
dictionary[ dictionary$Entry=="A0A0F7KJ89" , "Entry"] <- "A0A0S3AGH0"    # 92% identity according to blast (protein db)
dictionary[ dictionary$Entry=="A0A0S3AGH0" , "Taxon"] <- "Azoarcus sp."  # or, alternatively, it was assigned to "Aromatoleum"...
row.names(table_abs)[row.names(table_abs)=="A0A0F7KJ89"]<-"A0A0S3AGH0"

# A0A1H7JWE9 (supposed CspA of Nitrosomonas)
# NB: this one can be assigned either to Azoarcus or Thaeura (both >96.5% identity) but the abundance patterns mimic the DNA abundance of Thauera in 16S dataset
dictionary[ dictionary$Entry=="A0A1H7JWE9" , "Entry"] <- "A0A235ETV1"    # 92% identity according to blast (protein db)
dictionary[ dictionary$Entry=="A0A235ETV1" , "Taxon"] <- "Thauera propionica"  # or, alternatively, it was assigned to "Aromatoleum"...
row.names(table_abs)[row.names(table_abs)=="A0A1H7JWE9"]<-"A0A235ETV1"



### MOREOVER ...
# this one is a taxonomic synonymous (see LPSN)
dictionary[ dictionary$Taxon=="Aromatoleum aromaticum" , "Taxon"] <- "Azoarcus sp."    # 92% identity according to blast (protein db)




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
                                  grepl("Poly (3-hydroxyalkan", dictionary$Description , fixed=T)  |
                                  grepl("poly(3-hydroxyalkan", dictionary$Description, fixed=T )  |
                                  grepl("Poly(3-hydroxybutyrate)", dictionary$Description )  |
                                  grepl("Poly (3-hydroxybutyrate)", dictionary$Description )  |
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




################ *** DEFINING THE STARVATION REGULATORY GENES *** ###########################

starvation_genes<-dictionary [ grepl("[S-s]tarvation", dictionary$Description ) |    
                                 grepl("[S-s]stationary", dictionary$Description ) | 
                                 grepl("[S-s]tringent response", dictionary$Description ) | 
                                 grepl("[F-f]amine", dictionary$Description ) | 
                                 grepl("stationary phase", dictionary$Description ) | 
                                 grepl("deficiency", dictionary$Description ) | 
                                 grepl("ribosome silencing factor", dictionary$Description ) |    #  increased when nutrients are low, doi: 10.1016/j.str.2015.07.014
                                 grepl("^rsfS ", dictionary$Gene ) |    #  increased when nutrients are low, doi: 10.1016/j.str.2015.07.014
                                 grepl("[H-h]ybernation", dictionary$Description ) |    # ribosome Hybernation during stress
                                 grepl("ElaB", dictionary$Description ) |    # ElaB --> ribosome Hybernation during stress/stationary phase
                                 grepl("^elaB", dictionary$Gene ) |    # ElaB --> ribosome Hybernation during stress/stationary phase
                                 grepl("[S-s]tress", dictionary$Description ) | 
                                 grepl(" Usp", dictionary$Description ) |  # involved in environmental stress and also in nutrient limitation https://www.sciencedirect.com/science/article/pii/S0981942824010271
                                 grepl(" Usp$", dictionary$Description ) |  # involved in environmental stress and also in nutrient limitation https://www.sciencedirect.com/science/article/pii/S0981942824010271
                                 grepl("^usp", dictionary$Gene ) |   # # involved in environmental stress and also in nutrient limitation https://www.sciencedirect.com/science/article/pii/S0981942824010271
                                 grepl("Growth inhibitor", dictionary$Description ) |
                                 grepl("Inhibitor of Growth", dictionary$Description ) |
                                 grepl("inhibits growth", dictionary$Description ) |
                                 grepl("^cst", dictionary$Gene ) |  # literally "carbon starvation" operon
                                 grepl("^phoR", dictionary$Gene ) |  # membrane sensor of low phosphate, DOI: 10.1128/spectrum.02260-23
                                 grepl("PhoR", dictionary$Description ) |  # membrane sensor of low phosphate, DOI: 10.1128/spectrum.02260-23
                                 grepl("^phoB", dictionary$Gene ) |  # transcr regulator linked to low phosphate, DOI: 10.1128/spectrum.02260-23
                                 grepl("PhoB", dictionary$Description ) |  # transcr regulator linked to low phosphate, DOI: 10.1128/spectrum.02260-23
                                 grepl("PhoH", dictionary$Description ) |  # also involved in P starvation, see UniProt
                                 grepl("^phoH", dictionary$Gene ) |  # also involved in P starvation, see UniProt
                                 grepl("^ntrB", dictionary$Gene ) |  # Nitrogen starvation, see fig 1 of https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("^ntrC", dictionary$Gene ) |  # Nitrogen starvation, see fig 1 of https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("NtrB", dictionary$Description ) |  # Nitrogen starvation, see fig 1 of https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("NtrC", dictionary$Description ) |  # Nitrogen starvation, see fig 1 of https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("[S-s]igma 38", dictionary$Description ) |  # sigma 38 (S), starvation sigma factor, see doi: 10.1128/jb.178.2.470-476.1996
                                 grepl("Sigma-38", dictionary$Description , fixed=T ) |  # sigma 38 (S), starvation sigma factor, see doi: 10.1128/jb.178.2.470-476.1996
                                 grepl("sigma-38", dictionary$Description , fixed=T ) |  # sigma 38 (S), starvation sigma factor, see doi: 10.1128/jb.178.2.470-476.1996
                                 grepl("RpoS", dictionary$Description ) |  # sigma 38 (S), starvation sigma factor, see doi: 10.1128/jb.178.2.470-476.1996
                                 grepl("^rpoS", dictionary$Gene ) |  # sigma 38, starvation sigma factor, see doi: 10.1128/jb.178.2.470-476.1996
                                 grepl("^arcZ", dictionary$Gene ) |  # ArcZ activates the transcription of Sigma 38, see doi.org/10.1074/jbc.REV119.005593
                                 grepl("^rprA", dictionary$Gene ) |  # RprA too activates the transcription of Sigma 38, see doi.org/10.1074/jbc.REV119.005593
                                 grepl("^csbD", dictionary$Gene ) |  # Stress response protein induced by P limitation, see https://www.uniprot.org/uniprotkb/P70964/entry
                                 grepl("CsbD", dictionary$Description ) |  # Stress response protein induced by P limitation, see https://www.uniprot.org/uniprotkb/P70964/entry
                                 grepl("MFS transp", dictionary$Description ) |  # Energy efficient "Major Facilitator Superfamily" transporters, expressed when the nutrients are low, see doi: 10.1128/Spectrum.01846-21
                                 grepl("^rssB", dictionary$Gene ) |  # Adapter which helps ClpX to degrade Sigma 38, see https://doi.org/10.1099/mic.0.001195
                                 grepl("^csp", dictionary$Gene ) |  # the cold shock proteins are induced either by cold or nutrient starvation, especially cspD # doi: 10.3389/fmicb.2016.01151
                                 grepl("[C-c]old ", dictionary$Description ) |  # the cold shock proteins are induced either by cold or nutrient starvation, # doi: 10.3389/fmicb.2016.01151
                                 grepl("Cold-", dictionary$Description , fixed=T) |  # the cold shock proteins are induced either by cold or nutrient starvation, # doi: 10.3389/fmicb.2016.01151
                                 grepl("cold-", dictionary$Description , fixed=T) |  # the cold shock proteins are induced either by cold or nutrient starvation, # doi: 10.3389/fmicb.2016.01151
                                 grepl("^hsp", dictionary$Gene) |  # also the heat shock are involved in nutrient depletion response, see  doi: 10.1371/journal.pone.0165980 and DOI: 10.4315/0362-028x-66.11.2045
                                 grepl("[H-h]eat ", dictionary$Description) |  # also the heat shock are involved in nutrient depletion response, see  doi: 10.1371/journal.pone.0165980 and DOI: 10.4315/0362-028x-66.11.2045
                                 grepl("heat-", dictionary$Description , fixed=T) |  # also the heat shock are involved in nutrient depletion response, see  doi: 10.1371/journal.pone.0165980 and DOI: 10.4315/0362-028x-66.11.2045
                                 grepl("Heat-", dictionary$Description , fixed=T) |  # also the heat shock are involved in nutrient depletion response, see  doi: 10.1371/journal.pone.0165980 and DOI: 10.4315/0362-028x-66.11.2045
                                 grepl("Heat-", dictionary$Description , fixed=T) |  # also the heat shock are involved in nutrient depletion response, see  doi: 10.1371/journal.pone.0165980 and DOI: 10.4315/0362-028x-66.11.2045
                                 grepl("^ibp", dictionary$Gene) |  # ibp are small heat shock
                                 grepl("^grpE", dictionary$Gene) |  # heat shock induced by starvation , https://doi.org/10.1016/j.micres.2019.126309
                                 grepl(" GrpE", dictionary$Gene) |  # heat shock induced by starvation , https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("^hdeA", dictionary$Gene) |  # heat shock induced by starvation , https://doi.org/10.1016/j.micres.2019.126309
                                 grepl(" HdeA", dictionary$Description) |  # heat shock induced by starvation
                                 # The pattern HSP includes directly proteins such as Hsp70 and Hsp90... also these heat shocks are involved with environmental stress including starvation! See DOI:10.3390/ijms26020528 (bacteria) and DOI:10.3389/fphys.2021.753914 (insecta)
                                 grepl(" Hsp", dictionary$Description) |  # see above
                                 grepl(" HSP", dictionary$Description) |  # see above
                                 grepl("^Hsp", dictionary$Description) |  # see above
                                 grepl("^HSP", dictionary$Description) |  # see above
                                 grepl("GroEL", dictionary$Description) |  #  Hsp60 system
                                 grepl("GroES", dictionary$Description) |  # Hsp60 system
                                 # Moreover, new role was found for DnaK in regulating sigma factor abundance during starvation, see https://digitalcommons.unl.edu/dissertations/AAI9929226/
                                 grepl("^dnaK", dictionary$Gene) |  #  it is Hsp70
                                 grepl("DnaK", dictionary$Description) |  # it is Hsp70
                                 grepl("^dnaJ", dictionary$Gene) |  # it is Hsp40 , for additional reference see https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("DnaJ", dictionary$Description) |  # it is Hsp40 , for additional reference see https://doi.org/10.1016/j.micres.2019.126309
                                 grepl("^hrcA", dictionary$Gene) |  # hrcA is the repressor of Hsp60 transcription, see doi: 10.1128/jb.180.19.5129-5134.1998
                                 grepl("Dps ", dictionary$Description) |  # DNA-binding protein from starved cells, Dps, is a universally conserved prokaryotic ferritin that, in many species, also binds DNA , see doi: 10.1128/jb.00036-22 
                                 grepl("Dps$", dictionary$Description) |  # DNA-binding protein from starved cells, Dps, is a universally conserved prokaryotic ferritin that, in many species, also binds DNA , see doi: 10.1128/jb.00036-22 
                                 grepl("^dps", dictionary$Gene ) |  # DNA-binding protein from starved cells, Dps, is a universally conserved prokaryotic ferritin that, in many species, also binds DNA , see doi: 10.1128/jb.00036-22  
                                 grepl("CodY", dictionary$Description ) |   # metabolic regulator activated by high GTP level, it signals aa depletion, regulated by RSH, see https://doi.org/10.1128/iai.01439-09 
                                 grepl("^codY", dictionary$Gene ) |   # metabolic regulator activated by high GTP level, it signals aa depletion, regulated by RSH, see https://doi.org/10.1128/iai.01439-09 
                                 grepl("^sspA", dictionary$Gene ) |   # starvation protein , https://doi.org/10.1111/mmi.14996
                                 grepl("SspA", dictionary$Description ) |   # starvation protein , https://doi.org/10.1111/mmi.14996
                                 grepl("^hns", dictionary$Gene ) |  # gene of H-NS, against Fis activity, see https://www.uniprot.org/uniprotkb/P0ACF8/entry
                                 grepl("H-NS", dictionary$Description , fixed=T) | # see above
                                 grepl("histone-like nucleoid-structuring prot", dictionary$Description , fixed=T) | # see above
                                 #  grepl("lrp", dictionary$Gene) | # collaborates with H-NS, see DOI: 10.1111/j.1365-2958.2005.04873.x
                                 grepl("Lrp ", dictionary$Description , fixed=T) | # collaborates with H-NS, see DOI: 10.1111/j.1365-2958.2005.04873.x
                                 grepl("Lrp$", dictionary$Description) | # collaborates with H-NS, see DOI: 10.1111/j.1365-2958.2005.04873.x
                                 # (p)ppGpp , second messenger that is universally conserved among bacteria, triggered by nutrient starvation, see DOI: 10.1016/j.molcel.2020.08.005
                                 grepl("ppGpp", dictionary$Description ) | # it should be "(p)ppGpp" but the first p is omitted to avoid pattern misinterpretations 
                                 grepl("guanosine tetraphosph", dictionary$Description ) | # it should be "(p)ppGpp" but the first p is omitted to avoid pattern misinterpretations 
                                 grepl("guanosine pentaphosph", dictionary$Description ) | # it should be "(p)ppGpp" but the first p is omitted to avoid pattern misinterpretations 
                                 grepl("guanosine 5'-diphosphate 3'-diphosph", dictionary$Description ) | # it should be "(p)ppGpp" but the first p is omitted to avoid pattern misinterpretations 
                                 grepl("^relX", dictionary$Gene ) | # gene involved in ppGpp amount regulation , see https://pubmed.ncbi.nlm.nih.gov/342913/
                                 grepl("^relA", dictionary$Gene ) | # gene involved in ppGpp amount regulation , see https://pubmed.ncbi.nlm.nih.gov/342913/
                                 grepl("^relP", dictionary$Gene ) | # domain homology with (p)ppGpp synthetase, see https://doi.org/10.1128/iai.01439-09
                                 grepl("^relQ", dictionary$Gene ) | # domain homology with (p)ppGpp synthetase, see https://doi.org/10.1128/iai.01439-09
                                 grepl("RelX", dictionary$Description ) | 
                                 grepl("RelA", dictionary$Description ) |
                                 grepl("GTP pyrophosphokinase", dictionary$Description ) | # relA coded protein, see the gene ref above
                                 grepl("RelP", dictionary$Description ) | 
                                 grepl("RelQ", dictionary$Description ) | 
                                 grepl("^rsh", dictionary$Gene ) | # (p)ppGpp synthase, homolog of relA, see https://doi.org/10.1128/iai.01439-09
                                 grepl(" RSH", dictionary$Description ) |
                                 grepl("RSH$", dictionary$Description ) |
                                 grepl("^gppA", dictionary$Gene ) |  #  PppGpp 5'-phosphohydrolase and exopolyphosphatase, see https://www.uniprot.org/uniprotkb/Q74CD3/entry
                                 # *** RSH is responsible for (p)ppGpp synthesis in response to amino acid deprivation (lack of Leu/Val) or mupirocin treatment (see the ref above)
                                 grepl("^spoT", dictionary$Gene ) | # (p)ppGpp synthase, homolog of relA and rsh, see https://doi.org/10.1128/iai.01439-09
                                 grepl("SpoT", dictionary$Description) | # (p)ppGpp synthase, homolog of relA and rsh, see https://doi.org/10.1128/iai.01439-09
                                 grepl("Guanosine-3',5'-bis(DiP) 3'-pyrophosphohydrolas", dictionary$Gene , fixed=T ) | # description of spoT
                                 # *** From "PMC7356644" : In E. coli, synthesis of (p)ppGpp is carried out by RelA, during amino acid starvation, 
                                 # *** [...] and the bifunctional SpoT protein during other stress treatments [...]  SpoT is also responsible for degradation of (p)ppGpp .
                                 grepl("^dksA", dictionary$Gene ) |   # transcr fact, binds to RNApol augmentig ppGpp induced transcr initiation, see  www.ncbi.nlm.nih.gov/gene/944850 and https://journals.asm.org/doi/10.1128/iai.01439-09
                                 grepl("DksA", dictionary$Description ) |    # transcr fact, binds to RNApol augmentig ppGpp induced transcr initiation, see  www.ncbi.nlm.nih.gov/gene/944850 and https://journals.asm.org/doi/10.1128/iai.01439-09
                                 # c-di-GMP is a universal prokaryote intra-cellular messager signaling stress, e.g. see DOI: 10.1016/j.micres.2023.127302
                                 grepl("di-GMP|di-GMP|diGMP|[D-d]iguanylate", dictionary$Description) |  # see above
                                 grepl("cyclic diguanylic", dictionary$Description) |  # see above
                                 grepl("NrnA|NrnB|NrnC", dictionary$Description) |  # Oligoribonuclease which degrades pGpG (result of c-di-GMP degradation), which otherwise would inhibit c-di-GMP leading to its accumulation , see DOI:10.1002/2211-5463.13389
                                 grepl("^dgcB", dictionary$Gene) |  # diguanylate cyclase
                                 grepl("^yddV", dictionary$Gene) |  # diguanylate cyclase studied in E coli k-12, see doi:10.1002/2211-5463.13389
                                 grepl("YddV", dictionary$Description) |  # diguanylate cyclase studied in E coli k-12, see doi:10.1002/2211-5463.13389
                                 grepl("^ydeH", dictionary$Gene) |  # diguanylate cyclase studied in E coli k-12, see doi:10.1002/2211-5463.13389
                                 grepl("YdeH", dictionary$Description) |  # diguanylate cyclase studied in E coli k-12, see doi:10.1002/2211-5463.13389
                                 # Clp proteases are ENERGY DEPENDENT (ATP) and also involved in starvation REGULATION, as "Strains lacking ClpP failed to increase degradation of starvation proteins when glucose was added" see DOI: 10.1128/jb.175.1.53-63.1993
                                 grepl("^clp" ,dictionary$Gene) |
                                 grepl("Caseinolytic proteas" ,dictionary$Description) |   # which means "Clp"
                                 grepl("Caseinolytic peptidase" ,dictionary$Description) |   # which means "Clp"
                                 grepl("Clp " ,dictionary$Description) |
                                 grepl("Clp$" ,dictionary$Description) |
                                 grepl(" Clp " ,dictionary$Description) ,
]

# The σS-controlled starvation response is triggered when the concentration of the carbon source in the media falls 
# below 10−7 M, and is absent in both nutrient-excess (e.g. early and mid-exponential growth phases) and 
# nutrient-limited (e.g. early stationary and prolonged stationary phases) conditions  https://doi.org/10.1099/mic.0.001195

# In particular, Clp degrades the sigma S (involved in signaling stress and starvation) and Carbon stavation related proteins, then allowing to change behaviour when needed, see DOI:10.1046/j.1365-2958.1999.01357.x and DOI: 10.1128/jb.178.2.470-476.1996
# Moreover, c-di-GMP should also increases ClpP, see DOI: 10.1016/j.csbr.2024.100023

# [...] protease called ClpX to degrade RpoS [92]. An adaptor protein called RssB is essential to identify RpoS and help
# in the function of the protease. RssB is the limiting factor for the degradation of RpoS, and its intracellular level is tightly
# controlled with the help of three anti-adaptor proteins – IraP, IraM and IraD [94, 95]. These anti-adaptors are induced in response
# to environmental perturbations like DNA damage and cold shock, primarily via the (p)ppGpp-mediated stringent response  https://doi.org/10.1099/mic.0.001195


starvation_genes <- starvation_genes[!grepl("Adineta|Rotaria|Diploscapter|Carchesium|Stentor|Nippostrongylus",starvation_genes$Taxon) , ]
starvation_genes <- starvation_genes[!grepl("TipF",starvation_genes$Description) , ]
starvation_genes <- starvation_genes[!grepl("T6SS",starvation_genes$Description) , ]
starvation_genes <- starvation_genes[!grepl("mitochondrial",starvation_genes$Description) , ]




############## SAVING THE FINAL OBJECT TO ANALYZE ... ##############

save( file= "Tabelle_RNA_PHA_after_filters.RData",
      list = c("table_abs", "table_gene", "mapping", "unmapped_count", 
               "dictionary", "metadata", "remaining_message", "proof_filter",
               "only_PHAgenes_IDs","PHA_IDs","starvation_genes"
               ) 
      )



######## EXTRA: COMPUTING THE """RAREFACTION""" CURVE OF RNAs ################
# NB: This is a very slow step on standard computers

options(scipen=100)
library("vegan")

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as generated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of the original samples (and in the same order!!!)
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    # check<-"red"
    check<-"black"
    #the slope is estimated (average) at the last b rarefaction points
    differ<- length(v)-t
    if(differ<0){differ<-0} # to avoid the error with negative values
    slope<-mean(diff(v[ (differ):length(v) ])/dx)
    if(slope<lim) {
      # check<-"blue"
      check<-"black"
      sat = sat+1
    }
    text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],
         col=check,pch=16,cex=cex,labels=names[i])
  }
}



table_temp <- table_abs
colnames(table_temp)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")
colnames(table_temp)<- gsub("2_81th day", "2\n81th day ", colnames(table_temp))
png(file="FASTQ_check/Rarefaction_curve_RNA.png",width=2100,height=1650, res=300)
r<-rarecurve(t(as(table_temp,"matrix")), 
             col="darkgray",
             step=5000,label=F, ylab = "Transcripts", xlab= "Reads amount")
evalslopes(r,colnames(table_temp),lim=0.001,cex=0.6)
dev.off()
rm(r)


# Again, but only the most abundant bacteria
only_these_ones <- dictionary[grepl("Thauera|Paracoccus", dictionary$Taxon) , ]
table_temp <- table_abs[only_these_ones$Entry , ]
table_temp <- cbind.data.frame(table_temp,Descr=only_these_ones$Description)
table_temp <- aggregate( . ~ Descr , table_temp, FUN=sum )
table_temp$Descr<-NULL 
colnames(table_temp)<- paste0(metadata$Reactor, "_", metadata$Experiment_day)
png(file="FASTQ_check/Rarefaction_curve_FUNCTIONS_OnlyThaueraParacoccus.png",width=2100,height=1650, res=300)
r<-rarecurve(t(as(table_temp,"matrix")), 
             col="darkgray",
             step=10000,label=F, ylab = "Transcripts", xlab= "Reads amount")
evalslopes(r,colnames(table_temp),lim=0.001,cex=0.6)
dev.off()
rm(r)




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
