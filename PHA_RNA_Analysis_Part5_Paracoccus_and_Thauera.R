################ PREPARING THE ENVIRONMENT #################

{
  library("ggplot2")
  library("ggvenn")
  library("Hmisc")
  library("reshape2")
  library("ggh4x")
  library("ecodist")
  suppressPackageStartupMessages( library("vegan") ) # it has a verbose loading!
}

options(scipen=100)

suppressWarnings( dir.create("RNA_PHA_Results") )
dir.create("RNA_PHA_Results/Most_abundant_RNAs")

load(file="Tabelle_RNA_PHA_after_filters.RData")



#### COMPUTING THE PROPORTIONS

if(! "proof_filter" %in% ls()){
  stop("\n Wait! Did you perform the filtering step??? \n\n")
}

# Proportions computed after the filters and with the unmapped
prop_with_unmap <- apply(table_abs, MARGIN = 2, function(x) x/sum(x))
prop_with_unmap <- prop_with_unmap * 100

# Proportions computed after the filters BUT without the unmapped (for statistical tests or alike)
prop_with_NO_unmap <- apply(table_abs[! row.names(table_abs) %in% "UNMAPPED", ], MARGIN = 2, function(x) x/sum(x))
prop_with_NO_unmap <- prop_with_NO_unmap * 100

# Gene proportions
prop_GENE_with_unmap <- apply(table_gene, MARGIN = 2, function(x) x/sum(x))
prop_GENE_with_unmap <- as.data.frame(prop_GENE_with_unmap) * 100




#################### ADJUSTING THE DESCRIPTIONS ######################

infos <- dictionary$Description
{
  infos<- gsub("ATPase YchF","ATPase",infos, fixed=T)
  infos<-gsub("coenzyme A","CoA",infos,fixed=T)
  infos<-gsub("-coA ","-CoA",infos,fixed=T)
  # Clp proteins
  infos<- gsub("ATP-dependent","ATP-dep",infos)
  infos<- gsub("^Clp ","ATP-dep Clp ",infos)
  infos<- gsub(" prot Cl"," Cl",infos)
  infos<- gsub("receiver modulated","",infos)
  infos<- gsub("phosphodiesterase","PDE",infos)
  infos<- gsub("diguanylate cyclase","diG-cyclase",infos)
  infos<- gsub("Clp protease N-term domain-containing prot","ATP-dep protease Clp N-term domain-containing",infos)
  infos<- gsub("Clp protease ","protease ",infos)
  infos<- gsub("ATP-dep protease adapter ClpS", "ATP-dep protease ClpS", infos, fixed=T)  
  infos<- gsub("Universal stress protein Usp protein ", "Universal stress protein ", infos, fixed=T)  
  infos<- gsub("Endopeptidase Clp regulatory (ClpX)", "ATP-dep protease ClpX", infos, fixed=T)  
  infos<- gsub("Endopeptidase Clp regulatory ClpX", "ATP-dep protease ClpX", infos, fixed=T)  
  infos<- gsub(", clpA", " ClpA", infos, fixed=T)
  # Various related to starvation or stationary phase ...
  infos<- gsub("PAS domain-containing prot","PAS domain",infos)
  infos<- gsub("PAC/Chase","PAC",infos , fixed=T)
  infos<- gsub("PAC and GAF","PAC-GAF",infos , fixed=T)
  infos<- gsub("H-NS family nucleoid-associated regulatory$", "H-NS histone family", infos)  
  infos<- gsub("Nucleoid-structuring H-NS$", "H-NS histone family", infos)  
  infos<- gsub("H-NS histone$", "H-NS histone family", infos)  
  infos<- gsub("DNA-bind H-NS$", "H-NS histone family", infos)
  infos<- gsub("Bifunctional (p)ppGpp synthetase/guanosine-3,5-bis(diP) 3-pyrophosphohydrolase", "(P)ppGpp bifunctional synthase/hydrolase", infos, fixed = T)
  infos<- gsub("Bifunctional (p)ppGpp synthetase/guanosine-3',5'-bis(diP) 3'-pyrophosphohydrolase", "(P)ppGpp bifunctional synthase/hydrolase", infos, fixed = T)
  infos<- gsub("Bifunctional (P)ppGpp synthetase/guanosine-3,5-bis(DiP) 3-pyrophosphohydrolase", "(P)ppGpp bifunctional synthase/hydrolase", infos, fixed = T)
  infos<- gsub("GTP pyrophosphokinase, (P)ppGpp synthetase$", "(P)ppGpp synthetase", infos)
  infos<- gsub("GTP pyrophosphokinase rsh$", "(P)ppGpp synthetase", infos) # pyrophosphokinase = synthetase, see https://www.uniprot.org/uniprotkb/P0AG20/entry, 
  infos<- gsub("GTP pyrophosphokinase$", "(P)ppGpp synthetase", infos)  # pyrophosphokinase = synthetase, see https://www.uniprot.org/uniprotkb/P0AG20/entry
  infos<- gsub("TraR/DksA family transcript regulator", "RNA pol-bind transcription factor DksA", infos)
  infos<- gsub("RNA pol-bind DksA", "RNA pol-bind transcription factor DksA", infos)
  infos<- gsub("DNA-bind transcript regulator, Lrp family", "DNA-bind Lrp family transcript regulator", infos)
  infos<- gsub("Transcript regulator, AsnC/Lrp family", "DNA-bind Lrp family transcript regulator", infos)
  infos<- gsub("Universal stress protein Usp ", "Universal stress protein ", infos, fixed=T)
  infos<- gsub("Carbon starvation CstA.*", "Carbon starvation CstA", infos)  
  infos<- gsub("Carbon starvation A.*", "Carbon starvation CstA", infos)  
  infos<- gsub("Bifunctional oligoribonuclease and PAP phosphatase NrnA", "Bifunctional oligoribonuclease/PAP phosphatase NrnA", infos)  
  infos<- gsub("ATP-dep ATP-dep protease Clp N-term", "ATP-dep protease Clp N-term", infos)  
  infos<- gsub(", MazF antagonist", "", infos)  
  infos<- gsub("diguanylate cyclase","diG-cyclase",infos)
  infos<- gsub("sensor-containing diguanylate cyclase","sensor + diguanylate cyclase",infos)
  # Heat and Cold shocks
  infos<- gsub("Heat-", "Heat ", infos, fixed=T)
  infos<- gsub(" like prot", "-like", infos, fixed=T)
  infos<- gsub("HSP", "Hsp", infos, fixed=T)
  infos<- gsub("Molecular chaper ", "", infos, fixed=T)
  infos<- gsub("Chaperone DnaJ", "Hsp chaperone DnaJ", infos, fixed=T)
  infos<- gsub("^Hsp20 ", "Heat shock Hsp20 ", infos)
  infos<- gsub("^u", "U", infos)
  infos<- gsub("Chaperonin GroEL ", "Hsp60-system chaperonin GroEL ", infos, fixed=T)
  infos<- gsub("Co-chaperonin GroES ", "Hsp60-system cochaperonin GroES ", infos, fixed=T)
  infos<- gsub("Co-chaper GroES", "Hsp60-system cochaperonin GroES", infos, fixed=T)
  infos<- gsub("Co-chaper GroES", "Hsp60-system cochaperonin GroES", infos)
  infos<- gsub("^GroES$", "Hsp60-system cochaperonin GroES ", infos)
  # infos<- gsub("GroES-like", "Hsp60-system cochaperonin GroES ", infos)
  infos<- gsub(".*GroES L", "Hsp60-system cochaperonin GroES ", infos)
  infos<- gsub("^Chaperonin GroEL", "Hsp60-system chaperonin GroEL", infos)
  infos<- gsub("Chaperonin (groEL)", "Hsp60-system chaperonin GroEL", infos, fixed=T)
  infos<- gsub(".*chaperonin.*GroES", "Hsp60-system cochaperonin GroES ", infos)  
  infos<- gsub("^GroEL$", "Hsp60-system chaperonin GroEL", infos)
  infos<- gsub("Chaperone DnaK ", "Hsp70-system chaperone DnaK ", infos, fixed=T)
  infos<- gsub("Heat shock 70.*", "Hsp70-system chaperone DnaK", infos)  
  infos<- gsub(".*Heat shock 70.*", "Hsp70-system chaperone DnaK", infos)  
  infos<- gsub(".*70 kDa heat*", "Hsp70-system chaperone DnaK", infos)  
  infos<- gsub("18 kDa heat shock", "Heat shock Hsp18", infos)  
  infos<- gsub("Heat shock Hsp20.*", "Heat shock Hsp20", infos)  
  infos<- gsub("Ribosome-associated heat shock Hsp15", "Heat shock Hsp15", infos)
  infos<- gsub("GroES L", "GroES", infos)
  infos<- gsub("[C-c]old-shock", "Cold shock", infos)  
  infos<- gsub("[H-h]eat-shock", "Heat shock", infos) 
  infos<- gsub("[C-c]old shock prot, DNA bind$", "Cold shock family", infos)  
  infos<- gsub("[C-c]old-shock prot, DNA bind$", "Cold shock family", infos)  
  infos<- gsub("[C-c]old shock DNA-bind", "Cold shock", infos)  
  infos<- gsub("[C-c]old-shock DNA-bind", "Cold shock", infos)  
  # infos<- gsub("Heat inducible transcription repressor HrcA", "Hsp GroEL transcriptional repressor", infos, fixed=T)
  infos<- gsub("Heat inducible transcription repressor HrcA", "GroEL's transcriptional repressor", infos, fixed=T)
  infos<- gsub("Protein GrpE", "Hsp70-system cochaperone GrpE", infos, fixed=T)
  infos<- gsub("Protease HtpX homolog (htpX)", "HSP protease HtpX homolog", infos, fixed=T)  # yes, it is a HSP itself
  # General PHA names
  infos<-gsub("Poly(R)-hydroxyalkanoic acid","PHA",infos, fixed=T)
  infos<-gsub("poly(R)-hydroxyalkanoic acid","PHA",infos, fixed=T)
  infos<-gsub("PolyR-hydroxyalkanoic acid","PHA",infos, fixed=T)
  infos<- gsub("[P-p]oly-hydroxyalkanoic acid", "PHA", infos)
  infos<-gsub("[P-p]olyhydroxyalkanoic acid","PHA",infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("[P-p]olyhydroxyalkanoate", "PHA", infos)
  infos<- gsub("Poly-beta-hydroxyalkanoate","PHA",infos, fixed=T)
  infos<- gsub("Poly(3-hydroxyalkanoate)", "PHA", infos, fixed = T)
  infos<- gsub("poly(3-hydroxyalkanoate)", "PHA", infos, fixed = T)
  # PHB
  infos<- gsub("[P-p]oly-beta-hydroxybutyrate", "PHB", infos)
  infos<- gsub("Poly-beta-hydroxybutyrate", "PHB", infos)
  infos<- gsub("Poly(3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("poly(3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("Poly (3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("poly (3-hydroxybutyrate)", "PHB", infos, fixed = T)
  infos<- gsub("Polyhydroxybutyrate", "PHB", infos, fixed = T)
  # Phasin related
  infos<-gsub("Phasin fam prot","Phasin family",infos,fixed=T)
  infos<-gsub("[P-]phasin fam$","Phasin family",infos)
  infos<- gsub("Phasin prot","Phasin",infos)
  infos<- gsub("prot (PHASIN)","prot, phasin",infos, fixed=T)
  infos<- gsub("TIGR01841 fam phasin","Phasin family",infos)
  infos<- gsub("granule regulatory prot","phasin",infos , fixed=T)
  # PHA PHB Synthases
  infos<- gsub("PHA/PHB synth fam","PHA/PHB synthase",infos, fixed=T)
  infos<- gsub("3-hydroxyalkanoate synthetase","PHA/PHB synthase",infos, fixed=T) # see https://www.uniprot.org/uniprotkb/A0A133XN55/entry , search for the "InterPro" link
  infos<- gsub("PHA synth$","PHA synthase",infos)
  infos<- gsub("PHB synth$","PHB synthase",infos)
  infos<- gsub("PHA pol domain.*","PHA synthase",infos)
  infos<- gsub("PHA pol subunit.*","PHA synthase",infos)
  infos<- gsub("PHB pol family.*","PHB synthase",infos)
  infos<- gsub("PHB pol fam$","PHB synthase",infos)
  infos<- gsub("PHB pol domain.*","PHB synthase",infos)
  infos<- gsub("PHB pol subunit.*","PHB synthase",infos)
  infos<- gsub("PHB pol$","PHB synthase",infos)
  infos<- gsub("PHB pol N-term domain-containing","PHB synthase",infos)
  infos<- gsub("PhaM fam PHA granule multifunctional regulatory","PHB synthase activator",infos) #see doi: 10.1128/AEM.02935-13
  # Phosphate starvation
  infos<- gsub("Two-component system P regulon response regulator PhoB","P regulon transcript regulator PhoB",infos)
  infos<- gsub("P regulon transcript regulatory prot PhoB","P regulon transcript regulator PhoB",infos)
  infos<- gsub("P regulon sensor prot PhoR","P regulon sensor histidine kinase PhoR",infos)
  # Other enzymes
  infos<- gsub("Acetyl-CoA acetyltransf-1","Acetyl-CoA acetyltransferase",infos, fixed=T)
  infos<- gsub("Acetyl-CoA acetyltransfs","Acetyl-CoA acetyltransferase",infos, fixed=T)
  infos<- gsub("C-acetyltransferase","acetyltransferase",infos, fixed=T)
  infos<- gsub("AcetoAcetoacetyl","Acetoacetoacetyl",infos, fixed=T)
  infos<- gsub("(R)-specific enoyl-CoA hydratase","3-hydroxybutyryl-CoA dehydratase",infos, fixed=T)
  infos<- gsub("^[C-c]lass I ","",infos) # phaC
  infos<- gsub("synth, class I","synthase",infos, fixed=T) 
  infos<- gsub("synth class I","synthase",infos, fixed=T) 
  infos<- gsub("3-ketoacyl-CoA thiolase [fadN-fadA-fadE operon]","3-ketoacyl-CoA thiolase",infos, fixed=T) 
  infos<- gsub("[S-s]teroid 3-ketoacyl-CoA thiolase","3-ketoacyl-CoA thiolase",infos) 
  infos<- gsub("Beta-ketothiolase BktB","Beta-ketothiolase",infos) 
  infos<- gsub("Acyl-CoA dehydrogenaseFadE9","Acyl-CoA dehydrogenase",infos) 
  infos<- gsub("Short chain enoyl-CoA hydratase /3-hydroxyacyl-CoA dehydrogenase","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase/ enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T)  # fad J
  infos<- gsub("NAD-bind 3-hydroxyacyl-CoA dehydrogenase","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase NAD-bind.*","3-hydroxyacyl-CoA dehydrogenase",infos) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase NAD bind.*","3-hydroxyacyl-CoA dehydrogenase",infos) 
  infos<- gsub("(3R)-3-hydroxyacyl-CoA","3-hydroxyacyl-CoA",infos, fixed = T) 
  infos<- gsub("3-hydroxyacyl-CoA dehydrogenase type-2","3-hydroxyacyl-CoA dehydrogenase",infos, fixed = T) 
  infos<- gsub("(2S)-methylsuccinyl-CoA dehydro","Methylsuccinyl-CoA dehydro",infos, fixed = T) 
  infos<- gsub("DMT fam prot","DMT fam transp",infos, fixed = T) 
  # focus on depolimerases
  infos<- gsub("^[E-e]sterase/P","P",infos) 
  infos<- gsub("^[E-e]sterase P","P",infos) 
  infos<- gsub("^Esterase, P","P",infos) 
  infos<- gsub("depolimerase fam esterase$","depolimerase",infos) 
  infos<- gsub("depolimerase fam esterase ","depolimerase",infos) 
  infos<- gsub("depolimerase fam$","depolimerase",infos) 
  infos<- gsub("PHB depolimerase family","PHB depolimerase",infos) 
  infos<- gsub("depolimerase, intracellular","depolimerase",infos) 
  infos<- gsub("depol, intracellular","depolimerase",infos)
  # Removing sporadic extras
  infos<-gsub("alpha","α",infos,fixed=T)
  infos<- gsub(" prot$", "", infos)
  infos<- gsub("subunit","",infos)
  infos<- gsub("^[P-p]robable ", "", infos)  
  infos<- gsub("[P-p]utative ", "", infos)  
  infos<- gsub(", putative", "", infos, fixed=T)  
  infos<- gsub("'", "", infos, fixed=T)  
  infos<- gsub(" $", "", infos)  
  infos<- Hmisc::capitalize(infos)
}

updated_dictionary <- dictionary
updated_dictionary$Description <- infos



### few descriptions are still not unique for the same gene, thus I'm solving them manually according to their Description and/or gene name..
{
  # CLP
  updated_dictionary[ updated_dictionary$Entry %in% c("A0A037ZLN2","A0A1F3NAW7","UPI00048CEB86","UPI0009627D53") , "Description" ] <- "ATP-dep protease ClpP"
  updated_dictionary[ updated_dictionary$Entry %in% c("A0A7S6PV13","A0A4P5XDK1","A0A0D3MKN6","A0A077KDY8") , "Description" ] <- "ATP-dep protease ClpC"
  updated_dictionary[ grepl("ClpB",updated_dictionary$Description) , "Description" ] <- "ATP-dep protease ClpB"
  updated_dictionary[ grepl("^clpB",updated_dictionary$Gene) , "Description" ] <- "ATP-dep protease ClpB"
  updated_dictionary[ grepl("ClpA",updated_dictionary$Description) , "Description" ] <- "ATP-dep protease ClpA"
  updated_dictionary[ grepl("^clpA",updated_dictionary$Gene) , "Description" ] <- "ATP-dep protease ClpA"
  updated_dictionary[ grepl("ClpS",updated_dictionary$Description) , "Description" ] <- "ATP-dep protease adaptor ClpS"
  updated_dictionary[ grepl("^clpS",updated_dictionary$Description) , "Description" ] <- "ATP-dep protease adaptor ClpS"
  updated_dictionary[ grepl("ClpXP protease specificity-enhancing factor SspB",updated_dictionary$Description) , "Description" ] <- "ClpXP protease specificity-enhancing factor"
  updated_dictionary[ grepl("^clpP",updated_dictionary$Gene) , "Description" ] <- "ATP-dep protease ClpP"
  updated_dictionary[ grepl("^clpX",updated_dictionary$Gene) , "Description" ] <- "ATP-dep protease ClpX"
  updated_dictionary[ grepl("CplX",updated_dictionary$Description) , "Description" ] <- "ATP-dep protease ClpX"
  # HSP 
  updated_dictionary[ updated_dictionary$Entry %in% c("A0A0F2RJC5","UPI00041397C6") , "Description"] <- "Hsp60-system chaperonin GroEL"
  updated_dictionary[ updated_dictionary$Entry %in% c("A0A0B1RY98","A0A0K6HNJ5","A0A7J5V3I3","UPI00036A7D0A") , "Description"] <- "Hsp70-system chaperone DnaK"
  updated_dictionary[ updated_dictionary$Gene %in% c("dnak","dnaK") , "Description"] <- "Hsp70-system chaperone DnaK"
  updated_dictionary[ grepl("Hsp70",updated_dictionary$Description) , "Description"] <- "Hsp70-system chaperone DnaK"
  updated_dictionary[ grepl("^hsp70",updated_dictionary$Gene) , "Description"] <- "Hsp70-system chaperone DnaK"
  updated_dictionary[ grepl("^HSP70",updated_dictionary$Gene) , "Description"] <- "Hsp70-system chaperone DnaK"
  updated_dictionary[ grepl("Dna[J-j]",updated_dictionary$Description) , "Description"] <- "Hsp chaperone DnaJ"
  updated_dictionary[ grepl("^dnaJ",updated_dictionary$Description) , "Description"] <- "Hsp chaperone DnaJ"
  updated_dictionary[ grepl("dnaK",updated_dictionary$Gene) , "Description"] <- "Hsp70-system chaperone DnaK"
  updated_dictionary[ grepl("IbpA",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp20 family molecular chaper IbpA"
  updated_dictionary[ grepl("^ibpA",updated_dictionary$Gene) , "Description"] <- "Heat shock Hsp20 family molecular chaper IbpA"
  updated_dictionary[ grepl("Hsp20",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp20 family molecular chaper IbpA"
  updated_dictionary[ grepl(".*Hsp33.*",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp33 HslO"
  updated_dictionary[ grepl(".*Hsp90*",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp90"
  updated_dictionary[ grepl(".*hsp-90*",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp90"
  updated_dictionary[ grepl(".*hsp-90*",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp90"
  updated_dictionary[ grepl("Heat shock 90",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp90"
  updated_dictionary[ grepl("Heat shock 90",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp90"
  updated_dictionary[ grepl("Co-chaper GroES",updated_dictionary$Description) , "Description"] <- "Heat shock Hsp90"
  updated_dictionary[ grepl("Alcohol dehydrog GroES",updated_dictionary$Description) , "Description"] <- "Hsp90 GroES containing alcohol dehydrogenase"
  updated_dictionary[ grepl("Chaperonin GroEL",updated_dictionary$Description) , "Description"] <- "Hsp60-system chaperonin GroEL"
  updated_dictionary[ grepl("Molecular chaper GroEL",updated_dictionary$Description) , "Description"] <- "Hsp60-system chaperonin GroEL"
  updated_dictionary[ grepl("^hsp82",updated_dictionary$Gene) , "Description"] <- "Heat shock Hsp82"
  # c di GMP
  which_ones <- grepl("[D-d]iguanylate|diG-|DiG-cyclase|[D-d]iG cyclase|c-di-GMP|[C-c]yclic-di-GMP-|Cyclic di-GMP",updated_dictionary$Description) & grepl("cyclase",updated_dictionary$Description) & grepl("PDE|phosphodiester",updated_dictionary$Description)
  updated_dictionary[ which_ones , "Description"] <- "Diguanylate (c-di-GMP) cyclase/PDE"
  which_ones <- grepl("[D-d]iguanylate|diG-|DiG-cyclase|[D-d]iG cyclase|c-di-GMP|[C-c]yclic-di-GMP-|Cyclic di-GMP",updated_dictionary$Description) & grepl("cyclase",updated_dictionary$Description) & !grepl("PDE|phosphodiester",updated_dictionary$Description)
  updated_dictionary[ which_ones , "Description"] <- "Diguanylate (c-di-GMP) cyclase"
  which_ones <- grepl("[D-d]iguanylate|c-di-GMP|[C-c]yclic-di-GMP-|Cyclic di-GMP",updated_dictionary$Description) & grepl("PDE|phosphodiester",updated_dictionary$Description) & !grepl("cyclase",updated_dictionary$Description)
  updated_dictionary[ which_ones , "Description"] <- "Diguanylate (c-di-GMP) PDE"
  # Cold Shock
  updated_dictionary[ grepl("CspA",updated_dictionary$Description) , "Description"] <- "Cold shock CspA"
  updated_dictionary[ grepl("^cspA",updated_dictionary$Gene) , "Description"] <- "Cold shock CspA"
  updated_dictionary[ grepl("^cspC",updated_dictionary$Gene) , "Description"] <- "Cold shock CspC"
  # DPS
  updated_dictionary[ grepl("Dps|DPS",updated_dictionary$Description) & grepl("bind",updated_dictionary$Description) , "Description"] <- "DPS (DNA protection during starvation)"
  updated_dictionary[ grepl("^dps",updated_dictionary$Gene) , "Description"] <- "DPS (DNA protection during starvation)"
  updated_dictionary[ grepl("DNA starvation/stationary phase protection",updated_dictionary$Description) , "Description"] <- "DPS (DNA protection during starvation)"
  # Usp
  updated_dictionary[ grepl("UspA",updated_dictionary$Description) , "Description"] <- "Universal stress protein UspA"
  updated_dictionary[ grepl("^uspF",updated_dictionary$Gene) , "Description"] <- "Universal stress protein UspF"
  updated_dictionary[ grepl("UspE",updated_dictionary$Description) , "Description"] <- "Universal stress protein UspE"
  # P or C specific starvation
  updated_dictionary[ grepl("CstA",updated_dictionary$Description) , "Description"] <- "Carbon starvation CstA fam"
  updated_dictionary[ grepl("^Carbon starvation prot A$",updated_dictionary$Description) , "Description"] <- "Carbon starvation CstA fam"
  updated_dictionary[ grepl("PhoH|PsiH",updated_dictionary$Description) , "Description"] <- "P-starvation-inducible PsiH family" # They are the same https://pubmed.ncbi.nlm.nih.gov/8444794/
  updated_dictionary[ grepl("P-starvation-inducible E",updated_dictionary$Description) , "Description"] <- "P-starvation-inducible PsiE"
  updated_dictionary[ grepl("P-starvation-inducible PsiE",updated_dictionary$Description) , "Description"] <- "P-starvation-inducible PsiE family"
  # ppGpp
  updated_dictionary[ grepl("SpoT",updated_dictionary$Description) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase"
  updated_dictionary[ grepl("^spoT",updated_dictionary$Gene) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase"
  updated_dictionary[ grepl("RelA| rsh",updated_dictionary$Description) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase" # yes, double function also for relA
  updated_dictionary[ grepl("^relA",updated_dictionary$Gene) , "Description"] <- "(P)ppGpp bifunctional synthase/hydrolase" # yes, double function also for relA
  updated_dictionary[ grepl("GTP pyrophosphokinase",updated_dictionary$Gene) , "Description"] <- "(P)ppGpp synthetase" # yes, double function also for relA
  # DksA
  updated_dictionary[ grepl("^dksA",updated_dictionary$Gene) , "Description"] <- "RNA pol-bind transcription factor DksA"
  # Pha/phB
  updated_dictionary[ grepl("[P-p]haC",updated_dictionary$Description) , "Description"] <- "PHA synthase"
  updated_dictionary[ grepl("[P-p]haC",updated_dictionary$Description) , "Description"] <- "PHA synthase"
  updated_dictionary[ grepl("^[P-p]haC",updated_dictionary$Gene) , "Description"] <- "PHA synthase"
  updated_dictionary[ grepl("[P-p]hbB",updated_dictionary$Description) , "Description"] <- "Acetoacetyl-CoA reductase"
  updated_dictionary[ grepl("^phbB",updated_dictionary$Gene) , "Description"] <- "Acetoacetyl-CoA reductase"
  updated_dictionary[ grepl("^phaJ",updated_dictionary$Gene) , "Description"] <- "3-hydroxybutyryl-CoA dehydratase"
  updated_dictionary[ grepl("^phaR",updated_dictionary$Gene) , "Description"] <- "PHA synthesis repressor"
  updated_dictionary[ grepl("PhaR",updated_dictionary$Description) , "Description"] <- "PHA synthesis repressor"
  updated_dictionary[ grepl("PhbR",updated_dictionary$Description) , "Description"] <- "PHA synthesis repressor"
  # Phasins
  updated_dictionary[ grepl("^phaP",updated_dictionary$Gene) , "Description"] <- "Phasin"
  updated_dictionary[ grepl("^phbP",updated_dictionary$Gene) , "Description"] <- "Phasin"
  updated_dictionary[ grepl("[P-p]haP",updated_dictionary$Description) , "Description"] <- "Phasin"
  updated_dictionary[ grepl("PHA granule-associated",updated_dictionary$Description) , "Description"] <- "Phasin"
  updated_dictionary[ grepl("PHA-granule associated",updated_dictionary$Description) , "Description"] <- "Phasin"
  updated_dictionary[ grepl("Granule-associated prot",updated_dictionary$Description) , "Description"] <- "Phasin"
  # Metabolism
  updated_dictionary[ grepl("^fadD",updated_dictionary$Gene) , "Description"] <- "Long-chain acyl-CoA synthetase"
  updated_dictionary[ grepl("^acsA",updated_dictionary$Gene) , "Description"] <- "Acetyl-CoA synthetase"
  updated_dictionary[ grepl("AcsA",updated_dictionary$Description) , "Description"] <- "Acetyl-CoA synthetase"
  updated_dictionary[ grepl("^paaC",updated_dictionary$Gene) , "Description"] <- "3-hydroxyacyl-CoA dehydrogenase"
  updated_dictionary[ grepl("PaaC",updated_dictionary$Description) , "Description"] <- "3-hydroxyacyl-CoA dehydrogenase"
  updated_dictionary[ grepl("PaaH",updated_dictionary$Description) , "Description"] <- "3-hydroxyacyl-CoA dehydrogenase" # yes, the same as PaaC, see Uniprot descriptions
}

# View(updated_dictionary[!duplicated(updated_dictionary$Description) , ])




########### MOST ABUNDANT RNA OF PARACOCCUS * GLOMMING DESCRIPTIONS * ###############

# table<-table_abs
table <- as.data.frame( prop_with_unmap[ updated_dictionary$Entry[grepl("Paracoccus|Paracoccaceae",dictionary$Taxon)&
                                                                    !grepl("phage|virus|virales",updated_dictionary$Taxon )] , ] ) 
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

### Aggregating using the description
selected_infos<-updated_dictionary[updated_dictionary$Entry %in% row.names(table) , ]
row.names(selected_infos)<- selected_infos$Entry
selected_infos<- selected_infos[ row.names(table) ,  ]
selected_infos$Description[selected_infos$Description=="Uncharact"]<-paste( selected_infos$Description[selected_infos$Description=="Uncharact"], selected_infos$Entry[selected_infos$Description=="Uncharact"] )
# infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<-selected_infos$Description
infos<- gsub("Putative ","", infos)
infos<- gsub("Probable ","", infos)
table$Protein <- infos
table<-aggregate( .~Protein, table, FUN=sum )  # Protein = both function and taxon (species)
table <- table[!table$Protein%in%c("Hypothetical prot","Hypothetical","hypothetical"), ]  # they are of no use in this context ...
table <- table[!table$Protein%in%c("Phage tail"), ]  
top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", colnames(table)!="Protein" ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]

# plotting ...
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# aesthetic modifications
# table_to_plot$Protein <- gsub("Paracoccus","P.", table_to_plot$Protein)
table_to_plot$Protein <- gsub("dehydrog","dehydrogenase", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("transp","transporter", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("transf","transferase", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("Uncharact","Uncharacterized", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub(" prot$","", table_to_plot$Protein)
table_to_plot$Protein <- gsub("bind$","binding protein", table_to_plot$Protein)
table_to_plot$Protein <- gsub("fam","family", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("reduct","reductase", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("CsbD family","CsbD family (stress response)", table_to_plot$Protein, fixed = T)  # https://www.ncbi.nlm.nih.gov/Structure/cdd/cl22912
table_to_plot$Protein <- gsub("membrane","membrane protein", table_to_plot$Protein, fixed = T)  # https://www.ncbi.nlm.nih.gov/Structure/cdd/cl22912
table_to_plot$Protein <- Hmisc::capitalize(table_to_plot$Protein)

custom_colors<-c("chartreuse","deeppink", "firebrick3",
                 "orange","darkmagenta", 
                 "darkblue","grey65", "cadetblue",
                 "slateblue3","black","cyan","brown4",
                 "greenyellow", "darkslategray3",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue",
                  "chartreuse2","yellow4","yellow", "chocolate4","coral2","darkgreen",
                 "slateblue1", "orangered","springgreen3")
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9.5) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=custom_colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.4, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.065, "cm"),
        legend.key.width = unit(0.34, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        # legend.title = element_text ( size = 9 ),
        strip.text.x=element_text(size=9,colour="black", 
                                  lineheight = 1,
                                  margin = margin(2,3,2,3, "pt")
        ),
        strip.background = element_rect(color = "black", linewidth = 0.45),
        legend.text = element_text ( size = 7.7 ),
        legend.margin = margin(-18.5,0,2,-6),
        legend.position="bottom",
        plot.margin = margin(3,1,0.5,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       # title = "Most abundant mRNA of Paracoccus spp."
  )
ggsave("RNA_PHA_Results/Most_abundant_RNAs/ParacoccusFunctions_Most_abundant.png", width= 5.45, height = 4.65, dpi=300)




########### MOST ABUNDANT RNA OF THAUERA * GLOMMING DESCRIPTIONS * ###############

# table<-table_abs
table <- as.data.frame( prop_with_unmap[ updated_dictionary$Entry[grepl("Thauera",updated_dictionary$Taxon)&
                                                                    !grepl("phage|virus|virales",updated_dictionary$Taxon )] , ] ) 
colnames(table)<- paste0(metadata$Reactor, "_", metadata$Experiment_day, "th day")

### Aggregating using the description
selected_infos<-updated_dictionary[updated_dictionary$Entry %in% row.names(table) , ]
row.names(selected_infos)<- selected_infos$Entry
selected_infos<- selected_infos[ row.names(table) ,  ]
selected_infos$Description[selected_infos$Description=="Uncharact"]<-paste( selected_infos$Description[selected_infos$Description=="Uncharact"], selected_infos$Entry[selected_infos$Description=="Uncharact"] )
# infos<-paste0( selected_infos$Description," (", selected_infos$Taxon,")" )
infos<-selected_infos$Description
infos<- gsub("Putative ","", infos)
infos<- gsub("Probable ","", infos)
infos<- gsub("PEP-CTERM sorting","PEP-CTERM prot-sorting", infos, fixed=T) # otherwise same description aggregated in two different observations!
table$Protein <- infos
table <- table[!table$Protein%in%c("Hypothetical prot","Hypothetical","hypothetical"), ]  # they are of no use in this context ...
table <- table[!table$Protein%in%c("Phage tail"), ]  
table<-aggregate( .~Protein, table, FUN=sum )  # Protein = both function and taxon (species)
table <- table[table$Protein!="Hypothetical prot", ]  # they are of no use in this context ...
top <- names(sort(rowSums(table[!row.names(table) %in% "UNMAPPED", colnames(table)!="Protein" ]), decreasing=TRUE))[1:30]
table_top<- table[top, ]

# plotting ...
only_them<-table_top
only_them[ , colnames(only_them)!="Protein"] <- apply( only_them[ , colnames(only_them)!="Protein"],
                                                       MARGIN = 2, function(x) x/sum(x) )
only_them[ , colnames(only_them)!="Protein"] <- only_them[ ,colnames(only_them)!="Protein"] * 100
table_to_plot<-melt(only_them, id.vars = "Protein")
table_to_plot$Reactor <- gsub("_.*", "", table_to_plot$variable)
table_to_plot$variable <- gsub(".*_", "", table_to_plot$variable)
# aesthetic modifications
table_to_plot$Protein <- gsub(", methanol/ethanol fam",", meth/eth fam", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("family outer membrane","prot on outer membrane", table_to_plot$Protein, fixed = T)
# table_to_plot$Protein <- gsub("Thauera","P.", table_to_plot$Protein)
# table_to_plot$Protein <- gsub("P. sp.","Thauera sp.", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("dehydrog","dehydrogenase", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("transp","transporter", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("transf","transferase", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("Uncharact","Uncharacterized", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub(" prot$","", table_to_plot$Protein)
table_to_plot$Protein <- gsub("attachment$","attachment protein", table_to_plot$Protein)
table_to_plot$Protein <- gsub("Transmembrane$","Transmembrane protein", table_to_plot$Protein)
table_to_plot$Protein <- gsub("bind$","binding protein", table_to_plot$Protein)
table_to_plot$Protein <- gsub("fam$","family", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("reduct","reductase", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- gsub("  "," ", table_to_plot$Protein, fixed = T)
table_to_plot$Protein <- Hmisc::capitalize(table_to_plot$Protein)

custom_colors<-c("chartreuse","grey65","deeppink",
                 "orange","darkmagenta", 
                 "navyblue","black","cyan", "cadetblue",
                 "slateblue3","gold2","brown4",
                 "greenyellow", "springgreen3",
                 "deepskyblue2","darkgreen","orange3",
                 "darkorchid", "lightblue1" , "violet", "blue1",
                 "yellow", "chartreuse3","yellow4", "chocolate4","coral2","darkgreen",
                 "slateblue1", "red","darkslategray3")
ggplot(data=table_to_plot, aes(x=variable, y=value, fill=Protein)) +
  geom_bar(stat="identity", position="stack", linewidth=0.5) +
  theme_classic(base_size =9.5) +
  facet_grid2( . ~ variable+Reactor,
               scales = "free",
               space="free",
               switch = "x", # This places the grid at the bottom (as x axis)
               strip = strip_nested(size="constant"))+
  scale_fill_manual(values=custom_colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(angle=0, vjust=0.4, hjust=0.5, size=6.5),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.075, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        # legend.title = element_text ( size = 9 ),
        strip.text.x=element_text(size=9,colour="black", 
                                  lineheight = 1,
                                  margin = margin(2,3,2,3, "pt")
        ),
        strip.background = element_rect(color = "black", linewidth = 0.45),
        legend.text = element_text ( size = 7.7 ),
        legend.margin = margin(-18.5,0,2,-5.5),
        legend.position="bottom",
        plot.margin = margin(3,1,0.5,1)
  ) +
  scale_x_discrete (expand = c(0.01,0) ) +
  scale_y_continuous(expand = c(0.01,0)) +  #  bars closer to the x axis
  guides(fill=guide_legend(nrow=15)) +
  labs(x="", y="Percentages in samples\n(excluding other RNAs)",
       fill="",
       # title = "Most abundant mRNA of Thauera spp."
  )
ggsave("RNA_PHA_Results/Most_abundant_RNAs/ThaueraFunctions_Most_abundant.png", width= 5.45, height = 4.65, dpi=300)

