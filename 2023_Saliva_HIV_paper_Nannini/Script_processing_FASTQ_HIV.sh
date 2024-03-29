# operating system: Ubuntu 22.04 (Linux)
# required: conda, qiime2-2022.8, picrust-2.5, lefse-1.1.2, R-4.2


mkdir QIIME

######################## PROCESSING RAW FILES ########################

##### FASTQ CHECKSUM


md5sum raw_FASTQ_dir/* | sed 's+raw_FASTQ_dir/++g' | sed 's/ /\t/' > QIIME/checksum_FASTQ.tsv

dir raw_FASTQ_dir/* | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt


##### OPENING QIIME ENVIRONMENT

conda activate qiime2-2022.8

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_dir --input-format CasavaOneEightLanelessPerSampleDirFmt --output-path demux-paired-end.qza


##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 8 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

rm Log_di_cutadapt.txt

cat Log_reads_with_primers.txt

##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 275 --p-trunc-len-r 168 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 7

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

rm denoising-stats.qza demux-paired-end.qza
 

#### BAYESIAN CLASSIFICATION (using Scikit-learn retrained on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Scrivania/Tools/Database/SILVA_138_V3V4_BayesClassifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv


#### DE NOVO TREE ALIGNMENT

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza



#################### PICRUST2-LEFSE ##############################

#### PREPARING FILES FOR PICRUST 

mkdir PICRUST2_LEFSE

qiime tools export --input-path rep-seqs.qza --output-path PICRUST2_LEFSE # it will export "dna-sequences.fasta"
qiime tools export --input-path table.qza --output-path PICRUST2_LEFSE # it will export "feature-table.biom"

conda deactivate

cd PICRUST2_LEFSE


##### PICRUST2

conda activate picrust2

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 4 --verbose -e 0.0 -t sepp 
# the last two argument after verbose are for SEPP algorithm ("e" needed to avoid and error, see https://github.com/gavinmdouglas/q2-picrust2/issues/13)

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m metacyc -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

conda deactivate 


########## RENAMING FEATURES FOR LEFSE PLOT

R # open R

a <- read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
colnames(a)[1:4]

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]
colnames(a)<-gsub(".","-", colnames(a), fixed=T)

ID <-read.table(file="../../metadata_stool_saliva.csv", sep=",", header=T)
row.names(ID)<- ID$FASTQ_name # the actual names of the samples
original_length<-dim(ID)[1]

ID <- ID[colnames(a), ] # reordering
identical(dim(ID)[1], original_length)    # TRUE

colnames(a) <- ID$Sample_ID

a<-rbind.data.frame(ID$Condition,ID$Responders,a)
a<-cbind.data.frame(c("Condition","Responders",Descriptions$pathway),a)
colnames(a)[colnames(a)=="c(\"Condition\", \"Responders\", Descriptions$pathway)"]<-"Sample"
head(a, n=3)

original_a <- a   # to re-subset later (first H vs HIV, then again but with Responders yes vs no)
rm(a)


### HEALTHY vs HIV

a <- original_a[ -2 , original_a[1, ] != "Treatment" ] # removing the treatment (2o row) to focus on Healthy vs HIV
head(a,n=3)

Saliva<- ID[ ID$Sample_Type=="Saliva" & ID$Condition!="Treatment", "Sample_ID" ]

write.table(a[, c("Sample",Saliva)] ,file="Output_for_LefSe_Healthy_vs_HIV_SALIVA.tsv",row.names = F,quote = F, sep="\t")

rm(a, Saliva)


### RESPONDERS YES vs NOT RESPONDERS

a <- original_a[ , original_a[1, ] %in% c("Condition","Treatment") ] # removing the conditions (no healthy) but the samples column is the demarked thrugh  "condition" in this row
a <- a[ -1, ]    # the first row are the conditions, removing them

a[1,] <- gsub("1","Responders",as.character(a[1,]))   # the new first row are the responding conditions --> changing numbers with name
a[1,] <- gsub("0","Not responders",as.character(a[1,]))
head(a, n=3)

Saliva<- ID[ ID$Sample_Type=="Saliva" & ID$Condition=="Treatment", "Sample_ID" ]

write.table(a[, c("Sample",Saliva)] ,file="Output_for_LefSe_Responders_vs_Not_SALIVA.tsv",row.names = F,quote = F, sep="\t")

rm(a, Saliva)



quit(save="no")



######### LEFSE

conda activate lefse

# groups is second row (then -c2) and normalization to 1 milion (see Galaxy recommendations)


### for saliva HEALTHY_vs_HIV

lefse_format_input.py 'Output_for_LefSe_Healthy_vs_HIV_SALIVA.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000
lefse_run.py formatted_table_for_lefse.in -y 1 Result_LEFSE_Healthy_vs_HIV_SALIVA.res # y is for multiclass analysis

rm formatted_table_for_lefse.in


### for saliva RESPONDERS

lefse_format_input.py 'Output_for_LefSe_Responders_vs_Not_SALIVA.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 
lefse_run.py formatted_table_for_lefse.in -y 1 Result_LEFSE_Responders_vs_Not_SALIVA.res # y is for multiclass analysis

rm formatted_table_for_lefse.in


conda deactivate

cd ..

exit
