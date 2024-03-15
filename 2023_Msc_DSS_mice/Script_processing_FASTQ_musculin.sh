# requirements: linux bash, qiime, picrust2, lefse, R, raw FASTQ in "FASTQ_iniziali", Metadata.xlsx, RDP classifier


############################ PROCESSING THE RAW FASTQ FILES ########################

conda activate qiime2-2021.4

mkdir QIIME
cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../FASTQ_iniziali --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

rm Log_di_cutadapt.txt

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # del totale!

echo -e "\n to see quality put that demux.qzv in this site https://view.qiime2.org/ \n"

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 283 --p-trunc-len-r 168 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza
 
### bayesian classification (simil RDP retrained on V3-V4 region of SILVA 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Scrivania/Tools/Databases/SILVA_138_V3V4_BayesClassifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv


conda deactivate

# the files are now ready for R (see the relative stand alone R script of this analysis)


##################### PREPARING (SUBSETTING) FILES FOR PICRUST 2-LEFSE #####################


###### open R environment to select and copy

R

setwd("FASTQ_iniziali")

Initial<-length(dir())
Elenco<-as.data.frame(dir())
head(Elenco$'dir()')
Elenco$Code<-substr( Elenco$'dir()', start=16, stop = 17 )
Elenco$Code<-gsub("_","",Elenco$Code, fixed=T)
Elenco$Code
Elenco<-Elenco[Elenco$'dir()'!="Script per rimuovere specifici campioni.R", ]

Sample_Data <- readxl::read_excel("../Metadata.xlsx")
Subset<-subset(Sample_Data, Time=="T0")

dir.create("../Only_T0")
Elenco_1<-Elenco[Elenco$Code %in% Subset$FASTQ_ID,]
file.copy(from=Elenco_1$`dir()`, to=paste("../Only_T0",Elenco_1$`dir()`,sep="/") )

dir.create("../Only_T1")
Elenco_2<-Elenco[! Elenco$Code %in% Subset$FASTQ_ID,]
file.copy(from=Elenco_2$`dir()`, to=paste("../Only_T1",Elenco_2$`dir()`,sep="/") )

quit(save = "no") # Close R back to normal bash


####### open qiime again to create the files of corrispective ASV (same settings as whole data)


mkdir PICRUST2_LEFSE

conda activate qiime2-2021.4

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Only_T0 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 283 --p-trunc-len-r 168 --o-table table_T0.qza --o-representative-sequences rep-seqs_T0.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

rm trimmed-seqs.qza demux-paired-end.qza denoising-stats.qza

qiime tools export --input-path rep-seqs_T0.qza --output-path PICRUST2_LEFSE
qiime tools export --input-path table_T0.qza --output-path PICRUST2_LEFSE

mv PICRUST2_LEFSE/dna-sequences.fasta PICRUST2_LEFSE/dna-sequences_T0.fasta
mv PICRUST2_LEFSE/feature-table.biom PICRUST2_LEFSE/feature-table_T0.biom

rm rep-seqs_T0.qza table_T0.qza

rm Only_T0 -R

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Only_T1 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 283 --p-trunc-len-r 168 --o-table table_T1.qza --o-representative-sequences rep-seqs_T1.qza --o-denoising-stats denoising-stats.qza --p-n-threads 8

rm trimmed-seqs.qza demux-paired-end.qza denoising-stats.qza

qiime tools export --input-path rep-seqs_T1.qza --output-path PICRUST2_LEFSE
qiime tools export --input-path table_T1.qza --output-path PICRUST2_LEFSE

mv PICRUST2_LEFSE/dna-sequences.fasta PICRUST2_LEFSE/dna-sequences_T1.fasta
mv PICRUST2_LEFSE/feature-table.biom PICRUST2_LEFSE/feature-table_T1.biom

rm rep-seqs_T1.qza table_T1.qza

rm Only_T1 -R

conda deactivate



################################### PICRUST2-LEFSE #############################

cd PICRUST2_LEFSE

conda activate picrust2 ######### opne PICRUST2

picrust2_pipeline.py -s dna-sequences_T0.fasta -i feature-table_T0.biom -o picrust2_T0 -p 7 --verbose -e 0.0 -t sepp 

# the last two argument are for SEPP algorithm ("e" needed to avoid an error, see https://github-wiki-see.page/m/picrust/picrust2/wiki/# PICRUSt2-Tutorial-%28v2.3.0-beta%29#add-functional-descriptions)

add_descriptions.py -i picrust2_T0/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2_T0/pathways_out/path_abun_unstrat_descrip.tsv.gz


picrust2_pipeline.py -s dna-sequences_T1.fasta -i feature-table_T1.biom -o picrust2_T1 -p 7 --verbose -e 0.0 -t sepp 

add_descriptions.py -i picrust2_T1/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2_T1/pathways_out/path_abun_unstrat_descrip.tsv.gz

conda deactivate


R ########### open R Environment to prepare files for LEFSE (plot)

### for T0
a <- read.delim("picrust2_T0/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
if (!require('stringr')) {install.packages('stringr')}
colnames(a)[stringr::str_starts(colnames(a),pattern = "ID1535")]<-substring(colnames(a)[stringr::str_starts(colnames(a),pattern = "ID1535")],first = 16)

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]
Metadata <- as.data.frame(readxl::read_excel("../Metadata.xlsx"))
Metadata<-Metadata[Metadata$Time=="T0",]
Metadata$FASTQ_ID<-as.character(Metadata$FASTQ_ID)
rownames(Metadata)<-Metadata$FASTQ_ID
Metadata<-Metadata[colnames(a),]
identical(Metadata$FASTQ_ID,colnames(a))
a<-rbind.data.frame(Metadata$Condition,a)
head(a, n=2)
colnames(a)<-paste("M",colnames(a),sep="_")
a<-cbind.data.frame(c("Condition",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("Condition", Descriptions$pathway)']<-"Sample"
head(a)

write.table(a,file="Output_for_LefSe_T0.tsv",row.names = F,quote = F, sep="\t")

rm(a, Descriptions,Metadata)

### for T1
a <- read.delim("picrust2_T1/pathways_out/path_abun_unstrat_descrip.tsv.gz") 
if (!require('stringr')) {install.packages('stringr')}
colnames(a)[stringr::str_starts(colnames(a),pattern = "ID1535")]<-substring(colnames(a)[stringr::str_starts(colnames(a),pattern = "ID1535")],first = 16)

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]
Metadata <- as.data.frame(readxl::read_excel("../Metadata.xlsx"))
Metadata<-Metadata[Metadata$Time=="T1",]
Metadata$FASTQ_ID<-as.character(Metadata$FASTQ_ID)
rownames(Metadata)<-Metadata$FASTQ_ID
Metadata<-Metadata[colnames(a),]
identical(Metadata$FASTQ_ID,colnames(a))
a<-rbind.data.frame(Metadata$Condition,a)
head(a, n=2)
colnames(a)<-paste("M",colnames(a),sep="_")
a<-cbind.data.frame(c("Condition",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("Condition", Descriptions$pathway)']<-"Sample"
head(a)

write.table(a,file="Output_for_LefSe_T1.tsv",row.names = F,quote = F, sep="\t")

quit(save="no")


conda activate lefse  ######## now opening LEFSE

mkdir LEFSE
cd LEFSE

lefse_format_input.py ../'Output_for_LefSe_T0.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 # groups is second row, normalization to 1 milion (see Galaxy recommendation)

lefse_run.py formatted_table_for_lefse.in Risultato_WT_vs_KO_T0.res

rm formatted_table_for_lefse.in

lefse_format_input.py ../'Output_for_LefSe_T1.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000

lefse_run.py formatted_table_for_lefse.in Risultato_WT_vs_KO_T1.res

rm formatted_table_for_lefse.in

conda deactivate

# for plotting the results see the stand alone R script of this analysis
