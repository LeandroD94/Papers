######################## PROCESSING RAW FILES ########################

mkdir QIIME

##### FASTQ CHECKSUM

md5sum raw_FASTQ_dir/* | sed 's+raw_FASTQ_dir/++g' | sed 's/ /\t/' > QIIME/checksum_FASTQ.tsv

dir raw_FASTQ_dir/* | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt


##### OPENING QIIME ENVIRONMENT

conda activate qiime2-2021.4

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza

##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

rm Log_di_cutadapt.txt

##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 5 --p-trim-left-r 0 --p-trunc-len-f 274 --p-trunc-len-r 174 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 6

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza


#####################################################

 
#### TAXONOMIC ASSIGNMENT
 
 
### BAYESIAN CLASSIFICATION (using retrained RDP on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier ~/Scrivania/Tools/Databases/SILVA_138_V3V4_BayesClassifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv


#### DE NOVO TREE ALIGNMENT

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza


#### PREPARING FILES FOR PICRUST 

mkdir PICRUST2_LEFSE

qiime tools export --input-path rep-seqs.qza --output-path PICRUST2_LEFSE # it will export "dna-sequences.fasta"
qiime tools export --input-path table.qza --output-path PICRUST2_LEFSE # it will export "feature-table.biom"

conda deactivate


################################### PICRUST2-LEFSE #############################

cd PICRUST2_LEFSE

##### PICRUST2

conda activate picrust2

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 4 --verbose -e 0.0 -t sepp 

# the last two argument after verbose are for SEPP algorithm ("e" needed to avoid and error, see https://github-wiki-see.page/m/picrust/picrust2/wiki/# PICRUSt2-Tutorial-%28v2.3.0-beta%29#add-functional-descriptions)

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

conda deactivate


########## RENAMING FEATURES FOR LEFSE PLOT

R # open R

a <- read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 

colnames(a)[1:4]

if (!require('stringr')) {install.packages('stringr')}

temp<-substring(colnames(a)[stringr::str_starts(colnames(a),pattern = "X")],first = 10)
temp

colnames(a)[stringr::str_starts(colnames(a),pattern = "X")]<-temp

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]

Metadata <- as.data.frame(read.csv("../../Metadata.csv")) # import metadata file to regroups
head(Metadata)
Metadata$FASTQ_ID<-as.character(Metadata$FASTQ_ID)
rownames(Metadata)<-Metadata$FASTQ_ID
Metadata<-Metadata[colnames(a),]
if(identical(Metadata$FASTQ_ID,colnames(a))){a<-rbind.data.frame(Metadata$Mutation_Treatment,a)} else {cat ("\nsomething gone wrong...\n\n")}
a<-cbind.data.frame(c("Condition",Descriptions$pathway),a)
colnames(a)[colnames(a)=='c("Condition", Descriptions$pathway)']<-"Sample"
head(a,n=3)

write.table(a,file="Output_for_LefSe.tsv",row.names = F,quote = F, sep="\t")

rm(a, Descriptions, Metadata)

quit(save="no")


##### LEFSE

conda activate lefse

lefse_format_input.py 'Output_for_LefSe.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 # groups is second row, normalization to 1 milion (see Galaxy recommendation)

lefse_run.py formatted_table_for_lefse.in -y 1 Result_LEFSE.res # y is for multiclass analysis

rm formatted_table_for_lefse.in

conda deactivate

cd ..
