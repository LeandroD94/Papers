# operating system: Ubuntu 22.04 (Linux)
# required: conda, qiime2-2021.4, picrust-2.5, lefse-1.1.2, R-4.1


########################################### PROCESSING RAW FILES FOR ANALYSIS THROUGH R ####################################################

mkdir QIIME
chmod +rwx QIIME

##### FASTQ CHECKSUM

md5sum raw_EVERY_FASTQ_dir/* | sed 's+raw_EVERY_FASTQ_dir/++g' | sed 's/  / \t/' > QIIME/checksum_FASTQ.tsv

dir raw_EVERY_FASTQ_dir/ | sed '/R2_001/d' >R1.txt     ### to create the two columns reporting the FASTQ pair
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt


##### OPENING QIIME ENVIRONMENT

conda activate qiime2-2021.4

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_EVERY_FASTQ_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza


##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --p-minimum-length 100 --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt > Log_reads_with_primers.txt

cat Log_reads_with_primers.txt

##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 261 --p-trunc-len-r 184 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 7


qiime tools export --input-path denoising-stats.qza --output-path ./                ### to create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > Percentuals_Initial_and_processed_reads.tsv
cat Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv
rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv


######################### TAXONOMIC CLASSIFICATION ##################################

CLASSIFICATION THROUGH GLOBAL ALIGNMENT
qiime feature-classifier classify-consensus-vsearch --i-query rep-seqs.qza --i-reference-reads ~/Scrivania/Tools/Database/silva-138-99-seqs.qza --i-reference-taxonomy ~/Scrivania/Tools/Database/silva-138-99-tax.qza --p-perc-identity 0.99 --p-threads 8 --o-classification taxonomy_vsearch.qza --o-search.results vsearch_hits.qza --verbose 2> ASV_matched_in_SILVA.txt 



################# COMPUTING THE PHYLOGENETIC TREE DE NOVO #############################

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza





####################################################### PICRUST 2 v2.50 (STOOL) ############################################################


##### OPENING QIIME ENVIRONMENT TO PROCESS STOOL FASTQ ONLY

conda activate qiime2-2021.4

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Raw_FASTQ_Stool --input-format CasavaOneEightLanelessPerSampleDirFmt --output-path demux-paired-end_stool.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end_stool.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs_stool.qza --p-discard-untrimmed --p-minimum-length 100 --verbose

##### DENOISING, QUALITY LENGTH TRIMMING and MERGING (Same settings as the original one)
qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs_stool.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 261 --p-trunc-len-r 184 --o-table table_stool.qza --o-representative-sequences rep-seqs_stool.qza --o-denoising-stats denoising-stats_stool.qza --p-n-threads 7

rm demux-paired-end_stool.qza trimmed-seqs_stool.qza denoising-stats_stool.qza


qiime tools export --input-path rep-seqs_stool.qza --output-path PICRUST2_LEFSE_stool # it will export "dna-sequences.fasta"
qiime tools export --input-path table_stool.qza --output-path PICRUST2_LEFSE_stool # it will export "feature-table.biom"

rm rep-seqs_stool.qza  table_stool.qza




######## STARTING PICRUST2 STEP (Stool)

conda deactivate

conda activate picrust2

cd PICRUST2_LEFSE_stool # directory created through exporting

picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2 -p 4 --verbose -e 0.0 -t sepp 

# the last two argument after verbose are for SEPP algorithm ("e" is needed to avoid and error, see https://github.com/gavinmdouglas/q2-picrust2/issues/13)

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

conda deactivate

###### RENAMING FEATURES FOR LEFSE PLOT

R # open R

a <- read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 

colnames(a)[1:4]

if (!require('stringr')) {install.packages('stringr')}

# cleaning the fastq names in order to match them with the metadata column
temp<-gsub(".16S","",colnames(a), fixed=T)
temp<-substring(temp[stringr::str_starts(temp,pattern = "ID")],first=7, last=18)
temp<-gsub("_S[1-9]","",temp, fixed=F)
temp<-gsub("_S","",temp, fixed=T)
temp<-gsub("_","",temp, fixed=T)
temp<-gsub(".FAS.","-",temp, fixed=T)
temp<-gsub(".","",temp, fixed=T)

temp

colnames(a)[stringr::str_starts(colnames(a),pattern = "ID")]<-temp

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]

Metadata <- as.data.frame(read.table("../../metadata_R.tsv", sep="\t", header=T)) # import metadata file to regroups
head(Metadata)
Metadata<-Metadata[Metadata$Sampling_Type=="F",]
Metadata$FASTQ_code<-as.character(Metadata$FASTQ_code)
rownames(Metadata)<-Metadata$FASTQ_code
Metadata<-Metadata[colnames(a),]

if(identical(Metadata$FASTQ_code,colnames(a))) {a<-rbind.data.frame(Metadata$Tumor,a)} else {cat ("\nsomething gone wrong...\n\n")}
a<-cbind.data.frame(c("Condition",Descriptions$pathway),a) ### it adds the first (temp) column to the data set, which each element identify the corrispondent row
colnames(a)[colnames(a)=='c("Condition", Descriptions$pathway)']<-"Sample"
head(a,n=3)


write.table(a,file="Output_for_LefSe.tsv",row.names = F,quote = F, sep="\t")

rm(a, Descriptions, Metadata)

quit(save="no")


########## LEFSE (Stool)

conda activate lefse

lefse_format_input.py 'Output_for_LefSe.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 # groups is second row, normalization to 1 milion (see Galaxy recommendation)

lefse_run.py formatted_table_for_lefse.in Result_LEFSE.res

rm formatted_table_for_lefse.in

conda deactivate



cd ..   # back to QIIME directory




####################################################### PICRUST 2 v2.50 (BIOPSY) ############################################################


##### OPENING QIIME ENVIRONMENT TO PROCESS BIOPSY FASTQ ONLY

conda activate qiime2-2021.4

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Raw_FASTQ_Biopsy --input-format CasavaOneEightLanelessPerSampleDirFmt --output-path demux-paired-end_biopsy.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end_biopsy.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs_biopsy.qza --p-discard-untrimmed --p-minimum-length 100 --verbose

##### DENOISING, QUALITY LENGTH TRIMMING and MERGING (Same settings as the original one)
qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs_biopsy.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 261 --p-trunc-len-r 184 --o-table table_biopsy.qza --o-representative-sequences rep-seqs_biopsy.qza --o-denoising-stats denoising-stats_biopsy.qza --p-n-threads 7

rm demux-paired-end_biopsy.qza trimmed-seqs_biopsy.qza denoising-stats_biopsy.qza


qiime tools export --input-path rep-seqs_biopsy.qza --output-path PICRUST2_LEFSE_biopsy # it will export "dna-sequences.fasta"
qiime tools export --input-path table_biopsy.qza --output-path PICRUST2_LEFSE_biopsy # it will export "feature-table.biom"

rm rep-seqs_biopsy.qza  table_biopsy.qza



##### DECONTAMINATION FROM HUMAN ASV (To get prediction based on the TRUE bacterial asv)

cd PICRUST2_LEFSE_biopsy # directory created through exporting

# to create the list of contaminant ASVs
qiime tools export --input-path ../assigned_taxa.qza --output-path ./   ### the ASV codes match between the two rep-seq files (checked)
cat taxonomy.tsv | grep "Unassigned" | cut -f 1 > human_ASV_biopsy.tsv

cat dna-sequences.fasta | grep -A1 --no-group-separator -f human_ASV_biopsy.tsv > to_delete.txt
cat dna-sequences.fasta | grep -v -f to_delete.txt > dna-sequences_decontaminated.fasta    ### NB: two steps because it's impossible to use -A1 and -v together (using just -v will join two consequent ASVs!)
rm to_delete.txt


# just to check the correctness of decontamination
A=$(grep -c '>' to_delete.txt)
B=$(grep -c '>' dna-sequences_decontaminated.fasta)
C=$(grep -c '>' dna-sequences.fasta)
if [ $C = $(($B + $A)) ]; then echo  -e "\n\n\n  The decontamination was successfull!   \n\n\n"; else echo -e "\n\n\n\n   Something gone wrong, check the removed reads...   \n\n\n" ; fi
unset A B C


biom convert -i feature-table.biom -o feature-table.tsv --to-tsv ### the package biom convert is included in qiime enviroment it self
cat feature-table.tsv | grep '# Constructed from biom file' -v | sed 's/#OTU ID/OTU_ID/' > feature-table2.tsv   ### cleaning this output from extra headers added by the command, NB: different output name!

cat feature-table2.tsv | grep -v -f human_ASV_biopsy.tsv > feature_table_decontaminated.tsv

rm feature-table.tsv feature-table2.tsv taxonomy.tsv to_delete.txt



######## STARTING PICRUST2 STEP (Biopsy)

conda deactivate

conda activate picrust2

picrust2_pipeline.py -s dna-sequences_decontaminated.fasta -i feature_table_decontaminated.tsv -o picrust2 -p 4 --verbose -e 0.0 -t sepp   ### PICRUST it self checks the matching of ASV between the files

# the last two argument after verbose are for SEPP algorithm ("e" is needed to avoid and error, see https://github.com/gavinmdouglas/q2-picrust2/issues/13)

add_descriptions.py -i picrust2/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz

conda deactivate

###### RENAMING FEATURES FOR LEFSE PLOT

R ### open R

a <- read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz") 

colnames(a)[1:4]

if (!require('stringr')) {install.packages('stringr')}

# cleaning the fastq names in order to match them with the metadata column
temp<-gsub(".16S","",colnames(a), fixed=T)
temp<-substring(temp[stringr::str_starts(temp,pattern = "ID")],first=7, last=18)
temp<-gsub("_S[1-9]","",temp, fixed=F)
temp<-gsub("_S","",temp, fixed=T)
temp<-gsub("_","",temp, fixed=T)
temp<-gsub(".FAS.","-",temp, fixed=T)
temp<-gsub(".","",temp, fixed=T)

temp

colnames(a)[stringr::str_starts(colnames(a),pattern = "ID")]<-temp

Descriptions<-a[,c("pathway","description")]
a<-a[ ,! colnames(a) %in% c("pathway","description")]

Metadata <- as.data.frame(read.table("../../metadata_R.tsv", sep="\t", header=T)) # import metadata file to regroups
head(Metadata)
Metadata<-Metadata[Metadata$Sampling_Type=="B",]
Metadata$FASTQ_code<-as.character(Metadata$FASTQ_code)
rownames(Metadata)<-Metadata$FASTQ_code
Metadata<-Metadata[colnames(a),]

if(identical(Metadata$FASTQ_code,colnames(a))) {a<-rbind.data.frame(Metadata$Tumor,a)} else {cat ("\nsomething gone wrong...\n\n")}
a<-cbind.data.frame(c("Condition",Descriptions$pathway),a) ### it adds the first (temp) column to the data set, which each element identify the corrispondent row
colnames(a)[colnames(a)=='c("Condition", Descriptions$pathway)']<-"Sample"
head(a,n=3)


write.table(a,file="Output_for_LefSe.tsv",row.names = F,quote = F, sep="\t")

rm(a, Descriptions, Metadata)

quit(save="no")


########## LEFSE (Biopsy)

conda activate lefse

lefse_format_input.py 'Output_for_LefSe.tsv' formatted_table_for_lefse.in -c 2 -u 1 -o 1000000 # groups is second row, normalization to 1 milion (see Galaxy recommendation)

lefse_run.py formatted_table_for_lefse.in Result_LEFSE.res

rm formatted_table_for_lefse.in

conda deactivate

cd ..   # back to QIIME directory
