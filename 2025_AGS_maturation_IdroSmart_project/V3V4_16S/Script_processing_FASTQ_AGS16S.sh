#!/usr/bin/env

# Operative System: Debian GNU/Linux 11 (bullseye)
# programs and envs required: conda, qiime2-amplicon-2024.2




#################### PREPARING REQUIRED VARIABLES ############################


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows loads the environment variable "into the script"

SILVA_trained_db=/media/matteo/SSD1/reference/SILVA_138_V3V4_BayesClassifier.qza
MIDAS_trained_db=/media/matteo/SSD1/reference/MIDAS_v5_3_Bayes_classifier_QIIME2.qza # optional

threads=32




########################### PROCESSING RAW FILES #############################

mkdir QIIME
chmod +rwx QIIME



##### FASTQ CHECKSUM

md5sum raw_FASTQ_dir/* | sed 's+raw_FASTQ_dir/++g' | sed 's/  / \t/' > QIIME/checksum_FASTQ.tsv

dir raw_FASTQ_dir/ | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt



##### OPENING QIIME ENVIRONMENT

conda activate qiime2-amplicon-2024.2

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza



##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNBGCWSCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores $threads --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --p-minimum-length 100 --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

cat Log_reads_with_primers.txt



##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

qiime tools view demux.qzv
#echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


# the reads are cutted allowing a possible overlap up to 20nt (see https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html)
qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 272 --p-trunc-len-r 170 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads $threads

qiime tools export --input-path denoising-stats.qza --output-path ./      ### to export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > DADA2_Initial_and_processed_reads.tsv
cat DADA2_Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv




######################### TAXONOMIC CLASSIFICATION ##################################

#### CLASSIFICATION THROUGH GLOBAL ALIGNMENT
#qiime feature-classifier classify-consensus-vsearch --i-query clean_rep-seqs.qza --i-reference-reads silva-138-99-seqs.qza --i-reference-taxonomy ~/Scrivania/Tools/Database/silva-138-99-tax.qza --p-perc-identity 0.99 --p-threads $threads --o-classification taxonomy_vsearch.qza --o-search-results vsearch_hits.qza --verbose 2> ASV_matched_in_SILVA.txt 
 
### BAYESIAN CLASSIFICATION (using Scikit-learn retrained on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier $SILVA_trained_db --i-reads rep-seqs.qza --o-classification taxonomy_SILVA.qza
qiime metadata tabulate --m-input-file taxonomy_SILVA.qza --o-visualization confidence_taxonomy_SILVA.qzv
# again but using MIDAS
qiime feature-classifier classify-sklearn --i-classifier $MIDAS_trained_db --i-reads rep-seqs.qza --o-classification taxonomy_MIDAS.qza
qiime metadata tabulate --m-input-file taxonomy_MIDAS.qza --o-visualization confidence_taxonomy_MIDAS.qzv



################# COMPUTING THE PHYLOGENETIC TREE DE NOVO #############################


#### DE NOVO TREE ALIGNMENT

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza




############# EXPORTING ORIGINAL NUMBER OF READS ############################

qiime demux summarize --i-data demux-paired-end.qza --o-visualization original_reads.qzv 
qiime tools export --input-path original_reads.qzv  --output-path Raw_Reads_Info
mv Raw_Reads_Info/per-sample-fastq-counts.tsv Original_number_of_reads_for_sample.tsv
chmod +rw Original_number_of_reads_for_sample.tsv
chmod +rw Raw_Reads_Info -R
rm Raw_Reads_Info -R


############# CLEANING THE FOLDER FROM QIIME TEMP FILES #####################

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv original_reads.qzv 





