# Operative System: Debian GNU/Linux 11 (bullseye)
# required: conda, qiime2-amplicon-2024.2

threads=25


########################## PROCESSING RAW FILES ############################

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

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores $threads --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --p-minimum-length 100 --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

cat Log_reads_with_primers.txt


##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

#qiime tools view demux.qzv
#echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 276 --p-trunc-len-r 196 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads $threads

qiime tools export --input-path denoising-stats.qza --output-path ./      ### to export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > DADA2_Initial_and_processed_reads.tsv
cat DADA2_Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv
rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv



######### CHECK OF HUMAN CONTAMINANTS WITH BOWTIE2 ...

# this part of the script is not reported, because there was NO HUMAN OFF-TARGETS CONTAMINANTS !

mv rep-seqs.qza clean_rep-seqs.qza




######################### TAXONOMIC CLASSIFICATION ##################################

SILVA_trained_db=~/SSD1/reference/SILVA_138_V3V4_BayesClassifier.qza

### BAYESIAN CLASSIFICATION (using Scikit-learn retrained on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier $SILVA_trained_db --i-reads clean_rep-seqs.qza --o-classification taxonomy.qza --p-n-jobs $threads

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization confidence_taxonomy.qzv




################# COMPUTING THE PHYLOGENETIC TREE DE NOVO #############################


#### DE NOVO TREE ALIGNMENT

qiime alignment mafft --i-sequences clean_rep-seqs.qza --o-alignment aligned-seqs.qza --p-n-threads $threads
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza



############# EXPORTING ORIGINAL NUMBER OF READS ############################

qiime tools export --input-path demux.qzv --output-path Raw_Reads_Info
mv Raw_Reads_Info/per-sample-fastq-counts.tsv Original_number_of_reads_for_sample.tsv
chmod +rw Original_number_of_reads_for_sample.tsv
chmod +rw Raw_Reads_Info -R
rm Raw_Reads_Info -R



conda deactivate   # finished! Continue with R analysis script ...
