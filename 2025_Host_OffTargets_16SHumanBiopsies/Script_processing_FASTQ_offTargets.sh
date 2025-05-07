#!/usr/bin/env

# version 14/01/2025


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh


############## PREPARING REQUIRED VARIABLES ###################

Bowtie2db=/media/matteo/SSD1/reference/Bwt2_indexed_Human_chm13_v2/chm13v2.0

SILVA_trained_db=/media/matteo/SSD1/reference/SILVA_138_V3V4_BayesClassifier.qza

Threads=60



### The following variable contains the works to process (NB: each omonimous folder contains the samples of its work!)
# those folders will be searched in the ACTUAL working dir

targets=$(ls raw_FASTQ_collection_in_folders/ )


################# PROCESSING RAW FILES ########################

mkdir QIIME
chmod +rwx QIIME
cd QIIME   # NB: operating in a new working dir now!

mkdir Files_to_check_quality

conda activate qiime2-amplicon-2024.2  # opening QIIME environment


### WORKING ON EACH BATCH ( IMPORT, PRIMERS CHECK AND REMOVING UNTRIMMED, QUALITY SUMMARY )

for x in $targets
do

echo -e "\n\n ... Working on $x project ... \n\n\n"

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_collection_in_folders/$x --input-format CasavaOneEightLanelessPerSampleDirFmt --output-path demux-paired-end_$x.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end_$x.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores $Threads --o-trimmed-sequences trimmed-seqs_$x.qza --p-discard-untrimmed --p-minimum-length 150 --verbose >Log_di_cutadapt_$x.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt_$x.txt >Log_reads_with_primers_$x.txt

cat Log_reads_with_primers_$x.txt

qiime demux summarize --i-data trimmed-seqs_$x.qza --o-visualization demux_$x.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"

done


# ordering the files
mv *qzv  Files_to_check_quality/
mv Log_* Files_to_check_quality/



#### ANNOTATING THE ORIGINAL NUMBERS OF READS

for x in $targets
do

qiime demux summarize --i-data demux-paired-end_$x.qza --o-visualization original_reads.qzv
qiime tools export --input-path original_reads.qzv  --output-path Raw_Reads_Info
mv Raw_Reads_Info/per-sample-fastq-counts.tsv Original_number_of_reads_$x.qza.tsv

chmod +rw Original_number_of_reads_$x.qza.tsv
rm original_reads.qzv
chmod +rw Raw_Reads_Info -R
rm Raw_Reads_Info -R

done


####### GENERATING THE ASVs AND TABLES THROUGH DADA2 ######

for x in $targets
do

echo "Found $x ..."


##### Setting the optimal cut positions for each batch ...

if [ $x = 'PRJEB46353' ] ; then F=260 ; R=184 ; fi

if [ $x = 'PRJEB57580' ] ; then F=255 ; R=215 ; fi

if [ $x = 'PRJNA398187' ] ; then F=251 ; R=192 ; fi

if [ $x = 'PRJNA743150' ] ; then F=255 ; R=188 ; fi

#if [ $x = 'GSE217490' ] ; then F=261 ; R=183 ; fi
if [ $x = 'PRJNA899104' ] ; then F=261 ; R=183 ; fi

if [ $x = 'PRJNA507548' ] ; then F=256 ; R=187 ; fi

if [ $x = 'PRJNA995580' ] ; then F=275 ; R=188  ; fi

if [ $x = 'PRJNA840857' ] ; then F=233 ; R=210  ; fi

if [ $x = 'PRJNA298957' ] ; then F=256 ; R=187  ; fi

if [ $x = 'PRJNA325650' ] ; then F=230 ; R=212  ; fi


##### Looping DADA2 over the batches

echo -e "\n ... DADA2: $x ... \n\n"

qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs_$x.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f $F --p-trunc-len-r $R --o-table table_$x.qza --o-representative-sequences rep-seqs_$x.qza --o-denoising-stats denoising-stats_$x.qza --p-n-threads $Threads

qiime tools export --input-path denoising-stats_$x.qza --output-path ./     # to create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > Percentuals_Initial_and_processed_reads_$x.tsv
cat Percentuals_Initial_and_processed_reads_$x.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9

chmod +rwx stats.tsv
rm denoising-stats_$x.qza trimmed-seqs_$x.qza demux-paired-end_$x.qza stats.tsv

done




################ REMOVING EUCARYOTIC ASV #######################

# re-ordering the folders
mkdir rep_seqs_with_offtargets
mv rep-seqs* rep_seqs_with_offtargets/
mkdir rep_seqs_decontam


for x in $targets
do

qiime tools export --input-path rep_seqs_with_offtargets/rep-seqs_$x.qza --output-path ./
mv dna-sequences.fasta dna-sequences_$x.fasta

done


conda deactivate   # temporarily, to use bowtie2...



### Bowtie2
for x in $targets
do

bowtie2 -x $Bowtie2db -U dna-sequences_$x.fasta --al rep_seqs_decontam/eucar_contaminant_ASV_$x.fasta -S rep_seqs_decontam/bowtie2_SAM_$x.txt -f -p $Threads 2> rep_seqs_decontam/percentuals_aligned_bowtie2_$x.txt

grep -Fwv -f rep_seqs_decontam/eucar_contaminant_ASV_$x.fasta dna-sequences_$x.fasta > rep_seqs_decontam/decontaminated_seq_$x.fasta

echo -e " Number of human contaminant ASV sequences:" $(grep ">" rep_seqs_decontam/eucar_contaminant_ASV_$x.fasta | wc -l ) "of" $(grep ">" rep_seqs_decontam/decontaminated_seq_$x.fasta | wc -l) "procaryotic ASVs" > rep_seqs_decontam/Percent_of_eucar_Contaminant_$x.txt

chmod +rwx dna-sequences_$x.fasta
rm dna-sequences_$x.fasta

done



conda activate qiime2-amplicon-2024.2

for x in $targets
do
qiime tools import --input-path rep_seqs_decontam/decontaminated_seq_$x.fasta --output-path rep_seqs_decontam/clean_rep-seqs_$x.qza --type FeatureData[Sequence]
rm rep_seqs_decontam/decontaminated_seq_$x.fasta
done



####################### TAXONOMIC CLASSIFICATION ###########################

for x in $targets
do

qiime feature-classifier classify-sklearn --i-classifier $SILVA_trained_db --i-reads rep_seqs_decontam/clean_rep-seqs_$x.qza --o-classification rep_seqs_decontam/taxonomy_$x.qza --p-n-jobs $Threads --p-read-orientation 'same'

qiime metadata tabulate --m-input-file rep_seqs_decontam/taxonomy_$x.qza --o-visualization rep_seqs_decontam/confidence_taxonomy_$x.qzv


### BAYESIAN CLASSIFICATION ALSO ON THE HOST CONTAMINATED OBJECT
qiime feature-classifier classify-sklearn --i-classifier $SILVA_trained_db --i-reads rep_seqs_with_offtargets/rep-seqs_$x.qza --o-classification rep_seqs_with_offtargets/contam_taxonomy_$x.qza --p-n-jobs $Threads  --p-read-orientation 'same'

qiime metadata tabulate --m-input-file rep_seqs_with_offtargets/contam_taxonomy_$x.qza --o-visualization rep_seqs_with_offtargets/confidence_taxonomy_$x.qzv

done


### NB: without '--p-read-orientation same' the classifier will study by it self the first 100 reads of each object and if they are humans then it will think that the reads are reverse/complementary (the default is 'auto') --> reverse orientation --> most of the read unassigned at kingdom or phylum level in same work (e.g. GEO work) !!!

echo "The presence of ASV caused problems during the taxonomic identification, see the notes in the scripts"  > ../Warning_Problem_ident_caused_by_human.txt




############### COMPUTING THE PHYLOGENETIC TREE DE NOVO ######################

for x in $targets
do

qiime alignment mafft --i-sequences rep_seqs_decontam/clean_rep-seqs_$x.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rep_seqs_decontam/rooted-tree_$x.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza


### FOR THE CONTAMINATED SEQS TOO (TO PERFORM FURTHER CHECKS)
qiime alignment mafft --i-sequences rep_seqs_with_offtargets/rep-seqs_$x.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rep_seqs_with_offtargets/rooted-tree_with_contam_$x.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza

done


cd ..
