# operating system: Ubuntu 22.04 (Linux)
# programs required: conda, qiime2-2022.8, bowtie-2.2.5

############################## PREPARATION STEPS ##################################

mkdir QIIME
chmod +rwx QIIME

BOWTIE2_DATABASE='~/Scrivania/Tools/Database'   # path where the bowtie2 databases are stored



##### FASTQ CHECKSUM

md5sum raw_FASTQ_dir/* | sed 's+raw_FASTQ_dir/++g' | sed 's/  / \t/' > QIIME/checksum_FASTQ.tsv

dir raw_FASTQ_dir/ | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt



###################### COUNTING THE HOST OFF TARGET SEQUENCES  ######################

# This step is performed now on raw data and not after on the ASV to ensure that the off-targets are not lost during the quality filter alongside OTHERS reads (counting them now allows a more precise information)


mkdir QIIME/Host_off_target_counts


conda activate bowtie2


##### HUMAN OFF TARGET 
for i in raw_FASTQ_dir/* ;
do

echo -e "\n\nProcessing the sample ${i/'raw_FASTQ_dir/'/} \n\n"

bowtie2 -x $BOWTIE2_DATABASE/INDEXED_GRCh38_analysis_set_BOWTIE2/GRCh38_noalt_as -U $i --al-gz QIIME/Host_off_target_counts/${i/'raw_FASTQ_dir/'/} -S QIIME/Host_off_target_counts/SAM_${i/'raw_FASTQ_dir/'/}.sam -p 4 2> QIIME/Host_off_target_counts/${i/'raw_FASTQ_dir/'/}_percentuals_aligned.txt

echo ${i/'raw_FASTQ_dir/'/} "," $(grep "overall" QIIME/Host_off_target_counts/${i/'raw_FASTQ_dir/'/}_percentuals_aligned.txt) >> QIIME/Host_off_target_counts/HUMAN_off_target_for_each_FASTQ.csv

done

#cat QIIME/Host_off_target_counts/HUMAN_off_target_for_each_FASTQ.csv



### MICE OFF TARGET (TWO SAMPLES)
for i in $(ls raw_FASTQ_dir/ | grep -e "LDG[5-6]")
do

echo -e "\n\nProcessing the sample $i \n\n"

rm QIIME/Host_off_target_counts/$i   # those are from the HUMAN alignment perfomed before (they are not human)

bowtie2 -x $BOWTIE2_DATABASE/INDEXED_mouse_GRCm39_BOWTIE2/GRCm39 -U raw_FASTQ_dir/$i --al-gz QIIME/Host_off_target_counts/$i -S QIIME/Host_off_target_counts/SAM_${i/'.fastq.gz'/'.sam'} -p 4 2> QIIME/Host_off_target_counts/${i/'.fastq.gz'/'_percentuals_aligned.txt'}

echo $i "," $( grep "overall"  QIIME/Host_off_target_counts/${i/'.fastq.gz'/'_percentuals_aligned.txt'} ) >> QIIME/Host_off_target_counts/MICE_off_target_for_each_FASTQ.csv

done

#cat QIIME/Host_off_target_counts/MICE_off_target_for_each_FASTQ.csv



conda deactivate



########################## PROCESSING THE RAW READS #############################


##### OPENING QIIME ENVIRONMENT

conda activate qiime2-2022.8

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../raw_FASTQ_dir --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza


##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores 7 --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --p-minimum-length 100 --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

cat Log_reads_with_primers.txt


##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 # of total pool!

echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 268 --p-trunc-len-r 175 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 7

qiime tools export --input-path denoising-stats.qza --output-path ./      ### to export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > DADA2_Initial_and_processed_reads.tsv
cat DADA2_Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv



####################### REMOVING EUCARYOTIC ASV (CONTAMINANTS) ############################

qiime tools export --input-path rep-seqs.qza --output-path ./ # it will export "dna-sequences.fasta"

mv rep-seqs.qza eucar_contaminants_ASV/rep-seqs_with_euc_contam.qza


conda deactivate

mkdir eucar_contaminants_ASV

conda activate bowtie2


####### FOR HUMAN DNA (BLOOD SAMPLES)
bowtie2 -x $BOWTIE2_DATABASE/INDEXED_GRCh38_analysis_set_BOWTIE2/GRCh38_noalt_as -U dna-sequences.fasta --al eucar_contaminant_ASV_1.fasta -S eucar_contaminants_ASV/bowtie2_HUMAN_SAM_output.txt -f -p 4 2> eucar_contaminants_ASV/percentuals_aligned_HUMAN_bowtie2.txt

####### FOR MICE DNA (TRAP SAMPLES)
bowtie2 -x $BOWTIE2_DATABASE/INDEXED_mouse_GRCm39_BOWTIE2/GRCm39 -U dna-sequences.fasta --al eucar_contaminant_ASV_2.fasta -S eucar_contaminants_ASV/bowtie2_MICE_SAM_output.txt -f -p 4 2> eucar_contaminants_ASV/percentuals_aligned_MICE_bowtie2.txt


cat eucar_contaminant_ASV_1.fasta eucar_contaminant_ASV_2.fasta > Every_OFF_TARGETS_in_ASVs.txt
grep -Fwv -f Every_OFF_TARGETS_in_ASVs.txt dna-sequences.fasta > decontaminated_seq.fasta  # removing both the off-targets type from the ASV

mv eucar_contaminant_ASV_1.fasta eucar_contaminants_ASV/HUMAN_OFF_TARGETS_in_ASVs.fasta
mv eucar_contaminant_ASV_2.fasta eucar_contaminants_ASV/MICE_OFF_TARGETS_in_ASVs.fasta
mv Every_OFF_TARGETS_in_ASVs.txt eucar_contaminants_ASV/Every_OFF_TARGETS_in_ASVs.txt


conda deactivate 
chmod +rwx dna-sequences.fasta
rm dna-sequences.fasta


conda activate qiime2-2022.8

qiime tools import --input-path decontaminated_seq.fasta --output-path clean_rep-seqs.qza --type FeatureData[Sequence]



############################ TAXONOMIC CLASSIFICATION #####################################


#### CLASSIFICATION THROUGH GLOBAL ALIGNMENT
qiime feature-classifier classify-consensus-vsearch --i-query clean_rep-seqs.qza --i-reference-reads ~/Scrivania/Tools/Database/silva-138-99-seqs.qza --i-reference-taxonomy ~/Scrivania/Tools/Database/silva-138-99-tax.qza --p-perc-identity 0.99 --p-threads 8 --o-classification taxonomy_vsearch.qza --o-search-results vsearch_hits.qza --verbose 2> ASVpercentage_matched_in_SILVA_using_VSearch.txt




##################### COMPUTING THE PHYLOGENETIC TREE DE NOVO ################################


#### DE NOVO TREE ALIGNMENT

qiime alignment mafft --i-sequences clean_rep-seqs.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza


#### FOR THE CONTAMINATED SEQS TOO (TO PERFORM FURTHER CHECKS)

qiime alignment mafft --i-sequences eucar_contaminants/rep-seqs_with_euc_contam.qza --o-alignment aligned-seqs.qza
qiime alignment mask --i-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree eucar_contaminants/rooted-tree_with_contam.qza

rm aligned-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza




###################### EXPORTING ORIGINAL NUMBER OF READS ###################################

qiime demux summarize --i-data demux-paired-end.qza --o-visualization original_reads.qzv 
qiime tools export --input-path original_reads.qzv  --output-path Raw_Reads_Info
mv Raw_Reads_Info/per-sample-fastq-counts.tsv Original_number_of_reads_for_sample.tsv
chmod +rw Original_number_of_reads_for_sample.tsv
chmod +rw Raw_Reads_Info -R
rm Raw_Reads_Info -R



################### CLEANING THE FOLDER FROM QIIME TEMP FILES ###########################

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv original_reads.qzv 

conda deactivate
